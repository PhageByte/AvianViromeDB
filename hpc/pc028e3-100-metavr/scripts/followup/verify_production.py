"""Recount the returned production tables and join exact IMG dataset provenance.

Linear scans over source metadata and matched pairs; sets deduplicate each counting
unit. No cache. All outputs belong to the consuming analysis step's --out directory.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_table(path, records, fields):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def ecosystem_group(value):
    labels = {v.strip() for v in value.split(",")}
    known = labels & {"Host-associated", "Environmental", "Engineered"}
    if len(known) == 1:
        return next(iter(known))
    return "Mixed ecosystems" if known else "Unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("production", "lzani", "query-flags", "sequence-index", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results = args.production / "results"
    summary = json.loads((results / "summary.json").read_text())
    query_rows = list(rows(results / "query-status.tsv.gz"))
    query = {r["query_id"]: r for r in query_rows}
    index = {r["sequence_id"]: r for r in rows(args.sequence_index) if r["catalog"] == "AVIAN"}
    flags = list(rows(args.query_flags))
    assert len(query_rows) == len(query) == len(index) == len(flags) == 61608
    assert len({r["vOTU_id"] for r in query_rows}) == 61608
    assert set(query) == set(index) == {r["query_id"] for r in flags}
    assert all(query[q]["vOTU_id"] == index[q]["vOTU_id"] and
               query[q]["length"] == index[q]["length"] for q in query)
    three = {r["query_id"] for r in flags if r["matched_any"].lower() == "true"}
    assert three == {q for q, r in query.items() if r["matched_three_catalogs"] == "True"}
    blocks = list(rows(results / "block-status.tsv"))
    assert len(blocks) == 235 and len({r["block_id"] for r in blocks}) == 235
    assert all(r["state"] == "done" for r in blocks)
    assert summary["complete"] and summary["block_states"] == {"done": 235}

    references = list(rows(results / "matched-references.tsv.gz"))
    ref = {r["reference_id"]: r for r in references}
    assert len(ref) == len(references) == summary["metavr_sequences"]
    wanted = {r["taxon_oid"] for r in references}
    source_fields = ["IMG Genome ID", "Study Name", "Genome Name / Sample Name", "Ecosystem",
                     "Ecosystem Category", "Ecosystem Type", "Ecosystem Subtype", "Specific Ecosystem"]
    metadata = {}
    for r in rows(args.production / "Source_dataset_metadata.tsv.gz"):
        key = r["IMG Genome ID"]
        if key in wanted:
            value = {k: r[k] for k in source_fields}
            assert key not in metadata or metadata[key] == value, key
            metadata[key] = value
    group = {key: ecosystem_group(metadata.get(key, {}).get("Ecosystem", "")) for key in wanted}
    grouped = defaultdict(lambda: defaultdict(set))
    dataset = defaultdict(lambda: defaultdict(set))
    quality = defaultdict(set)
    pair_keys, qrefs, rqueries = set(), defaultdict(set), defaultdict(set)
    totals = Counter()
    matched_query_flags = defaultdict(set)
    for r in rows(results / "matched-pairs.tsv.gz"):
        q, t, taxon = r["query_id"], r["reference_id"], r["taxon_oid"]
        assert (q, t) not in pair_keys
        pair_keys.add((q, t))
        assert q in query and t in ref and ref[t]["taxon_oid"] == taxon
        assert r["avian_votu"] == query[q]["vOTU_id"] and r["uvig"] == ref[t]["uvig"]
        assert r["votu"] == ref[t]["votu"] and r["matched"] == "True"
        assert math.isclose(float(r["ani_percent"]),
                            100 * int(r["nident_sum"]) / int(r["alignment_length_sum"]),
                            rel_tol=1e-12)
        assert float(r["ani_percent"]) >= 95 and float(r["shorter_af_percent"]) >= 85
        if r["query_length"] == r["reference_length"]:
            assert float(r["query_af_percent"]) >= 85 and float(r["reference_af_percent"]) >= 85
        qrefs[q].add(t)
        rqueries[t].add(q)
        totals["both_af"] += r["both_af_ge_85"] == "True"
        totals["short_reference"] += r["short_reference_in_query"] == "True"
        totals["stacking"] += float(r["hsp_overlap_fold_value"]) >= 10
        matched_query_flags[q].add(group[taxon])
        quality[r["reference_quality"]].add(q)
        for bucket in (grouped[group[taxon]], dataset[taxon]):
            for unit, value in (("queries", q), ("references", t), ("uvigs", r["uvig"]),
                                ("reference_votus", r["votu"])):
                bucket[unit].add(value)
    matched = set(qrefs)
    new_queries = matched - three
    assert len(pair_keys) == summary["matched_pairs"]
    assert len(matched) == summary["avian_queries"] == summary["avian_votus"]
    assert len({r["uvig"] for r in references}) == summary["metavr_uvigs"]
    assert len({r["votu"] for r in references}) == summary["metavr_votus"]
    assert len(three | matched) == summary["union_with_metavr_queries"]
    assert len(new_queries) == summary["new_queries_beyond_three_catalogs"]
    assert totals["both_af"] == summary["matched_pairs_both_af_ge_85"]
    assert totals["short_reference"] == summary["matched_pairs_short_reference_in_query"]
    assert totals["stacking"] == summary["matched_pairs_hsp_stacking_fold_ge_10"]
    for q, r in query.items():
        assert int(r["matched_references"]) == len(qrefs[q])
        assert (r["route_c_status"] == "matched") == (q in matched)
        assert (r["new_beyond_three_catalogs"] == "True") == (q in new_queries)
    assert all(int(r["matched_avian_queries"]) == len(rqueries[t]) for t, r in ref.items())

    def counts(bucket):
        return {"matched_queries": len(bucket["queries"]),
                "new_beyond_three_catalogs": len(bucket["queries"] - three),
                **{unit: len(bucket[unit]) for unit in ("references", "uvigs", "reference_votus")}}

    eco_rows = [{"ecosystem_group": name, **counts(bucket)} for name, bucket in sorted(grouped.items())]
    dataset_rows = [{"taxon_oid": key, "ecosystem_group": group[key], **counts(bucket),
                     **{k: metadata.get(key, {}).get(k, "") for k in source_fields[1:]}}
                    for key, bucket in sorted(dataset.items(), key=lambda item: -len(item[1]["references"]))]
    write_table(args.out / "ecosystem-summary.tsv", eco_rows, list(eco_rows[0]))
    write_table(args.out / "dataset-summary.tsv", dataset_rows, list(dataset_rows[0]))
    flag_rows = [{**r, "ecosystem_groups": ";".join(sorted(matched_query_flags[q]))}
                 for q, r in query.items()]
    write_table(args.out / "query-ecosystem-flags.tsv.gz", flag_rows, list(flag_rows[0]))
    ref_rows = [{**r, "ecosystem_group": group[r["taxon_oid"]],
                 **{k: metadata.get(r["taxon_oid"], {}).get(k, "") for k in source_fields[1:]}}
                for r in references]
    write_table(args.out / "matched-reference-sources.tsv.gz", ref_rows, list(ref_rows[0]))

    slice_tables = []
    for filename in ("route-a-pairs.tsv.gz", "production-route-c-pairs.tsv.gz"):
        selected = {(r["query_id"], r["reference_id"]): r for r in
                    rows(args.production / "validation-slice" / filename) if r["matched"] == "True"}
        slice_tables.append(selected)
    a, c = slice_tables
    assert len(a) == 2354 and len(c) == 2352 and c.keys() <= a.keys() and c.keys() <= pair_keys
    metric_fields = ["ani_percent", "query_af_percent", "reference_af_percent",
                     "shorter_af_percent", "both_af_ge_85"]
    assert all(a[key][k] == c[key][k] for key in c for k in metric_fields)
    lost = set(a) - set(c)
    assert all(q in three for q, _ in lost)
    lz_summary = json.loads((args.lzani / "lzani-summary.json").read_text())
    lz_counts = {}
    for block, expected in lz_summary["blocks"].items():
        path = args.lzani / (block + "-seed") / "matched-either.tsv.gz"
        items = list(rows(path)) if path.exists() else []
        blast = {(r["query_id"], r["reference_id"]) for r in items if r["blast_matched"] == "True"}
        lz = {(r["query_id"], r["reference_id"]) for r in items if r["lz_matched"] == "True"}
        assert blast <= pair_keys
        assert len(blast) == expected["megablast_units_in_scope"]["matched_pairs"]
        assert len(lz) == expected["lzani_units"]["matched_pairs"]
        assert len(blast & lz) == expected["shared_units"]["matched_pairs"]
        lz_counts[block] = {"megablast": len(blast), "lzani": len(lz), "shared": len(blast & lz)}
    checked_inputs = [Path(__file__).resolve(), args.query_flags, args.sequence_index,
                      args.production / "Source_dataset_metadata.tsv.gz",
                      results / "matched-pairs.tsv.gz", results / "query-status.tsv.gz",
                      args.production / "SHA256SUMS", args.lzani / "SHA256SUMS"]
    report = {"all_passed": True, "source_input_sha256": {str(p): sha(p) for p in checked_inputs},
              "catalog_representatives": len(query), "matched_pairs": len(pair_keys),
              "matched_queries": len(matched), "three_catalog_union": len(three),
              "union_with_metavr": len(three | matched), "new_beyond_three_catalogs": len(matched - three),
              "ecosystem_counts": eco_rows, "largest_datasets": dataset_rows[:10],
              "datasets_without_source_metadata": len(wanted - metadata.keys()),
              "queries_by_reference_quality": {k: {"matched": len(v), "new": len(v - three)}
                                             for k, v in quality.items()},
              "validation_slice": {"route_a": len(a), "route_c": len(c), "lost_pairs": sorted(lost),
                                   "shared_metrics_identical": True, "lost_queries_in_three_catalogs": True},
              "lzani_pilot": lz_counts,
              "counting_scope": "Queries are deduplicated within each ecosystem; ecosystem counts overlap. Mixed, Host-associated follows the explicit Host-associated label. Source metadata denotes sampling context, not viral host range.",
              "validation_scope": "Recomputed returned matched-pair, query, reference, validation-slice and LZ-ANI tables; raw alignments and full-library recall were not re-evaluated."}
    (args.out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("all_passed", "matched_pairs", "union_with_metavr",
                                           "ecosystem_counts", "datasets_without_source_metadata")}, indent=2))


if __name__ == "__main__":
    main()
