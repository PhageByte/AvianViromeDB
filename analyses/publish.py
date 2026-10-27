# /// script
# requires-python = ">=3.13"
# dependencies = ["pandas>=2.2,<3"]
# ///
"""Export a consistent local release from the curated catalog.

Run with `uv run --script publish.py --output-dir <step-owned directory>`.
Input databases and source files are read-only. The destination must be new;
it becomes visible after validation of the complete staged package.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import tarfile
import tempfile

import pandas as pd


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_database(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_rows(actual, expected, columns):
    """Compare full row multisets, retaining duplicate rows and missing values."""
    pd.testing.assert_frame_equal(
        actual[columns].sort_values(columns).reset_index(drop=True),
        expected[columns].sort_values(columns).reset_index(drop=True),
        check_dtype=False, check_exact=True,
    )


def read_catalog(conn):
    catalog = pd.read_sql_query("SELECT * FROM contig_annotation_filtered_vOTU", conn)
    require("original_representative_id" in catalog, "Rebuild representatives with step 10 first")
    members = pd.read_sql_query("SELECT * FROM contig_annotation_filtered", conn)
    require(not catalog.vOTU_id.duplicated().any(), "Duplicate vOTU identifiers")
    require(not catalog.contig_repid.duplicated().any(), "Duplicate representative identifiers")
    require(not catalog.original_representative_id.duplicated().any(), "Duplicate original representative identifiers")
    require(not members.contig_id.duplicated().any(), "Duplicate member identifiers")
    selected = members.sort_values(
        ["checkv_completeness", "contig_id"], ascending=[False, True], kind="stable"
    ).drop_duplicates("vOTU_id").set_index("vOTU_id")
    actual = catalog.set_index("vOTU_id").sort_index()
    selected = selected.sort_index()
    require(actual.index.equals(selected.index), "Catalog groups differ from selected members")
    require(actual.contig_repid.equals(selected.contig_id), "Representative selection differs")
    require(actual.original_representative_id.equals(selected.contig_repid), "Original mapping differs")
    fields = [c for c in actual if c in selected and c != "contig_repid" and not c.startswith("host_")]
    pd.testing.assert_frame_equal(actual[fields], selected[fields], check_dtype=False, check_exact=True)
    return catalog.sort_values("vOTU_id").reset_index(drop=True)


def extract_fasta(source, ids_path, output):
    # seqkit matches complete IDs by default; no regular-expression matching.
    subprocess.run(["seqkit", "grep", "-j", "2", "-f", str(ids_path),
                    str(source), "-o", str(output)], check=True)


def fasta_lengths(path):
    result = subprocess.run(["seqkit", "fx2tab", "-j", "2", "-n", "-i", "-l", str(path)],
                            check=True, capture_output=True, text=True)
    rows = [line.rsplit("\t", 1) for line in result.stdout.splitlines()]
    require(len({r[0] for r in rows}) == len(rows), "Duplicate exported FASTA identifiers")
    return {key: int(length) for key, length in rows}


def create_vOTU_fasta(source, output, catalog, scratch):
    ids = scratch / "representatives.txt"
    catalog.contig_repid.to_csv(ids, header=False, index=False)
    extract_fasta(source, ids, output)
    actual = fasta_lengths(output)
    expected = dict(zip(catalog.contig_repid, catalog.contig_length.astype(int)))
    require(actual == expected, "Exported sequence IDs or lengths differ from annotation")
    return len(actual)


def create_vOTU_protein(source, output, catalog, scratch):
    parents = set(catalog.contig_repid)
    expected = set()
    represented = set()
    ids = scratch / "proteins.txt"
    process = subprocess.Popen(["seqkit", "seq", "-j", "2", "-n", "-i", str(source)],
                               stdout=subprocess.PIPE, text=True)
    with process.stdout, ids.open("w") as handle:
        for line in process.stdout:
            identifier = line.strip()
            parent, separator, ordinal = identifier.rpartition("_")
            if separator and ordinal.isdigit() and parent in parents:
                require(identifier not in expected, "Duplicate source protein ID")
                expected.add(identifier)
                represented.add(parent)
                handle.write(identifier + "\n")
    require(process.wait() == 0, "Protein identifier scan failed")
    require(represented == parents, "A selected representative has no source protein records")
    extract_fasta(source, ids, output)
    require(set(fasta_lengths(output)) == expected, "Exported protein identifiers differ")
    return len(expected)


def filter_by_vOTU(df, catalog, conn):
    """Preserve cluster-member rows and identify each record's actual source."""
    members = pd.read_sql_query("SELECT contig_id, vOTU_id FROM contig_annotation", conn)
    selected = df.merge(members, on="contig_id", how="inner", validate="many_to_one")
    selected = selected.merge(
        catalog[["vOTU_id", "contig_repid", "original_representative_id"]],
        on="vOTU_id", how="inner", validate="many_to_one",
    ).rename(columns={"contig_id": "source_contig_id"})
    selected["is_representative"] = selected.source_contig_id.eq(selected.contig_repid).astype(int)
    keys = ["contig_repid", "vOTU_id", "original_representative_id", "source_contig_id", "is_representative"]
    return selected[keys + [c for c in selected if c not in keys]].sort_values("vOTU_id", kind="stable")


def write_table(frame, name, connection, output):
    frame.to_sql(name, connection, index=False, if_exists="fail")
    frame.to_csv(output, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})


def build_release(work_db, previous_db, sources, output):
    require(not output.exists(), f"Destination already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs = {"working_catalog": work_db, "previous_release": previous_db,
              "genomes": sources / "ccflt.fasta.gz", "proteins": sources / "proteins.faa"}
    input_hashes = {key: {"path": str(path.resolve()), "sha256": sha256(path)} for key, path in inputs.items()}
    with tempfile.TemporaryDirectory(prefix=".release-", dir=output.parent) as temporary:
        stage = Path(temporary) / "package"
        stage.mkdir()
        scratch = Path(temporary) / "scratch"
        scratch.mkdir()
        raw, previous = read_database(work_db), read_database(previous_db)
        published = sqlite3.connect(stage / "AvianViromeDB.sqlite")
        try:
            catalog = read_catalog(raw)
            old = pd.read_sql_query("SELECT * FROM vOTU_annotation", previous)
            check_rows(catalog, old, [c for c in old if c != "contig_repid"])
            original = catalog[["vOTU_id", "original_representative_id"]].rename(
                columns={"original_representative_id": "contig_repid"})
            check_rows(original, old, ["vOTU_id", "contig_repid"])
            write_table(catalog, "vOTU_annotation", published, stage / "vOTU_annotation.tsv.gz")
            counts = {"vOTU_annotation": len(catalog)}
            specs = [
                ("host_filtered", "virus_prokyriotic_hosts", "vOTU_prokyriotic_hosts.tsv.gz"),
                ("amg_dramv_filtered", "amg_dramv", "AMG_DRAMv.tsv.gz"),
                ("amg_vibrant_filtered", "amg_vibrant", "AMG_VIBRANT.tsv.gz"),
                ("anno_prot_ARG", "arg", "ARG.tsv.gz"),
            ]
            for source_table, target_table, filename in specs:
                data = pd.read_sql_query(f"SELECT * FROM {source_table}", raw)
                exported = filter_by_vOTU(data, catalog, raw)
                baseline = pd.read_sql_query(f"SELECT * FROM {target_table}", previous)
                check_rows(exported, baseline, [c for c in baseline if c != "contig_repid"])
                write_table(exported, target_table, published, stage / filename)
                counts[target_table] = len(exported)
            counts["genomes"] = create_vOTU_fasta(inputs["genomes"], stage / "vOTUs.fasta.gz", catalog, scratch)
            counts["proteins"] = create_vOTU_protein(inputs["proteins"], stage / "vOTU_proteins.faa.gz", catalog, scratch)
            require(published.execute("PRAGMA integrity_check").fetchall() == [("ok",)], "SQLite integrity failed")
            with tarfile.open(stage / "gene_annotations.tar.gz", "w:gz") as archive:
                for filename in ("AMG_DRAMv.tsv.gz", "AMG_VIBRANT.tsv.gz", "ARG.tsv.gz"):
                    archive.add(stage / filename, arcname=filename)
            for filename in ("AMG_DRAMv.tsv.gz", "AMG_VIBRANT.tsv.gz", "ARG.tsv.gz"):
                (stage / filename).unlink()
            description = [
                "AvianViromeDB: local corrected release\n",
                "vOTU_id retains the existing cluster identity. contig_repid identifies the exported member.\n",
                "original_representative_id records the representative used by the existing network.\n",
                "Representative selection: highest CheckV completeness, then contig_id in ascending order.\n",
                "host_* and vOTU_prevalence summarize bird sources across cluster members.\n",
                "vc_* fields retain the existing network classification of the original representative.\n",
                "Host and gene tables retain member-level records: source_contig_id is their actual source;\n",
                "is_representative=1 marks records originating from the exported representative.\n",
                "The protein FASTA contains proteins of exported representatives; member-level annotations\n",
                "can refer to proteins outside this representative-only file.\n\n",
            ]
            for name in ["vOTU_annotation"] + [item[1] for item in specs]:
                columns = [row[1] for row in published.execute(f"PRAGMA table_info({name})")]
                description.append(f"{name} ({counts[name]} rows): " + ", ".join(columns) + "\n")
            (stage / "db_description.txt").write_text("".join(description))
            manifest = {"status": "local_corrected_release", "inputs": input_hashes,
                        "publisher_sha256": sha256(Path(__file__)),
                        "dependency_lock_sha256": sha256(Path(str(__file__) + ".lock")),
                        "pandas_version": pd.__version__, "counts": counts,
                        "representatives_changed": int(catalog.contig_repid.ne(catalog.original_representative_id).sum()),
                        "length_mismatches": 0, "identifier_mismatches": 0,
                        "catalog_attributes_unchanged": True, "member_annotation_rows_unchanged": True}
        finally:
            published.close()
            previous.close()
            raw.close()
        # Check again after generation, before exposing a package from mixed inputs.
        for key, path in inputs.items():
            require(sha256(path) == input_hashes[key]["sha256"], f"Input changed during export: {key}")
        manifest["files"] = {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                             for path in sorted(stage.iterdir())}
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(output)
        print(json.dumps({"output": str(output), **counts, "length_mismatches": 0}))


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    require(output.is_relative_to(base / "data") and not output.is_relative_to(base / "data/00-raw"),
            "Output must be inside an analysis step, outside raw inputs")
    build_release(base / "pc028e3.sqlite", base / "AvianViromeDB.sqlite",
                  base / "data/00-raw/d99-db_publish", output)


if __name__ == "__main__":
    main()
