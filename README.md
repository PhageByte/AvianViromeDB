# AvianViromeDB

Code for **AvianViromeDB**, a catalog of the avian gut DNA virome assembled from 2,658 public source records, mostly from poultry.

![AvianViromeDB](www/logo.png)

The catalog contains 61,608 vOTUs of varying completeness. Across CVD, CMKMC, MGV and IMG/VR v5, 38,093 representative sequences (61.83%) have no match at ≥95% ANI over ≥85% of the shorter sequence. Gene functions and microbial hosts are predicted computationally.

- Website and catalog browser: [AvianViromeDB](https://phagebyte.github.io/avianviromedb)
- Archived data release: [October 29, 2026 (version 4) · 10.5281/zenodo.23768304](https://doi.org/10.5281/zenodo.23768304)
- Data releases: Zenodo, concept DOI [10.5281/zenodo.15660655](https://doi.org/10.5281/zenodo.15660655)
- Feedback and study suggestions: [GitHub Issues](https://github.com/PhageByte/AvianViromeDB/issues)

## Contents

| Path | Purpose |
|---|---|
| [R/utils.R](R/utils.R) | Parsers for upstream tool outputs, loaded with `devtools::load_all()` |
| [09-taxonomy.qmd](analyses/09-taxonomy.qmd) | Merges VITAP and geNomad taxonomy |
| [10-import.qmd](analyses/10-import.qmd) | Imports pipeline outputs and builds the working SQLite database |
| [12-vc3_hq.qmd](analyses/12-vc3_hq.qmd) | Summarizes the vConTACT3 gene-sharing network |
| [99-db_publish.qmd](analyses/99-db_publish.qmd) | Archived release files |
| [publish.py](analyses/publish.py) | Builds the release package |
| [100-reviewer-revision.qmd](analyses/100-reviewer-revision.qmd) | Supplementary tables, external catalog comparisons and the release package |
| [101-revision-figures.qmd](analyses/101-revision-figures.qmd) | Manuscript figures |
| [104-network-explorer.qmd](analyses/104-network-explorer.qmd) | Network explorer data |
| [www/](www/) | Website source |

## Running

Render each step from `analyses/` with Quarto, for example `quarto render 100-reviewer-revision.qmd`. Steps read from `analyses/data/00-raw/` or earlier steps and write to `analyses/data/<step>/`. Raw inputs and the working database are stored outside this repository; Zenodo hosts the processed catalog, and Table S1 of the article lists the source accessions. Each QMD lists its packages, tools and input files; `publish.py` runs with uv and Python 3.13 or later.

To serve the website locally, copy `catalog-summary.json` from step 100 into `www/` and the step 104 `web/` files into `www/network/data/`.

## Releases and maintenance

The Deng laboratory (Technical University of Munich and Helmholtz Munich) maintains AvianViromeDB and archives each release on Zenodo with a version DOI. Report issues or suggest public studies through the issue tracker.

## License

GPL-3.0-or-later; see [LICENSE.md](LICENSE.md).
