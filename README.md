# AvianViromeDB

Code for **AvianViromeDB**, a catalog of the avian gut DNA virome assembled from 2,658 public metagenomes, most of them from poultry.

![AvianViromeDB](www/logo.png)

The catalog contains 61,608 vOTUs of varying completeness from 2,658 public source records. Across CVD, CMKMC, MGV and IMG/VR v5, 38,093 representatives (61.83%) remain unmatched at ≥95% ANI over ≥85% of the shorter sequence. Gene functions and microbial hosts are computational predictions.

- Website and catalog browser: [AvianViromeDB](https://phagebyte.github.io/avianviromedb)
- Archived data release: [October 29, 2026 (version 4) · 10.5281/zenodo.23768304](https://doi.org/10.5281/zenodo.23768304)
- Data releases: Zenodo, concept DOI [10.5281/zenodo.15660655](https://doi.org/10.5281/zenodo.15660655)
- Feedback and study suggestions: [GitHub Issues](https://github.com/PhageByte/AvianViromeDB/issues)

## Contents

| Path | Purpose |
|---|---|
| [R/utils.R](R/utils.R) | Parsers for upstream tool outputs (VIBRANT, VirSorter2, CheckV, ANI clustering, vConTACT3, BACPHLIP, geNomad, VITAP and AMG tables), loaded with `devtools::load_all()` |
| [09-taxonomy.qmd](analyses/09-taxonomy.qmd) | Merges VITAP and geNomad taxonomy |
| [10-import.qmd](analyses/10-import.qmd) | Imports pipeline outputs, applies the catalog filters and builds the working SQLite database |
| [12-vc3_hq.qmd](analyses/12-vc3_hq.qmd) | Summarizes the vConTACT3 gene-sharing network |
| [99-db_publish.qmd](analyses/99-db_publish.qmd), [publish.py](analyses/publish.py) | Build the release package |
| [100-reviewer-revision.qmd](analyses/100-reviewer-revision.qmd) | Produces supplementary tables, external catalog comparisons, taxonomy corroboration, screening summaries and the release package |
| [101-revision-figures.qmd](analyses/101-revision-figures.qmd) | Draws figures from the step 100 outputs |
| [104-network-explorer.qmd](analyses/104-network-explorer.qmd) | Exports the node, module and adjacency files of the network explorer from the Fig. 6b layout |
| [verify_production.py](hpc/pc028e3-100-metavr/scripts/followup/verify_production.py) | Validates completed IMG/VR comparisons and summarizes matched reference sampling contexts for step 100 |
| [www/](www/) | Website source; the catalog browser reads `catalog-summary.json`, which step 100 generates, and the [network explorer](www/network.html) reads the files that step 104 writes |

## Running

The analysis steps follow a qproj layout and are rendered from `analyses/` with Quarto, for example `quarto render 100-reviewer-revision.qmd`. Each step reads its inputs from `analyses/data/00-raw/` or from earlier steps and writes to `analyses/data/<step>/`. The raw inputs and the working database `pc028e3.sqlite` are not part of this repository; the processed catalog is available from Zenodo, and Table S1 of the article lists the source sequencing accessions.

Each QMD declares its R packages and input files. The comparison scripts require Python 3.11 or newer with Biopython, NumPy, pandas, vpipe and taxMyPhage, plus BLAST+, Mash and the other executables named in the scripts. The included IMG/VR verification helper uses the Python standard library and reads completed comparison tables from `analyses/data/00-raw/hpc/`. Generated Slurm jobs inherit the submitting environment; activate these dependencies before submission and adapt resource requests to your cluster. Supplementary workbooks and source figures must be supplied at the paths named in the analysis steps.

Step 100 writes the browser metadata to `analyses/data/100-reviewer-revision/remaining-evidence/website/catalog-summary.json`. Copy that file into [www/](www/) and serve the directory over HTTP to use the browser. Step 104 writes the explorer files to `analyses/data/104-network-explorer/web/`; copy them into `www/network/data/` to use the network explorer. The browser metadata and the sequence downloads describe the same October 29, 2026 release.

## Releases and maintenance

The Deng laboratory (Technical University of Munich and Helmholtz Munich) maintains AvianViromeDB. Corrections and additions are released as needed; each release is archived separately on Zenodo and earlier versions remain citable. Release notes record the date, data scope and changes to identifiers or annotations. Public studies can be suggested for review and possible inclusion through the issue tracker.

## License

GPL-3.0-or-later; see [LICENSE.md](LICENSE.md).
