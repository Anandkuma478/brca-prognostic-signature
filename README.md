# Prognostic gene signature in breast cancer (TCGA, GEO and METABRIC)

Analysis code for "Identification and validation of a prognostic gene signature in breast cancer using TCGA and GEO data".

Two LASSO Cox models for overall survival were developed in SCAN-B (GSE96058), frozen, and tested without retraining in TCGA-BRCA and METABRIC.

## Data

Raw data are not included in this repository. They are available from:

| Dataset | Source |
|---|---|
| TCGA-BRCA expression, clinical and mutation data | Genomic Data Commons, https://portal.gdc.cancer.gov/projects/TCGA-BRCA |
| SCAN-B | Gene Expression Omnibus, accession GSE96058 |
| METABRIC | cBioPortal, study `brca_metabric` (https://datahub.assets.cbioportal.org/brca_metabric.tar.gz), extracted to `data/metabric/raw/` |
| Human Protein Atlas | https://www.proteinatlas.org/download/tsv/ (`cancer_data.tsv.zip`, `normal_ihc_data.tsv.zip`), extracted to `data/hpa/` |

## Running the analysis

Run every script from the project root, in the order below. Package versions are listed in `SOFTWARE_VERSIONS.txt`. Two conda environments are used: `brca_genomics` and `brca_modules`.

| Step | Script | Environment | Purpose |
|---|---|---|---|
| 1 | `01_download.R` | brca_genomics | Download TCGA-BRCA |
| 2 | `02_qc.R` | brca_genomics | Quality control and filtering |
| 3 | `03_de.R`, `03b_volcano_heatmap.R` | brca_genomics | Differential expression and figures |
| 4 | `04_survival_prep.R` | brca_genomics | TCGA-BRCA survival data |
| 5 | `06a_scanb_prep.R` | brca_genomics | SCAN-B preparation |
| 6 | `06b_scanb_signature.R`, `06c_signature_check.R` | brca_genomics | Model development, nested cross-validation, frozen models |
| 7 | `07_stability.R` | brca_genomics | Bootstrap stability |
| 8 | `05_signature.R`, `05b_signature.R` | brca_genomics | Exploratory TCGA-BRCA-only signatures |
| 9 | `08_wgcna.R` | brca_genomics | Co-expression modules |
| 10 | `10_string_hubs.R`, `11_enrichment.R` | brca_modules | Hub genes and pathway enrichment |
| 11 | `13_batch_check.R` | brca_genomics | Batch-effect check |
| 12 | `09_metabric_validation.R` | brca_genomics | METABRIC validation (primary analysis, run once on the frozen models) |
| 13 | `12_survival_evaluation.R` | brca_modules | Kaplan-Meier, time-dependent ROC, calibration, hub gene survival |
| 14 | `14_cindex_ci.R` | brca_genomics | C-index confidence intervals, proportional hazards check, TCGA-BRCA added value |
| 15 | `15_batch_correction_sensitivity.R` | brca_genomics | PCA before and after batch correction; sensitivity analysis |
| 16 | `16_immune_association.R` | brca_modules | Immune features |
| 17 | `17_tcga_maf_download.R`, `18_mutation_association.R` | brca_genomics | Somatic mutations |
| 18 | `19_protein_atlas.R` | brca_genomics | Human Protein Atlas summary |
| 19 | `20_published_and_random_benchmark.R` | brca_modules | Published and random signature benchmark |

## Outputs

`results/` holds the result tables and `figures/` holds the figures.
