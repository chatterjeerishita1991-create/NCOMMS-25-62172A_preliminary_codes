# TP53 TCGA figure reproduction

This repository is the minimal, from-source workflow for eight requested TP53
figures. All small, immutable analysis inputs are versioned under `data/input/`;
only downloaded expression data and generated results are ignored by Git.

## Included figure targets

| Figure | Output | Reproduction analysis |
| --- | --- | --- |
| 1A | `figures/Figure_1A_survival.pdf` | Pan-cancer TCGA-CDR overall survival, R175 vs R273 |
| 1B | `figures/Figure_1B_CD274.pdf` | CD274 VST boxplots in UCEC, STAD, LGG, COAD, and BRCA |
| 2G | `figures/Figure_2G_emt_enrichment.pdf` | OV R175H vs R273H, Hallmark EMT |
| 3E | `figures/Figure_3E_tnfa_nfkb_enrichment.pdf` | OV R175H vs R273H, Hallmark TNFA/NFKB |
| S2G | `figures/Figure_S2G_volcano.pdf` | OV R175H vs R273H volcano |
| S2H | `figures/Figure_S2H_emt_volcano.pdf` | Volcano with Hallmark EMT genes highlighted |
| S2I | `figures/Figure_S2I_kegg_oxphos_volcano.pdf` | Volcano with KEGG oxidative-phosphorylation genes highlighted |
| S3E | `figures/Figure_S3E_ifngamma_enrichment.pdf` | OV R175H vs R273H, Hallmark interferon-gamma response |

## Setup

Use R 4.3+ and install the packages declared in `DESCRIPTION`. Bioconductor
packages are `TCGAbiolinks`, `SummarizedExperiment`, `DESeq2`, and `fgsea`.
For example, install Bioconductor packages with `BiocManager::install()` and
the remaining packages with `install.packages()`.

Before running the ovarian workflow, add these exact MSigDB 7.5.1 GMT files to
`data/input/gene_sets/` and commit them to the repository:

- `h.all.v7.5.1.symbols.gmt`
- `c2.cp.kegg.v7.5.1.symbols.gmt`

The compact `tcga_cdr_survival.csv` already contains the TCGA-CDR patient ID,
overall-survival event, and overall-survival time fields needed for Figure 1A.
STAR count data are acquired directly from the GDC through TCGAbiolinks; no
preprocessed RDS objects are used.

Run the complete workflow from the repository root:

```sh
make all
```

Or run the stages explicitly:

```sh
Rscript scripts/01_fetch_expression.R
Rscript scripts/02_pan_cancer_figures.R
Rscript scripts/03_ovarian_de_gsea_figures.R
```

`01_fetch_expression.R` downloads only BRCA, COAD, LGG, STAD, UCEC, and OV.
Downloaded GDC files are cached in `data/cache/`; later runs reuse the prepared
cohort objects.

## Reproduction choices retained from the source analyses

- `data/input/tcga_group_labels.csv` is the sole source of cohort labels. Its
  patient barcode is the matching key for TCGA RNA samples; unlisted samples
  receive `None`. It has 445 pan-cancer labels (170 R175; 275 R273) and 25
  ovarian labels (15 R175H; 10 R273H). No MAF is read by this repository.
- Figure 1B includes only the five displayed cancer types and uses blind DESeq2
  VST, primary-tumor STAR unstranded counts, a two-sided unpaired Wilcoxon test,
  and unadjusted significance labels.
- Ovarian DE/GSEA uses the explicit R175H/R273H labels, primary tumors, and the
  DESeq2 R175H/R273H contrast. Its enrichment and volcano plots use local,
  namespaced copies of the archived CMCCTools plotting functions, including
  `EnhancedVolcano`, the signed natural-log p-value ranking, the CMCCTools
  pathway ordering, and seed-controlled null enrichment traces.
- Enrichment rankings are signed `-log10(pvalue)` values, with positive scores
  denoting higher expression in R175H. `results/` preserves the DE and GSEA
  tables and a session record after each stage.

The upstream data services and package releases can evolve, so regenerated
counts or numerical results may differ from the archived figures even when the
analysis definition is unchanged.
