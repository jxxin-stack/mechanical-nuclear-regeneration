# Mechanical Nuclear Regeneration

Analysis code for the manuscript **“Mechanical injury licenses adult nuclei to rebuild functional cells.”**

This repository contains only the custom scripts needed to reproduce the computational analyses and displayed figures. Raw sequencing data, large processed matrices and generated figures are not stored in Git.

## Analysis modules

| Script | Main outputs |
|---|---|
| `scripts/01_bulk_rna_analysis.R` | Bulk RNA expression panels, Fresh-to-NRC ranking, GSEA, blood-reference correlations and Fig. 3f heatmaps |
| `scripts/02_scrna_integration.R` | Fig. 4c Harmony-PCA/UMAP, Fig. 4d Spearman correlations and the monocyte/DC marker dot plot |
| `scripts/03_scrna_gsea.R` | Fig. 5b enrichment of NRC clusters relative to PBMC monocytes |
| `scripts/04_atac_peak_metrics.R` | Fig. 2c called-region counts and FRiP |
| `scripts/05_atac_regulatory_coupling.R` | Extended Data Fig. 3 accessibility changes and regulatory-element-linked RNA shifts |

The scripts start from processed matrices described in `data/README.md`. Public software such as MACS2, esATAC, Harmony, UMAP, fgsea and PECA is not redistributed here.

## Repository layout

```text
mechanical-nuclear-regeneration/
├── README.md
├── LICENSE
├── .gitignore
├── config/
├── data/
├── results/
└── scripts/
```

## Running the analyses

1. Obtain the public and study-derived inputs listed in `data/README.md`.
2. Place the processed files under `data/` using the documented names.
3. Install the packages listed in `config/R-packages.txt`.
4. Run the scripts in numerical order from the repository root:

```bash
Rscript scripts/01_bulk_rna_analysis.R
Rscript scripts/02_scrna_integration.R
Rscript scripts/03_scrna_gsea.R
Rscript scripts/04_atac_peak_metrics.R
Rscript scripts/05_atac_regulatory_coupling.R
```

Outputs are written to `results/`. Each script records `sessionInfo()` alongside its results.

## Reproducibility scope

- Bulk RNA analyses use the processed FPKM matrices specified in the manuscript; raw reads are not reprocessed.
- Single-cell analyses start from gene-by-cell UMI matrices and cell annotations deposited with the study.
- ATAC alignment, quality control and peak calling use the published esATAC and MACS2 workflows described in Methods. The custom downstream analyses start from the union-peak openness matrix, library QC table and state-specific regulatory networks.
- Dataset-source correction by Harmony is used only for joint visualization. Cell-type correlations and pathway rankings use uncorrected log-normalized expression.


