# Input data

Large data files are intentionally excluded from Git. Deposit the minimum processed inputs with the study data archive or in the same Zenodo record as the code.

## Bulk RNA (`data/bulk/`)

- `merged_fpkm_qn.tsv`: gene-level, median-quantile-normalized FPKM matrix used for the Fresh/NRC and ENCODE plots. The first columns must include `Symbol`; remaining columns are samples.
- `sample_metadata.tsv`: tab-separated columns `sample_id`, `state`, `donor`, `cell_type`, and `reference_panel`. Use `Fresh` and `NRC` in `state`. Technical libraries from the same donor must share a donor label.
- Optional separate Corces et al. and ENCODE matrices may be represented in the same matrix and distinguished by `reference_panel`.

The study RNA-seq data are available from GEO under GSE242346. The Corces et al. blood reference is available under GSE74246. ENCODE accessions are listed in the manuscript tables.

## Single-cell RNA (`data/scrna/`)

- `nrc_counts.rds`: sparse gene-by-cell raw UMI matrix.
- `nrc_metadata.tsv`: columns `cell_id`, `cluster`, and `qc_pass`; clusters are `C0` to `C4`.
- `pbmc_v31_counts.rds`: sparse gene-by-cell raw UMI matrix for the 10x 3-prime v3.1 PBMC reference.
- `pbmc_v31_metadata.tsv`: columns `cell_id`, `cell_type`, and `qc_pass` after the annotation and quality-control procedure described in Methods.

The PBMC reference is the 10x Genomics `SC3_v3_NextGem_DI_PBMC_10K` filtered-feature matrix. Processed cell annotations used in the manuscript should be deposited as source data.

## ATAC and regulatory networks (`data/atac/`)

- `library_qc.tsv`: columns `sample_id`, `state`, `donor`, `peak_count`, and `frip`.
- `union_peak_openness.tsv`: union-peak openness scores, with peak identifiers in the first column and ATAC libraries in subsequent columns.
- `sample_metadata.tsv`: columns `sample_id`, `state`, and `donor` matching the openness matrix.
- `blacklisted_peaks.txt`: union-peak identifiers excluded from analysis.
- `fresh_network.tsv` and `nrc_network.tsv`: state-specific PECA regulatory-element-to-target-gene links with the columns documented in the deposited source data.
- `rna_fpkm_qn.tsv`: the processed RNA matrix used for regulatory-element-linked target analysis.

Raw ATAC reads and the processed union-peak matrix must receive a persistent accession before publication.

## Gene sets

- `config/bulk_gene_sets.gmt`: gene sets displayed in the bulk ranked-expression analyses; included with the code.
- `ReactomePathways.gmt`: Reactome pathway definitions corresponding to the release used in the manuscript.
- `config/macrophage_signature.tsv`: the a priori macrophage-associated signature; included with the code.
