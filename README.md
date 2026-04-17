# EoE Single-Cell RNA-seq Analysis

Reproduction and extension of https://www.nature.com/articles/s41467-024-47647-0#Sec2 (Ding et al. 2024).

---

## Data

orig source: https://singlecell.broadinstitute.org/single_cell/study/SCP1242/eoe-eosinophilic-esophagitis#study-download

our processed data: https://drive.google.com/drive/folders/1JNF3krwtIS6jUbGPptXp2njbLXli91A_?usp=sharing

| File                     | Description                                                                                                                 |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------------- |
| `EoE.mtx`                | Raw count matrix (33,694 genes × 393,763 cells), Matrix Market format                                                       |
| `EoE_cell.tsv`           | Cell barcodes for raw matrix                                                                                                |
| `EoE_gene.tsv`           | Gene names for raw matrix                                                                                                   |
| `EoE_processed.mtx`      | Normalized/processed count matrix (same dimensions)                                                                         |
| `EoE_cell_processed.tsv` | Cell barcodes for processed matrix                                                                                          |
| `EoE_gene_processed.tsv` | Gene names for processed matrix                                                                                             |
| `EoE_meta.txt`           | Cell-level metadata: disease status (Active/Remission/Ctrl), donor ID, cell type annotation, sex, organ, etc.               |
| `EoE_coord_2d.txt`       | 2D UMAP coordinates for each cell                                                                                           |
| `eoe_sub.rda`            | Seurat object: stratified 20% subsample (~79k cells), used for all downstream analysis                                      |
| `de_results.rda`         | DE results: Active vs Ctrl — includes `cluster.marker.lr.log2`, `marker.up.filt`, `marker.down.filt`, `marker.filt.ambient` |
| `de_rem_ctrl_full.rda`   | DE results: Remission vs Ctrl — same structure as above in `res.rem.ctrl`                                                   |
| `de_active_rem_full.rda` | DE results: Active vs Remission — same structure in `res.active.rem`                                                        |

---

## Code

| File               | Description                                                                                                                                                                    |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `convert_eoe.R`    | construct eoe.rda from mtx tsv; subsampling.                                                                                                                                   |
| `de_analysis_re.R` | Main DE pipeline: Active vs Ctrl. Runs LR test per cell type, ambient RNA filter, gene selection. Saves `de_results.rda` `de_rem_ctrl_full.rda` and `de_active_rem_full.rda`   |
| `fig2a.R`          | Reproduces Figure 2a. Visualizes the percentage shifts of key clinical markers (e.g., Eosinophils, Th2 cells) across disease states using boxplots.                            |
| `fig2b.R`          | Reproduces Figure 2b (Volcano plot for correlation analysis, integrating scRNA-seq cell fractions with clinical Eos/HPF metadata via Spearman correlation and FDR correction). |
| `fig2c.R`          | Reproduces Figure 2c. (Regression scatter plot comparing Eos/HPF and CLR-transformed cell abundance).                                                                          |
| `fig2d.R`          | Reproduces Figure 2d. (PCA Biplot, Dual-panel)                                                                                                                                 |
| `fig5a_b.R`        | Reproduces Fig. 5a (pairwise log2FC scatter plots) and Fig. 5b (GO BP enrichment heatmap)                                                                                      |
| `figure5cde.R`     | Reproduces Fig. 5c (violin plots of key fibroblast genes), Fig. 5d (IL13RA2+ fibroblast fraction vs PC1), and Fig. 5e (Lasso regression)                                       |
| `script.py`     | Script for running cNMF for full cell types, epi subsets                                       |
| `cnmf.ipynb`     | notebook for processing the cnmf output, performing analysis                                      |


