rm(list = ls())
library(Seurat)
library(dplyr)

counts <- Read10X(data.dir = "data_orig/")

meta_raw <- read.delim("EoE_meta.txt", header = TRUE, stringsAsFactors = FALSE)
meta <- meta_raw[-1, ] 
rownames(meta) <- meta$NAME 

eoe <- CreateSeuratObject(counts = counts, project = "EoE_Atlas", meta.data = meta)

save(eoe, file = "eoe.rda")


# ====subsampling, keep 20% <
load('eoe.rda') # load eoe.rda (already created using mtx. gene/cells.tsv)
DefaultAssay(eoe) = 'RNA'
eoe = NormalizeData(eoe)
gc()
set.seed(42)
cells.keep = unlist(lapply(cell.type.lev, function(ct) {
  cells.ct = rownames(eoe@meta.data)[eoe@meta.data$cell_type_anno == ct]
  n.keep = max(50, round(length(cells.ct) * 0.2))  # at least 50
  sample(cells.ct, min(n.keep, length(cells.ct)))
}))
eoe.sub = subset(eoe, cells = cells.keep)
cat("Original cells:", ncol(eoe), "\n")
cat("After 20% subsample:", ncol(eoe.sub), "\n")
table(eoe.sub@meta.data$disease_status)
DefaultAssay(eoe.sub) = 'RNA'