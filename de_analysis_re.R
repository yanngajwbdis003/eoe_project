rm(list = ls())
library(Seurat)
library(dplyr)

# ====subsampling, keep 20% <
load('eoe.rda') # load eoe.rda
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
#===================


load('eoe_sub.rda')
cat("Loaded eoe.sub:", ncol(eoe.sub), "cells\n")

# ===== gene filter
gene.all         = rownames(eoe.sub)
mt.gene          = grep('^MT-', gene.all, value=TRUE)
ribosome.gene    = grep('^RPL|^RPS|^MRPS|^MRPL', gene.all, value=TRUE)
problematic.gene = c('MALAT1', 'EEF1A1', 'TPT1')
gene.to.remove   = c(mt.gene, ribosome.gene, problematic.gene)
features         = setdiff(gene.all, gene.to.remove)
cat("Features after filter:", length(features), "\n")

eoe.sub@meta.data$nFeature_RNA_log2 = log2(eoe.sub@meta.data$nFeature_RNA + 1)
Idents(eoe.sub) = eoe.sub@meta.data$cell_type_anno
cell.type.lev   = levels(factor(eoe.sub@meta.data$cell_type_anno))

# ===== DE =======
run_de = function(disease1, disease2) {
  sapply(cell.type.lev, function(z) {
    ident.1 = names(Idents(eoe.sub))[Idents(eoe.sub) == z &
                                       eoe.sub@meta.data$disease_status == disease1]
    ident.2 = names(Idents(eoe.sub))[Idents(eoe.sub) == z &
                                       eoe.sub@meta.data$disease_status == disease2]
    
    if (length(ident.1) < 10 | length(ident.2) < 10) {
      cat("Skipping", z, "| ", disease1, ":", length(ident.1),
          disease2, ":", length(ident.2), "\n")
      return(NULL)
    }
    cat("Running:", z, "|", disease1, "vs", disease2, "\n")
    
    FindMarkers(eoe.sub, ident.1, ident.2,
                test.use    = 'LR',
                latent.vars = 'nFeature_RNA_log2',
                features    = features)
  }, simplify=FALSE)
}

# ===== function for filter up/down ======
filter_de = function(de.list) {
  list(
    up = sapply(de.list, function(z) {
      if (is.null(z) | length(z) <= 0) return(NULL)
      z %>% tibble::rownames_to_column(var='gene') %>%
        filter(!(gene %in% gene.to.remove)) %>%
        filter(avg_log2FC > log2(1.5)) %>%
        filter(pct.1 > 0.2) %>%
        filter(p_val_adj < 0.01)
    }, simplify=FALSE),
    down = sapply(de.list, function(z) {
      if (is.null(z) | length(z) <= 0) return(NULL)
      z %>% tibble::rownames_to_column(var='gene') %>%
        filter(!(gene %in% gene.to.remove)) %>%
        filter(avg_log2FC < -log2(1.5)) %>%
        filter(pct.2 > 0.2) %>%
        filter(p_val_adj < 0.01)
    }, simplify=FALSE)
  )
}

# ===== FilterAmbient ======
FilterAmbient = function(case, gene, cell.type) {
  gene = intersect(gene, rownames(eoe.sub))
  if (length(gene) <= 0) return(NULL)
  
  id      = eoe.sub@meta.data$donor_id == case
  cell    = rownames(eoe.sub@meta.data)[id]
  eoe.tmp = subset(eoe.sub, cells=cell)
  
  if (sum(Idents(eoe.tmp) == cell.type) <= 10) return(NULL)
  
  FindMarkers(eoe.tmp, cell.type,
              setdiff(unique(Idents(eoe.tmp)), cell.type),
              test.use    = 'LR',
              latent.vars = 'nFeature_RNA_log2',
              features    = gene)
}

# ===== SelectDEGene
SelectDEGene = function(cell.type, de.list.up, de.list.down,
                        disease='Active', ctrl='Ctrl') {
  # up
  id   = eoe.sub@meta.data$disease_status == disease
  case = unique(eoe.sub@meta.data$donor_id[id])
  gene = de.list.up$gene
  
  out = sapply(case, function(z) {
    cat(z, sep='\n')
    FilterAmbient(z, gene=gene, cell.type=cell.type)
  }, simplify=FALSE)
  
  res = sapply(out, function(z) {
    if (is.null(z)) return(list(pos=character(0), neg=character(0)))
    list(pos = rownames(z)[z$avg_log2FC > 0 & z$p_val_adj < 0.01],
         neg = rownames(z)[z$avg_log2FC < 0])
  }, simplify=FALSE)
  
  pos = table(unlist(sapply(res, function(z) z$pos)))
  neg = table(unlist(sapply(res, function(z) z$neg)))
  g   = intersect(names(pos), names(neg))
  gene.up.filt = setdiff(names(pos), g[which(pos[g] - neg[g] < 0)])
  
  # down
  id   = eoe.sub@meta.data$disease_status %in% ctrl
  case = unique(eoe.sub@meta.data$donor_id[id])
  gene = de.list.down$gene
  
  out = sapply(case, function(z) {
    cat(z, sep='\n')
    FilterAmbient(z, gene=gene, cell.type=cell.type)
  }, simplify=FALSE)
  
  res = sapply(out, function(z) {
    if (is.null(z)) return(list(pos=character(0), neg=character(0)))
    list(pos = rownames(z)[z$avg_log2FC > 0 & z$p_val_adj < 0.01],
         neg = rownames(z)[z$avg_log2FC < 0])
  }, simplify=FALSE)
  
  pos = table(unlist(sapply(res, function(z) z$pos)))
  neg = table(unlist(sapply(res, function(z) z$neg)))
  g   = intersect(names(pos), names(neg))
  gene.down.filt = setdiff(names(pos), g[which(pos[g] - neg[g] < 0)])
  
  list(up=gene.up.filt, down=gene.down.filt)
}

# ===== full pipeline
run_full_pipeline = function(disease1, disease2) {
  cat("\n========================================\n")
  cat("DE:", disease1, "vs", disease2, "\n")
  cat("========================================\n")
  
  de = run_de(disease1, disease2)
  ran = names(which(!sapply(de, is.null)))
  cat("Completed DE for", length(ran), "cell types\n")
  
  filt = filter_de(de)
  cat("\nUP genes:\n")
  print(sapply(filt$up, function(z) if(is.null(z)) 0 else nrow(z)))
  cat("\nDOWN genes:\n")
  print(sapply(filt$down, function(z) if(is.null(z)) 0 else nrow(z)))
  
  cat("\nRunning ambient filter...\n")
  ambient = sapply(cell.type.lev, function(z) {
    cat("\nAmbient filter:", z, "\n")
    SelectDEGene(z,
                 de.list.up   = filt$up[[z]],
                 de.list.down = filt$down[[z]],
                 disease      = disease1,
                 ctrl         = disease2)
  }, simplify=FALSE)
  
  cat("\nAfter ambient - UP:\n")
  print(sapply(ambient, function(z) length(z$up)))
  cat("\nAfter ambient - DOWN:\n")
  print(sapply(ambient, function(z) length(z$down)))
  
  list(de=de, up=filt$up, down=filt$down, ambient=ambient)
}

# ===== Remission vs Ctrl
res.rem.ctrl = run_full_pipeline('Active', 'Ctrl')
save(res.rem.ctrl, file='de_active_ctrl_full.rda')
cat("Saved de_active_ctrl_full.rda\n")

# ===== Remission vs Ctrl
res.rem.ctrl = run_full_pipeline('Remission', 'Ctrl')
save(res.rem.ctrl, file='de_rem_ctrl_full.rda')
cat("Saved de_rem_ctrl_full.rda\n")

# ===== Active vs Remission
res.active.rem = run_full_pipeline('Active', 'Remission')
save(res.active.rem, file='de_active_rem_full.rda')
cat("Saved de_active_rem_full.rda\n")

cat("\nAll done.\n")