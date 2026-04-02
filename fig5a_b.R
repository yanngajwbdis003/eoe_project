
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

load('de_results.rda')          # Active vs Ctrl
load('de_rem_ctrl_full.rda')    # Remission vs Ctrl
load('de_active_rem_full.rda')  # Active vs Remission

# check cell type
cat("Active vs Ctrl:\n")
print(sapply(cluster.marker.lr.log2, function(z) if(is.null(z)) 0 else nrow(z)))

cat("\nRemission vs Ctrl:\n")
print(sapply(res.rem.ctrl$de, function(z) if(is.null(z)) 0 else nrow(z)))

cat("\nActive vs Remission:\n")
print(sapply(res.active.rem$de, function(z) if(is.null(z)) 0 else nrow(z)))

# ====== 5a
library(ggplot2)
library(dplyr)
library(patchwork)

# only take signficant DEs
filter_sig = function(de.list) {
  lapply(names(de.list), function(ct) {
    z = de.list[[ct]]
    if (is.null(z) || length(z) == 0 || nrow(z) == 0) return(NULL)
    z = z[z$p_val_adj < 0.05, ]  
    if (nrow(z) == 0) return(NULL)
    data.frame(gene=rownames(z), cell_type=ct, fc=z$avg_log2FC)
  }) %>% bind_rows()
}

fc.active.ctrl = filter_sig(cluster.marker.lr.log2) %>% rename(fc_active_ctrl=fc)
fc.rem.ctrl    = filter_sig(res.rem.ctrl$de)         %>% rename(fc_rem_ctrl=fc)
fc.active.rem  = filter_sig(res.active.rem$de)       %>% rename(fc_active_rem=fc)

# combine
df1 = inner_join(fc.active.ctrl, fc.active.rem, by=c('gene','cell_type'))  # top
df2 = inner_join(fc.active.ctrl, fc.rem.ctrl,   by=c('gene','cell_type'))  # mid
df3 = inner_join(fc.rem.ctrl,    fc.active.ctrl, by=c('gene','cell_type')) # bottom

scatter_panel = function(df, x, y, xlabel, ylabel) {
  ct  = cor.test(df[[x]], df[[y]])
  r   = round(ct$estimate, 2)
  p   = formatC(ct$p.value, format='e', digits=2)
  lm_fit    = lm(as.formula(paste(y, '~', x)), data=df)
  intercept = round(coef(lm_fit)[1], 3)
  slope     = round(coef(lm_fit)[2], 2)
  
  ggplot(df, aes_string(x=x, y=y)) +
    geom_point(alpha=0.15, size=0.3, color='grey40') +
    geom_smooth(method='lm', se=FALSE, color='black', linewidth=0.5) +
    geom_hline(yintercept=0, linetype='dashed', color='grey60') +
    geom_vline(xintercept=0, linetype='dashed', color='grey60') +
    annotate('text', x=Inf, y=Inf, hjust=1.1, vjust=1.5, size=2.8,
             label=paste0('R = ', r, ', p = ', p,
                          '\ny = ', intercept, ' + ', slope, 'x')) +
    labs(x=xlabel, y=ylabel) +
    theme_classic(base_size=10)
}

p1 = scatter_panel(df1, 'fc_active_ctrl', 'fc_active_rem',
                   'Active vs healthy', 'Active vs remission')
p2 = scatter_panel(df2, 'fc_active_ctrl', 'fc_rem_ctrl',
                   'Active vs healthy', 'Remission vs healthy')
p3 = scatter_panel(df3, 'fc_rem_ctrl',    'fc_active_ctrl',
                   'Remission vs healthy', 'Active vs healthy')

p_a = p1 / p2 / p3 +
  plot_annotation(caption='Log2(FC) of differentially expressed genes')
ggsave('fig5a.pdf', p_a, width=4, height=10)
cat("Saved fig5a.pdf\n")



# ======5b=======
library(ggplot2)
library(tidyr)
library(dplyr)

load('de_results.rda')
load('de_active_rem_full.rda')

genes.for.go = sapply(cell.type.lev, function(ct) {
  up.ctrl = marker.filt.ambient[[ct]]$up
  up.rem  = res.active.rem$ambient[[ct]]$up
  
  down.ctrl = marker.filt.ambient[[ct]]$down
  down.rem  = res.active.rem$ambient[[ct]]$down
  
  list(
    up   = intersect(up.ctrl,   up.rem),
    down = intersect(down.ctrl, down.rem)
  )
}, simplify=FALSE)

# check overlapping cells
cat("UP genes after intersection:\n")
print(sapply(genes.for.go, function(z) length(z$up)))
cat("DOWN genes after intersection:\n")
print(sapply(genes.for.go, function(z) length(z$down)))

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)

go.results.b = sapply(names(genes.combined), function(ct) {
  genes = genes.combined[[ct]]
  if (length(genes) < 5) return(NULL)
  cat("GO for:", ct, "\n")
  
  ego = enrichGO(
    gene          = genes,
    universe      = features,  
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff  = 0.05
  )
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>% filter(p.adjust < 0.05)
}, simplify=FALSE)

cat("\nGO terms per cell type:\n")
print(sapply(go.results.b, function(z) if(is.null(z)) 0 else nrow(z)))

# ===== heatmap =====
go.combined = bind_rows(lapply(names(go.results.b), function(ct) {
  df = go.results.b[[ct]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(cell_type     = ct,
           neg_log10_fdr = -log10(p.adjust),
           gene_count    = Count)
}))

top.terms = go.combined %>%
  group_by(Description) %>%
  summarise(max_fdr = max(neg_log10_fdr)) %>%
  arrange(desc(max_fdr)) %>%
  head(25) %>%
  pull(Description)

heat.df = go.combined %>%
  filter(Description %in% top.terms) %>%
  dplyr::select(Description, cell_type, neg_log10_fdr) %>%
  tidyr::complete(Description, cell_type, fill=list(neg_log10_fdr=0)) %>%
  tidyr::pivot_wider(names_from=cell_type, 
                     values_from=neg_log10_fdr, values_fill=0)

heat.mat = as.matrix(heat.df[, -1])
rownames(heat.mat) = heat.df$Description

count.df = go.combined %>%
  filter(Description %in% top.terms) %>%
  dplyr::select(Description, cell_type, gene_count) %>%
  tidyr::complete(Description, cell_type, fill=list(gene_count=0)) %>%
  tidyr::pivot_wider(names_from=cell_type,
                     values_from=gene_count, values_fill=0)

count.mat = as.matrix(count.df[, -1])
rownames(count.mat) = count.df$Description
count.mat = count.mat[rownames(heat.mat), colnames(heat.mat)]

pheatmap(heat.mat,
         display_numbers = count.mat,
         number_format   = '%d',
         number_color    = 'black',
         color           = colorRampPalette(c('white','#2166ac'))(50),
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 8,
         fontsize_number = 6,
         main            = 'GO BP enrichment (Active vs Ctrl & Remission)',
         filename        = 'fig5b.pdf')

cat("Saved fig5b.pdf\n")