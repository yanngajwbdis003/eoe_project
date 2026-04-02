library(Seurat)
library(ggplot2)
library(dplyr)
install.packages('ggrepel')
library(ggrepel)

load('eoe_sub.rda')
DefaultAssay(eoe.sub) = 'RNA'


#----fig5c----
# Fibroblast
eoe.fib = subset(eoe.sub, cell_type_anno == 'Fibroblast')

# target gene
genes.c = c('IL13RA2', 'IGF2', 'C3', 'CFD', 'IGFBP2', 'GAS6', 'C1S')

cat(genes.c %in% rownames(eoe.fib), "\n")

#  violin plot
VlnPlot(eoe.fib, 
        features = genes.c,
        group.by = 'disease_status',
        cols = c('Ctrl'='#7CAE00', 'Remission'='#00BFC4', 'Active'='#F8766D'),
        pt.size = 0,
        ncol = 4) 

ggsave('fig5c.pdf', width=14, height=6)

#---fig5d----
library(glmnet)
library(compositions)

load('eoe_sub.rda')
DefaultAssay(eoe.sub) = 'RNA'

# =====  IL13RA2+ fibroblast =====
eoe.fib = subset(eoe.sub, cell_type_anno == 'Fibroblast')
il13ra2_expr = FetchData(eoe.fib, vars='IL13RA2')
eoe.fib@meta.data$is_il13ra2_pos = il13ra2_expr$IL13RA2 > 0
fib_ratio = eoe.fib@meta.data %>%
  group_by(donor_id, disease_status) %>%
  summarise(
    n_fib     = n(),
    n_pos     = sum(is_il13ra2_pos),
    ratio_pos = n_pos / n_fib,
    .groups   = 'drop'
  )

cat("IL13RA2+ fibroblast ratio per donor:\n")
print(fib_ratio)


# ===== Panel d：CLR + PCA =====
library(glmnet)
library(compositions)  


# ===== IL13RA2+ fibroblast =====
eoe.fib = subset(eoe.sub, cell_type_anno == 'Fibroblast')
il13ra2_expr = FetchData(eoe.fib, vars='IL13RA2')
eoe.fib@meta.data$is_il13ra2_pos = il13ra2_expr$IL13RA2 > 0
fib_ratio = eoe.fib@meta.data %>%
  group_by(donor_id, disease_status) %>%
  summarise(
    n_fib     = n(),
    n_pos     = sum(is_il13ra2_pos),
    ratio_pos = n_pos / n_fib,
    .groups   = 'drop'
  )

cat("IL13RA2+ fibroblast ratio per donor:\n")
print(fib_ratio)


comp_mat = eoe.sub@meta.data %>%
  group_by(donor_id, cell_type_anno) %>%
  summarise(n=n(), .groups='drop') %>%
  tidyr::pivot_wider(names_from=cell_type_anno, 
                     values_from=n, values_fill=0) %>%
  tibble::column_to_rownames('donor_id')

clr_mat = as.data.frame(clr(comp_mat + 0.5))

# PCA
pca_res = prcomp(clr_mat, scale.=TRUE)
pc1 = pca_res$x[, 1]
plot_d = fib_ratio %>%
  left_join(data.frame(donor_id=names(pc1), PC1=pc1), by='donor_id')
disease_cols = c('Ctrl'='#7CAE00', 'Remission'='#00BFC4', 'Active'='#F8766D')
plot_list = lapply(c('Ctrl','Remission','Active'), function(ds) {
  df = plot_d %>% filter(disease_status == ds)
  
  ct = if(nrow(df) > 2) cor.test(df$PC1, df$ratio_pos) else NULL
  r_val = if(!is.null(ct)) round(ct$estimate, 2) else NA
  p_val = if(!is.null(ct)) round(ct$p.value, 2) else NA
  
  ggplot(df, aes(x=PC1, y=ratio_pos, label=donor_id)) +
    geom_point(color=disease_cols[ds], size=2) +
    geom_smooth(method='lm', se=FALSE, color='black', linewidth=0.5) +
    ggrepel::geom_text_repel(size=2.5, max.overlaps=20) +  
    annotate('text', x=Inf, y=Inf, hjust=1.1, vjust=1.5,
             label=paste0('R = ', r_val, '\nFDR = ', p_val), size=3) +
    labs(title=ds, x='PC1', y='Fraction of IL13RA2+ fibroblasts') +
    theme_classic(base_size=10)
})

library(patchwork)
p_d = wrap_plots(plot_list, nrow=1)
ggsave('fig5d.pdf', p_d, width=12, height=4)
cat("Saved fig5d.pdf\n")


# ===== Panel e：Lasso =====
comp_ratio = sweep(as.matrix(comp_mat), 1, rowSums(comp_mat), '/')
common_donors = intersect(rownames(comp_ratio), fib_ratio$donor_id)
X = comp_ratio[common_donors, ]
y = fib_ratio$ratio_pos[match(common_donors, fib_ratio$donor_id)]

# Lasso
set.seed(42)
cv_fit = cv.glmnet(X, y, alpha=1)
lasso_fit = glmnet(X, y, alpha=1, lambda=cv_fit$lambda.min)
coef_mat = coef(lasso_fit)
coef_df = data.frame(
  cell_type = rownames(coef_mat)[-1],
  coef      = as.vector(coef_mat)[-1]
) %>% arrange(coef)

cat("Lasso coefficients:\n")
print(coef_df %>% filter(coef != 0))
coef_df$index = seq_along(coef_df$cell_type)

p_e_left = ggplot(coef_df, aes(x=index, y=coef, label=cell_type)) +
  geom_point(size=2, color='black') +
  geom_hline(yintercept=0, linetype='dashed') +
  ggrepel::geom_text_repel(
    data = coef_df %>% filter(coef != 0), 
    size=2.5, max.overlaps=20
  ) +
  labs(x='Cell type index', y='Lasso linear regression coefficient') +
  theme_classic(base_size=10)

top_ct = coef_df %>% filter(coef != 0) %>% 
  arrange(desc(abs(coef))) %>% slice(1) %>% pull(cell_type)
cat("Top predictor:", top_ct, "\n")

x_top = X[, top_ct]
plot_e_right = data.frame(
  donor_id  = common_donors,
  x_ratio   = x_top,
  il13ra2   = y,
  disease_status = fib_ratio$disease_status[match(common_donors, fib_ratio$donor_id)]
  
)

p_e_right = ggplot(plot_e_right, 
                   aes(x=x_ratio, y=il13ra2, 
                       color=disease_status, 
                       label=donor_id)) +
  geom_point(size=2) +
  geom_smooth(method='lm', se=FALSE, color='black', linewidth=0.5) +
  ggrepel::geom_text_repel(size=2.5, max.overlaps=20) +
  scale_color_manual(values=c('Ctrl'='#7CAE00', 
                              'Remission'='#00BFC4', 
                              'Active'='#F8766D')) +
  labs(x='IL13RA2+ fibroblast ratios', y='Lasso prediction') +
  theme_classic(base_size=10)

p_e = p_e_left + p_e_right
ggsave('fig5e.pdf', p_e, width=10, height=4)
cat("Saved fig5e.pdf\n")