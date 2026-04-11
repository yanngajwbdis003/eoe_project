# --- 1. LOAD LIBRARIES ---
library(tidyverse)
library(data.table)
library(compositions) 
library(ggrepel)
library(patchwork)

# --- 2. LOAD DATA ---
metadata <- fread("EoE_meta.txt", fill = TRUE)

# --- 3. DATA PREPARATION & CLR TRANSFORMATION ---
count_data <- metadata %>%
  filter(NAME != "TYPE", cell_type_anno != "" & !is.na(cell_type_anno)) %>%
  group_by(biosample_id, disease_status, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop")

count_matrix <- count_data %>%
  select(biosample_id, cell_type_anno, Count) %>%
  pivot_wider(names_from = cell_type_anno, values_from = Count, values_fill = 0) %>%
  column_to_rownames("biosample_id")

clr_transformed <- clr(count_matrix + 1)

# --- 4. PCA & AXIS FLIPPING ---
pca_results <- prcomp(clr_transformed)

# Calculate Variance Explained
vars <- pca_results$sdev^2
vars_perc <- round(vars / sum(vars) * 100, 1)

# Prepare Sample Coordinates, FLIP PC1, and RENAME "Ctrl"
pca_samples <- as.data.frame(pca_results$x) %>%
  rownames_to_column("biosample_id") %>%
  mutate(PC1 = PC1 * -1) %>% 
  left_join(distinct(metadata, biosample_id, disease_status), by = "biosample_id") %>%
  # --- FIX: RENAME LABEL HERE ---
  mutate(disease_status = recode(disease_status, "Ctrl" = "Healthy"))

# Prepare Cell Type Loadings & FLIP PC1
pca_loadings <- as.data.frame(pca_results$rotation) %>%
  rownames_to_column("cell_type_anno") %>%
  mutate(PC1 = PC1 * -1) 

# --- 5. VISUALIZATION ---

# Plot A: Sample Composition (Left)
plot_left <- ggplot(pca_samples, aes(x = PC1, y = PC2, color = disease_status)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = biosample_id), size = 2.5, show.legend = FALSE, max.overlaps = 5) +
  scale_color_manual(values = c(
    "Healthy"   = "#99cc33", 
    "Remission" = "#6699cc", 
    "Active"    = "#ff6666"
  )) +
  labs(
    x = paste0("PC1 (", vars_perc[1], "%)"),
    y = paste0("PC2 (", vars_perc[2], "%)"),
    color = "Condition"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top")

# Plot B: Cell Type Loadings (Right)
plot_right <- ggplot(pca_loadings, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(color = "black", size = 1.5) +
  geom_text_repel(data = pca_loadings %>% 
                    mutate(dist = sqrt(PC1^2 + PC2^2)) %>% 
                    top_n(25, dist), 
                  aes(label = cell_type_anno), 
                  size = 2.8, fontface = "italic", box.padding = 0.5) +
  labs(
    x = paste0("PC1 (", vars_perc[1], "%)"),
    y = paste0("PC2 (", vars_perc[2], "%)")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Combine panels
final_fig2d <- plot_left + plot_right + plot_annotation(tag_levels = 'A')
print(final_fig2d)

