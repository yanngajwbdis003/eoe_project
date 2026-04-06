# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(tidyr)
library(compositions) # Now that it's installed
library(ggrepel)


metadata <- fread("EoE_meta.txt", fill = TRUE)

# --- 2. CREATE cell_counts (The Patient-Level Table) ---
# This aggregates the single-cell data into the percentages used in Figure 2a
cell_counts <- metadata %>%
  filter(NAME != "TYPE") %>%
  group_by(biosample_id, disease_status, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(TotalCells = sum(Count),
         Percentage = (Count / TotalCells) * 100) %>%
  ungroup()

# --- 3. PREPARE FOR PCA (Figure 2d) ---
# The paper uses a "sample-cell type count profile matrix" 
# after centered log-ratio (CLR) transformation

# --- FIX: Filter out empty strings before pivoting ---
count_matrix <- cell_counts %>%
  filter(cell_type_anno != "" & !is.na(cell_type_anno)) %>% # This stops the spec$.name error
  select(biosample_id, cell_type_anno, Count) %>%
  pivot_wider(names_from = cell_type_anno, values_from = Count, values_fill = 0) %>%
  tibble::column_to_rownames("biosample_id")

# --- CLR & PCA (Fig 2d) ---
clr_matrix <- clr(count_matrix + 1)
pca_results <- prcomp(clr_matrix)

pca_samples <- as.data.frame(pca_results$x) %>%
  tibble::rownames_to_column("biosample_id") %>%
  left_join(distinct(metadata, biosample_id, disease_status), by = "biosample_id")

# View the first few rows of your new PCA coordinates
head(pca_samples)


## PLOTTING
library(ggplot2)

# Create Figure 2d (Left Panel)
# --- 1. PREPARE THE DATA ---
# (Ensure you've run the PCA and 'pca_samples' code from before)

# --- 2. CREATE THE FINAL FIGURE 2D ---
ggplot(pca_samples, aes(x = PC1, y = PC2, color = disease_status)) +
  # Add the dotted crosshair lines first (so they are behind the points)
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  
  # Add the patient sample points
  geom_point(size = 3.5, alpha = 0.8) +
  
  # Apply the paper's colors
  scale_color_manual(values = c(
    "Ctrl" = "#4DAF4A",      # Green
    "Remission" = "#377EB8", # Blue
    "Active" = "#E41A1C"     # Red
  )) +
  
  # Match the paper's exact axis labels
  labs(
    title = "Figure 2d: Sample Composition Profiles",
    x = "PC1 (24.4%)",
    y = "PC2 (8.6%)",
    color = "Condition"
  ) +
  
  # Clean, professional theme
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "top"
  )


# 1. Extract the cell type loadings (how much each cell contributes to PC1 and PC2)
pca_loadings <- as.data.frame(pca_results$rotation) %>%
  tibble::rownames_to_column("cell_type_anno")

# 2. Plot Figure 2d (Right Panel)
ggplot(pca_loadings, aes(x = PC1, y = PC2)) +
  geom_point(color = "black", alpha = 0.5) +
  # Label the "Famous" cell types mentioned in the paper
  geom_text_repel(data = filter(pca_loadings, 
                                abs(PC1) > 0.15 | abs(PC2) > 0.15),
                  aes(label = cell_type_anno), size = 3) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Figure 2d (Right): Cell Type Loadings",
       x = "PC1 (Cell types that define 'Active' vs 'Ctrl')",
       y = "PC2")