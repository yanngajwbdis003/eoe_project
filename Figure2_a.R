# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

metadata <- fread("EoE_meta.txt", fill = TRUE)

# --- 2. AGGREGATE DATA TO PATIENT PROPORTIONS ---
cell_counts <- metadata %>%
  filter(NAME != "TYPE") %>%
  group_by(biosample_id, disease_status, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop") %>%
  
  complete(nesting(biosample_id, disease_status), 
           cell_type_anno, 
           fill = list(Count = 0)) %>%
  group_by(biosample_id) %>%
  mutate(TotalCells = sum(Count),
         Percentage = (Count / TotalCells) * 100) %>%
  ungroup()

# Force the factor levels to match the paper's progression (Healthy -> Remission -> Active)
cell_counts$disease_status <- factor(cell_counts$disease_status, 
                                     levels = c("Ctrl", "Remission", "Active"))

target_cells <- c(
  # --- Clinical Markers ---
  "Eosinophil", 
  "Mast cell (cycling)",
  
  # --- The "Invasion" (Increased in Active) ---
  "Macrophage (PLAC8)",       # This is the PLAC8+ macrophage
  "Macrophage (ALOX15)",      # This is the ALOX15+ macrophage
  "cDC2C (PRDM16)",           # This is the PRDM16+ DC
  "pDC",                      # Plasmacytoid DCs
  "Th2",                      # Major source of IL-4/IL-13
  "CD4+ T cell (cycling)",    # Cycling CD4+ T cells
  "Plasma cell (IgG+)",       # Specifically highlighted in the B cell section
  "Plasma cell (IgM+)",
  "Plasma cell (cycling)",
  
  # --- The "Attrition" (Decreased in Active) ---
  "Macrophage (SEPP1/FOLR2)", # This is the FOLR2+ macrophage
  "Apical cell",              # The top layer of the epithelium
  
  # --- Tissue Recovery Markers (The "Mismatch" ones) ---
  "CD4+ Trm (CCL5)", 
  "Treg", 
  "Th17",
  "Basal cell (cycling)"
)

# --- 4. GENERATE THE REPRODUCED FIGURE 2a ---
ggplot(cell_counts %>% filter(cell_type_anno %in% target_cells), 
       aes(x = disease_status, y = Percentage, fill = disease_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", size = 0.4) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5, color = "grey20") +
  # Faceting allows each cell type to have its own y-axis scale
  facet_wrap(~cell_type_anno, scales = "free_y", ncol = 3) + 
  
  # 1. Update the Legend labels
  scale_fill_manual(
    values = c("Ctrl" = "#4DAF4A", "Remission" = "#377EB8", "Active" = "#E41A1C"),
    labels = c("Ctrl" = "Healthy", "Remission" = "Remission", "Active" = "Active")
  ) +
  
  # 2. Update the X-axis labels
  scale_x_discrete(
    labels = c("Ctrl" = "Healthy", "Remission" = "Remission", "Active" = "Active")
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Figure 2a Reproduction: Cellular Composition Shifts",
    subtitle = "Comparing Clinical Success (Eosinophils) vs. Tissue Failure (T-cells)",
    y = "Percent of Total Cells (%)",
    x = "Disease Condition",
    fill = "Status" # Clean up the legend title
  )

