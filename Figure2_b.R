# --- 1. LOAD LIBRARIES ---
library(tidyverse)
library(ggrepel)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

# --- 2. LOAD & CLEAN CLINICAL DATA ---
clinical_data <- read_excel("supplementary_data_1.xlsx", sheet = "metadata")

clinical_clean <- clinical_data %>%
  mutate(
     ID = as.character(ID),
     distal_hpf = as.numeric(case_when(
    `Eos/HPF (distal)` == "No tissue" ~ NA_character_,
    `Eos/HPF (distal)` == "rare"      ~ "0.5",
    # This regex removes EVERYTHING except numbers and decimals 
    # (handles ">", "<", "-", or spaces)
    TRUE ~ str_replace_all(as.character(`Eos/HPF (distal)`), "[^0-9.]", "")
  ))) %>%
  select(ID, distal_hpf)

# --- 3. AGGREGATE DATA BY DONOR ---
metadata <- fread("EoE_meta.txt", fill = TRUE)
stats_data_2b <- metadata %>%
  filter(NAME != "TYPE", cell_type_anno != "") %>%
  mutate(donor_id = as.character(donor_id)) %>%
  group_by(donor_id, biosample_id, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(TotalCells = sum(Count),
         Percentage = (Count / TotalCells) * 100) %>%
  ungroup() %>%
  left_join(clinical_clean, by = c("donor_id" = "ID")) %>%
  filter(!is.na(distal_hpf))

# --- 4. CALCULATE CORRELATION AND FDR ---
stats_2b_final <- stats_data_2b %>%
  group_by(cell_type_anno) %>%
  summarise(
    # Use cor.test to get p-values directly
    test = list(cor.test(Percentage, distal_hpf, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(
    rho = map_dbl(test, "estimate"),
    p_val = map_dbl(test, "p.value"),
    fdr = p.adjust(p_val, method = "BH")
  ) 
  # Filter to show only the significant or nearly significant ones if desired
  # filter(fdr < 0.1)

# --- 5. DEFINE LINEAGES (CRITICAL FOR MATCHING TEXT) ---
stats_2b_final <- stats_2b_final %>%
  mutate(lineage = case_when(
    cell_type_anno %in% c("Eosinophil", "Mast cell", "Mast cell (cycling)") ~ "Granulocyte",
    # Myeloid (Teal) - includes DCs and Macrophages
    str_detect(cell_type_anno, "Macrophage|DC|Monocyte|Langerhans|cDC") ~ "Myeloid",
    # Lymphocyte (Red) - includes T, B, NK, and ILCs
    str_detect(cell_type_anno, "CD4|CD8|Th|Treg|B cell|Plasma|ILC|NK") ~ "Lymphocyte",
    cell_type_anno %in% c("Apical cell", "Basal cell (cycling)", "Suprabasal") ~ "Epithelial",
    TRUE ~ "Other"
  ))

# --- 6. GENERATE THE VOLCANO PLOT ---
# 1. Define labels
labels_to_show <- c("Eosinophil", "Th2", "cDC2C (PRDM16)", "pDC", 
                    "Mast cell (cycling)", "ILC2", "Langerhans cell", 
                    "Macrophage (PLAC8)", "Apical cell")

# 2. Build the plot
ggplot(stats_2b_final, aes(x = rho, y = fdr, color = lineage)) + 
  # 1. Thinner, more transparent points to see overlapping data
  geom_point(alpha = 0.6, size = 1.8) +
  
  # 2. Reverse scale handles the '1.0 at the top' logic more robustly
  scale_y_continuous(limits = c(1, 0),
                     breaks = seq(0, 1, 0.25)) +
  
  # 3. Match the specific color hex codes from the original figure
  scale_color_manual(values = c(
    "Myeloid"    = "#80CED7", # Teal-ish
    "Granulocyte" = "#FFB347", # Golden Orange
    "Lymphocyte"  = "#C3B1E1", # Soft Purple
    "Epithelial"  = "#A2B9D6", # Steel Blue
    "Other"       = "#D3D3D3"  # Light Grey
  )) +
  
  # 4. Text repel adjustments for better legibility
  geom_text_repel(data = filter(stats_2b_final, cell_type_anno %in% labels_to_show),
                  aes(label = cell_type_anno), 
                  color = "black", 
                  size = 3.2, 
                  box.padding = 0.6,
                  segment.size = 0.3,
                  segment.color = 'grey70',
                  fontface = "italic",
                  min.segment.length = 0) + 
  
  theme_classic() + 
  # 5. This is the "secret sauce" - forcing a square/portrait aspect ratio
  theme(
    aspect.ratio = 1.2,
    legend.position = "right",
    axis.title = element_text(size = 10)
  ) +
  labs(
    x = "Spearman correlation coefficients between the\nfraction of each cell type and # eosinophils/HPF",
    y = "FDR (correlation analysis)",
    color = NULL
  )

