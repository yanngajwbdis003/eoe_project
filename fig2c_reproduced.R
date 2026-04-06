# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl) 

# --- 2. LOAD & CLEAN CLINICAL EXCEL DATA ---
# This version prevents the "NAs introduced by coercion" warning
clinical_data <- read_excel("supp_data_1.xlsx", sheet = "metadata")

clinical_clean <- clinical_data %>%
  mutate(distal_hpf = as.numeric(case_when(
    `Eos/HPF (distal)` == "No tissue" ~ NA_character_,
    `Eos/HPF (distal)` == "rare"      ~ "0.5",
    `Eos/HPF (distal)` == "N/A"       ~ NA_character_,
    # Clean string characters first (remove ">"), then the result becomes numeric
    TRUE ~ gsub(">", "", as.character(`Eos/HPF (distal)`))
  ))) %>%
  select(ID, distal_hpf)

# --- 3. AGGREGATE CELL DATA ---
# Ensure donor_id is included so we can merge with the clinical 'ID'
cell_counts <- metadata %>%
  filter(NAME != "TYPE") %>%
  group_by(donor_id, biosample_id, disease_status, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(TotalCells = sum(Count),
         Percentage = (Count / TotalCells) * 100) %>%
  ungroup()

# --- 4. MERGE DATA ---
plot_data_2c <- cell_counts %>%
  left_join(clinical_clean, by = c("donor_id" = "ID")) %>%
  filter(!is.na(distal_hpf))

# Targets for Figure 2c
target_cells_2c <- c("Eosinophil", "Th2", "cDC2C (PRDM16)")

# --- 5. GENERATE FIGURE 2c ---
ggplot(plot_data_2c %>% filter(cell_type_anno %in% target_cells_2c), 
       aes(x = distal_hpf, y = Percentage)) +
  # Use actual points colored by condition
  geom_point(aes(color = disease_status), size = 2.5, alpha = 0.8) +
  # Add regression line
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.8) +
  # Add Spearman correlation (Matches the paper's math exactly)
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  # Formatting to match the paper
  facet_wrap(~cell_type_anno, scales = "free_y") +
  scale_color_manual(values = c("Ctrl" = "#4DAF4A", "Remission" = "#377EB8", "Active" = "#E41A1C")) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 11)
  ) +
  labs(
    title = "Figure 2c:Spearman Correlation with Clinical Eos/HPF Counts",
    x = "Eosinophils / HPF (distal)",
    y = "Abundance (%)",
    color = "Status"
  )