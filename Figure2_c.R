# --- 1. LOAD LIBRARIES ---
library(tidyverse)
library(ggpubr)
library(readxl)
library(data.table)

# --- 2. LOAD & CLEAN DATA ---
# Clinical Excel Data
clinical_data <- read_excel("supplementary_data_1.xlsx", sheet = "metadata")

# Cleaning clinical IDs and HPF counts
clinical_clean <- clinical_data %>%
  mutate(distal_hpf = as.numeric(case_when(
    `Eos/HPF (distal)` == "No tissue" ~ NA_character_,
    `Eos/HPF (distal)` == "rare"      ~ "0.5",
    TRUE ~ gsub("[^0-9.]", "", as.character(`Eos/HPF (distal)`))
  ))) %>%
  select(donor_id = ID, distal_hpf)

# Single Cell Metadata
# fread is faster for large .txt files
metadata <- fread("EoE_meta.txt", fill = TRUE)

# --- 3. AGGREGATE BY DONOR WITH CLR TRANSFORMATION & IMPUTATION ---
wide_counts <- metadata %>%
  filter(NAME != "TYPE", cell_type_anno != "") %>%
  group_by(donor_id, disease_status, cell_type_anno) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = cell_type_anno, values_from = Count, values_fill = 0)

# 3.2 Impute zeros (add 1) and calculate CLR
# We use a matrix for the math
count_matrix <- wide_counts %>% 
  select(-donor_id, -disease_status) %>% 
  as.matrix()

# Add pseudocount of 1 to handle zeros for the log math
imputed_matrix <- count_matrix + 1

# CLR = log(x / geometric_mean(x))
geom_mean_vector <- exp(rowMeans(log(imputed_matrix)))
clr_matrix <- log(imputed_matrix / geom_mean_vector)

# 3.3 Convert back to long format for ggplot and correlation
plot_data_transformed <- clr_matrix %>%
  as.data.frame() %>%
  bind_cols(wide_counts %>% select(donor_id, disease_status)) %>%
  pivot_longer(cols = -c(donor_id, disease_status), 
               names_to = "cell_type_anno", 
               values_to = "CLR_Value")

plot_data_2c <- plot_data_transformed %>%
  left_join(clinical_clean, by = "donor_id") %>%
  filter(!is.na(distal_hpf)) %>%
  # Update the label here
  mutate(disease_status = recode(disease_status, "Ctrl" = "Healthy"))

# --- 4. CALCULATE STATS FOR LABELS ---
target_cells <- c("Eosinophil", "Th2", "cDC2C (PRDM16)")

stats_labels <- plot_data_2c %>%
  filter(cell_type_anno %in% target_cells) %>%
  group_by(cell_type_anno) %>%
  summarise(
    # Use CLR_Value instead of Percentage to reflect compositional analysis
    res = list(cor.test(CLR_Value, distal_hpf, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(res, ~.$estimate),
    p = map_dbl(res, ~.$p.value),
    # This correctly performs the Benjamini-Hochberg FDR estimate
    fdr = p.adjust(p, method = "BH"),
    label = paste0("R = ", round(r, 2), "; FDR = ", signif(fdr, 2))
  )

# --- 5. FINAL PLOT ---
ggplot(plot_data_2c %>% filter(cell_type_anno %in% target_cells), 
       aes(x = distal_hpf, y = CLR_Value)) + 
  geom_point(aes(color = disease_status), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.6) +
  
  geom_text(data = stats_labels, aes(x = 0, y = Inf, label = label), 
            hjust = 0, vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  
  facet_wrap(~cell_type_anno, scales = "free_y") +
  
  scale_color_manual(values = c(
    "Healthy" = "#4DAF4A",
    "Remission" = "#377EB8", 
    "Active" = "#E41A1C"     
  )) +
  
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 12),
    legend.position = "bottom"
  ) +
  
  labs(
    x = "Eosinophils / HPF (distal)",
    y = "Abundance (CLR-transformed)", 
    color = "Status"
  )
