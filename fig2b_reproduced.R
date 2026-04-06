# --- 1. LOAD LIBRARIES ---
library(tidyverse)
library(ggrepel)
library(readxl)

# --- 2. LOAD & CLEAN CLINICAL DATA ---
clinical_data <- read_excel("supp_data_1.xlsx", sheet = "metadata")

clinical_clean <- clinical_data %>%
  mutate(distal_hpf = as.numeric(case_when(
    `Eos/HPF (distal)` == "No tissue" ~ NA_character_,
    `Eos/HPF (distal)` == "rare"      ~ "0.5",
    # This regex removes EVERYTHING except numbers and decimals 
    # (handles ">", "<", "-", or spaces)
    TRUE ~ str_replace_all(as.character(`Eos/HPF (distal)`), "[^0-9.]", "")
  ))) %>%
  select(ID, distal_hpf)

# --- 3. AGGREGATE DATA BY DONOR ---
stats_data_2b <- metadata %>%
  filter(NAME != "TYPE", cell_type_anno != "") %>%
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
    rho = cor(Percentage, distal_hpf, method = "spearman", use = "complete.obs"),
    p_val = if(n() >= 3 && sd(Percentage) > 0) {
      cor.test(Percentage, distal_hpf, method = "spearman", exact = FALSE)$p.value
    } else { NA_real_ }
  ) %>%
  filter(!is.na(p_val)) %>%
  mutate(
    fdr = p.adjust(p_val, method = "BH"),
    log_fdr = -log10(fdr)
  )

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
labels_to_show <- c("Eosinophil", "Th2", "cDC2C (PRDM16)", "pDC", 
                    "Mast cell (cycling)", "ILC2", "Langerhans cell", 
                    "Macrophage (PLAC8)", "Apical cell")

ggplot(stats_2b_final, aes(x = rho, y = log_fdr, color = lineage)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_vline(xintercept = c(-0.48, 0.48), linetype = "dashed", color = "blue", linewidth = 0.8) +
  scale_y_continuous(breaks = c(0, 1.0, 2.0, 3.0), labels = c("1.0", "0.1", "0.01", "0.001")) +
  geom_text_repel(data = filter(stats_2b_final, cell_type_anno %in% labels_to_show),
                  aes(label = cell_type_anno), color = "black", size = 3.5, fontface = "italic") +
  scale_color_manual(values = c("Granulocyte" = "#E7B800", "Lymphocyte" = "#FC4E07", 
                                "Myeloid" = "#00AFBB", "Epithelial" = "#2E9FDF", "Other" = "grey")) +
  theme_bw() +
  labs(x = "Spearman correlation coefficients (ρ)", y = "FDR (Significance)", color = "Lineage")