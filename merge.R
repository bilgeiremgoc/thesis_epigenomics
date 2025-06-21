library(readxl)
library(dplyr)

endo <- read_excel("endo/dmp_teze_uygun.xlsx")
ra <- read_excel("ra/dmp_teze_uygun.xlsx")
ss <- read_excel("SS/SS_teze_uygun.xlsx")
sle <- read_excel("sle/SLE_dmp_teze_uygun.xlsx")
uc <- read_excel("UC/UC_teze_uygun.xlsx")

endo_df <- dplyr::select(endo, Gene = UCSC_RefGene_Name, Endo_change = logFC)
ra_df   <- dplyr::select(ra, Gene = UCSC_RefGene_Name, RA_change = logFC)
ss_df   <- dplyr::select(ss, Gene = Gene, SS_change = deltaBeta)
sle_df  <- dplyr::select(sle, Gene = Gene, SLE_change = logFC)
uc_df   <- dplyr::select(uc, Gene = Gene, UC_change = logFC)


endo_df <- endo_df %>% group_by(Gene) %>% summarise(Endo_change = mean(Endo_change, na.rm = TRUE))
ra_df   <- ra_df   %>% group_by(Gene) %>% summarise(RA_change   = mean(RA_change, na.rm = TRUE))
ss_df   <- ss_df   %>% group_by(Gene) %>% summarise(SS_change   = mean(SS_change, na.rm = TRUE))
sle_df  <- sle_df  %>% group_by(Gene) %>% summarise(SLE_change  = mean(SLE_change, na.rm = TRUE))
uc_df   <- uc_df   %>% group_by(Gene) %>% summarise(UC_change   = mean(UC_change, na.rm = TRUE))

merged_epigenetic <- full_join(endo_df, ra_df, by = "Gene") %>%
  full_join(ss_df, by = "Gene") %>%
  full_join(sle_df, by = "Gene") %>%
  full_join(uc_df, by = "Gene")

merged_epigenetic_clean <- merged_epigenetic %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)

merged_epigenetic_clean$Gene <- sapply(strsplit(merged_epigenetic_clean$Gene, ";"), `[`, 1)

merged_epigenetic_clean[is.na(merged_epigenetic_clean)] <- 0


writexl::write_xlsx(merged_epigenetic_clean, path = "epigenetic_expression_matrix.xlsx")

ml_data <- ml_data[!duplicated(ml_data$Gene), ]
ml_data <- as.data.frame(ml_data)  # tibble → data.frame
rownames(ml_data) <- ml_data$Gene
ml_data <- ml_data[, -1]  

writexl::write_xlsx(ml_data, path = "ml_data.xlsx")


