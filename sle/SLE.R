
library(limma)

beta_file <- "C:/Users/bilge/Documents/epigenomic_analysis/sle/GSE59250_average_betas.txt"

beta_matrix <- read.delim(beta_file, row.names = 1, check.names = FALSE)

dim(beta_matrix)
head(colnames(beta_matrix))

group_labels <- ifelse(grepl("^SLE", colnames(beta_matrix)), "SLE", "Control")
group_labels <- factor(group_labels)

table(group_labels)

beta_to_m <- function(beta) log2(beta / (1 - beta))
m_matrix <- beta_to_m(beta_matrix)

design <- model.matrix(~0 + group_labels)
colnames(design) <- levels(group_labels)

fit <- lmFit(m_matrix, design)
contrast <- makeContrasts(SLEvsControl = SLE - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

dmp_all <- topTable(fit2, coef = "SLEvsControl", number = Inf, adjust.method = "BH")
head(dmp_all)


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

dmp_all$Gene <- ann450k[rownames(dmp_all), "UCSC_RefGene_Name"]
dmp_all$Region <- ann450k[rownames(dmp_all), "UCSC_RefGene_Group"]
dmp_all$Island <- ann450k[rownames(dmp_all), "Relation_to_Island"]
dmp_all$chr <- ann450k[rownames(dmp_all), "chr"]
dmp_all$pos <- ann450k[rownames(dmp_all), "pos"]
dmp_all$strand <- ann450k[rownames(dmp_all), "strand"]


dmp_sig <- subset(dmp_all, adj.P.Val < 0.01 & abs(logFC) > 1 & !is.na(Gene))
nrow(dmp_sig)

library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(dplyr)

dmp_teze_uygun <- dmp_sig %>%
  filter(adj.P.Val < 0.01, abs(logFC) > 1.5, !is.na(Gene)) %>%
  mutate(start = pos - 100,
         end = pos + 100,
         name = rownames(.)) %>%
  select(chr, start, end, name, strand)

seqs <- getSeq(Hsapiens,
               names = dmp_teze_uygun$chr,
               start = dmp_teze_uygun$start,
               end = dmp_teze_uygun$end,
               strand = dmp_teze_uygun$strand)

names(seqs) <- dmp_teze_uygun$name

writeXStringSet(seqs, filepath = "SLE_dmp_cpg_sequences.fasta")

writexl::write_xlsx(dmp_teze_uygun, "SLE_dmp_teze_uygun.xlsx")

