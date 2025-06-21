library(minfi)
library(limma)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

gset <- getGEO("GSE73950", GSEMatrix = TRUE)


idat_dir <- "C:/Users/bilge/Documents/epigenomic_analysis/endo/GSE73950"

rgSet <- read.metharray.exp(base = idat_dir)

qcReport(rgSet, pdf = "qc_report.pdf")


densityPlot(rgSet, main = "Raw Intensity Distribution", legend = FALSE)

sample_ids <- c(
  "GSM1906605", "GSM1906606", "GSM1906607", "GSM1906608", "GSM1906609", "GSM1906610", "GSM1906611",
  "GSM1906612", "GSM1906613", "GSM1906614", "GSM1906615", "GSM1906616", "GSM1906617", "GSM1906618",
  "GSM1906619", "GSM1906620", "GSM1906621", "GSM1906622", "GSM1906623", "GSM1906624", "GSM1906625",
  "GSM1906626", "GSM1906628", "GSM1906629", "GSM1906630", "GSM1906631", "GSM1906632", "GSM1906633",
  "GSM1906634", "GSM1906635", "GSM1906636", "GSM1906637", "GSM1906638", "GSM1906639", "GSM1906640",
  "GSM1906641", "GSM1906642", "GSM1906643", "GSM1906644", "GSM1906645", "GSM1906646", "GSM1906647",
  "GSM1906648"
)

sample_labels <- c(
  "E9_1_ES", "E12_1_ES", "E13_1_ES", "E11_1_ES", "E8_1_ES", "E10_1_ES", "E7_1_ES",
  "E14_1_MS", "E20_1_MS", "E21_1_MS", "E18_1_MS", "E17_1_MS", "E22_1_MS", "E15_1_MS_r1",
  "E15_1_MS_r2", "E19_1_MS_r1", "E19_1_MS_r2", "E23_1_LS", "E27_1_LS", "E29_1_LS", "E31_1_LS",
  "E28_1_LS", "E24_1_LS", "E26_1_LS", "E30_1_LS", "E25_1_LS", "H14_1_MS", "H15_1_MS", "H6_1_MS",
  "H12_1_MS", "H8_1_MS", "H11_1_MS", "H20_1_MS", "H13_1_MS", "H7_1_MS", "H21_1_MS", "H18_1_MS",
  "H9_1_MS", "H10_1_MS", "H16_1_MS", "H17_1_MS", "H22_1_MS", "H19_1_MS"
)

group <- ifelse(grepl("^E", sample_labels), "Endo", "Control")
group <- factor(group)
table(group)

mSet <- preprocessIllumina(rgSet) 

betaVals <- getBeta(mSet)
mVals <- getM(mSet)

design <- model.matrix(~group)
colnames(design)


library(limma)

fit <- lmFit(mVals, design)
fit2 <- eBayes(fit)

top_dmp <- topTable(fit2, coef = 2, number = Inf, adjust = "BH")
head(top_dmp)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

top_dmp$CpG <- rownames(top_dmp)

top_dmp_annotated <- merge(top_dmp, annotation, by.x = "CpG", by.y = "")

head(top_dmp_annotated)

top_dmp_df <- as.data.frame(top_dmp_annotated)

top_dmp_clean <- top_dmp_df %>%
  filter(!is.na(UCSC_RefGene_Name))

top_dmp_sig <- subset(top_dmp_clean, adj.P.Val < 0.05)

top_dmp_loose <- top_dmp_clean %>%
  filter(adj.P.Val < 0.1)

dmp_teze_uygun <- top_dmp_df %>%
  filter(P.Value < 0.001 & abs(logFC) > 1 & UCSC_RefGene_Name != "" & !is.na(UCSC_RefGene_Name))

nrow(dmp_teze_uygun)
head(dmp_teze_uygun)

writexl::write_xlsx(dmp_teze_uygun, "dmp_teze_uygun.xlsx")

first10 <- dmp_teze_uygun %>%
  select(CpG, logFC, P.Value, UCSC_RefGene_Name, Relation_to_Island) %>%
  arrange(P.Value) %>%
  head(10)

print(first10)


write.csv(first10$UCSC_RefGene_Name, "gen_isimleri_first10.csv")









BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

library(Biostrings)
library(dplyr)

bed_df <- dmp_teze_uygun %>%
  mutate(start = pos - 100,
         end = pos + 100,
         name = CpG) %>%
  select(chr, start, end, name, strand)

seqs <- getSeq(genome,
               names = bed_df$chr,
               start = bed_df$start,
               end = bed_df$end,
               strand = bed_df$strand)

names(seqs) <- bed_df$name
writeXStringSet(seqs, filepath = "dmp_cpg_sequences.fasta")
