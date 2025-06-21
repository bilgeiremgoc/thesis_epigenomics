library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

idat_dir <- "C:/Users/bilge/Documents/epigenomic_analysis/ra/E-MTAB-8085_"

rgSet <- read.metharray.exp(base = idat_dir)

densityPlot(rgSet, main = "Raw Intensity Distribution", legend = FALSE)

mSet <- preprocessIllumina(rgSet)
betaVals <- getBeta(mSet)
mVals <- getM(mSet)


# 450K CpG'lerle sınırlamak için
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
common_cpgs <- intersect(rownames(betaVals), anno_450k$Name)
beta_450k <- betaVals[common_cpgs, ]

targets <- data.frame(
  Sample_Name = c("RA6731", "RA5384", "RA5815", "RA6426", "RA6596",
                  "H5491", "H5531", "H6102", "H6430"),
  Slide = rep("9374342150", 9),
  Array = c("R02C02", "R03C01", "R05C01", "R06C02", "R06C01",
            "R02C01", "R03C02", "R01C02", "R04C02"),
  Group = c(rep("RA", 5), rep("Control", 4))
)


targets$Slide_Array <- c(
  "9374342150_R01C02",  # H6102
  "9374342150_R02C01",  # H5491
  "9374342150_R03C02",  # H5531
  "9374342150_R04C02",  # H6430
  "9374342150_R02C02",  # RA6731
  "9374342150_R03C01",  # RA5384
  "9374342150_R05C01",  # RA5815
  "9374342150_R06C02",  # RA6426
  "9374342150_R06C01"   # RA6596
)


common_samples <- intersect(colnames(mVals), targets$Slide_Array)

mVals_filtered <- mVals[, common_samples]
targets_filtered <- targets[targets$Slide_Array %in% common_samples, ]

targets_filtered <- targets_filtered[match(common_samples, targets_filtered$Slide_Array), ]

group <- factor(targets_filtered$Group)  # "RA" ve "Control"
design <- model.matrix(~group)

fit <- lmFit(mVals_filtered, design)
fit2 <- eBayes(fit)
top_dmp <- topTable(fit2, coef = 2, number = Inf, adjust = "BH")


anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
top_dmp$CpG <- rownames(top_dmp)

top_annotated <- merge(top_dmp, anno, by.x = "CpG", by.y = "Name")

top_dmp_df <- as.data.frame(top_annotated)

top_dmp_clean <- top_dmp_df %>%
  filter(!is.na(UCSC_RefGene_Name))


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


write.csv(first10$UCSC_RefGene_Name, "ra_gen_isimleri_first10.csv")

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





