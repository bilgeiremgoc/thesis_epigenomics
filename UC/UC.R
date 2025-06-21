library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

idat_dir <- "C:/Users/bilge/Documents/epigenomic_analysis/UC/GSE81211"

rgSet <- read.metharray.exp(base = idat_dir)

densityPlot(rgSet, main = "Raw Intensity Distribution", legend = FALSE)

baseDir <- "C:/Users/bilge/Documents/epigenomic_analysis/UC/GSE81211"  
RGset <- read.metharray.exp(base = baseDir, recursive = TRUE)
RGset

gsm_ids <- sapply(strsplit(colnames(RGset), "_"), `[`, 1)

group_mapping <- data.frame(
  GSM = c(
    "GSM2144712", "GSM2144713", "GSM2144714",
    "GSM2144715", "GSM2144716", "GSM2144717",
    "GSM2144718", "GSM2144719", "GSM2144720",
    "GSM2144721", "GSM2144722", "GSM2144723"
  ),
  Group = c(
    rep("Normal", 3),
    rep("UC", 8),
    "HCT116"
  )
)

group_vector <- group_mapping$Group[match(gsm_ids, group_mapping$GSM)]

pData(RGset)$Group <- factor(group_vector)

RGset <- RGset[, pData(RGset)$Group != "HCT116"]
pData(RGset)$Group <- droplevels(pData(RGset)$Group)

table(pData(RGset)$Group)

GRset <- preprocessQuantile(RGset)

mVals <- getM(GRset)

group <- pData(GRset)$Group
group <- relevel(group, ref = "Normal") 

design <- model.matrix(~group)
colnames(design)

fit <- lmFit(mVals, design)
fit <- eBayes(fit)

dmp_all <- topTable(fit, coef = "groupUC", number = Inf, adjust.method = "BH")
head(dmp_all)


ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

dmp_all$Gene <- ann450k[rownames(dmp_all), "UCSC_RefGene_Name"]
dmp_all$Region <- ann450k[rownames(dmp_all), "UCSC_RefGene_Group"]
dmp_all$Island <- ann450k[rownames(dmp_all), "Relation_to_Island"]

annot_sub <- ann450k[, c("chr", "pos", "strand")]

dmp_all$chr    <- annot_sub[rownames(dmp_all), "chr"]
dmp_all$pos    <- annot_sub[rownames(dmp_all), "pos"]
dmp_all$strand <- annot_sub[rownames(dmp_all), "strand"]

dmp_sig <- subset(dmp_all, 
                  adj.P.Val < 0.05 &
                    abs(logFC) > 1 &
                    !is.na(Gene))
dmp_teze_uygun <- dmp_sig %>%
  filter(P.Value < 0.001 & abs(logFC) > 1 & Gene != "" & !is.na(Gene))


dmp_sig_clean <- subset(dmp_all, 
                  adj.P.Val < 0.01 &
                    abs(logFC) > 1.5 &
                    !is.na(Gene))

top10_clean <- head(dmp_sig_clean[, c("logFC", "P.Value", "adj.P.Val", "Gene", "Region", "Island")], 10)
write.csv(top10_clean, "top10_uc_clean_cpg.csv", row.names = TRUE)


writexl::write_xlsx(dmp_sig_clean, "UC_teze_uygun.xlsx")



library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(dplyr)

bed_df <- dmp_sig_clean %>%
  mutate(start = pos - 100,
         end = pos + 100,
         name = rownames(.)) %>%
  select(chr, start, end, name, strand)


seqs <- getSeq(Hsapiens,
               names = bed_df$chr,
               start = bed_df$start,
               end = bed_df$end,
               strand = bed_df$strand)

names(seqs) <- bed_df$name

writeXStringSet(seqs, filepath = "UC_dmp_cpg_sequences.fasta")










