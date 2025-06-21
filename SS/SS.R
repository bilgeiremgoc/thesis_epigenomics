dmp_pss <- read.delim("C:/Users/bilge/Documents/epigenomic_analysis/SS/GSE146116_DMP_ALL_SS.VS.HC_BHadjust.txt.gz", header = TRUE)

str(dmp_pss)
head(dmp_pss)
colnames(dmp_pss)

dmp_pss_clean <- subset(dmp_pss,
                        adj.P.Val < 0.01 &
                          abs(deltaBeta) > 0.3)

nrow(dmp_pss_clean) 


library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ann_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

dmp_pss_clean$Gene <- ann_epic[dmp_pss_clean$probeID, "UCSC_RefGene_Name"]
dmp_pss_clean$Region <- ann_epic[dmp_pss_clean$probeID, "UCSC_RefGene_Group"]
dmp_pss_clean$Island <- ann_epic[dmp_pss_clean$probeID, "Relation_to_Island"]
dmp_pss_clean$chr <- ann_epic[dmp_pss_clean$probeID, "chr"]
dmp_pss_clean$pos <- ann_epic[dmp_pss_clean$probeID, "pos"]
dmp_pss_clean$strand <- ann_epic[dmp_pss_clean$probeID, "strand"]


dmp_sig <- subset(dmp_pss_clean, 
                  adj.P.Val < 0.05 &
                    abs(deltaBeta) > 0.3 &
                    !is.na(Gene))

dmp_teze_uygun <- dmp_sig %>%
  filter(adj.P.Val < 0.01 & abs(deltaBeta) > 0.3 & Gene != "" & !is.na(Gene))


top10_clean <- head(dmp_teze_uygun[, c("deltaBeta", "P.Value", "adj.P.Val", "Gene", "Region", "Island")], 10)

write.csv(top10_clean, "SS_top10_pss_clean_cpg.csv", row.names = TRUE)


writexl::write_xlsx(dmp_teze_uygun, "SS_teze_uygun.xlsx")


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(dplyr)

bed_df <- dmp_teze_uygun %>%
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

writeXStringSet(seqs, filepath = "SS_dmp_cpg_sequences.fasta")









