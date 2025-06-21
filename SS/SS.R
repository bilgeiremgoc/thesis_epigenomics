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











df <- read_excel("SS/SS_teze_uygun.xlsx")

gene_column <- df$Gene

genes <- unlist(strsplit(gene_column, ";"))
genes <- unique(trimws(genes))
genes <- genes[genes != ""]

gene_entrez <- bitr(genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

go_results <- enrichGO(gene = gene_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP", # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

head(go_results)
dotplot(go_results, showCategory = 15)

go_results_df <- as.data.frame(go_results)

writexl::write_xlsx(go_results_df, "SS_GO_results.xlsx")


kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID,
                           organism = 'hsa',
                           pvalueCutoff = 0.05)
dotplot(kegg_results, showCategory = 15)

KEGG_results_df <- as.data.frame(kegg_results)

writexl::write_xlsx(KEGG_results_df, "SS_kegg_results.xlsx")



reactome_results <- enrichPathway(gene = gene_entrez$ENTREZID,
                                  organism = "human",
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)

dotplot(reactome_results, showCategory = 20)


reactom_results_df <- as.data.frame(reactome_results)

writexl::write_xlsx(reactom_results_df, "SS_reactome_results.xlsx")







