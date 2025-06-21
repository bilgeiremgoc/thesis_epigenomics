BiocManager::install("ChIPseeker")
library(ChIPseeker)

BiocManager::install("GenomicRanges") 
library(GenomicRanges)

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") # İnsan genomu (hg38) gen anotasyonu için
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(dplyr)
library(readr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(BSgenome.Hsapiens.UCSC.hg38) # Genom dizilerini almak için
library(RCy3)

bed_file_path <- "C:/Users/bilge/Documents/epigenomic_analysis/ENCFF218WQM.bed"
endometriosis_data_df <- read_tsv(bed_file_path)

head(endometriosis_data_df)
str(endometriosis_data_df)
summary(endometriosis_data_df)


endometriosis_data_df_cleaned <- endometriosis_data_df %>%
  filter(!is.na(start) & !is.na(end))

message("Orijinal satır sayısı: ", nrow(endometriosis_data_df))
message("Temizlenmiş satır sayısı: ", nrow(endometriosis_data_df_cleaned))

endometriosis_elements_gr <- GRanges(
  seqnames = endometriosis_data_df_cleaned$`#chr`,
  ranges = IRanges(start = endometriosis_data_df_cleaned$start, end = endometriosis_data_df_cleaned$end),
  strand = "*",
  name = endometriosis_data_df_cleaned$name,
  class = endometriosis_data_df_cleaned$class,
  TargetGene = endometriosis_data_df_cleaned$TargetGene,
  TargetGeneEnsemblID = endometriosis_data_df_cleaned$TargetGeneEnsemblID,
  TargetGeneTSS = endometriosis_data_df_cleaned$TargetGeneTSS,
  CellType = endometriosis_data_df_cleaned$CellType,
  Score = endometriosis_data_df_cleaned$Score,
  DistanceToTSS = endometriosis_data_df_cleaned$DistanceToTSS,
  H3K27ac = endometriosis_data_df_cleaned$H3K27ac,
  Open = endometriosis_data_df_cleaned$Open,
  Cofactor = endometriosis_data_df_cleaned$Cofactor,
  Activity = endometriosis_data_df_cleaned$Activity,
  HiC_Contacts = endometriosis_data_df_cleaned$HiC_Contacts,
  HiC_FoldChange = endometriosis_data_df_cleaned$HiC_FoldChange
)

print(endometriosis_elements_gr)


enhancer_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, endometriosis_elements_gr)

print(head(enhancer_sequences))
message("Çekilen toplam dizi sayısı: ", length(enhancer_sequences))




# Enhancer bölgelerini genlere göre anote etme
# Gen anotasyonu için TxDb objesini kullanıyoruz
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# annotatePeak fonksiyonu, GRanges objelerini gen anotasyonu ile eşleştirir
peak_annotation <- annotatePeak(
  endometriosis_elements_gr,
  tssRegion = c(-3000, 3000), 
  TxDb = txdb,
  level = "gene"
)


# Anotasyon özetini görme
plotAnnoPie(peak_annotation)
plotAnnoBar(peak_annotation) 

# Anotasyon sonuçlarını data frame olarak al
annotation_df <- as.data.frame(peak_annotation)
head(annotation_df)

BiocManager::install(c("JASPAR2024", "TFBSTools", "universalmotif"))
BiocManager::install("motifmatchr")
BiocManager::install("JASPAR2022")

library(JASPAR2024) 
library(TFBSTools)
library(universalmotif) 
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)

opts <- list(species = 9606, all_versions = FALSE)

pfm_list <- getMatrixSet(JASPAR2022, opts)

motif_matches <- matchMotifs(pwms = pfm_list,
                             subject = endometriosis_elements_gr,
                             genome = BSgenome.Hsapiens.UCSC.hg38)

motif_summary <- motifMatches(motif_matches)


motif_matrix <- motifMatches(motif_matches)  # binary matrix (TRUE / FALSE)

# Her motifin toplam kaç bölgeyle eşleştiğini sayıyoruz
motif_counts <- colSums(motif_matrix)

# Data frame'e çeviriyoruz
motif_summary_df <- data.frame(
  motif_id = names(motif_counts),
  match_count = as.integer(motif_counts)
)

motif_info <- lapply(pfm_list, function(x) data.frame(
  motif_id = ID(x),
  tf_name = name(x)
))

motif_info_df <- do.call(rbind, motif_info)


motif_summary_annotated <- merge(motif_summary_df, motif_info_df, by = "motif_id")

# En çok eşleşen ilk 10 TF
top_motifs <- motif_summary_annotated[order(-motif_summary_annotated$match_count), ][1:10, ]
print(top_motifs)



# En çok eşleşen ilk 10 motif
top_motifs_50 <- motif_summary_annotated[order(-motif_summary_annotated$match_count), ][1:50, ]


writexl::write_xlsx(top_motifs_100, "top100_endo_motifs.xlsx")



ggplot(top_motifs_50, aes(x = reorder(tf_name, -match_count), y = match_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Motifs", x = "Motif", y = "Match Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))



promoters <- getPromoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream=1000, downstream=200)
motif_matches_promoters <- matchMotifs(pfm_list, promoters, genome = BSgenome.Hsapiens.UCSC.hg38)

motif_summary_annotated$normalized_count <- motif_summary_annotated$match_count / max(motif_summary_annotated$match_count)

ggplot(motif_summary_annotated, aes(x = reorder(tf_name, -normalized_count), y = normalized_count, fill = normalized_count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Normalized Motif Frequency", x = "Motif", y = "Normalized Match Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  theme(legend.position = "none")


# İlgili motifleri filtrele (örneğin en fazla eşleşenler)
motif_summary_annotated_filtered <- motif_summary_annotated[motif_summary_annotated$match_count > 50, ]

# Eşleşen motiflerin bar grafiği
ggplot(motif_summary_annotated_filtered, aes(x = reorder(tf_name, -match_count), y = match_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Motif Frequency (Filtered)", x = "Motif", y = "Match Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# Gene Ontology ve KEGG Pathway analizi yapmak için
library(clusterProfiler)
motif_genes <- motif_summary_annotated$TargetGene
go_results <- enrichGO(motif_genes, OrgDb = org.Hs.eg.db, ont = "BP")  # Biological Process
kegg_results <- enrichKEGG(motif_genes, organism = 'hsa')



promoters <- getPromoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 1000, downstream = 200)
motif_matches_promoters <- matchMotifs(pfm_list, promoters, genome = BSgenome.Hsapiens.UCSC.hg38)
