if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma", "minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

beta_raw <- read.table("GSE87621_Endomtriosis_450K_raw_signal_matrix.txt.gz", 
                       header = TRUE, row.names = 1, sep = "\t")

meth <- beta_raw[, grep("Methylated_Signal$", colnames(beta_raw))]
unmeth <- beta_raw[, grep("Unmethylated_Signal$", colnames(beta_raw))]

# 3. Beta değerlerinin hesaplanması (metilasyon oranı)
beta <- meth / (meth + unmeth + 100)  # +100 küçük sayı stabilite için

colnames(beta) <- gsub("_Methylated_Signal", "", colnames(beta))

group <- factor(c(rep("OESC", 4), rep("CESC", 5)))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

library(limma)

beta_adj <- pmin(pmax(beta, 0.001), 0.999)
Mval <- log2(beta_adj / (1 - beta_adj))

# Modelleme
fit <- lmFit(Mval, design)
contrast.matrix <- makeContrasts(OESC - CESC, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Sonuçları al
top_dmp <- topTable(fit2, number = Inf, adjust = "BH")
head(top_dmp)
significant_probes <- top_dmp$ID[top_dmp$adj.P.Val < 0.05]





library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

library(dplyr)
library(tidyr)

# 1. Anlamlı probe ID’ler
significant_probes <- rownames(top_dmp)[top_dmp$adj.P.Val < 0.05]

# 2. Ortak probe ID’ler ile anotasyon
common_probes <- intersect(significant_probes, rownames(anno))
annotated_probes <- anno[common_probes, ]

# 3. DataFrame’e dönüştür ve probe ID’yi sütun yap
df <- as.data.frame(annotated_probes)
df$Name <- rownames(annotated_probes)

# 4. Gen isimlerini çıkar
genes <- df %>%
  filter(UCSC_RefGene_Name != "") %>%
  separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  distinct(UCSC_RefGene_Name)

head(genes)



writexl::write_xlsx(df, "dna_methylation-endo_results.xlsx")





gene_symbols <- genes$UCSC_RefGene_Name  

gene_symbols <- gene_symbols[gene_symbols != "NA"]

library(clusterProfiler)
library(org.Hs.eg.db)

entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

head(entrez_ids)

# KEGG yolu zenginleştirme analizi
kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

# Sonuçları tablo olarak görelim
head(kegg_enrich)

kegg_enrich_df <- as.data.frame(kegg_enrich)


writexl::write_xlsx(kegg_enrich_df, "epi-endo-kegg-pathway.xlsx")

# GO zenginleştirme analizi (BP - biyolojik süreç)
go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)

head(go_enrich)

go_enrich_df <- as.data.frame(go_enrich)

library(enrichplot)

barplot(kegg_enrich, showCategory=10, title="KEGG Pathway Enrichment")

dotplot(go_enrich, showCategory=10, title="GO Biological Process Enrichment")


library(pheatmap)

# Anlamlı CpG'lerin ID'leri
sig_cpg <- rownames(top_dmp)[top_dmp$adj.P.Val < 0.05]

# Beta değerlerinden sadece anlamlı CpG’leri seç
beta_sig <- beta[sig_cpg, ]

# Örnek isimlerini sadeleştir (grup isimleri)
sample_annotation <- data.frame(Group = group)
rownames(sample_annotation) <- colnames(beta_sig)

pheatmap(beta_sig, 
         annotation_col = sample_annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Heatmap of Significant CpG Methylation")




table(df$UCSC_RefGene_Group)
barplot(table(df$UCSC_RefGene_Group), las=2, main="CpG Region Distribution")


barplot(table(df$Relation_to_Island), las=2, main="CpG Island Distribution")


pca <- prcomp(t(beta_sig))
plot(pca$x[,1:2], col=group, pch=19, main="PCA of Beta Values")
legend("center", legend=levels(group), col=1:length(group), pch=19)


library(EnhancedVolcano)
EnhancedVolcano(top_dmp,
                lab = rownames(top_dmp),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Volcano Plot of CpG Sites',
                pCutoff = 0.05)


emapplot(go_enrich, showCategory = 20)
cnetplot(go_enrich, showCategory = 10)

genes_freq <- table(genes$UCSC_RefGene_Name)
sort(genes_freq, decreasing = TRUE)[1:10]


# CpG başına gen ismini ekleyelim
df$GENE <- rownames(df)
cpg_gene_map <- df %>%
  dplyr::filter(UCSC_RefGene_Name != "") %>%
  tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  dplyr::select(GENE, UCSC_RefGene_Name)

# Her gen için en anlamlı CpG'yi seçelim
top_cpg_per_gene <- top_dmp %>%
  rownames_to_column("GENE") %>%
  inner_join(cpg_gene_map, by = "GENE") %>%
  group_by(UCSC_RefGene_Name) %>%
  slice_min(order_by = adj.P.Val, n = 1)

# Örnek: hub genlerin listesi
hub_genes <- 


# Epigenetik analizle kesişim:
epigenetik_genler <- unique(top_cpg_per_gene$UCSC_RefGene_Name)
common_genes <- intersect(hub_genes, epigenetik_genler)




