
setRepositories(ind=1:8)  # Adjust the range depending on how many repos are available

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("vsn")

### Loaded Libraries
library(DESeq2)           
library(ReportingTools)
library(Glimma)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggpmisc)
library(multcomp)
library(tximport)
library(tximportData)
library(readr)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)
library(variancePartition)
library(EnhancedVolcano)
library(fgsea)
library(data.table)
#install with BiocManager::install("reactome.db")
library(reactome.db)
library(vsn)

#############################################################################################################
### make tx2gene from gtf
### use version 35? release 101
TxDb.2 <- makeTxDbFromGFF(file = "Homo_sapiens.GRCh38.101.gtf")       ### use this
###
k.2 <- keys(TxDb.2, keytype = "TXNAME")
###
tx2gene.2 <- AnnotationDbi::select(TxDb.2, k.2, "GENEID", "TXNAME")
gtf <- rtracklayer::import('Homo_sapiens.GRCh38.101.gtf')
gtf_df=as.data.frame(gtf)
gtf_df=subset(gtf_df, select=c("gene_id", "gene_name")) %>% distinct(gene_id, .keep_all = TRUE)
rm(gtf)

#############################################################################################################
### Read in salmon *.quant files
samples=read.delim("samples.txt", stringsAsFactors = TRUE)
### insert metrics into sample file
### 
files=file.path("../results/ur_rna4_results/ur_rna4_salmon/", samples$sample_id_2, "quant.sf")
file.exists(files)
names(files) <- paste0(samples$sample_id_2)
###
txi.or_rna.2 <- tximport(files, type = "salmon", tx2gene = tx2gene.2, importer = read_tsv, ignoreTxVersion=TRUE)
#test=as.data.frame(txi.or_rna.2)
###
### Create dds file
dds_or_rna <- DESeqDataSetFromTximport(txi.or_rna.2, colData = samples, design = ~ Treatment)
dds_or_rna <- DESeqDataSetFromTximport(txi.or_rna.2, colData = samples, design = ~ Treatment + mapped_passed_pct + properly.paired_passed_pct)

# Filtering poor genes first ################################################################################
smallestGroupSize <- 9
keep <- rowSums(counts(dds_or_rna) >= 10) >= smallestGroupSize
dds_or_rna.filt <- dds_or_rna[keep,]

#############################################################################################################
### QC
### histogram of read depth 
dds_or_rna.filt <- estimateSizeFactors(dds_or_rna.filt)
count_1=counts(dds_or_rna.filt)
hist(as.matrix(count_1), col="blue", border="white", breaks=20000, xlim=c(0,2000), main="Counts per gene", xlab="Counts (truncated axis)", ylab="Number of genes", las=1, cex.axis=0.7)
### boxplot per sample
boxplot(count_1, main = "Boxplot of Counts Per Sample", las = 2,  col = "lightblue",  ylab = "Count", xlab = "Sample")    


#############################################################################################################
### Variance Partition  
isexpr.2 <- rowSums(fpm(dds_or_rna.filt) > 1) >= 0.5 * ncol(dds_or_rna.filt)
###
quantLog.2 <- log2(fpm(dds_or_rna.filt)[isexpr.2, ] + 1)
form <- ~ (1 | Treatment) + mapped_passed_pct + properly.paired_passed_pct
varPart.2 <- fitExtractVarPartModel(quantLog.2, form, samples)
varPart.2 <- sortCols(varPart.2)
  variance=as.data.frame(varPart.2)
### plot
plotVarPart(varPart.2)


################################################################################
### PCA 
vsd <- vst(dds_or_rna.filt, blind=FALSE)
plotPCA(vsd, intgroup=c("Treatment"))+ geom_label_repel(aes(label = name))

### Extract Loading


ntd <- normTransform(dds_or_rna.filt)
meanSdPlot(assay(ntd))

### HeatMap of reads.
select <- order(rowMeans(counts(dds_or_rna.filt,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_or_rna.filt)[,c("Treatment")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)


################################################################################
### Gene Expression
dds_or_rna <- DESeqDataSetFromTximport(txi.or_rna.2, colData = samples, design = ~ Treatment)
smallestGroupSize <- 9
keep <- rowSums(counts(dds_or_rna) >= 10) >= smallestGroupSize
dds_or_rna.filt <- dds_or_rna[keep,]

### GE
dds_or_rna <- DESeq(dds_or_rna.filt)
### results for each contrast
res_B_W <- results(dds_or_rna, contrast=c("Treatment", "B", "W"))
res_D_W <- results(dds_or_rna, contrast=c("Treatment", "D", "W"))
res_G_W <- results(dds_or_rna, contrast=c("Treatment", "G", "W"))
#comparison of G vs D
res_G_D <- results(dds_or_rna, contrast=c("Treatment", "G", "D"))
#added DvsPB
res_D_B <- results(dds_or_rna, contrast =c("Treatment", "D", "B"))
### Merge Data Frames
res_B_W_df=as.data.frame(res_B_W) %>% rename_with(~ paste0("res_B_W_", .)) %>% mutate(gene_id = rownames(.))
res_D_W_df=as.data.frame(res_D_W) %>% rename_with(~ paste0("res_D_W_", .)) %>% mutate(gene_id = rownames(.))
res_G_W_df=as.data.frame(res_G_W) %>% rename_with(~ paste0("res_G_W_", .)) %>% mutate(gene_id = rownames(.))
# need add G vs D/D vs B
res_G_D_df=as.data.frame(res_G_D) %>% rename_with(~ paste0("res_G_D_", .)) %>% mutate(gene_id = rownames(.))
res_D_B_df=as.data.frame(res_D_B) %>% rename_with(~ paste0("res_D_B_",.)) %>% mutate(gene_id = rownames(.))
### Then into one
res_merged_df <- res_B_W_df %>% full_join(res_D_W_df, by = "gene_id") %>% full_join(res_G_W_df, by = "gene_id") %>% full_join(res_G_D_df, by = "gene_id") %>% full_join(res_D_B_df, by = "gene_id") %>% full_join(gtf_df, by ="gene_id")%>% select(gene_id, gene_name, everything())
# Add additional for G vs D
### export 
write.csv(res_merged_df, 'ue_results_rnaseq.csv')



################################################################################
#### fgsea with hallmarks dataset 
### https://biostatsquid.com/fgsea-tutorial-gsea/
### https://www.gsea-msigdb.org/gsea/index.jsp
### Run set comparison then make plot
### B versus W (pBABE vs WT)
res_B_W_gea=as.data.frame(res_B_W)
res_B_W_gea= res_B_W_gea %>% rownames_to_column("gene_id")
res_B_W_gea <- res_B_W_gea[order(res_B_W_gea$pvalue),]
res_B_W_gea=merge(res_B_W_gea, gtf_df, by ="gene_id")
res_B_W_gea <- res_B_W_gea %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% group_by(gene_name) %>% summarize(stat=mean(stat))
ranks <- deframe(res_B_W_gea)
ranks <- sort(ranks, decreasing = TRUE)
plot(ranks)
### set pathways
pathways.hallmark <- gmtPathways("gmt/h.all.v2023.2.Hs.symbols.gmt")
pathways.hallmark %>% head() %>% lapply(head)
### Run
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, scoreType = 'std',  nproc = 1)
sum(fgseaRes[, padj < 0.05])
head(fgseaRes[order(pval), ])
### output 
fwrite(fgseaRes, file="ur_rnaseq_BvW_fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

### (G184R vs WT)
### G versus W
res_G_W_gea=as.data.frame(res_G_W)
res_G_W_gea= res_G_W_gea %>% rownames_to_column("gene_id")
res_G_W_gea <- res_G_W_gea[order(res_G_W_gea$pvalue),]
res_G_W_gea=merge(res_G_W_gea, gtf_df, by ="gene_id")
res_G_W_gea <- res_G_W_gea %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% group_by(gene_name) %>% summarize(stat=mean(stat))
ranks <- deframe(res_G_W_gea)
ranks <- sort(ranks, decreasing = TRUE)
plot(ranks)
### set pathways
pathways.hallmark <- gmtPathways("gmt/h.all.v2023.2.Hs.symbols.gmt")
pathways.hallmark %>% head() %>% lapply(head)
### Run
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, scoreType = 'std',  nproc = 1)
sum(fgseaRes[, padj < 0.05])
head(fgseaRes[order(pval), ])
### output 
fwrite(fgseaRes, file="ur_rnaseq_GvW_fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


### D versus W (R113W/G184R vs WT)
res_D_W_gea=as.data.frame(res_D_W)
res_D_W_gea= res_D_W_gea %>% rownames_to_column("gene_id")
res_D_W_gea <- res_D_W_gea[order(res_D_W_gea$pvalue),]
res_D_W_gea=merge(res_D_W_gea, gtf_df, by ="gene_id")
res_D_W_gea <- res_D_W_gea %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% group_by(gene_name) %>% summarize(stat=mean(stat))
ranks <- deframe(res_D_W_gea)
ranks <- sort(ranks, decreasing = TRUE)
plot(ranks)
### set pathways
pathways.hallmark <- gmtPathways("gmt/h.all.v2023.2.Hs.symbols.gmt")
pathways.hallmark %>% head() %>% lapply(head)
### Run
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, scoreType = 'std',  nproc = 1)
sum(fgseaRes[, padj < 0.05])
head(fgseaRes[order(pval), ])
### output 
fwrite(fgseaRes, file="ur_rnaseq_DvW_fgseaRes.txt", sep="\t", sep2=c("", " ", ""))
# make these (fgseaRes) into dataframes for each then run the plots to change the titles


#FGSEA plot

# Clean pathway names
fgseaRes$pathway <- gsub("^HALLMARK_", "", fgseaRes$pathway)
fgseaRes$pathway <- gsub("_", " ", fgseaRes$pathway)

# Plot with larger font for pathway labels
ggplot(fgseaRes, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "#56B4E9")) +
  coord_flip() +
  labs(
    title = "Normalized Enrichment Scores (NES) of Pathways G184R vs R113W/G184R",
    x = "Pathway", y = "NES"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6),   # Increase font size of pathway labels
    axis.text.x = element_text(size = 10),   # Optional: increase x-axis font size
    axis.title = element_text(size = 10),    # Optional: increase axis titles
    plot.title = element_text(size = 10)  # Optional: title
  )

######################### Top 10 pathways (change name accordingly for each plot)
# Filter significant pathways
fgsea_sig <- fgseaRes[fgseaRes$padj < 0.05, ]

# Sort by NES
fgsea_sorted <- fgsea_sig[order(fgsea_sig$NES, decreasing = TRUE), ]

# Top 5 upregulated
top_up <- head(fgsea_sorted, 5)

# Top 5 downregulated
top_down <- tail(fgsea_sorted, 5)

# Combine
top_pathways <- rbind(top_up, top_down)

# Clean pathway names
top_pathways$pathway <- gsub("^HALLMARK_", "", top_pathways$pathway)
top_pathways$pathway <- gsub("_", " ", top_pathways$pathway)

# Plot
ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#E69F00", "FALSE" = "#56B4E9")) +
  scale_y_continuous(limits = c(-2.5, 3.5)) +  # Fixed NES axis limits
  coord_flip() +
  labs(
    title = "Top 10 Upregulated and Downregulated Pathways pBABE vs WT (GSEA) ",
    x = "Pathway", y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10)
  )





