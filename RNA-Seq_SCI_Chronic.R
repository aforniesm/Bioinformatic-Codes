load("C:\\Users\\aforn\\OneDrive\\Downloads\\SCI_Chronic.Rdata")
library(dplyr)
# Counts - liver dataset
library(SummarizedExperiment)
counts <- assay(rse_gene,1)
lib.size <- colSums(counts)
df <- rowData(rse_gene)
geneLengths <- df$bp_length
#Check names
genes.ok <- intersect(as.character(rownames(counts)),
                      as.character(rownames(df)))
geneAnnot <- df[genes.ok,]
counts.ok <- counts[genes.ok,]
identical(rownames(geneAnnot), rownames(counts.ok))

# Conditions
condition <- rse_gene$time
ids.early <- grep(paste("Uninjured", condition))
ids.late <- grep(paste("Chronic", condition))

colData(rse_gene)$GROUP <- rep(NA, ncol(rse_gene))
colData(rse_gene)$GROUP[ids.early] <- "early"
colData(rse_gene)$GROUP[ids.late] <- "late"

# Pheno data and remove NA
pheno_data <- colData(rse_gene)
colData <- as.data.frame(pheno_data) %>% filter (!is.na(GROUP))
f <-rownames(colData)
counts.ok <- counts.ok[,f]
dim(counts.ok)
dim(colData)

########################## DESeq 2 ####################################
library(DESeq2)
# Differential expression analysis
countData <- as.matrix(counts.ok)
colData <- as.data.frame(colData)
designFormula <- as.formula(~ GROUP)
#create a DESeq dataset object from the count matrix and the colData
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = designFormula)
# DE ANALYSIS
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
DEresults <- results(dds, contrast = c("GROUP", 'late', 'early'))
DEresults <- DEresults[order(DEresults$pvalue),]
DEresults <- na.omit(DEresults)
# MA plot - Normalization is necessary
library(edgeR)
maPlot(counts.ok[,1], counts.ok[,2], pch=19, cex=.5, ylim=c(-8,8),
       allCol="darkgray", lowess=TRUE,
       xlab=expression( A == log[2] (sqrt(Sample1 %.% Sample3)) ),
       ylab=expression(M == log[2](Sample1/Sample3)))
grid(col="black")
title("Raw data")
DESeq2::plotMA(object = dds, ylim = c(-5, 5), main = "DESeq2
Normalization")
# RLE plots - Normalization is necessary
library(EDASeq)
plotRLE(DESeq2::counts(dds, normalized = FALSE),
        outline=FALSE, ylim=c(-2, 2),
        col = as.numeric(colData$GROUP),
        main = 'Raw Counts',
        xaxt = "n")
mtext("Individuals", side = 1, line = 2)
plotRLE(DESeq2::counts(dds, normalized = TRUE),
        outline=FALSE, ylim=c(-2, 2),
        col = as.numeric(colData$GROUP),
        main = 'Normalized Counts (DESeq2)',
        xaxt = "n")
mtext("Individuals", side = 1, line = 2)
library(pheatmap)
library(ggfortify)
# PCA & HEATMAP
vsd <- vst(dds, blind = FALSE)
NormByDESeq2 <- assay(vsd)
V <- apply(NormByDESeq2, 1, var)
selectedGenes <- names(V[order(V, decreasing = TRUE)][1:50])
annotation_col <- data.frame(GROUP = colData$GROUP)
rownames(annotation_col) <- colnames(NormByDESeq2)
pheatmap(NormByDESeq2[selectedGenes, ], scale = 'column', show_rownames =
           FALSE, show_colnames = FALSE, annotation_col = annotation_col)
M <- t(NormByDESeq2[selectedGenes, ])
M <- log2(M + 1)
pcaResults <- prcomp(M)
autoplot(pcaResults, data = as.data.frame(colData(dds)), colour =
           'xml_gender')
autoplot(pcaResults, data = as.data.frame(colData(dds)), colour =
           'GROUP')
library(GWASTools)
# QQ-Plot
GWASTools::qqPlot(DEresults$pvalue) # Inflation of p-vals - Remove unwanted variability
################### Remove unwanted variability ###############
library(sva)
set.seed(162)
mod1 <- model.matrix( ~ GROUP, data=colData)
mod0 <- model.matrix( ~ 1, data=colData)
countData <- countData[ rowSums(countData) > 10, ]
sv <- svaseq(countData, mod1, mod0)
# Number of estimated surrogate variables
sv$n.sv

# Re run DESeq2
temp <- DataFrame(colData, sv$sv)
colData(dds) <- temp
# update the design
colData(dds)$GROUP <- as.factor(colData(dds)$GROUP)
for (i in 1:10) {
  colname <- colnames(colData(dds))[i]
  colData(dds)[[colname]] <- as.factor(colData(dds)[[colname]])
}
design(dds) <- as.formula( ~ GROUP + V1 + V2)
# re-run the analysis
dds <- dds[rowSums(DESeq2::counts(dds)) > 10]
dds2 <- DESeq(dds, betaPrior = T)
res <- results(dds2, contrast = c("GROUP", 'late', 'early'))
# QQ-PLOT
GWASTools::qqPlot(res$pvalue)
# Volcano Plot
genes <- rownames(res)
genes <- sub("\\..*","", genes)
rownames(res) <- sub("\\..*","", rownames(res))
library(org.Hs.eg.db)
library(AnnotationDbi)
gene_symbol <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL",
                      keytype = "ENSEMBL")

gene_symbol <- unname(gene_symbol)
res$gene_symbol <- gene_symbol
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = res$gene_symbol,
                x = 'log2FoldChange',
                y = 'pvalue')

# P VALUE DISTRIBUTION
library(ggplot2)
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) +
  geom_histogram(bins = 100)
############## ENRICHMENT ANALYSIS ################
res <- res[order(res$pvalue),]
# DEG
mask <- res$padj < 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) >
  log2(2)
deGenes <- rownames(res[mask, ])
# Gene universe
geneUniverse <- rownames(res[!is.na(res$pvalue), ])
# Clean
deGenes_cleaned <- sub("\\..*","", deGenes)
geneUniverse_clean <- sub("\\..*","", geneUniverse)
# List
library(org.Hs.eg.db)
deGenes <- unlist(mget(deGenes_cleaned, envir=org.Hs.egENSEMBL2EG,
                       ifnotfound = NA))
geneUniverse <- unlist(mget(geneUniverse_clean,
                            envir=org.Hs.egENSEMBL2EG,
                            ifnotfound = NA))
# Cluster profiler - enrichKEGG
library(clusterProfiler)
ans.kegg <- enrichKEGG(gene = deGenes, organism = "hsa",
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)

library(enrichplot)
p2 <- dotplot(ans.kegg, showCategory=20) + ggtitle("KEGG enrichment
analysis")
p2

library(DOSE)
disease_enrichment <- enrichDGN(deGenes)
barplot(disease_enrichment, showCategory = 15) + ggtitle("Disease
enrichment analysis")

# Drug enrichment - DGIdb (interactions load into R environment)
gene_names <- mapIds(org.Hs.eg.db, keys = deGenes, column = "SYMBOL",
                     keytype = "ENTREZID", multiVals = "first")

drug2gene <- interactions[, c("drug_concept_id", "gene_name")]
drug2name <- interactions[, c("drug_concept_id", "drug_name")]
gene_names <- na.omit(gene_names)
ans.comp <- enricher(gene_names, TERM2GENE=drug2gene,
                     TERM2NAME=drug2name)
barplot(ans.comp, showCategory = 10) + ggtitle("Drug enrichment
analysis")

# String DB
library(STRINGdb)
deGenes <- as.data.frame(deGenes)
topdeGenes <- as.data.frame(deGenes[1:30,])
colnames(topdeGenes) <- "Gene_ID"
topdeGenes$ENTREZID <-topdeGenes$Gene_ID
string_db <- STRINGdb$new(version = "11.5", species = 9606,
                          score_threshold = 400)
mapped_genes <- string_db$map(topdeGenes, "ENTREZID", removeUnmappedRows
                              = TRUE)
network <- string_db$get_interactions(mapped_genes$STRING_id)
string_db$plot_network(mapped_genes$STRING_id)