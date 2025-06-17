# Packages
# install.packages("SNPassoc")
# BiocManager::install("snpStats")
# BiocManager::install("SNPRelate")
# install.packages("remotes")
# remotes::install_github("bioc/biomaRt")
# BiocManager::install("AnnotationFilter")
# BiocManager::install("ensembldb")
# BiocManager::install("GenomicRanges")
# BiocManager::install("biomaRt")
# install.packages("PredictABEL")
# install.packages("MASS")
# BiocManager::install("VariantAnnotation", ask = FALSE, force = TRUE)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(readxl)
library(SNPassoc)
library(snpStats)
library(SNPRelate)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(remotes)
library(biomaRt)
library(GenomicRanges)
library(PredictABEL)
library(MASS)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(caret)
library(pROC)

# LOAD DATA
pheno <- read.delim("phenotypes.tsv")
plink_data <- read.plink(bed = "snps.bed",
                         fam = "snps.fam",
                         bim = "snps.bim")
excel <- read_excel("data_dictionary.xlsx")
geno <- plink_data$genotypes
annotation <- plink_data$map
family <- plink_data$fam

# CHECK ROWNAMES
rownames(pheno) <- pheno$subject_id
identical(rownames(pheno), rownames(geno))

ids <- intersect(rownames(pheno), rownames(geno))
geno <- geno[ids,]
pheno <- pheno[ids,]
identical(rownames(geno), rownames(pheno))

# CASE-CONTROL - Remove controls with other circulatory disorders
pheno$casecont <- 0
pheno$casecont[pheno$circulatory_disorders=="haemorrhoids"] <- 1
pheno$casecont[pheno$circulatory_disorders%in%c("varicose_veins_of_lower_
extremities","cardiomyopathy")] <- NA
table(pheno$casecont)

# SNP QC
info.snps <- col.summary(geno)
controls <- pheno$casecont ==0 & !is.na(pheno$casecont)
geno.controls <- geno[controls,]
info.controls <- col.summary(geno.controls)
use <- info.snps$Call.rate > 0.95 &
  info.snps$MAF > 0.05 &
  abs(info.controls$z.HWE < 3.3)
mask.snps <- use & !is.na(use)
# Filter
geno.qc.snps <- geno[ , mask.snps]
annotation <- annotation[mask.snps, ]
# SNPs remove for each reason
paste("SNPs removed for Call rate:", sum(info.snps$Call.rate < 0.95))

paste("SNPs removed for MAF:",sum(info.snps$MAF < 0.05, na.rm=TRUE))

paste("SNPs removed for HWE:",sum(abs(info.controls$z.HWE > 3.3),
                                  na.rm=TRUE))

# The total number of SNPs do not pass QC
paste("Total SNPs removed:",sum(!mask.snps))

# INDIVIDUALS QC
info.indv <- row.summary(geno.qc.snps)
geno.X <- geno.qc.snps[,annotation$chromosome=="23" &
                         !is.na(annotation$chromosome)]
info.X <- row.summary(geno.X)
mycol <- ifelse(pheno$sex=="female", "pink", "lightblue")
plot(info.X$Heterozygosity, col=mycol,
     pch=16, xlab="Individuals",
     ylab="Heterozygosity in chromosome X")
legend("topright", c("Females", "Males"), col=c("pink", "lightblue"),
       pch=16)
# Remove sex discrepancies
sex.discrep <- (pheno$sex=="male" & info.X$Heterozygosity > 0.2) |
  (pheno$sex=="female" & info.X$Heterozygosity < 0.2)
# F statistic filtering
mean_het <- mean(info.indv$Heterozygosity)
sd_het <- sd(info.indv$Heterozygosity)
lower_bound <- mean_het - 3 * sd_het
upper_bound <- mean_het + 3 * sd_het
# Transform PLINK data into GDS format
snpgdsBED2GDS("snps.bed",
              "snps.fam",
              "snps.bim",
              out="GDS")
genofile <- snpgdsOpen("GDS")
# Prune SNPs for IBD analysis
set.seed(12345)
snps.qc <- colnames(geno.qc.snps)
snp.prune <- snpgdsLDpruning(genofile, ld.threshold = 0.2,
                             snp.id = snps.qc)
snps.ibd <- unlist(snp.prune, use.names=FALSE)
# IBD coefficients
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    snp.id = snps.ibd,
                    num.thread = 4)
ibd.kin <- snpgdsIBDSelection(ibd)
head(ibd.kin)

# Individuals with higher than expected relatedness are considered with
kinship score > 0.1
ibd.kin.thres <- subset(ibd.kin, kinship > 0.1)
head(ibd.kin.thres)

ids.rel <- SNPassoc::related(ibd.kin.thres)
# Remove Individuals
use <- info.indv$Call.rate > 0.95 &
  info.indv$Heterozygosity > lower_bound &
  info.indv$Heterozygosity < upper_bound &
  !sex.discrep &
  !rownames(info.indv) %in% ids.rel
mask.indiv <- use & !is.na(use)
# Filter
geno.qc <- geno.qc.snps[mask.indiv, ]
pheno.qc <- pheno[mask.indiv, ]
# Checking
identical(rownames(pheno.qc), rownames(geno.qc))

# Count of individual remove for each method
paste("Individuals removed for Call rate:", sum(info.indv$Call.rate <
                                                  0.95))
het1 <- sum(info.indv$Heterozygosity < lower_bound)
het2 <- sum(info.indv$Heterozygosity > upper_bound)
rm_het <- het1+het2
paste("Individuals removed for Heterozigosity:", rm_het)

paste("Individuals removed for Sex Discrepancies:",sum(sex.discrep))

paste("Individuals removed for Relatedness:", length(ids.rel))

paste("Total Individuals removed:",sum(!mask.indiv))

table(pheno$casecont)

# PCA
pca <- snpgdsPCA(genofile, sample.id = rownames(geno.qc),
                 snp.id = snps.ibd,
                 num.thread=4)
with(pca, plot(eigenvect[,1], eigenvect[,2],
               xlab="1st Principal Component",
               ylab="2nd Principal Component",
               main = "Ancestry Plot",
               pch=20, bg="gray90", cex=0.8, col ="gray10"))
# table(pheno.qc$ethnicity)
# Save Principal components
pheno.qc <- data.frame(pheno.qc, pca$eigenvect[, 1:5])
closefn.gds(genofile)
# GWAS
# Remove NA in case-control column from both pheno.qc and geno.qc
pheno.qc <- pheno.qc[!is.na(pheno.qc$casecont),]
ids <- intersect(rownames(pheno.qc), rownames(geno.qc))
geno.qc <- geno.qc[ids, ]
identical(rownames(pheno.qc), rownames(geno.qc))

# GWAS without adjusting
res <- single.snp.tests(pheno.qc$casecont, data=pheno.qc,
                        snp.data=geno.qc)
# QQ-plot
chi2 <- chi.squared(res, df=1)
qq.chisq(chi2)

# Adjusted for ethnicity, X1 and X2
res.adj <- snp.rhs.tests(casecont ~ ethnicity + X1 + X2 + sex +
                           age_recruitment + bmi, data = pheno.qc, snp.data = geno.qc)
# MANHATTAN PLOT
nCHR <- length(unique(annotation$chromosome))
annotation$position_cum <- NA
s <- 0
nbp <- c()
for (i in unique(annotation$chromosome)){
  nbp[i] <- max(annotation[annotation$chromosome == i,]$position)
  annotation[annotation$chromosome == i,"position_cum"] <-
    annotation[annotation$chromosome == i,"position"] + s
  s <- s + nbp[i]
}
axisdf <- annotation %>%
  group_by(chromosome) %>%
  summarize(center=(max(position_cum) + min(position_cum))/2)
axisdf$chromosome <- as.factor(axisdf$chromosome)
annotation$chromosome <- as.factor(annotation$chromosome)
# Create the required data frame
pvals <- data.frame(SNP=annotation$snp.name,
                    CHR=annotation$chromosome,
                    BP=annotation$position,
                    P=p.value(res.adj),
                    BPcum=annotation$position_cum,
                    LOG_P=-log10(p.value(res.adj)))
# missing data is not allowed
pvals <- subset(pvals, !is.na(CHR) & !is.na(P))
# Significance threshold - Adjusted p-value
sig <- 0.05/length(colnames(geno)) # No significant SNPs
new_threshold <- 2e-5
# Colors
mypalette <- c("grey40", "grey70")
manhattanPlot1 <- ggplot(pvals, aes(x = BPcum, y = LOG_P, color = CHR)) +
  geom_point() +
  scale_x_continuous(label = axisdf$chromosome, breaks = axisdf$center) +
  labs(title = "Manhattan Plot", x = "Chromosome", y = "- Log10(Pvalue)") +
  theme(axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)), plot.title =
          element_text(hjust = 0.5, size = 14)) + geom_hline(yintercept = -
                                                               log10(sig), linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme(legend.position = "none") + scale_colour_manual(values =
                                                          rep(mypalette, length(unique(annotation$chromosome)))) +
  geom_label_repel(data=pvals[pvals$P<2e-5,], aes(label=as.factor(SNP),
                                                  alpha=0.7), size=3.5, force=1.3) + geom_hline(yintercept = -
                                                                                                  log10(new_threshold), linetype = "dashed", color = "red", linewidth = 1)
manhattanPlot1 + theme(axis.text.x = element_text(size = 5,hjust = 1,
                                                  angle = 60))
# Top p-values
topPvals <- subset(pvals, P<2e-5)
nrow(topPvals)

topSNPs <- as.character(topPvals$SNP)
# subset top SNPs
geno.topSNPs <- geno.qc[, topSNPs]
geno.topSNPs
write.SnpMatrix(geno.topSNPs, file="topSNPs.txt")
# For locus zoom
write.table(topPvals, file = "topPvals.txt", sep = "\t", quote = F,
            row.names = F)
# SNPassoc preparation and analysis
top_p <- read.delim("topSNPs.txt", sep="")
top_p <- cbind(top_p,pheno.qc)
ii <- grep("^rs", names(top_p))
top.s <- setupSNP(top_p, colSNPs = ii,
                  name.genotypes=c(2,1,0))
# Max-stat
maxstat(top.s$casecont, top.s$rs8851)

# In the vast majority of SNPs the best model is the additive one
# Run WGassociation
wg_asso <- WGassociation(casecont, top.s)
# ODDS RATIO
OR_add <- odds(wg_asso, model = "log-additive")
OR_dom <- odds(wg_asso, model = "dominant")
OR_rec <- odds(wg_asso, model = "recessive")
# ADD OR TO THE TABLE
topPvals <- topPvals[,-5]
topPvals <- topPvals[,-5]
topPvals$OR_additive <- OR_add$OR
topPvals$OR_dominant <- OR_dom$OR
topPvals$OR_recessive <- OR_rec$OR
# Gene annotation
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
snpsList <- topPvals$SNP
snpsList <- as.vector((snpsList))
snpInfo <- getBM(c("refsnp_id", "chr_name", "chrom_start",
                   "minor_allele","minor_allele_freq", "ensembl_gene_name",
                   "associated_gene"), filters = c("snp_filter"), values = snpsList, mart =
                   snpmart)
snpInfo <- snpInfo[-nrow(snpInfo), ]
snpInfo[snpInfo == ""] <- NA
# Merge snpInfo with the OR values
or_columns <- grep("^OR", names(topPvals), value = TRUE)
topPvals_subset <- topPvals[, c("SNP", or_columns)]
snpInfo <- merge(snpInfo, topPvals_subset, by.x = "refsnp_id", by.y =
                   "SNP", all.x = TRUE)
# Complete missing minor allels with annotation data
annotation_sel <- annotation [snpInfo$refsnp_id,]
snpInfo <- merge(snpInfo, annotation_sel[, c("snp.name", "allele.1")],
                 by.x = "refsnp_id", by.y = "snp.name", all.x = TRUE)
snpInfo$minor_allele <- ifelse(is.na(snpInfo$minor_allele),
                               snpInfo$allele.1, snpInfo$minor_allele)
snpInfo$allele.1 <- NULL
# Annotate genes in other way
snpsInfo.gr <- makeGRangesFromDataFrame(snpInfo,
                                        start.field="chrom_start",keep.extra.columns = TRUE,seqnames.field =
                                          "chr_name",end.field = "chrom_start")
seqnames(snpsInfo.gr)

# USCS manner
seqlevelsStyle(snpsInfo.gr) <- "UCSC"
seqlevels(snpsInfo.gr, pruning.mode = "coarse") <- paste0("chr",
                                                          c(1:22,"X","Y"))
coding <- locateVariants(snpsInfo.gr,
                         TxDb.Hsapiens.UCSC.hg38.knownGene,
                         CodingVariants())

allvar <- locateVariants(snpsInfo.gr,
                         TxDb.Hsapiens.UCSC.hg38.knownGene,
                         AllVariants())
keys <- allvar$GENEID
genes <- AnnotationDbi::select(org.Hs.eg.db,
                               columns=c("GENENAME","ENSEMBL", "SYMBOL"),
                               key = keys, keytype = "ENTREZID")
# Complete summary table with more gene symbols that biomaRt doesn't
provide
for (i in seq_len(nrow(snpInfo))) {
  if (is.na(snpInfo$associated_gene[i])) {
    match_row <- which(genes$ENSEMBL == snpInfo$ensembl_gene_name[i])
    if (length(match_row) > 0) {
      snpInfo$associated_gene[i] <- genes$SYMBOL[match_row[1]]
    }
  }
}
snpInfo <- snpInfo %>%
  mutate(
    ensembl_gene_name = ifelse(is.na(ensembl_gene_name), "intergenic
region", ensembl_gene_name),
    associated_gene = ifelse(is.na(associated_gene), "intergenic region",
                             associated_gene)
  )
snpInfo <- snpInfo %>%
  mutate(
    minor_allele_freq = ifelse(is.na(minor_allele_freq), "not provided by
biomaRt", minor_allele_freq)
  )
# Discrepancies between genomic position: Annotation data vs BiomaRt
# We conserve the Annotation Position
snpInfo <- merge(snpInfo, annotation_sel, by.x = "refsnp_id", by.y =
                   "snp.name", all.x = TRUE)
snpInfo$chrom_start <- snpInfo$position
# Now snpInfo is ready for constructing the summary table
# Genetic score
ans <- WGassociation(casecont, top.s, model="log-add")
sel <- labels(top.s)[additive(ans)<0.1]
top.sel <- top.s[,sel]
top.sel <- data.frame(lapply(top.sel, additive))
dd.end <- data.frame(casecontrol=top.s$casecont, top.sel)
# Complete cases and AIC
dd.end.complete <- dd.end[complete.cases(dd.end),]
mod <- stepAIC(glm(casecontrol ~ .,
                   dd.end.complete,family="binomial"),method="forward", trace=0)
# Create the genetic score
snps.score <- names(coef(mod))[-1]
snps.score

pos <- which(names(dd.end.complete)%in%snps.score)
score <- riskScore(mod, data=dd.end.complete,
                   cGenPreds=pos,
                   Type="unweighted")
table(score)

# Association with our disease
mod.lin <- glm(casecontrol~score, dd.end.complete,family="binomial")
exp(coef(mod.lin)[2])

predrisk <- predRisk(mod.lin, dd.end.complete)
plotROC(data=dd.end.complete, cOutcome=1,
        predrisk = predrisk)

# 10-fold cross-validation procedure
folds <- createFolds(dd.end.complete$casecontrol, k = 10, list = TRUE,
                     returnTrain = FALSE)
auc_scores <- c()
dd.end.complete$GeneticScore <- riskScore(mod, data = dd.end.complete,
                                          cGenPreds = pos,
                                          Type = "unweighted")
for (i in seq_along(folds)) {
  test_indices <- folds[[i]]
  train_data <- dd.end.complete[-test_indices, ]
  test_data <- dd.end.complete[test_indices, ]
  
  model <- glm(casecontrol ~ GeneticScore, data = train_data, family =
                 "binomial")
  
  predictions <- predict(model, newdata = test_data, type = "response")
  
  auc <- roc(test_data$casecontrol, predictions)$auc
  auc_scores <- c(auc_scores, auc)
}
mean_auc <- mean(auc_scores)
sd_auc <- sd(auc_scores)
paste("Mean AUC:", round(mean_auc, 3))

paste("Standard Deviation of AUC:", round(sd_auc, 3))

