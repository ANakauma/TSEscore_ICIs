# Author:  J. Alberto Nakauma Gonzalez
# Date:   28-06-2022
# Function: This script generates figures from RNA-Seq data for DR-176 mUC/CPCT-02 cohort
# e-mail: j.nakaumagonzalez@erasmusmc.nl

# Clean everything before starting a new project and set working directory (default R script path)---------------------------
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(plyr)
library(dplyr)
library(DESeq2)
library(ConsensusClusterPlus)
library(consensusMIBC)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(eulerr)
library(dndscv)
library(rtracklayer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(reshape2)
library(tcR)

# Import sample information -----------------------------------------------

# path to data and output dir
pathHMF <- "/HMF/DR176/postHMF/"
odir <- "/output/figures/"

# Import sample information
load(paste0(pathHMF,"RData/DR176.MetaData.RData"))
load(paste0(pathHMF, "RData/DESeq2.DR176.WT.RespNonResp.RData"))
load(paste0(pathHMF,"RData/data.Cohort.RData"))
load(paste0(pathHMF,"RData/results.Cohort.RData"))
load(paste0(pathHMF, "RData/immuneCFResults.RData"))

# keep only samples with RNA-seq data
DR176.MetaData <- DR176.MetaData %>% dplyr::filter(sampleId %in% colnames(immuneCFResults))


# define colors for different classes
colorConsensusClass = c("Basal/Squamous" = "#65AF45", "Stroma-rich" = "#A9A9A9", "Luminal unstable" = "#EBE377",
                        "Luminal papillary" = "#E85252", "Luminal non-specified" = "#D88735", "Neuroendocrine-like" = "#6877F7")
colorRNAsubtypes <- c("Non-specified" = "#ca3b3b", "Luminal-b" = "#E69F00",
                      "Luminal-a" = "#56B4E9", "Stroma-rich" = "#A9A9A9", "Basal/Squamous" = "#65AF45")
color_immunoClusters <- c("Cluster 1" = "#7dd49b", "Cluster 3" = "#ff4d5b", "Cluster 2" = "#ab6a20", "Cluster 4" = "#96A5AF")


# Helper functions --------------------------------------------------------

annotationTheme <- function(){
  theme(legend.position = 'bottom', axis.ticks = element_blank(), axis.title.y = element_text(size = 8), axis.text.x = element_blank(), text=element_text(size=8, family='Helvetica'),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(NULL),
        plot.margin = ggplot2::margin(t = 0, 0,0,0, 'cm'),
        legend.margin = ggplot2::margin(t=-1, r=0, b=.5, l=0, unit='cm')
  )
}

# Sorting of oncoplot. ----------------------------------------------------

memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}


#=====================================  RNA subtypes ====================================
#===============================================================================================

# load RNA genes from Nakauma-Gonzalez, Rijnders et. al, 2021
mUC_subtypesGenes <- readxl::read_xlsx("/databases/Genes_mUC_RNAsubtypes.xlsx")

# normalize counts with Variance Stabilizing Transformation
DESeq2.DR176.WT.vst <- DESeq2::vst(DESeq2.DR176.WT[rowSums(DESeq2::counts(DESeq2.DR176.WT, normalized = T)) > 0,], blind = T)

#get normalized count Matrix
normalizedCountMatrix_all <- SummarizedExperiment::assay(DESeq2.DR176.WT.vst)

# filter and keep only relevant genes for mUC RNA subtypes
normMatrixRNA_clean <- normalizedCountMatrix_all[rownames(normalizedCountMatrix_all) %in% mUC_subtypesGenes$geneName, ]

# Normalize data by centering around the median each gene value across samples (data already log2 transformed)
normMatrixRNA_clean  <- normMatrixRNA_clean - rowMedians(normMatrixRNA_clean)

# get data frame from expression data
normMatrixRNA_clean <- as.data.frame(normMatrixRNA_clean)
normMatrixRNA_clean$geneName <- rownames(normMatrixRNA_clean)
rownames(normMatrixRNA_clean) <- NULL
normMatrixRNA_clean <- normMatrixRNA_clean %>% reshape2::melt(id = "geneName", value.name = "Expression")


# Add phenotypic label of genes
normMatrixRNA_clean <- normMatrixRNA_clean %>%
  dplyr::left_join(dplyr::select(mUC_subtypesGenes, geneName, geneRNAsubtype_mUCDR31 = RNAcluster), by = "geneName")

# Calculate score for gene signatures for Up and Down regulated genes
normMatrixRNA_clean <- normMatrixRNA_clean %>% dplyr::group_by(geneRNAsubtype_mUCDR31, variable) %>%
  dplyr::mutate(meanExp = mean(Expression)) %>% dplyr::ungroup()

# order genenames
normMatrixRNA_clean <- normMatrixRNA_clean %>%
  dplyr::arrange(factor(geneRNAsubtype_mUCDR31, levels = c("Basal/Squamous", "Stroma-rich",
                                          "Luminal-a", "Luminal-b", "Non-specified")))

# get final scores for each phenotype
normMatrixRNA_clean <- normMatrixRNA_clean %>%
  dplyr::select(variable, geneRNAsubtype_mUCDR31, meanExp) %>% dplyr::distinct() %>%
  dplyr::select(sample = variable, geneRNAsubtype_mUCDR31, subtypeScore = meanExp) %>% dplyr::distinct()


# Based on Signature score for each mUC subtype, we can define transcriptomic subtype per sample 
normMatrixRNA_clean_LumA <- normMatrixRNA_clean %>% 
  dplyr::filter(geneRNAsubtype_mUCDR31 == "Luminal-a") %>% dplyr::select(sample, LumA = subtypeScore)
normMatrixRNA_clean_LumB <- normMatrixRNA_clean %>%
  dplyr::filter(geneRNAsubtype_mUCDR31 == "Luminal-b") %>% dplyr::select(sample, LumB = subtypeScore)
normMatrixRNA_clean_BaSq <-  normMatrixRNA_clean %>%
  dplyr::filter(geneRNAsubtype_mUCDR31 == "Basal/Squamous") %>% dplyr::select(sample, BaSq = subtypeScore)
normMatrixRNA_clean_Stroma <-  normMatrixRNA_clean %>%
  dplyr::filter(geneRNAsubtype_mUCDR31 == "Stroma-rich") %>% dplyr::select(sample, Stroma = subtypeScore)
normMatrixRNA_clean_NonSpe <- normMatrixRNA_clean %>%
  dplyr::filter(geneRNAsubtype_mUCDR31 == "Non-specified") %>% dplyr::select(sample, NonSpe = subtypeScore)

# get all scores per sample
normMatrixRNA_clean_mUCSubtypes <- normMatrixRNA_clean_LumA %>% dplyr::full_join(normMatrixRNA_clean_LumB, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_clean_BaSq, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_clean_Stroma, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_clean_NonSpe, by = "sample")

# Assign each sample it mUC subtype
results_RNAmUCsubtype <- normMatrixRNA_clean %>% dplyr::arrange(-subtypeScore) %>%
  dplyr::distinct(sample, geneRNAsubtype_mUCDR31) %>% dplyr::distinct(sample, .keep_all = TRUE) %>%
  dplyr::full_join(normMatrixRNA_clean_mUCSubtypes, by = "sample")

# order samples
colnames(results_RNAmUCsubtype)[1:2] <- c("sampleId", "RNA_mUCsubtype")
results_RNAmUCsubtype$RNA_mUCsubtype <- factor(results_RNAmUCsubtype$RNA_mUCsubtype,
                                             levels = c("Luminal-a", "Luminal-b", "Basal/Squamous",
                                                        "Stroma-rich", "Non-specified"))


# plot cluster number labels
plot_RNAsubtype_DR176 <- ggplot(results_RNAmUCsubtype, aes(sampleId, y = "Transcriptomic subtype", fill = RNA_mUCsubtype)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T, alpha=0.75) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = colorRNAsubtypes,
                    guide = guide_legend(title = 'Transcriptomic subtype', title.position = 'top',
                                         title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_RNAsubtype_DR176 <- cowplot::get_legend(
  plot_RNAsubtype_DR176 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_RNAsubtype_DR176 <- plot_RNAsubtype_DR176 + theme(legend.position = 'none')




#=========================  DGE between responders and non-responders ==========================
#===============================================================================================
# Import annotations
genesDB <- rtracklayer::import.gff('/annotation/hg19/GENCODE/gencode.v35lift37.annotation.gtf')
genesDB <- tibble::as_tibble(mcols(genesDB))
geneInfo <- genesDB %>% dplyr::distinct(gene_id, gene_name)


# Find diff expressed genes
DESeq2.DR176.WTDiff <- DESeq2::results(DESeq2.DR176.WT, contrast = c("Outcome6M", "non-responder", "responder"), pAdjustMethod = 'BH')

# get differential expressed genes
DESeq2.DR176.WTDiff_result <- tibble::as_tibble(DESeq2::results(DESeq2.DR176.WT, contrast = c("Outcome6M", "non-responder", "responder"), pAdjustMethod = 'BH', tidy = T))

# discard low expressed genes
DESeq2.DR176.WTDiff_result <- DESeq2.DR176.WTDiff_result %>% dplyr::filter(baseMean > 10)

# Find max and min values to use in Volcano plot
max(DESeq2.DR176.WTDiff_result$log2FoldChange)
min(DESeq2.DR176.WTDiff_result$log2FoldChange)
min(DESeq2.DR176.WTDiff_result$padj)

# define min padj and log2FC for colors
minPadj <- 0.01
minLog2FC <- 1
DESeq2.DR176.WTDiff_result <- DESeq2.DR176.WTDiff_result %>%
  dplyr::mutate(expUpOrDown = ifelse(log2FoldChange > minLog2FC & padj < minPadj, "Up",
                                     ifelse(log2FoldChange < -minLog2FC & padj < minPadj, "Down", "Stable"))) %>%
  dplyr::mutate(labelGene = ifelse(expUpOrDown == "Stable", NA, row))
DESeq2.DR176.WTDiff_result <- DESeq2.DR176.WTDiff_result %>% dplyr::filter(!is.na(padj))

# volcanoplot
volcanoPlot_RespNonResp <- ggplot(data = DESeq2.DR176.WTDiff_result, aes(x = log2FoldChange, y = -log10(padj), colour=expUpOrDown, label = labelGene)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("Down" = "blue", "Stable" = "grey", "Up" = "red")) +
  geom_label_repel(aes(label = labelGene),
                   box.padding   = 0.05, 
                   label.padding = 0.05,
                   point.padding = 0.05,
                   segment.color = 'grey50',
                   max.overlaps = Inf,
                   fill = NA,
                   label.size = NA,
                   size = 1.5) +
  xlim(c(-5.5, 6.5)) +
  ylim(c(-0, 6.5)) +
  geom_vline(xintercept=c(-minLog2FC, minLog2FC),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(minPadj),lty=4,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adj p-value)",
       title= "Responder/Non-responder")  +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8), 
    legend.position="none", 
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey80', linetype = 'dotted'),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 6),
    text = element_text(size = 8, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    strip.background = element_rect(colour = 'grey20', fill = 'white'),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    legend.title = element_blank())


# Get a list of log2FC
# get differential expressed genes
MIN_ABS_LOG2FC <- 1
MIN_PVALUE <- 0.1
MINPATHWAYS <- 10
DESeq2.DR176.WTDiff_result <- tibble::as_tibble(DESeq2::results(DESeq2.DR176.WT, contrast = c("Outcome6M", "non-responder", "responder"), pAdjustMethod = 'BH', tidy = T))
DiffExpGenes <- DESeq2.DR176.WTDiff_result %>% dplyr::filter(padj <= MIN_PVALUE, abs(log2FoldChange) > MIN_ABS_LOG2FC, baseMean > 10)

# get name of genes 
DiffExpGenes <- dplyr::left_join(DiffExpGenes, geneInfo, by = c("row" = "gene_name")) %>%
  dplyr::distinct(row, .keep_all = TRUE) %>% dplyr::mutate(gene_id = gsub("\\..*", "", gene_id))

# get list of genes with log2FC
geneDiffExpList <- DiffExpGenes$log2FoldChange
names(geneDiffExpList) <- DiffExpGenes$gene_id

# Change ENSEMBL to ENTREZID.
geneIDs <- clusterProfiler::bitr(names(geneDiffExpList), fromType = "ENSEMBL", toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = FALSE)
names(geneDiffExpList) <- geneIDs[!duplicated(geneIDs$ENSEMBL), "ENTREZID"]
#geneDiffExpList <- sort(geneDiffExpList[!is.na(names(geneDiffExpList))], decreasing = TRUE)
geneDiffExpList <- geneDiffExpList[!is.na(names(geneDiffExpList))]

# Perform pathway enrichment for up-regulated genes
resultPathwayEnrich_Up_complete <- ReactomePA::enrichPathway(gene = names(geneDiffExpList)[geneDiffExpList > MIN_ABS_LOG2FC], maxGSSize = 500,
                                                             pvalueCutoff=0.05, readable=T)


# get results of enrichment and add extra information
resultPathwayEnrich_Up_complete <- resultPathwayEnrich_Up_complete@result %>% 
  dplyr::mutate(pathwayIs = "Up-regulated", log10_pAdj = -log10(p.adjust))

# List of top enriched pathways
resultPathwayEnrich_Up <- resultPathwayEnrich_Up_complete[resultPathwayEnrich_Up_complete$p.adjust < 0.05, ]
resultPathwayEnrich_Up <- resultPathwayEnrich_Up[0:min(nrow(resultPathwayEnrich_Up), MINPATHWAYS), ]


# Perform pathway enrichment for down-regulated genes
resultPathwayEnrich_Down_complete <- ReactomePA::enrichPathway(gene = names(geneDiffExpList)[geneDiffExpList < -MIN_ABS_LOG2FC], maxGSSize = 500,
                                                               pvalueCutoff=0.05, readable=T)

# get results of erichment and add extra information
resultPathwayEnrich_Down_complete <- resultPathwayEnrich_Down_complete@result %>%
  dplyr::mutate(pathwayIs = "Down-regulated", log10_pAdj = -log10(p.adjust))

# List of top enriched pathways
resultPathwayEnrich_Down <- resultPathwayEnrich_Down_complete[resultPathwayEnrich_Down_complete$p.adjust < 0.05, ]
resultPathwayEnrich_Down <- resultPathwayEnrich_Down[0:min(nrow(resultPathwayEnrich_Down), MINPATHWAYS), ] 

resultPathwayEnrich <- rbind(resultPathwayEnrich_Up, resultPathwayEnrich_Down)

# Order data
resultPathwayEnrich$Description <- factor(resultPathwayEnrich$Description, levels = rev(resultPathwayEnrich$Description) )
resultPathwayEnrich$pathwayIs <- factor(resultPathwayEnrich$pathwayIs, levels = c("Up-regulated", "Down-regulated") )

# Plot results
plotPathwayEnrich_RespNonResp <- ggplot(resultPathwayEnrich, aes(x = Description, y = log10_pAdj, fill = pathwayIs)) +
  geom_bar(stat = 'identity', col = 'grey20') +
  scale_y_continuous(expand = expansion(0,0)) +
  labs(x = NULL, y = "-log10(pAdj)") +
  geom_hline(yintercept=-log10(0.05),lty=4,col="black",lwd=0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Up-regulated" = "brown", "Down-regulated" = "dodgerblue3"),
                    guide = guide_legend(title = NULL, title.position = 'top',
                                         title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5))  +
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# export Figure 3
pdf(paste0(odir,"RNAexp_Volcano_NonRespVsResp.pdf"), width = 8, height = 2.5)
cowplot::plot_grid(volcanoPlot_RespNonResp,
                   plotPathwayEnrich_RespNonResp,
  nrow = 1, rel_widths = c(0.35, 0.65), rel_heights = c(1, 0.8))#, align = 'h', axis = 'tblr')

dev.off()






#================================  gene signature scores =======================================
#===============================================================================================


# Prepare normalized counts ---------------------------------------

# Normalize with DESeq2 with Variance Stabilizing Transformation
DESeq2.DR176.WT.vst <- DESeq2::vst(DESeq2.DR176.WT[rowSums(DESeq2::counts(DESeq2.DR176.WT, normalized = T)) > 10,], blind = T)

# Get sample annotations and clean names ---------------------------------------
colAnno <-  as.data.frame(SummarizedExperiment::colData(DESeq2.DR176.WT.vst)[, c("RESP_subject", "primaryTumorSubLocation", "gender", "tumorPurity",
                                                           "Biopsy site annotation", "BOR", "Outcome6M")])

#Change name from gender to sex
colnames(colAnno)[3] <- "sex"

# Add New names to Biopsy site to reduce variables
dplyr::count(colAnno, Biopsy.site.annotation, sort = T)
colAnno$sampleId <- rownames(colAnno)
colAnno <- colAnno %>% dplyr::mutate(Biopsy.site.annotation = ifelse(Biopsy.site.annotation %in%  c("Liver", "Lymph node", "Soft tissue"), Biopsy.site.annotation, "Other"))
rownames(colAnno) <- colAnno$sampleId

# Response to treatment
colAnno <- colAnno %>% dplyr::left_join(dplyr::select(results.Cohort$mutationalBurden, sample, Genome.TMB, tmbStatus), by = c("sampleId" = "sample")) %>%
  dplyr::mutate(Outcome6M = ifelse(Outcome6M == "responder", "Responder", "Non-responder")) %>%
  dplyr::mutate(Outcome6M = factor(Outcome6M, levels = c("Responder", "Non-responder"))) %>%
  dplyr::mutate(tmbStatus = ifelse(tmbStatus == 'High TMB (≥10)', 'High TMB (>10)',
                                   ifelse(tmbStatus == 'Medium TMB (≥5-10)', 'Medium TMB (5-10)', 'Low TMB (0-5)'))) %>%
  dplyr::mutate(tmbStatus = factor(tmbStatus, levels = c('High TMB (>10)', 'Medium TMB (5-10)', "Low TMB (0-5)")))

#get normalized count Matrix
normalizedCountMatrix <- SummarizedExperiment::assay(DESeq2.DR176.WT.vst)


# Get list of all gene signatures
gene_signatures <- readxl::read_xlsx("/databases/GeneSignatures_complete.xlsx",
                                     sheet = 1)
gene_signatures <- gene_signatures %>% dplyr::filter(Gene %in% rownames(normalizedCountMatrix))

# Get expression of markers
matrixRNA_sigExpression <- normalizedCountMatrix[rownames(normalizedCountMatrix) %in% unique(gene_signatures$Gene), ]

# Center gene expression around the median across samples, data is already log2 transformed
matrixRNA_sigExpression  <- matrixRNA_sigExpression - rowMedians(matrixRNA_sigExpression)

# Add column with gene names
matrixRNA_sigExpression <- as.data.frame(matrixRNA_sigExpression) %>% dplyr::mutate(geneName = rownames(matrixRNA_sigExpression))

# Melt results
signatureScoreMarkers <- reshape2::melt(matrixRNA_sigExpression, id.vars = "geneName")

# Add marker type information
signatureScoreMarkers <- gene_signatures %>% dplyr::right_join(signatureScoreMarkers, by = c("Gene" = "geneName")) %>%
  dplyr::full_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M), by = c("variable" = "sampleId"))

# get signature score
signatureScoreMarkers <- signatureScoreMarkers %>% dplyr::mutate(signatureID = Subcategory) %>%
  dplyr::group_by(signatureID, variable) %>%
  dplyr::mutate(sigScore = mean(value)) %>% dplyr::ungroup() %>%
  dplyr::distinct(signatureID, variable, .keep_all = TRUE)

# Get signature scores to use for ROC and AUC (later)
signatureScoreMarkers_all <- signatureScoreMarkers

# Calculate variability of each signature across the entire cohort
signatureScoreMarkers_sd <- signatureScoreMarkers %>% dplyr::group_by(signatureID) %>%
  dplyr::mutate(sd_signature = sd(sigScore)) %>% dplyr::ungroup() %>%
  dplyr::arrange(-sd_signature) %>% dplyr::distinct(signatureID, sd_signature, MainCategory)

# get order data by rows (genes)
signatureScoreMarkers_tmp <- reshape2::dcast(dplyr::select(signatureScoreMarkers, variable, signatureID, sigScore), signatureID ~ variable)
rownames(signatureScoreMarkers_tmp) <- signatureScoreMarkers_tmp$signatureID
signatureScoreMarkers_tmp$signatureID <- NULL
distRows <- dist(signatureScoreMarkers_tmp, method = "euclidean") # distance matrix
hclustRows <- hclust(distRows, method="complete")
plot(hclustRows) # display dendogram



# export hclust Figure S4
pdf(paste0(odir,"hclust_signatures.pdf"), width = 11, height = 8)
plot(hclustRows)
dev.off()



# get the order of signature clusters from dendogram
orderSignatureNames <- hclustRows$labels[c(hclustRows$order)]
orderSignatureNames <- c(orderSignatureNames[1:5], orderSignatureNames[13:36], orderSignatureNames[6:12])

# Select signatures that contribute most to T-cell and stromal compartments (based on plot(hclustRows) and later on higuest AUC)
selectedSignaturesGlobalSig <- c("CAF", "EMT/stroma core genes", "Fibroblasts", "Stromal signature", "TBRS", 
                                 "IFN gamma", "tGE8", "T cell signature", "Immune gene signature", "T cell inflamed GEP", "Chemoattractants", "Cytotoxic CD8 T cell")


# finishing getting the order order of signatures considering the main category and Response
signatureScoreMarkers <- signatureScoreMarkers %>%
  dplyr::arrange(factor(signatureID, levels = orderSignatureNames)) %>%
  dplyr::arrange(factor(MainCategory, levels = c("Stromal resident cells", #"Oncogenic signaling", "Antigen presentation and TLS",
                                                 "Immune cells (non T cell)", "T cells")))
orderSignatureNames <- unique(signatureScoreMarkers$signatureID)



# Calculate global signature scores for all immune cells and stroma
signatureScoreMarkers_tmp <- signatureScoreMarkers %>%
  dplyr::mutate(Immune_Stroma = ifelse(MainCategory != "Stromal resident cells",
                                                    "Immune cells", MainCategory)) %>%
  dplyr::group_by(Immune_Stroma, variable) %>%
  dplyr::mutate(globalSigScore_Immune_Stroma = mean(sigScore)) %>% dplyr::ungroup()

# Calculate global signature scores (for all Immune cells (non-T cells) and for selected Stroma and T-cells)
signatureScoreMarkers_tmp <- signatureScoreMarkers_tmp %>%
  dplyr::filter(signatureID %in% selectedSignaturesGlobalSig | MainCategory == "Immune cells (non T cell)") %>%
  dplyr::group_by(MainCategory, variable) %>%
  dplyr::mutate(globalSigScore_mainCategories = mean(sigScore)) %>% dplyr::ungroup()
  
# calculate (get) global signature score
signatureScoreMarkers_globalSigScore_stroma <-  signatureScoreMarkers_tmp %>% dplyr::filter(MainCategory == 'Stromal resident cells') %>%
   dplyr::select(variable, globalSigScore_stroma = globalSigScore_mainCategories) %>%
   dplyr::distinct(variable, globalSigScore_stroma)
signatureScoreMarkers_globalSigScore_Tcells <-  signatureScoreMarkers_tmp %>% dplyr::filter(MainCategory == 'T cells') %>%
   dplyr::select(variable, globalSigScore_Tcells = globalSigScore_mainCategories) %>%
   dplyr::distinct(variable, globalSigScore_Tcells)
signatureScoreMarkers_globalSigScore_immuneCells <-  signatureScoreMarkers_tmp %>% dplyr::filter(MainCategory == 'Immune cells (non T cell)') %>%
  dplyr::select(variable, globalSigScore_immuneCells = globalSigScore_mainCategories) %>%
  dplyr::distinct(variable, globalSigScore_immuneCells)
signatureScoreMarkers_globalSigScore_ALLimmuneCells <-  signatureScoreMarkers_tmp %>% dplyr::filter(Immune_Stroma == 'Immune cells') %>%
  dplyr::select(variable, globalSigScore_ALLimmuneCells = globalSigScore_Immune_Stroma) %>%
  dplyr::distinct(variable, globalSigScore_ALLimmuneCells)

signatureScoreMarkers_globalSigScore <- signatureScoreMarkers_globalSigScore_stroma %>%
  dplyr::full_join(signatureScoreMarkers_globalSigScore_Tcells, by = "variable") %>%
  dplyr::full_join(signatureScoreMarkers_globalSigScore_immuneCells, by = "variable") %>%
  dplyr::full_join(signatureScoreMarkers_globalSigScore_ALLimmuneCells, by = 'variable') %>%
  dplyr::mutate(signatureImbalance_TcellsStroma = globalSigScore_Tcells - globalSigScore_stroma,
                signatureImbalance_ImmuneCellsStroma = globalSigScore_immuneCells - globalSigScore_stroma,
                signatureImbalance_ALLimmuneCellsStroma = globalSigScore_ALLimmuneCells - globalSigScore_stroma) %>%
  dplyr::mutate(signatureImbalance_tmp_TcellsStroma = ifelse(signatureImbalance_TcellsStroma > 1, 1, 
                                                ifelse(signatureImbalance_TcellsStroma < -1, -1, signatureImbalance_TcellsStroma)),
                signatureImbalance_tmp_ImmuneCellsStroma = ifelse(signatureImbalance_ImmuneCellsStroma > 1, 1, 
                                                             ifelse(signatureImbalance_ImmuneCellsStroma < -1, -1, signatureImbalance_ImmuneCellsStroma)),
                signatureImbalance_tmp_ALLimmuneCellsStroma = ifelse(signatureImbalance_ALLimmuneCellsStroma > 1, 1,
                                                                     ifelse(signatureImbalance_ALLimmuneCellsStroma < -1, -1, signatureImbalance_ALLimmuneCellsStroma))) %>%
  dplyr::full_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M), by = c("variable" = "sampleId"))




# #### Consensus clustering for all signatures
# # ++++++++++++++++++++++++++++++ Calculate Clusters
# Cluster with consensus cluster
clusterDir <- paste0(odir, "RNA_immuneCluster")
ConsensusMatrix_immuneSig <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(signatureScoreMarkers_tmp), pItem=0.9, pFeature=1,
                                                              maxK = 10,reps = 1000, title = clusterDir, plot = "pdf",
                                                              distance = "pearson", finalLinkage="ward.D2") #, finalLinkage="ward.D2"
icl = ConsensusClusterPlus::calcICL(ConsensusMatrix_immuneSig, title=clusterDir, plot="pdf")

# save results
save(ConsensusMatrix_immuneSig, file = paste0(pathHMF, "RData/ConsensusMatrixRNA_immuneSig.RData"))

# load data
load(paste0(pathHMF, "RData/ConsensusMatrixRNA_immuneSig.RData"))

#Choose num of clusters and Get Consensus Matrix
numClusters <- 4
results_RNAConsensusClass <- as.data.frame(ConsensusMatrix_immuneSig[[numClusters]]$consensusClass)
colnames(results_RNAConsensusClass)[1] <- "Cluster"
results_RNAConsensusClass$sample <- rownames(results_RNAConsensusClass)


# Name clusters based on expression of RNAsubtype genes (STS) (see cluster results = consensus.pdf))
results_RNAConsensusClass <- results_RNAConsensusClass %>%
  dplyr::mutate(immunoCluster = ifelse(Cluster == 1, "Cluster 3",
                                        ifelse(Cluster == 2, "Cluster 1", "Cluster 2")))
results_RNAConsensusClass$immunoCluster <- factor(results_RNAConsensusClass$immunoCluster,
                                               levels = c("Cluster 1", "Cluster 2", "Cluster 3"))

results_RNAConsensusClass$Cluster <- NULL

# Save MSig DBS Clusters
resultsRNA_immuneClusters <- results_RNAConsensusClass
# ++++++++++++++++++++++++++++++ Finished -- Calculate Clusters


# order of data based on the consensus clustering
numClusters <- 4
orderSampleNames <- names(ConsensusMatrix_immuneSig[[numClusters]][["consensusClass"]][ConsensusMatrix_immuneSig[[numClusters]][["consensusTree"]][["order"]]])
# get order of samples (clusters should be in order, cluster 3 goes to the end)
resultsRNA_immuneClusters <- resultsRNA_immuneClusters %>%
  dplyr::arrange(factor(sample, levels = orderSampleNames)) %>%
  dplyr::arrange(immunoCluster)
orderSampleNames <- as.character(resultsRNA_immuneClusters$sample)


# order samples and signatures
signatureScoreMarkers$variable <- factor(signatureScoreMarkers$variable, levels = orderSampleNames)
signatureScoreMarkers <- signatureScoreMarkers %>%
  dplyr::mutate(signatureID = factor(signatureID, levels = orderSignatureNames))



# Set max.min values to -2 or 2
minExp <- 2
signatureScoreMarkers  <- signatureScoreMarkers %>%
  dplyr::mutate(sigScore = ifelse(sigScore > minExp, minExp, ifelse(sigScore < -minExp, -minExp, sigScore)))

# plot markers score
plot_sigScoreRNAcluster <- ggplot(signatureScoreMarkers, aes(x = variable, y = signatureID, fill = sigScore)) +
  geom_tile(size = 0.25, na.rm = F, alpha = 1, height = 1, width = 1) +
  scale_fill_gradient2(low = "#3288BD", mid = 'white', high = '#9E0142',
                       breaks=c(-2, -1, 0, 1, 2),labels=c("<-2", -1, 0, 1, ">2"),
                       guide = guide_colorbar(title = 'Signature score',
                                              title.position = 'top', title.hjust = 0.5, barwidth = 5, barheight = 0.5)) +
  annotationTheme() +
  labs(x = NULL, y = 'Gene signature') +
  theme(
    legend.position = 'bottom',
    text=element_text(size=8, family='Helvetica'),
    axis.text=element_text(size=10),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_sigScoreRNAcluster <- cowplot::get_legend(
  plot_sigScoreRNAcluster + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_sigScoreRNAcluster <- plot_sigScoreRNAcluster + theme(legend.position = 'none')


signatureScoreMarkers$MainCategory <- factor(signatureScoreMarkers$MainCategory,
                                               levels = c('T cells', 'Immune cells (non T cell)', 'Stromal resident cells'))

# plot Signature category
plot.sigCategory <- ggplot(dplyr::distinct(signatureScoreMarkers, signatureID, MainCategory),
                           aes(x = "Signature category", y = signatureID, fill = MainCategory)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = "Gene signature", x = NULL) +
  annotationTheme() +
  scale_fill_manual('Signature category',
                    values = c('T cells' = '#4CA64C', 'Immune cells (non T cell)' = '#B33636', 'Stromal resident cells' = '#306a91'),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.sigCategory <- cowplot::get_legend(
  plot.sigCategory + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.sigCategory <- plot.sigCategory + theme(legend.position = 'none')



# order samples
signatureScoreMarkers_globalSigScore$variable <- factor(signatureScoreMarkers_globalSigScore$variable, levels = orderSampleNames)

# plot diff between T cell signatures and stroma signatures
plot_sigImbalaceTcellsStroma <- ggplot(signatureScoreMarkers_globalSigScore, aes(x = variable, y = "T cell-to-Stroma Enrichment (TSE) score", fill = signatureImbalance_tmp_TcellsStroma)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = '#FF0000', mid = 'white', high = '#006D2C',
                       breaks = c(-1, -0.5, 0, 0.5, 1), labels = c('<-1', '-0.5', '0', '0.5', '>1'),
                      guide = guide_colorbar(title = 'Stromal-to-T cell score\n(Stromal-to-Immune cell score)',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_sigImbalaceTcellsStroma <- cowplot::get_legend(
  plot_sigImbalaceTcellsStroma + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_sigImbalaceTcellsStroma <- plot_sigImbalaceTcellsStroma + theme(legend.position = 'none')


# signature imbalance category
cutOffImbalance <- 0.5
signatureScoreMarkers_globalSigScore <- signatureScoreMarkers_globalSigScore %>%
  dplyr::mutate(sigImbalanceCategory = ifelse(signatureImbalance_TcellsStroma > cutOffImbalance, 'Positive',
                                              ifelse(signatureImbalance_TcellsStroma < -cutOffImbalance, 'Negative', 'Neutral'))) %>%
  dplyr::mutate(sigImbalanceCategory_Resp = ifelse(sigImbalanceCategory != "Negative", ifelse(Outcome6M == "responder", "Positive+Neutral (responder)", "Positive+Neutral (non-responder)"), "Negative"),
                sigImbalanceCategory = factor(sigImbalanceCategory, levels = c('Positive', 'Neutral', 'Negative'))) %>%
  dplyr::mutate(sigImbalanceCategory_Resp = factor(sigImbalanceCategory_Resp, levels = c("Positive+Neutral (responder)", "Positive+Neutral (non-responder)", "Negative"))) %>%
  dplyr::mutate(sigImbalanceCategory_ImmuneCells = ifelse(signatureImbalance_ImmuneCellsStroma > cutOffImbalance, 'Poitive',
                                                         ifelse(signatureImbalance_ImmuneCellsStroma < -cutOffImbalance, 'Negative', 'Neutral')),
                sigImbalanceCategory_ALLImmuneCells = ifelse(signatureImbalance_ALLimmuneCellsStroma > cutOffImbalance, 'Positive',
                                                         ifelse(signatureImbalance_ALLimmuneCellsStroma < -cutOffImbalance, 'Negative', 'Neutral')))

# order data
signatureScoreMarkers_globalSigScore <- signatureScoreMarkers_globalSigScore %>%
  arrange(factor(variable, levels = levels(variable)))

# plot signature imbalance category
plot_sigImbalaceTCellsStroma_category <- ggplot(signatureScoreMarkers_globalSigScore, aes(x = variable, y = "TSE Category", fill = sigImbalanceCategory)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C'),
                       guide = guide_legend(title = 'STS Category',
                                              title.position = 'top', title.hjust = 0, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +

  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_sigImbalaceTCellsStroma_category <- cowplot::get_legend(
  plot_sigImbalaceTCellsStroma_category + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_sigImbalaceTCellsStroma_category <- plot_sigImbalaceTCellsStroma_category + theme(legend.position = 'none')


cowplot::plot_grid(plot_sigScoreRNAcluster,
                   plot_sigImbalaceTcellsStroma,
                   plot_sigImbalaceTCellsStroma_category,
                   ncol = 1, rel_heights = c(1, 0.2, 0.2),
                   align = 'vh', axis = 'tblr')


# plot dendogram (modify manually in pdf as cluster of negative imbalance should go last)
dend <- as.dendrogram(ConsensusMatrix_immuneSig[[numClusters]]$consensusTree ) %>%
  dendextend::set("branches_k_color", k = 3, value = as.character(c(color_immunoClusters[2], color_immunoClusters[1], color_immunoClusters[3]))) %>% dendextend::set("branches_lwd", 0.7) %>%
  dendextend::set("labels", "")
ggd1 <- dendextend::as.ggdend(dend)
plot_dendro <- ggplot(ggd1, horiz = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



#================================== Expression of Selected Genes ===============================
#===============================================================================================
# Get only PDCD1 and CD274
normMatrixRNA_CD274 <- normalizedCountMatrix[rownames(normalizedCountMatrix) %in% c("PDCD1", "CD274"),]
normMatrixRNA_CD274  <- normMatrixRNA_CD274 - rowMedians(normMatrixRNA_CD274)

# Get data frame with APOBEC expression separately
normMatrixRNA_CD274 <- as.data.frame(t(normMatrixRNA_CD274))
normMatrixRNA_CD274$sample <- rownames(normMatrixRNA_CD274)

# melt data
normMatrixRNA_CD274 <- reshape2::melt(normMatrixRNA_CD274, id.vars = c("sample"), variable.name = "gene", value.name="Expression")

# get max Values
normMatrixRNA_CD274 <- normMatrixRNA_CD274 %>% dplyr::mutate(Expression_tmp = ifelse(Expression > 1, 1, 
                                                                            ifelse(Expression < -1, -1, Expression)))

# order samples
normMatrixRNA_CD274$sample <- factor(normMatrixRNA_CD274$sample, levels = orderSampleNames)
normMatrixRNA_CD274$gene <- factor(normMatrixRNA_CD274$gene, levels = c("PDCD1", "CD274"))

# plot
plot_CD274exp <- ggplot(normMatrixRNA_CD274, aes(x = sample, y = gene, fill = Expression_tmp)) +
  geom_tile(size = 0.25, na.rm = F, alpha = 1, height = 1, width = 1) +
  scale_fill_gradient2(low = "green", mid = 'black', high = 'red',
                       breaks=c(-1, 0, 1),labels=c("<-1", 0, ">1"),
                       guide = guide_colorbar(title = 'Normalized expression\n(mean centered)',
                                              title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5)) +
  annotationTheme() +
  labs(x = NULL, y = "Gene expression") +
  theme(
    legend.position = 'bottom',
    text=element_text(size=8, family='Helvetica'),
    axis.title.y = element_text(size = 10, angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_CD274exp <- cowplot::get_legend(
  plot_CD274exp + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_CD274exp <- plot_CD274exp + theme(legend.position = 'none')




# protein PD-L1 expression
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(PDL1exp = ifelse(PDL1 == "pos", "pos", 'non-pos'),
                                                   sampleId = factor(sampleId, orderSampleNames))
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(PDL1exp = factor(PDL1exp, levels = c('pos', 'non-pos')))

# plot
plot_PDL1exp <- ggplot(DR176.MetaData, aes(x = sampleId, y = "PD-L1 protein expression", fill = PDL1exp)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = F) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('pos' = 'black', 'non-pos' = 'white', na.value = "grey60"),
                    labels = c('pos' = 'Positive', 'non-pos' = 'Negative', 'na.value' = "NA"),
                    guide = guide_legend(title = 'PD-L1 CPS',
                                         title.position = 'top', title.hjust = 0, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_PDL1exp <- cowplot::get_legend(
  plot_PDL1exp + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_PDL1exp <- plot_PDL1exp + theme(legend.position = 'none')



# compare expression between responders and non responders
normMatrixRNA_CD274 <- normMatrixRNA_CD274 %>% dplyr::left_join(signatureScoreMarkers_globalSigScore, by = c("sample" = "variable")) %>%
  dplyr::left_join(resultsRNA_immuneClusters, by = 'sample')
normMatrixRNA_CD274_stats_RespNonResp <-  normMatrixRNA_CD274 %>% dplyr::group_by(gene) %>%
  dplyr::summarize(p.value = wilcox.test(Expression ~ Outcome6M)$p.value) %>% dplyr::ungroup()
normMatrixRNA_CD274_stats_RespNonResp$pAdj <- c(p.adjust(normMatrixRNA_CD274_stats_RespNonResp$p.value, method = "BH"))

normMatrixRNA_CD274$Outcome6M <- factor(normMatrixRNA_CD274$Outcome6M, levels = c("responder", "non-responder"))
normMatrixRNA_CD274$gene <- factor(normMatrixRNA_CD274$gene, levels = c('CD274', 'PDCD1'))

# Plot comparison of each cell type per RNA cluster
plot_CD274_stats_outcome6M <- ggplot(normMatrixRNA_CD274, aes(x = Outcome6M, y = Expression, fill = Outcome6M)) +
  geom_boxplot(notch = F, width = .33, alpha = 1, outlier.shape = NA, size = 0.25) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y.., 2))) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=normMatrixRNA_CD274_stats_RespNonResp,
            aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 2,
                label=paste0("Adj. p = ", format(round(pAdj, 3)))),
            size = 3, inherit.aes = FALSE) +
  #expand_limits(y=0) +
  labs(x = "Response to treatment", y = 'RNA gene expression') +
  scale_fill_manual(guide = guide_legend(title = "Response to treatment", title.position = 'top', title.hjust = 0.5,
                                         ncol = 2, keywidth = 0.75, keyheight = 0.75), name = NULL,
                    values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder')) +
  theme(
    legend.position = 'bottom',
    text=element_text(size=10, family='Helvetica'),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# compare expression between responders and non-responders based on signature imbalance
normMatrixRNA_CD274_stats_RespNonResp <-  normMatrixRNA_CD274 %>%
  dplyr::filter(sigImbalanceCategory != "Neutral") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(p.value = wilcox.test(Expression ~ sigImbalanceCategory)$p.value) %>% dplyr::ungroup()
normMatrixRNA_CD274_stats_RespNonResp$pAdj <- c(p.adjust(normMatrixRNA_CD274_stats_RespNonResp$p.value, method = "BH"))
normMatrixRNA_CD274_stats_RespNonResp <- normMatrixRNA_CD274_stats_RespNonResp %>%
  dplyr::mutate(pAdj_label = ifelse(pAdj < 0.001, 'Adj. p < 0.001', paste0("Adj. p = ", format(round(pAdj, 3))) ))

# Plot comparison of each cell type per sig. imbalance category
plot_CD274_stats_SigImbalance <- ggplot(normMatrixRNA_CD274, aes(x = sigImbalanceCategory, y = Expression, fill = sigImbalanceCategory)) +
  geom_boxplot(notch = F, width = .33, alpha = 1, outlier.shape = NA, size = 0.25) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y.., 2))) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=normMatrixRNA_CD274_stats_RespNonResp,
            aes(x = -Inf, y = Inf, hjust = -0.3, vjust = 2,
                label=pAdj_label),
            size = 3, inherit.aes = FALSE) +
  labs(x = "Signature imbalance category", y = 'RNA gene expression') +
  scale_fill_manual(guide = guide_legend(title = "Signature imbalance category", title.position = 'top', title.hjust = 0.5,
                                         ncol = 2, keywidth = 0.75, keyheight = 0.75), name = NULL,
                    values =  c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C')) +
  theme(
    legend.position = 'bottom',
    text=element_text(size=10, family='Helvetica'),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# export Figure S6
pdf(paste0(odir,"CD274expression_stats.pdf"), width = 9, height = 2.5)
cowplot::plot_grid(plot_CD274_stats_outcome6M, plot_CD274_stats_SigImbalance,
                   rel_widths = c(0.45, 0.55),  ncol = 2, align = 'h', axis = 'tblr')
dev.off()



# Continue with plotting
# ----------------------------------------------------------------------------------------------------------------

# order data
DR176.MetaData$sampleId <- factor(DR176.MetaData$sampleId, levels = orderSampleNames)
DR176.MetaData$Outcome6M <- factor(DR176.MetaData$Outcome6M, levels = c("responder", "non-responder"))

# plot Responder labels
plot.Responders <- ggplot(DR176.MetaData, aes(sampleId, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = 'Response to treatment', title.position = 'top',
                                         title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.Responders <- cowplot::get_legend(
  plot.Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.Responders <- plot.Responders + theme(legend.position = 'none')


# Plot TMB
#order data
TMBresults <- results.Cohort$mutationalBurden %>% dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::mutate(tmbStatus = ifelse(tmbStatus == 'High TMB (≥10)', 'High', 'Low' )) %>%
  dplyr::mutate(tmbStatus = factor(tmbStatus, levels = c('High', 'Low')),
                sample = factor(sample, levels = orderSampleNames))

# plot TMB > 10
plot.TMB <- ggplot(TMBresults, aes(sample, y = 'Tumor mutational burden (TMB)', fill = tmbStatus)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('High' = '#8E1B37', "Low" = '#7A9BCF'),
                    guide = guide_legend(title = "TMB", title.position = 'top', title.hjust = 0,
                                         nrow = 2, keywidth = 0.5, keyheight = 0.5)) + 
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.TMB <- cowplot::get_legend(
  plot.TMB + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.TMB <- plot.TMB + theme(legend.position = 'none')


# plot APOBEC enrichment analysis
# order samples

APOBEC_enrichResults <- data.Cohort$APOBEC_enrich %>% dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSampleNames))

# plot APOBEC mutagenesis low, medium, high, and non
plot.APOBECenrich <- ggplot(APOBEC_enrichResults, aes(sample, y = 'APOBEC mutagenesis', fill = APOBEC_mutagenesis)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  scale_fill_manual(values = c('No' = 'white', 'Low' = '#F5BFAC', 'Medium' = '#EC7F59', 'High' = '#BB0303'),
                    guide = guide_legend(title = 'APOBEC mutagenesis',
                                         title.position = 'top',
                                         title.hjust = 0.5, ncol = 3,
                                         keywidth = 0.5, keyheight = 0.5)) +
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.APOBECenrich <- cowplot::get_legend(
  plot.APOBECenrich + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.APOBECenrich <- plot.APOBECenrich + theme(legend.position = 'none')

# plot tumor purity
plot_tumorPurity <- ggplot(DR176.MetaData, aes(x = sampleId, y = "Tumor purity", fill = tumorPurity)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = 'white', mid = '#1d5c8f', high = '#1d5c8f', midpoint = 0.7,
                       guide = guide_colorbar(title = 'Tumor purity',
                                              title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_tumorPurity <- cowplot::get_legend(
  plot_tumorPurity + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_tumorPurity <- plot_tumorPurity + theme(legend.position = 'none')


# ---- plot mUC subtypes
load(paste0(pathHMF, "RData/results_RNAmUCsubtype.RData"))

results_RNAmUCsubtype <- results_RNAmUCsubtype %>%
  dplyr::mutate(sampleId = factor(sampleId, level = orderSampleNames))



#brewer.pal(n = 4, name = "Dark2")
# plot cluster number labels
plot_RNASubtype_DR176 <- ggplot(results_RNAmUCsubtype, aes(sampleId, y = "Transcriptomic mUC subtype", fill = RNA_mUCsubtype)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T, alpha=1) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = colorRNAsubtypes,
                    guide = guide_legend(title = 'Transcriptomic mUC subtype', title.position = 'top',
                                         title.hjust = 0.5, nrow = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_RNASubtype_DR176 <- cowplot::get_legend(
  plot_RNASubtype_DR176 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_RNASubtype_DR176 <- plot_RNASubtype_DR176 + theme(legend.position = 'none')


# order samples
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(sampleId = factor(sampleId, levels = orderSampleNames))

# plot patients with pre-treatment
plot.preTreatment <- ggplot(DR176.MetaData, aes(sampleId, y = 'Systemic pretreatment', fill = hasSystemicPreTreatment)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('No' = 'white', 'Yes' = 'black'),
                    guide = guide_legend(title = 'Systemic pretreatment',
                                         title.position = 'top',
                                         title.hjust = 0, ncol = 3,
                                         keywidth = 0.5, keyheight = 0.5)) + 
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.preTreatment <- cowplot::get_legend(
  plot.preTreatment + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.preTreatment <- plot.preTreatment + theme(legend.position = 'none')



DR176.MetaData <- DR176.MetaData %>%
  dplyr::mutate(primaryCancerSubtype = ifelse(primaryTumorSubLocation == "Pyelum", "UTUC", primaryTumorSubLocation)) %>%
  dplyr::mutate(primaryCancerSubtype = factor(primaryCancerSubtype, levels = c("Bladder", "UTUC")))

plot.cancerSubtype <- ggplot(DR176.MetaData, aes(sampleId, y = 'Primary tumor location', fill = primaryCancerSubtype)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + annotationTheme() + 
  scale_fill_manual(values = c('grey60','#ffe596', '#b3e6ff', '#41406d', '#fffffd'),
                    labels = c('Bladder', 'UTUC', 'Unknown'),
                    guide = guide_legend(title = "Primary tumor location", title.position = 'top', title.hjust = 0.5,
                                         nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.cancerSubtype <- cowplot::get_legend(
  plot.cancerSubtype + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.cancerSubtype <- plot.cancerSubtype + theme(legend.position = 'none')





# Plot TCR diversity (previously analyzed => see below)
load(paste0(pathHMF, 'RData/data.MiXCR.RData'))

# count number of unique clonotypes and hyperexpanded
plotData.numClones <- ldply(data.MiXCR, data.frame, .id = "sample") %>%
  dplyr::group_by(sample) %>% dplyr::mutate(total = sum(Read.count)) %>% dplyr::ungroup() %>%
  dplyr::filter(total >= 100) %>%
  dplyr::mutate(cloneCategory = ifelse(Read.count == 1, "Rare (1 clonotype)",
                                       ifelse(Read.proportion < 0.01, "Small (>1 clonotype - 1%)",
                                              ifelse(Read.proportion < 0.1, "Large (1%-10%)", "Hyperexpanded (>10%)")))) %>%
  dplyr::group_by(sample, cloneCategory) %>% dplyr::mutate(rel_cloneCategory = sum(Read.proportion)) %>%
  dplyr::ungroup() %>% dplyr::distinct(sample, cloneCategory, .keep_all = TRUE) %>%
  dplyr::right_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M),
                    by = c("sample" = "sampleId"))

# Evaluate the diversity of clones by the ecological diversity index.
plotData.diversity <- reshape2::melt(data.frame(Diversity = sapply(data.MiXCR, function (x) tcR::diversity(x$Read.count)), sample = names(data.MiXCR)))
plotData.diversity <- plotData.diversity %>%
  dplyr::right_join(dplyr::distinct(plotData.numClones, sample, .keep_all = TRUE), by = "sample") %>%
  dplyr::mutate(value = ifelse(is.na(total), NA, value))
plotData.diversity <- plotData.diversity %>% dplyr::mutate(diversity_tmp = ifelse(value > 100, 100, value),
                                                           sample = factor(sample, levels = orderSampleNames),
                                                           Outcome6M = factor(Outcome6M, levels = c("responder", "non-responder")))


# plot diversity stats
plot_diversity <- ggplot(plotData.diversity, aes(x = sample, y = "TCR diversity", fill = diversity_tmp)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = F) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = 'white', mid = 'white', high = '#fa0a8a', midpoint = 10, limits = c(0,100),
                      breaks = c(0, 50, 100),
                      labels = c(0, 50, ">100"),
                      guide = guide_colorbar(title = 'Diversity index (TCR)',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_diversity <- cowplot::get_legend(
  plot_diversity + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_diversity <- plot_diversity + theme(legend.position = 'none')



# add extra label for biopsy site (only include liver and lymph node)
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(biopsySite_tmp = ifelse(`Biopsy site annotation` %in% c("Lymph node", "Liver"), `Biopsy site annotation`, "Other"))

# plot biopsy location
plot.biopsySite <- ggplot(DR176.MetaData, aes(sampleId, y = 'Biopsy site', fill = biopsySite_tmp)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('cornflowerblue', 'goldenrod2', 'grey90'),
                    guide = guide_legend(title = 'Biopsy site',
                                         title.position = 'top',
                                         title.hjust = 0, ncol = 3,
                                         keywidth = 0.5, keyheight = 0.5)) +
  
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# extract legends
legend_plot.biopsySite <- cowplot::get_legend(
  plot.biopsySite + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.biopsySite <- plot.biopsySite + theme(legend.position = 'none')


# plot sex of patient
plot.sex <- ggplot(TMBresults, aes(sample, y = 'Female', fill = gender)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('black', 'white')) + 
  annotationTheme() + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# Export plots transcriptomic overview ---------------------------------------------------

PAR_rel_heights <- c(rep(0.03, 10), 0.06, 0.1, 1, 0.03, 0.03, 0.03, 0.03)

samplePlot <- cowplot::plot_grid(
  plot.Responders,
  plot.TMB,
  plot.APOBECenrich,
  plot.sex,
  plot_RNASubtype_DR176,
  plot.biopsySite,
  plot.cancerSubtype,
  plot_tumorPurity,
  plot.preTreatment,
  plot_PDL1exp,
  plot_CD274exp,
  plot_dendro,
  plot_sigScoreRNAcluster,
  plot_sigImbalaceTcellsStroma,
  plot_sigImbalaceImmuneCellsStroma,
  plot_sigImbalaceTCellsStroma_category,
  
  plot_diversity,
  ncol = 1, align = 'v',
  axis = 'tblr',
  rel_heights = PAR_rel_heights
)

# Plot sample-level overview.
samplePlot_2 <- cowplot::plot_grid(
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
  plot.sigCategory,
  NULL, NULL, NULL, NULL, NULL,
  ncol = 1, align = 'v',
  axis = 'tblr',
  rel_heights = PAR_rel_heights
)

legend_samplePlot <- cowplot::plot_grid(
  legend_plot.Responders, legend_plot_RNASubtype_DR176, legend_plot_CD274exp, legend_plot_sigImbalaceTcellsStroma, legend_plot_tumorPurity,
  
  
  legend_plot.APOBECenrich, legend_plot.cancerSubtype, legend_plot_sigScoreRNAcluster, legend_plot_sigImbalaceTCellsStroma_category, legend_plot.sigCategory,
  
  
  legend_plot.TMB, legend_plot_PDL1exp, legend_plot.preTreatment, legend_plot_diversity, legend_plot.biopsySite,
  nrow = 3, align = 'v', rel_heights = c(0.5, 0.5, 0.5),
  rel_widths = c(1, 0.5, 0.5, 0.5, 0.5))


# Export results (dendogram has to be manually modified as cluster 3 (negative STS) has be last)
# Figure 4
pdf(paste0(odir,"RNA_overviewClusters.pdf"),width = 9, height = 10)#, width = 14, height = 21)
cowplot::plot_grid(samplePlot, samplePlot_2, legend_samplePlot,
                   ncol = 2,
                   rel_heights = c(1, 0.2),
                   rel_widths = c(1, 0.025),
                   align = 'vh', axis = 'tblr')
dev.off()





#============================  Immune Infiltration (cell fraction) =============================
#===============================================================================================
load(paste0(pathHMF, "RData/globalSigScore.RData"))

# First, order samples based on Signature imbalance and
signatureScoreMarkers_globalSigScore <- signatureScoreMarkers_globalSigScore %>%
  dplyr::arrange(-signatureImbalance_TcellsStroma) %>%
  dplyr::arrange(factor(Outcome6M, levels = c("responder", "non-responder"))) %>%
  dplyr::arrange(factor(sigImbalanceCategory, levels = c("Positive", "Neutral", "Negative")))
  
orderSample_sigImbalance <- as.character(signatureScoreMarkers_globalSigScore$variable)


 # plot Sig imbalance category
# order data
signatureScoreMarkers_globalSigScore$variable <- factor(signatureScoreMarkers_globalSigScore$variable, levels = orderSample_sigImbalance)
# plot signature imbalance category
plot_sigImbalaceTCellsStroma_category <- ggplot(signatureScoreMarkers_globalSigScore, aes(x = variable, y = "STS category", fill = sigImbalanceCategory)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C'),
                    guide = guide_legend(title = 'TSE category',
                                         title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_sigImbalaceTCellsStroma_category <- cowplot::get_legend(
  plot_sigImbalaceTCellsStroma_category + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_sigImbalaceTCellsStroma_category <- plot_sigImbalaceTCellsStroma_category + theme(legend.position = 'none')




#  Immune cell fraction analysis
# Melt results
immuneResults_molten <- reshape2::melt(immuneCFResults, id.vars = "cell_type")

# Delete uncharacterized cell
immuneResults_molten <- dplyr::filter(immuneResults_molten, cell_type != "uncharacterized cell")

#Add cluster information
immuneResults_molten <- resultsRNA_immuneClusters %>% dplyr::full_join(immuneResults_molten, by = c("sample" = "variable"))


# Order samples to match heat map of RNA TopGenes clusters
immuneResults_molten$sample <- factor(immuneResults_molten$sample, levels = orderSample_sigImbalance)

plot_immuneCells <- ggplot(immuneResults_molten, aes(sample, y = value, fill = cell_type)) +
  geom_bar(stat='identity', colour = 'black', width = 1, size = 0.2) +
  scale_y_continuous(expand = expand_scale(0,0)) +
  labs(y = "Immune cell\nfraction", x = NULL) +
  scale_fill_brewer(palette = 'Paired', guide = guide_legend(title = 'Immune cell type', title.position = 'top',
                                                             title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme(
    legend.position = 'right',
    text=element_text(size=8, family='Helvetica'),
    axis.text=element_text(size=10),
    axis.title.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# extract legends
legend_plot_immuneCells <- cowplot::get_legend(
  plot_immuneCells + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_immuneCells <- plot_immuneCells + theme(legend.position = 'none')



# Get immune cell fraction stats for responder vs Non-responder

#  Perform Kruskal-Wallis test
immuneResults_molten_2 <- immuneResults_molten %>% dplyr::group_by(sample) %>%
  dplyr::mutate(totalImmuneCell = sum(value), cell_type = "Total immune cell") %>%
  dplyr::mutate(value = totalImmuneCell, totalImmuneCell = NULL) %>%
  dplyr::ungroup() %>% dplyr::distinct()
immuneResults_molten <- rbind(immuneResults_molten, immuneResults_molten_2)
immuneResults_molten <- immuneResults_molten %>%
  dplyr::left_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M), by = c("sample" = "sampleId")) %>%
  dplyr::full_join(dplyr::select(signatureScoreMarkers_globalSigScore, variable, sigImbalanceCategory), by = c('sample' = 'variable'))

# pValue for sig imbalance group
immuneResults_molten_statResults_sigImbalanceGroup <- dplyr::filter(immuneResults_molten, sigImbalanceCategory != "Neutral") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(p.value = wilcox.test(value ~ sigImbalanceCategory)$p.value) %>% dplyr::ungroup()
immuneResults_molten_statResults_sigImbalanceGroup$pAdj <- c(p.adjust(immuneResults_molten_statResults_sigImbalanceGroup$p.value[1:10], method = "BH"),
                                                             immuneResults_molten_statResults_sigImbalanceGroup$p.value[11])
# pValue for resopnder vs non-responder
immuneResults_molten_statResults_respNonresp <- immuneResults_molten %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(p.value = wilcox.test(value ~ Outcome6M)$p.value) %>% dplyr::ungroup()
immuneResults_molten_statResults_respNonresp$pAdj <- c(p.adjust(immuneResults_molten_statResults_respNonresp$p.value[1:10], method = "BH"),
                                                       immuneResults_molten_statResults_respNonresp$p.value[11])

# pValue for cluster 1 vs cluster 2
immuneResults_molten_statResults_sigCluster <- dplyr::filter(immuneResults_molten, immunoCluster != "Cluster 2") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(p.value = wilcox.test(value ~ immunoCluster)$p.value) %>% dplyr::ungroup()
immuneResults_molten_statResults_sigCluster$pAdj <- c(p.adjust(immuneResults_molten_statResults_sigCluster$p.value[1:10], method = "BH"),
                                                      immuneResults_molten_statResults_sigCluster$p.value[11]) 

# order outcome
immuneResults_molten <- immuneResults_molten %>%
  dplyr::mutate(Outcome6M = factor(Outcome6M, levels = c("responder", "non-responder")))


# Plot comparison of each cell type comparing Responder vs Non-responder
plot_immuneCells_stats_outcome6M <- ggplot(immuneResults_molten, aes(x = Outcome6M, y = value, fill = Outcome6M)) +
  geom_boxplot(notch = F, width = .33, alpha = .9, outlier.shape = NA, size = 0.25) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  facet_wrap(~cell_type, scales = "free", nrow = 3) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=immuneResults_molten_statResults_respNonresp,
            aes(x = -Inf, y = Inf, hjust = -0.4, vjust = 2,
                label=paste0("Adj. p = ", format(round(pAdj, 2)))),
            size = 2.5, inherit.aes = FALSE) +
  expand_limits(y=0) +
  labs(x = NULL, y = 'Immune cell fraction') +
  scale_fill_manual(guide = guide_legend(title = "Outcome", title.position = 'top', title.hjust = 0.5,
                                         ncol = 1, keywidth = 0.75, keyheight = 0.75), name = NULL,
                    values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder')) +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 9),
    text=element_text(size=9, family='Helvetica'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )

# Plot comparison of each cell type per Sig. imbalance category
plot_immuneCells_stats_sigImbalanceGroup <- ggplot(immuneResults_molten, aes(x = sigImbalanceCategory, y = value, fill = sigImbalanceCategory)) +
  geom_boxplot(notch = F, width = .5, alpha = 1, outlier.shape = NA, size = 0.25) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  facet_wrap(~cell_type, scales = "free", nrow = 3) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=immuneResults_molten_statResults_sigImbalanceGroup,
            aes(x = -Inf, y = Inf, hjust = -0.4, vjust = 2,
                label=paste0("Adj. p = ", format(round(pAdj, 2)))),
            size = 2.5, inherit.aes = FALSE) +
  expand_limits(y=0) +
  labs(x = NULL, y = 'Immune cell fraction') +
  scale_fill_manual(guide = guide_legend(title = "Sig. imbalance\ncategory", title.position = 'top',
                                         title.hjust = 0.5, ncol = 1, keywidth = 0.75, keyheight = 0.75),
                    values =  c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C')) +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 9),
    text=element_text(size=9, family='Helvetica'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )




# Export Figure S8
pdf(paste0(odir,"ImmunoCellsStats_respNonresp.pdf"), width = 7, height = 3)
plot_immuneCells_stats_outcome6M
dev.off()

pdf(paste0(odir,"ImmunoCellsStats_sigImbalance.pdf"), width = 7, height = 3)
plot_immuneCells_stats_sigImbalanceGroup
dev.off()






# Import MiXCR files ---------------------------------------------------------
library(tcR)
# Fix integer problem.
lapply(list.files(paste0(pathHMF, 'RNAseq/MiXCR/'), pattern = 'output.txt', full.names = T), function(x){
  data.MiXCR <- readr::read_delim(x, delim = '\t')
  data.MiXCR$cloneCount <- as.integer(data.MiXCR$cloneCount)
  write.table(data.MiXCR, gsub('.txt', '_fix.txt', x), row.names = F, quote = F, sep = '\t')
})


# Import with tcR.
data.MiXCR <- lapply(list.files(paste0(pathHMF, 'RNAseq/MiXCR/'), pattern = 'output_fix.txt', full.names = T), tcR::parse.mixcr)

# Add sample names.
names(data.MiXCR) <- gsub('_output.*', '', basename(list.files(paste0(pathHMF, 'RNAseq/MiXCR/'), pattern = '_fix.txt', full.names = T)))






# Analyze MiXCR -----------------------------------------------------------

# count number of unique clonotypes and hyperexpanded
plotData.numClones <- ldply(data.MiXCR, data.frame, .id = "sample") %>%
  dplyr::group_by(sample) %>% dplyr::mutate(total = sum(Read.count)) %>% dplyr::ungroup() %>%
  dplyr::filter(total >= 100) %>%
  dplyr::mutate(cloneCategory = ifelse(Read.count == 1, "Rare (1 clonotype)",
                                       ifelse(Read.proportion < 0.01, "Small (>1 clonotype - 1%)",
                                              ifelse(Read.proportion < 0.1, "Large (1%-10%)", "Hyperexpanded (>10%)")))) %>%
  dplyr::group_by(sample, cloneCategory) %>% dplyr::mutate(rel_cloneCategory = sum(Read.proportion)) %>%
  dplyr::ungroup() %>% dplyr::distinct(sample, cloneCategory, .keep_all = TRUE) %>%
  dplyr::right_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M),
                    by = c("sample" = "sampleId"))

# plot distribution of TCR repertoire
# order samples
plotData.numClones$sample <- factor(plotData.numClones$sample, levels = orderSample_sigImbalance)
plotData.numClones$cloneCategory <- factor(plotData.numClones$cloneCategory, levels = c("Hyperexpanded (>10%)", "Large (1%-10%)",
                                                                                        "Small (>1 clonotype - 1%)", "Rare (1 clonotype)"))

# Barplot of clone sizes per Sample
plot_cloneSize <- ggplot(plotData.numClones, aes(x = sample, y = rel_cloneCategory*total, fill = cloneCategory)) +
  coord_trans(y = 'sqrt') +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2, width = 1) +
  scale_y_continuous(breaks = c(0, 1000, 10000, 25000, 50000, 75000), expand = expand_scale(0,0), labels = scales::comma) +
  scale_x_discrete(expand = expand_scale(0,0)) +
  scale_fill_manual('Clonotypes', values = c('Rare (1 clonotype)' = 'orchid', 'Large (1%-10%)' = '#D29958',
                                             "Small (>1 clonotype - 1%)" = '#677EAA', 'Hyperexpanded (>10%)' = "black"),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  labs(x = NULL, y = 'Clonotype size\n(TCR repertoire)') +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey60', linetype = 'longdash'),
    panel.grid.minor.y = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )

# extract legends
legend_plot_cloneSize <- cowplot::get_legend(
  plot_cloneSize + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_cloneSize <- plot_cloneSize + theme(legend.position = 'none')


# Proportion of clones per Sample
plot_cloneSizeProp <- ggplot(plotData.numClones, aes(x = sample, y = rel_cloneCategory, fill = cloneCategory)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = expand_scale(0,0)) +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2, width = 1) +
  scale_fill_manual('TCR clonotypes', values = c('Rare (1 clonotype)' = 'orchid', 'Large (1%-10%)' = '#D29958',
                                                 "Small (>1 clonotype - 1%)" = '#677EAA', 'Hyperexpanded (>10%)' = "black"),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  labs(x = NULL, y = 'Rel. abundance\n(TCR repertoire)') +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )

# extract legends
legend_plot_cloneSizeProp <- cowplot::get_legend(
  plot_cloneSizeProp + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_cloneSizeProp <- plot_cloneSizeProp + theme(legend.position = 'none')



# evaluate clonotype size and stats
plotData.numClones_stats <- plotData.numClones %>% dplyr::filter(!is.na(rel_cloneCategory)) %>%
  tidyr::complete(sample, cloneCategory) %>% dplyr::mutate(Outcome6M = NULL) %>%
  dplyr::mutate(rel_cloneCategory = ifelse(is.na(rel_cloneCategory), 0, rel_cloneCategory)) %>%
  dplyr::left_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M), by = c('sample' = 'sampleId')) %>%
  dplyr::left_join(dplyr::select(signatureScoreMarkers_globalSigScore, variable, sigImbalanceCategory), by = c('sample' = 'variable')) %>%
  dplyr::left_join(resultsRNA_immuneClusters, by = 'sample') %>%
  dplyr::mutate(sample = factor(sample, levels = orderSampleNames))
  
# get pValues for Responder vs Non-responder
plotData.numClones_pValue <- plotData.numClones_stats %>%
  dplyr::group_by(cloneCategory) %>%
  dplyr::summarize(p.value = wilcox.test(rel_cloneCategory ~ Outcome6M)$p.value) %>% dplyr::ungroup()
plotData.numClones_pValue$pAdj <- p.adjust(plotData.numClones_pValue$p.value, method = "BH")


# plot clonotype size stats
plot_clonotypeSize_stat_RespNonResp <- ggplot(plotData.numClones_stats, aes(x = Outcome6M, y = rel_cloneCategory, fill = Outcome6M)) +
  stat_boxplot(geom ='errorbar', width = 0) +
  geom_boxplot(width = .33, notch = F, outlier.shape = NA, alpha = 1) +
  ylim(c(0, 1)) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y.., 2))) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  facet_wrap(~cloneCategory, nrow = 1) +
  geom_text(data=plotData.numClones_pValue,
            aes(x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
                label=paste0("Adj. p = ", format(round(pAdj, 2)))),
            size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = 'Response', title.position = 'top',
                                         title.hjust = 0.5, ncol = 2, keywidth = 0.75, keyheight = 0.75)) +
  labs(x = 'Response to treatment', y = 'Relative clone size\n(TCR-Vb Repertoire)') + 
  # Theme
  theme(legend.position = 'bottom', 
        axis.text.x = element_blank(),
        text=element_text(size=10, family='Helvetica'),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        strip.background = element_rect(colour = 'grey20', fill = 'white'),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  )

# p-values for sig imbalance category
plotData.numClones_pValue <- plotData.numClones_stats %>%
  dplyr::filter(sigImbalanceCategory!= 'Neutral') %>%
  dplyr::group_by(cloneCategory) %>%
  dplyr::summarize(p.value = wilcox.test(rel_cloneCategory ~ sigImbalanceCategory)$p.value) %>% dplyr::ungroup()
plotData.numClones_pValue$pAdj <- p.adjust(plotData.numClones_pValue$p.value, method = "BH")

# plot clonotype size stats
plot_clonotypeSize_stat_sigImbalance <- ggplot(plotData.numClones_stats, aes(x = sigImbalanceCategory, y = rel_cloneCategory, fill = sigImbalanceCategory)) +
  stat_boxplot(geom ='errorbar', width = 0) +
  geom_boxplot(width = 0.5, notch = F, outlier.shape = NA, alpha = 1) +
  ylim(c(0, 1)) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y.., 2))) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, cex = .25) +
  facet_wrap(~cloneCategory, nrow = 1) +
  geom_text(data=plotData.numClones_pValue,
            aes(x = -Inf, y = Inf, hjust = -0.2, vjust = 2,
                label=paste0("Adj. p = ", format(round(pAdj, 3)))),
            size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values =  c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C'),
                    guide = guide_legend(title = 'STS category', title.position = 'top',
                                         title.hjust = 0.5, ncol = 3, keywidth = 0.75, keyheight = 0.75)) +
  labs(x = 'T cell-to-Stroma Enrichment score', y = 'Relative clone size\n(TCR-Vb Repertoire)') + 
  # Theme
  theme(legend.position = 'bottom', 
        axis.text.x = element_blank(),
        text=element_text(size=10, family='Helvetica'),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        strip.background = element_rect(colour = 'grey20', fill = 'white'),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  )



# export Figure S9a
pdf(paste0(odir,"TCRCloneSize_stats.pdf"), width = 5, height = 4.75*1.5)
  cowplot::plot_grid(plot_clonotypeSize_stat_RespNonResp,
                     plot_clonotypeSize_stat_sigImbalance,
                     ncol = 1, align = 'hv', axis = 'tblr')
dev.off()



# Evaluate the diversity of clones by the ecological diversity index.
plotData.diversity <- reshape2::melt(data.frame(Diversity = sapply(data.MiXCR, function (x) tcR::diversity(x$Read.count)), sample = names(data.MiXCR)))
plotData.diversity <- plotData.diversity %>%
  dplyr::right_join(dplyr::distinct(plotData.numClones, sample, .keep_all = TRUE), by = "sample") %>%
  dplyr::mutate(value = ifelse(is.na(total), NA, value))
plotData.diversity <- plotData.diversity %>% dplyr::mutate(diversity_tmp = ifelse(value > 100, 100, value),
                                                           sample = factor(sample, levels = orderSample_sigImbalance),
                                                           Outcome6M = factor(Outcome6M, levels = c("responder", "non-responder")))


# plot diversity stats
plot_diversity <- ggplot(plotData.diversity, aes(x = sample, y = "TCR diversity", fill = diversity_tmp)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = F) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient(low = 'white', high = '#fa0a8a', limits = c(0,100),
                      breaks = c(0, 50, 100),
                      labels = c(0, 50, ">100"),
                      guide = guide_colorbar(title = 'Diversity index (TCR)',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_diversity <- cowplot::get_legend(
  plot_diversity + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_diversity <- plot_diversity + theme(legend.position = 'none')

# plot diversity stats
plot_diversity_stat_RespNonResp <- ggplot(plotData.diversity, aes(x = Outcome6M, y = value, fill = Outcome6M)) +
  stat_boxplot(geom ='errorbar', width = 0) +
  geom_boxplot(width = 0.33, notch = F, outlier.shape = NA, alpha = 1) +
  ylim(c(0, 165)) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y..))) +
  ggbeeswarm::geom_quasirandom(alpha = 0.33, cex = .25) +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = 'Response', title.position = 'top',
                                         title.hjust = 0.5, ncol = 1, keywidth = 0.75, keyheight = 0.75)) +
  labs(x = 'Response to treatment', y = 'Diversity index\n(TCR-Vb Repertoire)') + 
  # Theme
  theme(legend.position = 'bottom', 
        axis.text.x = element_blank(),
        text=element_text(size=10, family='Helvetica'),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  ) + 
  ggpubr::stat_compare_means(method = "wilcox.test",  comparisons = list(c("responder", "non-responder")),
                             size = 3, paired = F, label.x = 1.5, label = "p.format")


# diversity index between sig Imbalance category
plotData.diversity <- plotData.diversity %>%
  dplyr::full_join(dplyr::select(signatureScoreMarkers_globalSigScore, variable, sigImbalanceCategory), by = c('sample' = 'variable')) %>%
  dplyr::left_join(resultsRNA_immuneClusters, by = 'sample')


# plot diversity stats
plot_diversity_statSigImbalance <- ggplot(plotData.diversity, aes(x = sigImbalanceCategory, y = value, fill = sigImbalanceCategory)) +
  stat_boxplot(geom ='errorbar', width = 0) +
  geom_boxplot(width = .5, notch = F, outlier.shape = NA, alpha = 1) +
  ylim(c(0, 165)) +
  stat_summary(fun.y=median, colour="black", geom="text", size = 3, show.legend = FALSE, vjust=-1, angle = 90, hjust = .5, aes( label=round(..y..))) +
  ggbeeswarm::geom_quasirandom(cex = .25, alpha = .5) +
  scale_fill_manual(values =  c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C'),
                    guide = guide_legend(title = 'STS category', title.position = 'top',
                                         title.hjust = 0.5, ncol = 2, keywidth = 0.75, keyheight = 0.75)) +
  labs(x = 'Stromal-to-T cell\nscore', y = 'Diversity index\n(TCR-Vb Repertoire)') + 
  # Theme
  theme(legend.position = 'bottom', 
        axis.text.x = element_blank(),
        text=element_text(size=10, family='Helvetica'),
        panel.grid.major = element_line(colour = 'grey80', linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = 'white', colour = NA),
        panel.border = element_rect(fill = NA, colour = 'grey20')
  ) +
  ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(c("Negative", "Positive")),
                             size = 3, paired = F, label.x = 2.5)


# export Figure S9b
pdf(paste0(odir,"TCRdiversity_stats.pdf"), width = 5.2, height = 2.5)
cowplot::plot_grid(plot_diversity_stat_RespNonResp, plot_diversity_statSigImbalance,
                  ncol = 2, align = 'h', axis = 'tblr')
dev.off()






# order data
DR176.MetaData$sampleId <- factor(DR176.MetaData$sampleId, levels = orderSample_sigImbalance)
DR176.MetaData$Outcome6M <- factor(DR176.MetaData$Outcome6M, levels = c("responder", "non-responder"))

# plot Responder labels
plot.Responders <- ggplot(DR176.MetaData, aes(sampleId, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = 'Response to treatment', title.position = 'top',
                                         title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.Responders <- cowplot::get_legend(
  plot.Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.Responders <- plot.Responders + theme(legend.position = 'none')



# export immune analysis stratified by sig. imbalance

# export figure S7
pdf(paste0(odir,"ImmuneAndTCR_SigImbalance.pdf"), width = 9, height = 5)
cowplot::plot_grid(
  cowplot::plot_grid(plot_sigImbalaceTCellsStroma_category,
                     plot.Responders,
                   plot_immuneCells,
                   plot_cloneSize,
                   plot_cloneSizeProp,
                   plot_diversity,
                   ncol = 1, rel_heights = c(0.2, 0.2, 1.2, 1, 1, 0.2), align = 'v', axis = 'tblr'),
  cowplot::plot_grid(legend_plot_sigImbalaceTCellsStroma_category,
                     legend_plot.Responders,
                     legend_plot_immuneCells,
                     legend_plot_cloneSize,
                     legend_plot_diversity, ncol = 1, rel_heights = c(1, 0.5, 0.5, 1, 1), align = 'v', axis = 'tblr'),
  ncol = 2, rel_widths = c(1, 0.3))

dev.off()




#============ plot oncoplot and stats (positive vs negative TSE score) ===========================
#=================================================================================================


# First, order samples based on Signature imbalance and
signatureScoreMarkers_globalSigScore <- signatureScoreMarkers_globalSigScore %>%
  dplyr::arrange(-signatureImbalance_TcellsStroma) %>%
  dplyr::arrange(factor(Outcome6M, levels = c("responder", "non-responder"))) %>%
  dplyr::arrange(factor(sigImbalanceCategory, levels = c("Positive", "Neutral", "Negative")))

orderSample_sigImbalance <- as.character(signatureScoreMarkers_globalSigScore$variable)


# use the same code from mutLandscape analysis
# Melt the Gistic peaks
gisticPeaks <- reshape2::melt(data.frame(mcols(results.Cohort$GISTIC2$fullCohort$gisticNarrowPeaksWithAnno)), id.vars = c('Unique.Name', 'Descriptor', 'q.values', 'nGenes.GENCODE', 'overlapGenes.Final'))

#Change some col names
names(gisticPeaks)[names(gisticPeaks) == "nGenes.GENCODE"] <- "nGenes"
names(gisticPeaks)[names(gisticPeaks) == "overlapGenes.Final"] <- "overlappingGenes"

# keep samples
gisticPeaks$variable <- droplevels(gisticPeaks$variable)

# Only keep peaks with q < 0.05
q.keep <- 0.05
gisticPeaks <- gisticPeaks %>% dplyr::filter(q.values <= q.keep)

# Generate identifier. nGenes = 6
gisticPeaks <- gisticPeaks %>% dplyr::mutate(identifier = sprintf('%s (%s) - %s',
                                                                  Descriptor,
                                                                  ifelse(nGenes >= 6, sprintf('%s genes', nGenes), overlappingGenes),
                                                                  Unique.Name
))

# Add copy-number category
gisticPeaks <- gisticPeaks %>% dplyr::mutate(category = ifelse(value == 0, 'Neutral', ifelse(value == 1, 'Shallow', 'Deep')),
                                             category = paste(category, ifelse(grepl('Amplifi', identifier), 'Gain', 'Deletion')))


# Only keep deep gain / deletions
gisticPeaks <- droplevels(gisticPeaks)
gisticPeaks <- gisticPeaks %>% dplyr::filter(grepl('Deep', category))


#Get clean names of Genes
gisticPeaks$overlappingGenes <- gsub('\\\  \\(.*', '', gisticPeaks$overlappingGenes)
gisticPeaks$overlappingGenes <- gsub('\\\ \\(.*', '', gisticPeaks$overlappingGenes)

# Genes significant from GISTIC (check peaks with multiple genes, separate Add important genes to the list)
geneListGISTIC <- unique(gisticPeaks$overlappingGenes)
geneListGISTIC <- c(geneListGISTIC, 'ARID1A', "PDE4DIP ", "PPARG", "PDE4D", 'SOX4', "MTRES1", "AHR",
                    "KAT6A", "GATA3", "PAMR1", "CCND1", "MDM2", "ERBB2", "IKZF3", "LRP1B", "CCSER1", 
                    "BEND3", "PTEN", "AHR", "CDKN2A", "CDKN2B", "RB1", "WWOX", "EGR3", "LEPROTL1",
                    "ARID1B", "FOXO3", "PTEN", "HRAS", "TWIST2", "TSC2", "CREBBP", "FHIT")

# List of Known driver genes
geneListKnownDrivers <- c("TP53", "KDM6A", "PIK3CA", "KMT2D", "ARID1A", "FGFR3", "RB1",
                          "KMT2C", "CREBBP", "ERBB2", "STAG2", "ELF3", "EP300", "ERCC2",
                          "FAT1", "ATM", "ERBB3", "BIRC6", "CDKN1A", "KMT2A", "TSC1", "HRAS", "FBXW7", "BRCA2",
                          # the following are included from the mUC paper (n=116)
                          "RBM10", "ZFP36L1", "RHOA", "RXRA", "CNTNAP5", "RARG", "MGP", "FGF19", "NECTIN4")

#Get all significant genes and those from the known driver list
dataOncoplot <- results.Cohort$combinedReport %>% dplyr::filter(dNdS == "Significant" | SYMBOL %in% geneListKnownDrivers |
                                                                  SYMBOL %in% geneListGISTIC)


# Only keep coding mutations (incl. SVs) or CNA.
dataOncoplot <- dataOncoplot %>% dplyr::filter(!is.na(Consequence.Mut) | !is.na(Consequence.SV) | grepl('Deep', Consequence.CNA))

# Replace synonymous variants, neutral and shallow CN for NA to not show in plot
dataOncoplot <- dataOncoplot %>% dplyr::mutate(Consequence.CNA = ifelse(Consequence.CNA %in% c('Neutral', 'Deletion', 'Amplification'), NA, Consequence.CNA))

# rename deep amplification/deletion to just amplification and deletion
dataOncoplot <- dataOncoplot %>% dplyr::mutate(Consequence.CNA = gsub("Deep ", "", Consequence.CNA),
                                               Consequence.SV = gsub("Va", "va", Consequence.SV))

# Sometime results repeated for a sample (amplification and deletion in the same patient)
# discard these duplicate samples
dataOncoplot <- dataOncoplot %>% dplyr::distinct(sample, SYMBOL, .keep_all = TRUE)

# keep only included samples
dataOncoplot <- dataOncoplot %>% dplyr::filter(sample %in% DR176.MetaData$sampleId)

# number of samples
nSamples <- length(unique(DR176.MetaData$sampleId))

# Add axis names showing number of unique mutated samples.
dataOncoplot <- dataOncoplot %>% dplyr::group_by(ENSEMBL) %>% dplyr::add_tally() %>%
  dplyr::mutate(axisName = sprintf('%s (%s%%)', SYMBOL, round(100*n/nSamples))) %>% dplyr::ungroup()
dataOncoplot$sample <- factor(dataOncoplot$sample, level = DR176.MetaData$sampleId)
dataOncoplot <- dataOncoplot %>% dplyr::filter(n/nSamples > 0.05)
dataOncoplot <- dataOncoplot %>% tidyr::complete(sample, axisName)


# Color of the mutations.
colorMuts <- c('Deletion' = 'mediumblue',
               'Amplification' = 'red',
               'Inframe insertion' = '#CD7213', # Bronze
               'Inframe deletion' = '#02755c',
               'Frameshift variant' = 'purple',
               'Missense variant' = '#000000', # Black
               'Multiple coding mutations' = 'yellow3', # Yellow
               'Splice region variant' ='#F9A603', # Orange
               'Stop gained' = 'seagreen3', # Green
               'Splice acceptor variant' = 'yellow',
               'Start lost' = '#F0EB53', # Yellow
               'Splice donor variant' = 'thistle2', #Grey
               'Structural variant' = 'red' # red

)

dataOncoplot$isMutant <- ifelse(!is.na(dataOncoplot$ENSEMBL), 1, 0)

#count number of times a gene appears and sort from high to low
nGenesCounts <- dplyr::distinct(dataOncoplot, SYMBOL, .keep_all = T) %>% dplyr::filter(!(is.na(SYMBOL)))
nGenesCounts <- nGenesCounts[order(nGenesCounts$n),]

# Sort on mutually exclusiveness. APOBEC samples
memoData <- reshape2::dcast(dataOncoplot, axisName~sample, value.var = 'isMutant')
rownames(memoData) <- memoData$axisName; memoData$axisName <- NULL
memoData[is.na(memoData)] <- 0
memoData <- memoSort(memoData)

# Order axis
dataOncoplot$axisName <- factor(dataOncoplot$axisName, levels = rev(rownames(memoData)))

# Order samples
dataOncoplot$sample <- factor(dataOncoplot$sample, levels = orderSample_sigImbalance)


# Order by p-value Responder vs Non-resopnder

# prepare data for Fisher exact test
countPosNegImabalnce <- dplyr::select(signatureScoreMarkers_globalSigScore, variable, sigImbalanceCategory) %>%
  dplyr::add_count(sigImbalanceCategory, name = "nSizeSigImbCategory")
dataOncoplot_test <- dataOncoplot %>%
  dplyr::full_join(countPosNegImabalnce, by = c("sample" = "variable")) %>%
  dplyr::group_by(sigImbalanceCategory, axisName) %>%
  dplyr::mutate(nMutsPerImbCat = sum(isMutant)) %>%
  dplyr::ungroup(sigImbalanceCategory, axisName) %>%
  dplyr::distinct(sigImbalanceCategory, axisName, .keep_all = TRUE)
dataOncoplot_test_Positive <- dplyr::filter(dataOncoplot_test, sigImbalanceCategory == "Positive") %>% dplyr::select(axisName, groupSize1 = nSizeSigImbCategory, nGenesCount1 = nMutsPerImbCat)
dataOncoplot_test_Negative <- dplyr::filter(dataOncoplot_test, sigImbalanceCategory == "Negative") %>% dplyr::select(axisName, groupSize2 = nSizeSigImbCategory, nGenesCount2 = nMutsPerImbCat)
dataOncoplot_test <- dplyr::full_join(dataOncoplot_test_Positive, dataOncoplot_test_Negative, by = c("axisName"))


# Calculate Pvalues
dataOncoplot_test$pValue <- 1

for (iPatient in c(1:nrow(dataOncoplot_test))) {
  challenge.df = matrix(c(dataOncoplot_test$nGenesCount1[iPatient], dataOncoplot_test$nGenesCount2[iPatient],
                          dataOncoplot_test$groupSize1[iPatient]-dataOncoplot_test$nGenesCount1[iPatient], dataOncoplot_test$groupSize2[iPatient]-dataOncoplot_test$nGenesCount2[iPatient]), nrow = 2)
  
  dataOncoplot_test$pValue[iPatient] <- fisher.test(challenge.df)$p.value
}

# Adjust pValues for multiple testing
dataOncoplot_test <- dataOncoplot_test %>% dplyr::mutate(pAdj = p.adjust(pValue, method = "BH")) %>% dplyr::ungroup()

# order from low to high pValue
dataOncoplot_test <- dataOncoplot_test %>% dplyr::arrange(pValue)

# order data
dataOncoplot$axisName <- factor(dataOncoplot$axisName, levels = as.character(rev(dataOncoplot_test$axisName)))




#Plot oncoplot
plot_oncoResponderCohort <- ggplot(dataOncoplot, aes(x = sample, y = axisName)) +
  geom_tile(aes(fill = Consequence.CNA), size = 0.25, na.rm = T, alpha = 0.5, height = .9, width = .9) +
  geom_tile(aes(fill = Consequence.Mut), size = 0.25, na.rm = T, alpha = 1, height = .5, width = .5) +
  geom_tile(aes(fill = Consequence.SV), size = 0.25, na.rm = T, alpha = 1, height = .2, width = .9) +
  # Colors of mutations.
  scale_fill_manual(values = colorMuts, drop = F, na.translate = FALSE) +
  theme_minimal() +
  scale_x_discrete(drop = T) +
  scale_y_discrete(na.translate = FALSE) +
  labs(x = NULL, y = "Genes") +
  # Legend settings.
  guides(fill = guide_legend(title = 'Mutational Category',
                             title.position = 'top',
                             title.hjust = 0.5,
                             ncol = 1,
                             keywidth = 0.5,
                             keyheight = 0.5)) +
  # Change theme.
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 12),
    text=element_text(size=8, family='Helvetica'),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    #axis.text.x=element_text(angle = 45, hjust = 1),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_oncoResponderCohort <- cowplot::get_legend(
  plot_oncoResponderCohort + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_oncoResponderCohort <- plot_oncoResponderCohort + theme(legend.position = 'none')


# Plot P values for Oncoplot
dataOncoplot_test$axisName <- factor(dataOncoplot_test$axisName, levels = rev(dataOncoplot_test$axisName))
dataOncoplot_test <- dataOncoplot_test %>%
  dplyr::mutate(pValue_tmp = ifelse(pValue < 0.01, 0.01, pValue))

# plot p-value
plot_pValues_oncoplot <- ggplot(dataOncoplot_test, aes(x = "P-value", y = axisName, fill = -log10(pValue_tmp))) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = "white", mid = "#FFFFBF", high = '#9E0142', midpoint = 1,
                       breaks=c(0, 1, 2),labels=c("1", "0.1", "<0.01"),
                       limits = c(0, 2),
                       guide = guide_colorbar(title = 'P-value', reverse = TRUE,
                                              title.position = 'top', title.hjust = 0.5, barwidth = 3.5, barheight = 0.5)) +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        #axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5),
        axis.text.y = element_blank(),
        #plot.margin = unit(c(0, -1, -1, 0), "cm"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# extract legends
legend_plot_pValues_oncoplot <- cowplot::get_legend(
  plot_pValues_oncoplot + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_pValues_oncoplot <- plot_pValues_oncoplot + theme(legend.position = 'none')



# order data
DR176.MetaData$sampleId <- factor(DR176.MetaData$sampleId, levels = orderSample_sigImbalance)

# plot Responder labels
plot.Responders <- ggplot(DR176.MetaData, aes(sampleId, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = 'Response to treatment', title.position = 'top',
                                         title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.Responders <- cowplot::get_legend(
  plot.Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.Responders <- plot.Responders + theme(legend.position = 'none')



# order data
signatureScoreMarkers_globalSigScore$variable <- factor(signatureScoreMarkers_globalSigScore$variable, levels = orderSample_sigImbalance)

# plot signature imbalance category
plot_sigImbalaceTCellsStroma_category <- ggplot(signatureScoreMarkers_globalSigScore, aes(x = variable, y = "TSE category", fill = sigImbalanceCategory)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('Negative' = '#FF0000', 'Neutral' = '#FBE577', 'Positive' = '#006D2C'),
                    guide = guide_legend(title = 'TSE category',
                                         title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_sigImbalaceTCellsStroma_category <- cowplot::get_legend(
  plot_sigImbalaceTCellsStroma_category + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_sigImbalaceTCellsStroma_category <- plot_sigImbalaceTCellsStroma_category + theme(legend.position = 'none')





#++++++++++++++++++++++++++++++++++++ Plot Fusion genes
Fusions_data <- data.Cohort$fusionsLINX

# only keep samples with RNA
Fusions_data <- Fusions_data %>% dplyr::filter(sample %in% DR176.MetaData$sampleId)

# count frequency of genes appearing in fusion events
Fusions_data_summary <- rbind(dplyr::select(Fusions_data, sample, geneName = GeneStart),
                              dplyr::select(Fusions_data, sample, geneName = GeneEnd)) %>%
  dplyr::distinct() %>% dplyr::add_count(geneName, sort = TRUE)

# keep only fusions above 10% (n >= 4) and known fusion genes 
Fusions_data_summary <- Fusions_data_summary %>% dplyr::filter(n >= 4 | geneName %in% c('FGFR3')) %>%
  dplyr::filter(geneName %in% c('AHR', 'ARID1A', 'FGFR3'))

# number of samples
nSamples <- length(unique(DR176.MetaData$sampleId))

# Add axis names showing number of unique mutated samples.
Fusions_data_summary <- Fusions_data_summary %>%
  dplyr::mutate(axisName = sprintf('%s (%s%%)', geneName, round(100*n/nSamples)))
Fusions_data_summary <- Fusions_data_summary %>% tidyr::complete(sample, axisName)

# complete data set and order fusion genes
Fusions_data_summary$isMutant <- ifelse(!is.na(Fusions_data_summary$geneName), 'yes', 'no')
Fusions_data_summary <- Fusions_data_summary[order(Fusions_data_summary$n), ]

# Order samples and Axis
Fusions_data_summary$sample <- factor(Fusions_data_summary$sample, levels = orderSample_sigImbalance)
Fusions_data_summary$axisName <- factor(Fusions_data_summary$axisName, levels = unique(Fusions_data_summary$axisName))

#Plot oncoplot
plot_fusionGenes <- ggplot(Fusions_data_summary, aes(x = sample, y = axisName)) +
  geom_tile(aes(fill = isMutant), size = 0.25, na.rm = T, alpha = 1, height = .5, width = .5) +
  # Colors of mutations.
  scale_fill_manual(values = c("white", 'black'), drop = F, na.translate = FALSE) +
  theme_minimal() +
  scale_x_discrete(drop = T) +
  scale_y_discrete(na.translate = FALSE) +
  labs(x = NULL, y = "Fusion genes") +
  # Legend settings.
  guides(fill = guide_legend(title = 'Mutational Category',
                             title.position = 'top',
                             title.hjust = 0.5,
                             ncol = 1,
                             keywidth = 0.5,
                             keyheight = 0.5)) +
  # Change theme.
  theme(
    legend.position = 'none',
    axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5),
    text=element_text(size=10, family='Helvetica'),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))




#++++++++++++++++++++++++++++++++++++ Plot hotspot mutations
SNV_data_tmp <- data.Cohort$somaticVariants[data.Cohort$somaticVariants$mutType == "SNV"]
SNV_data <- data.frame(sample = SNV_data_tmp$sample,
                       chr = as.character(GenomeInfoDb::seqnames(SNV_data_tmp)),
                       pos = SNV_data_tmp@ranges@start,
                       TNC = SNV_data_tmp$TNC,
                       ref = as.character(SNV_data_tmp@ref),
                       alt = as.character(SNV_data_tmp@alt),
                       biotype = SNV_data_tmp@elementMetadata@listData[["ANN.BIOTYPE"]],
                       consequence = SNV_data_tmp@elementMetadata@listData[["ANN.Consequence"]],
                       impact = SNV_data_tmp@elementMetadata@listData[["ANN.IMPACT"]],
                       geneName = SNV_data_tmp@elementMetadata@listData[["ANN.SYMBOL"]],
                       HGVSp = SNV_data_tmp@elementMetadata@listData[["ANN.HGVSp"]]) %>%
  dplyr::mutate(ref = as.character(ref), alt = as.character(alt),
                snv_id = paste0(chr, ":", pos))

# keep only samples with RNA
SNV_data <- SNV_data %>% dplyr::filter(sample %in% DR176.MetaData$sampleId)
SNV_data <- droplevels(SNV_data)

# Keep only genes with hotspot mutations (from whole cohort) Find HotSpot SNVs
hotspotMutList <- c("chr21:18984827", #BTG3
                    "chr10:115511590", "chr10:115511593", #PLEKHS1
                    "chr5:1295228", "chr5:1295250", #TERT
                    "chr8:29952919", "chr8:29952921", #LEPROTL1
                    "chr6:142706206", "chr6:142706209", #ADGRG6
                    "chr10:96162368", "chr10:96162370", #"TBC1D12",
                    "chrX:151686628")
SNV_data_summary <- SNV_data %>% dplyr::filter(snv_id %in% hotspotMutList) %>%
  dplyr::group_by(snv_id) %>% dplyr::add_count(name = "totalHotspots", sort = TRUE) %>%
  dplyr::filter(totalHotspots>1) %>% dplyr::ungroup()


# keep only protein coding and intergenic variants
SNV_data_summary <- SNV_data_summary %>% dplyr::filter(biotype == 'protein_coding' | consequence == 'intergenic_variant')

# Change ID of hotspot mutations to include gene name
SNV_data_summary <- SNV_data_summary %>% dplyr::mutate(snv_id_aux = ifelse(is.na(HGVSp), snv_id, gsub(".*:", "", HGVSp))) %>%
  dplyr::mutate(snv_id_aux = ifelse(is.na(geneName), snv_id_aux, paste0(geneName, ", ", snv_id_aux)))

# Get % of hotspot muts in the cohort
SNV_data_summary <- SNV_data_summary %>%
  dplyr::mutate(axisName = paste0(snv_id_aux, " (", round(100 * totalHotspots/length(unique(DR176.MetaData$sampleId))), "%)"))

# order axis of hotspots
HotSpot_order <- rev(unique(SNV_data_summary$axisName))
SNV_data_summary$axisName <- factor(SNV_data_summary$axisName, levels = HotSpot_order)

# order samples
SNV_data_summary$sample <- factor(SNV_data_summary$sample, levels = orderSample_sigImbalance)


# Define consequence of hotspot
SNV_data_summary <- SNV_data_summary %>%
  dplyr::mutate(consequence = ifelse(consequence == "5_prime_UTR_variant", "5' UTR",
                                     ifelse(consequence == "intron_variant", "Intron variant",
                                            ifelse(consequence == "missense_variant", "Missense variant",
                                                   ifelse(consequence == "upstream_gene_variant", "Upstream gene variant",
                                                          ifelse(consequence == "downstream_gene_variant", "Downstream gene variant", "Intergenic variant"))))))
SNV_data_summary$consequence <- factor(SNV_data_summary$consequence, levels = c("Upstream gene variant", "Downstream gene variant", "5' UTR", "Intron variant", "Missense variant", "Intergenic variant"))

# Colors for Consequence of hotpot
# Color of the mutations.
colorHotSpot <- c(#'Missense variant' = '#000000', # Black
                  "Upstream gene variant" = "#4E96CB", # blue
                  "Downstream gene variant" = "#D4774C", # red
                  "Intron variant" = "#D9CE86", # yellow
                  "5' UTR"  = "#8FC8A3", # green
                  "Intergenic variant" = "#C8DCDE") # grey


plot_HotSpot_Cohort <- ggplot(SNV_data_summary, aes(x = sample, y = axisName, fill = consequence)) +
  geom_tile(size = 0.25, na.rm = F, alpha = 1, height = .5, width = .5) +
  # Legend settings.
  scale_fill_manual(name = "Hotspot consequence", values = colorHotSpot, drop = T, na.translate = F) +
  guides(fill = guide_legend(title = "Hotspot consequence",
                             title.position = 'top', title.hjust = 0.5, ncol = 1,
                             keywidth = 0.5, keyheight = 0.5)) +
  annotationTheme() +
  scale_x_discrete(drop = F) +
  scale_y_discrete(drop = T) +
  labs(x = NULL, y = "Hotspot mutations") +
  # Change theme.
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_HotSpot_Cohort <- cowplot::get_legend(
  plot_HotSpot_Cohort + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_HotSpot_Cohort <- plot_HotSpot_Cohort + theme(legend.position = 'none')





# export oncoplot
samplePlot <- cowplot::plot_grid(
  plot_sigImbalaceTCellsStroma_category,
  plot.Responders,

  
  plot_oncoResponderCohort,
  plot_fusionGenes,
  plot_HotSpot_Cohort,
  ncol = 1, align = 'v',
  axis = 'tblr',
  rel_heights = c(0.03, 0.03, 1, 0.08, 0.25)
)


# p values
samplePlot_pValue <- cowplot::plot_grid(NULL,
                                              NULL,
                                              plot_pValues_oncoplot,
                                              NULL,
                                              NULL,
                                              ncol = 1,
                                              align = 'v',
                                              axis = 'tblr',
                                              rel_heights = c(0.03, 0.03, 1, 0.07, 0.2)
)


legend_samplePlot <- cowplot::plot_grid(legend_plot_sigImbalaceTCellsStroma_category,
                     legend_plot.Responders,

                     legend_plot_oncoResponderCohort,
                     legend_plot_pValues_oncoplot,
                     legend_plot_HotSpot_Cohort,
                     ncol = 1, align = 'h', rel_heights = c(1, 0.5, 1.5, 0.5, 0.5))


# export Figure S5
pdf(paste0(odir,"oncoplotSigImbalance.pdf"),width = 9, height = 9)#, width = 14, height = 21)
cowplot::plot_grid(samplePlot, samplePlot_pValue, legend_samplePlot,
                   ncol = 3, rel_widths = c(0.8, 0.02, 0.2), align = 'h', axis = 'tblr')
dev.off()






#================================= Plot Figure S3c - Pathway analysis + DDR genes ============================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# load gene pathways
pathwayGenes <- xlsx::read.xlsx(file = "/databases/GeneAndPathways_TCGA.xlsx", sheetIndex = 2)
DDRgenes <- readxl::read_xlsx("/databases/DDRgenes.xlsx")

# One gene changed name from the list
pathwayGenes <- pathwayGenes %>% dplyr::mutate(Gene = ifelse(Gene == "NOV", "CCN3", as.character(Gene)))

# add DDR genes as extra pathway
DDRgenes <- DDRgenes %>% dplyr::mutate(Pathway = "DDR") %>% dplyr::select(Pathway, geneName)
pathwayGenes <- pathwayGenes %>% dplyr::select(Pathway, geneName = Gene) %>%
  rbind(DDRgenes) %>% dplyr::group_by(Pathway) %>% dplyr::add_count(name = 'sizePathway')

#Get all Genes for pathway analysis from combined report
dataMutPathway <- results.Cohort$combinedReport %>% dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::filter(SYMBOL %in% pathwayGenes$geneName)
  

# check which genes did not appear, maybe they have another name or no mutations
subset(pathwayGenes, !(pathwayGenes$geneName %in% dataMutPathway$SYMBOL))

# Only keep coding mutations or deep CNA.
dataMutPathway <- dataMutPathway %>% dplyr::filter(!is.na(Consequence.Mut) | grepl('Deep', Consequence.CNA))

# Replace synonymous variants, neutral and shallow CN for NA to not count them
dataMutPathway <- dataMutPathway %>% dplyr::mutate(Consequence.CNA = ifelse(Consequence.CNA %in% c('Neutral', 'Deletion', 'Amplification'), NA, Consequence.CNA))

# rename deep amplification/deletion to just amplification and deletion
dataMutPathway <- dataMutPathway %>% dplyr::mutate(Consequence.CNA = gsub("Deep ", "", Consequence.CNA),
                                                   Consequence.SV = gsub("Va", "va", Consequence.SV))

# number of samples
nSamples <- length(unique(DR176.MetaData$sampleId))

# count number of mutations per Pathway
dataMutPathway <- dataMutPathway %>% dplyr::left_join(pathwayGenes, by = c('SYMBOL' = 'geneName')) %>%
  dplyr::group_by(Pathway, sample) %>% dplyr::count() %>%
  dplyr::mutate(Pathway = factor(Pathway, levels = unique(pathwayGenes$Pathway)),
                sample = factor(sample, levels = unique(DR176.MetaData$sampleId))) %>%
  tidyr::complete(Pathway, sample) %>% dplyr::distinct() %>%
  dplyr::left_join(dplyr::distinct(pathwayGenes, Pathway, sizePathway), by = c('Pathway')) %>%
  dplyr::mutate(propMutations = ifelse(is.na(n), 0, n/sizePathway),
                nMutations = ifelse(is.na(n), 0, n),
                isMutant = ifelse(is.na(n), 0, 1))

# include a threshold
dataMutPathway <- dataMutPathway %>% dplyr::mutate(propMutations = ifelse(propMutations > 0.2, 0.21, propMutations))
dataMutPathway <- dataMutPathway %>% dplyr::mutate(sample = factor(sample, levels = orderSample_sigImbalance))

# order pathways
orderPathways <- unique((dataMutPathway %>% dplyr::group_by(Pathway) %>% dplyr::mutate(nMutantSamples = sum(isMutant)) %>%
                           dplyr::ungroup() %>% dplyr::arrange(-nMutantSamples))$Pathway)
dataMutPathway <- dataMutPathway %>% dplyr::mutate(Pathway = factor(Pathway, levels = rev(orderPathways)))

# order again (after looking at the first result)
dataMutPathway <- dataMutPathway %>% dplyr::mutate(Pathway = factor(Pathway, levels = rev(c("p53", "Cell cycle", "RTK-RAS", "Notch", "Hippo", "WNT", "PI3K", "DDR", "JAK-STAT", "Myc", "Nrf2", "TGF-beta"))))


#Plot oncoplot
plot_pathwayAlterations <- ggplot(dataMutPathway, aes(x = sample, y = Pathway)) +
  geom_tile(aes(fill = propMutations), size = 0.25, na.rm = T, alpha = 1, height = .9, width = .9) +
  scale_fill_gradient(low = 'white', high = 'red',
                      breaks=c(0, 0.1, 0.2),labels=c("0.0", 0.1, ">0.2"),
                      limits = c(0, 0.2), na.value = "#bd0a04",
                      guide = guide_colorbar(title = 'Proportion of\npathway mutated',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 5, barheight = 0.5))  + 
  theme_minimal() +
  scale_x_discrete(drop = F) +
  scale_y_discrete(na.translate = FALSE) +
  labs(x = NULL, y = "Pathway") +
  # Change theme.
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    #axis.text.x=element_text(angle = 45, hjust = 1),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_pathwayAlterations <- cowplot::get_legend(
  plot_pathwayAlterations + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_pathwayAlterations <- plot_pathwayAlterations + theme(legend.position = 'none')



############# Plot pathway activities #########################
# from Nakauma Gonzalez et al, 2022
TGFb_genes <- c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
                'RFLNB', "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A",
                "TAGLN", "TGFBI", "TNS1", "TPM1")
CellCycle_genes <- c("MKI67", "CCNE1", "BUB1", "BUB1B", "CCNB2", "CDC25C", "CDK2", "MCM4", "MCM6", "MCM2")
WNT_genes <- c("EFNB3", "MYC", "TCF12", "VEGFA")
Notch_genes <- c("HES1", "HES5", "HEY1")
PI3K_genes <- c("AGRP", "BCL2L11", "BCL6", "BNIP3", "BTG1", "CAT", "CAV1", "CCND1", "CCND2", "CCNG2",
                "CDKN1A", "CDKN1B", "ESR1", "FASLG", "FBXO32", "GADD45A", "INSR", "MXI1", "NOS3", "PCK1",
                "POMC", "PPARGC1A", "PRDX3", "RBL2", "SOD2", "TNFSF10")
Hippo_genes <- c("TAZ", "YAP1")
p53_genes <- c("CDKN1A", "RRM2B", "GDF15", "SUSD6", "BTG2", "DDB2", "GADD45A", "PLK3", "TIGAR", "RPS27L",
               "TNFRSF10B", "TRIAP1", "ZMAT3", "BAX", "BLOC1S2", "PGF", "POLH", "PPM1D", "PSTPIP2", "SULF2",
               "XPC")
NRF2_genes <- c("GCLM", "NQO1", "PHGDH", "PSAT1", "SHMT2")
MYC_genes <- c("TFAP4", "BMP7", "CCNB1", "CCND2", "CCNE1", "CDC25A", "CDK4", "CDT1", "E2F1",
               "GATA4", "HMGA1", "HSP90AA1", "JAG2", "CDCA7", "LDHA", "MCL1", "NDUFAF2", "MTA1", "MYCT1",
               "NPM1", "ODC1", "SPP1", "PIN1", "PTMA", "PRDX3", "PRMT5", "DNPH1", "TFRC", "EMP1",
               "PMEL", "C1QBP") # "H19" is not included because it is lncRNA
RTK_RAS_genes <- c("SPRY2", "SPRY4", "ETV4", "ETV5", "DUSP4", "DUSP6", "CCND1", "EPHA2", "EPHA4")
JAK_STAT <- c("IRGM", "ISG15", # All cell types
              "GATA3", "FCER2", "THY1", "NFIL3", # T/B cell
              "ARG1", "RETNLB", "CLEC7A", "CHIA", # macrophages
              "OSM", "BCL2L1", "CISH", "PIM1", # hematopoietic linaege
              "SOCS2", "GRB10") # Mammary epithelia

listPahwayTargetGenes <- data.frame(gene = c(TGFb_genes, CellCycle_genes, WNT_genes, Notch_genes, PI3K_genes,
                                             Hippo_genes, p53_genes, NRF2_genes, MYC_genes, RTK_RAS_genes, JAK_STAT),
                                    pathway = c(rep("TGF-beta", length(TGFb_genes)), rep("Cell cycle", length(CellCycle_genes)),
                                                rep("WNT", length(WNT_genes)), rep("Notch", length(Notch_genes)), 
                                                rep("PI3K", length(PI3K_genes)), rep("Hippo", length(Hippo_genes)),
                                                rep("p53", length(p53_genes)), rep("Nrf2", length(NRF2_genes)),
                                                rep("Myc", length(MYC_genes)), rep("RTK-RAS", length(RTK_RAS_genes)),
                                                rep("JAK-STAT", length(JAK_STAT)) ))


length(rownames(normalizedCountMatrix[rownames(normalizedCountMatrix) %in% listPahwayTargetGenes$gene, ]))
# Get pathway genes and center around the mean
signatureScore_Pathways <- normalizedCountMatrix[rownames(normalizedCountMatrix) %in% listPahwayTargetGenes$gene, ]
signatureScore_Pathways  <- signatureScore_Pathways - rowMedians(signatureScore_Pathways)

# Reshape data
signatureScore_Pathways <- as.data.frame(signatureScore_Pathways)
signatureScore_Pathways$gene <- rownames(signatureScore_Pathways)
rownames(signatureScore_Pathways) <- NULL
signatureScore_Pathways <- reshape2::melt(signatureScore_Pathways,
                                          id.vars = "gene", variable.name = "sample", value.name = "normExp")

# Add pathway ID to genes
signatureScore_Pathways <- signatureScore_Pathways %>%
  dplyr::left_join(listPahwayTargetGenes, by = "gene")

# Get score (mean expression per patient per pathway)
signatureScore_Pathways <- signatureScore_Pathways %>%
  dplyr::group_by(sample, pathway) %>% dplyr::mutate(scorePathway = mean(normExp)) %>%
  dplyr::ungroup() %>% dplyr::distinct(sample, pathway, .keep_all = TRUE)

# Add extra information such as if it is mutated or not
signatureScore_Pathways <- signatureScore_Pathways %>%
  dplyr::left_join(dplyr::select(dataMutPathway, sample, Pathway, isMutant), by = c("sample", "pathway" = "Pathway"))

#  Perform Kruskal-Wallis test
pValue_KruskalTest <- signatureScore_Pathways %>% dplyr::group_by(pathway) %>%
  dplyr::summarize(p.value = wilcox.test(scorePathway ~ isMutant)$p.value) %>% dplyr::ungroup()
pValue_KruskalTest$pAdj <- p.adjust(pValue_KruskalTest$p.value, method = "BH")
pValue_KruskalTest <- pValue_KruskalTest %>% #dplyr::mutate(pAdj_label = format(pAdj, scientific = TRUE, digits = 3)) %>%
  dplyr::mutate(pAdj_label = ifelse(pAdj <= 0.01, format(round(pAdj, 3)), round(pAdj, 2)))

plot_PathwayRNAScores <- ggplot(signatureScore_Pathways, aes(x = as.factor(isMutant), y = scorePathway, fill = as.factor(isMutant))) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey60") +
  geom_boxplot(notch = F, alpha = .66, outlier.shape = NA, width = .5) +
  geom_point(aes(x = as.factor(isMutant), fill = as.factor(isMutant)), alpha = 0.33, cex = .75,
             shape = 21, position = position_jitterdodge()) +
  facet_wrap(~pathway, nrow = 2) +
  ylim(c(-1.5, 1.5)) +
  ylab('Pathway activity (RNA)') +
  
  scale_fill_manual(na.translate=FALSE,
                    label = c("Non-mutant", "Mutant"),
                    guide = guide_legend(title = 'Pathway', title.position = 'top',
                                         title.hjust = 0.5, ncol = 1, keywidth = 1, keyheight = 1),
                    name = NULL, values = c("white", "red")) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=pValue_KruskalTest, aes(x = 1.5, y = 1.3,
                                         label = pAdj_label), size = 2.5, inherit.aes = FALSE) +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 8.2),
    text=element_text(size=8, family='Helvetica'),
    panel.spacing=unit(0, "lines"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    #axis.text.x = element_text(size = 9, angle=30, hjust=1, vjust=1.0),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )



# export Figure S3c
pdf(paste0(odir,"PathwayScore_MutNonMuts.pdf"), width = 6, height = 1.8)
plot_PathwayRNAScores
dev.off()




#======================================= plot ROC curves with TMB and imbalance ====================================
#=================================================================================================
# plot ROC curve for responders on TMB-high 
data_CohortROC <- results.Cohort$mutationalBurden %>% 
  dplyr::left_join(data.Cohort$APOBEC_enrich, by = "sample") %>%
  dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::left_join(dplyr::select(DR176.MetaData, sampleId, `PFS days`, `PFS status`, `OS days`, `Survival status`), by = c('sample' = 'sampleId')) %>%
  dplyr::full_join(signatureScoreMarkers_globalSigScore, by = c('sample' = 'variable')) %>%
  dplyr::mutate(Outcome6M = Outcome6M_tmp) %>%
  dplyr::mutate(OutcomeResponse = ifelse(Outcome6M == "Responder", 1, 0),
                isAPOBEChigh = ifelse(APOBEC_mutagenesis == "High", 1, 0),
                isTMBhigh = ifelse(tmbStatus == "High TMB (≥10)", 1, 0),
                isPositiveImbalance = ifelse(signatureImbalance > 0, 1, 0)) %>%
  dplyr::mutate(isWeakImbalanceAndHighTMB = ifelse(sigImbalanceCategory == 'Neutral' & isTMBhigh == 1, 'Positive',
                                                   ifelse(sigImbalanceCategory == 'Neutral' & isTMBhigh != 1, 'Negative', as.character(sigImbalanceCategory)))) %>%
  dplyr::mutate(isWeakImbalanceAndHighAPOBEC = ifelse(sigImbalanceCategory == 'Neutral' & isAPOBEChigh == 1, 'Positive',
                                                   ifelse(sigImbalanceCategory == 'Neutral' & isAPOBEChigh != 1, 'Negative', as.character(sigImbalanceCategory)))) %>%
  dplyr::mutate(isWeakImbalanceAndHighTMBHighAPOBEC = ifelse(sigImbalanceCategory == 'Neutral' & isAPOBEChigh == 1 & isTMBhigh == 1, 'Positive',
                                                      ifelse(sigImbalanceCategory == 'Neutral' & !(isAPOBEChigh == 1 & isTMBhigh == 1), 'Negative', as.character(sigImbalanceCategory))))


# sort data from low TMB to high
data_CohortROC <- data_CohortROC %>%
  dplyr::arrange(signatureImbalance)

# calculate glm fit for imbalance T cell - Stroma
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance, family = binomial)
ROC_sigImbalance <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# Multivariate logistic regression
summary(stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$globalSigScore_Tcells+data_CohortROC$globalSigScore_stroma, family = binomial))

# calculate glm fit for imbalance T cell - Stroma and TMB
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance+data_CohortROC$Genome.TMB, family = binomial)
ROC_sigImbalance_TMB <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)


# calculate glm fit for imbalance T cell - Stroma and APOBEC mutagenesis
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance+data_CohortROC$foldEnrichment, family = binomial)
ROC_sigImbalance_APOBEC <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# calculate glm fit for imbalance T cell - Stroma, TMB and APOBEC mutagenesis
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance+data_CohortROC$Genome.TMB+data_CohortROC$foldEnrichment, family = binomial)
ROC_sigImbalance_TMBAPOBEC <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)


# test sig between two ROC curves
pValueAUC_TMB <- round((pROC::roc.test(ROC_sigImbalance, ROC_sigImbalance_TMB))$p.value, 2)
pValueAUC_APOBEC <- round((pROC::roc.test(ROC_sigImbalance, ROC_sigImbalance_APOBEC))$p.value, 2)
pValueAUC_TMBAPOBEC <- round((pROC::roc.test(ROC_sigImbalance, ROC_sigImbalance_TMBAPOBEC))$p.value, 2)


# Multiple curves:
plot_ROC_tmb_APOBECmut <- pROC::ggroc(list(ROC_TMB,
                                           ROC_APOBECmut,
                                           ROC_sigImbalance,
                                           ROC_sigImbalance_TMB,
                                           ROC_sigImbalance_APOBEC,
                                           ROC_sigImbalance_TMBAPOBEC),
                                      legacy.axes=TRUE) +
  labs(y = 'True positive rate', x = 'False positive rate') +
  #geom_line(size=1) +
  scale_color_manual(values = c(alpha("#86a38b", 0.75),
                                alpha("#ff94fd", 0.75),
                                alpha("#1762ad", 0.75),
                                alpha("#f0a630", 0.75),
                                alpha("#a7cef2", 0.75),
                                alpha('#008837', 0.75)),
                     labels = c('1' = 'TMB',
                                '2' = 'APOBEC',
                                '3' = 'TSE',
                                '4' = 'TSE + TMB',
                                '5' = 'TSE + APOBEC',
                                '6' = 'TSE + TMB + APOBEC'),
                     guide = guide_legend(title = 'ROC curve',
                                          title.position = 'top',
                                          title.hjust = 0, nrow = 2,
                                          keywidth = 0.5, keyheight = 0.5,
                                          override.aes = list(size=3))) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", size = 0.2) +
  annotate(geom="text", x=0.5, y=0.45, size=2,
           label = paste0("AUC = ", format(round(ROC_TMB$auc, 2))), color="#86a38b", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.375, size=2,
           label = paste0("AUC = ", format(round(ROC_APOBECmut$auc, 2))), color="#ff94fd", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.3, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance$auc, 2))), color="#1762ad", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.225, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance_TMB$auc, 2))), color="#f0a630", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.15, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance_APOBEC$auc, 2))), color="#a7cef2", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.075, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance_TMBAPOBEC$auc, 2))), color="#008837", hjust = 0) +
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 9),
    text=element_text(size=9, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )



# calculate glm fit for imbalance Immune cells - Stroma
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance_ImmuneCellsStroma, family = binomial)
ROC_sigImbalance_ImmuneCells <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# calculate glm fit for imbalance Immune cells - Stroma
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$signatureImbalance_ALLimmuneCellsStroma, family = binomial)
ROC_sigImbalance_ALLImmuneCells <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# test sig between two ROC curves
pValueAUC_ImmuneCells <- round((pROC::roc.test(ROC_sigImbalance, ROC_sigImbalance_ImmuneCells))$p.value, 2)
pValueAUC_ALLImmuneCells <- round((pROC::roc.test(ROC_sigImbalance, ROC_sigImbalance_ALLImmuneCells))$p.value, 2)


# Multiple curves:
plot_ROC_tmb_APOBECmut_stroma <- pROC::ggroc(list(ROC_sigImbalance,
                                           ROC_sigImbalance_ImmuneCells,
                                           ROC_sigImbalance_ALLImmuneCells),
                                      legacy.axes=TRUE) +
  labs(y = 'True positive rate', x = 'False positive rate') +
  scale_color_manual(values = c(alpha("#1762ad", 0.75), alpha("#f0a630", 0.75), alpha("#008837", 0.75), alpha('gray60', 0.75)),
                     labels = c('1' = 'STS',
                                '2' = 'Stroma - Immune cells (non T cells)',
                                '3' = 'Stroma - All immune cells'),
                     guide = guide_legend(title = 'ROC curve',
                                          title.position = 'top',
                                          title.hjust = 0, ncol = 1,
                                          keywidth = 0.5, keyheight = 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", size = 0.2) +
  annotate(geom="text", x=0.5, y=0.3, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance$auc, 2))), color="#1762ad", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.225, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance_ImmuneCells$auc, 2)),
                          ", p = ", pValueAUC_TMB), color="#f0a630", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.15, size=2,
           label = paste0("AUC = ", format(round(ROC_sigImbalance_ALLImmuneCells$auc, 2)),
                          ", p = ", pValueAUC_APOBEC), color="#008837", hjust = 0) +
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 9),
    text=element_text(size=9, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )



# plot Response rate
data_CohortROC_responseRate <- dplyr::count(data_CohortROC, Outcome6M, imbalanceCategory = sigImbalanceCategory) %>%
  tidyr::complete(Outcome6M, imbalanceCategory) %>%
  dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
  dplyr::group_by(imbalanceCategory) %>% dplyr::mutate(sizePerGroup = sum(n)) %>% dplyr::ungroup() %>%
  dplyr::mutate(responderRate = round(100*n/sizePerGroup)) %>%
  dplyr::mutate(imbalanceCategory = factor(imbalanceCategory, levels = c('Positive', 'Neutral', 'Negative')))


# calculate fisher exact test
nResponderGroup_Positive <- data_CohortROC_responseRate %>% dplyr::filter(imbalanceCategory == "Positive" &
                                                                    Outcome6M == "Responder")
nResponderGroup_Negative <- data_CohortROC_responseRate %>% dplyr::filter(imbalanceCategory == "Negative" &
                                                                    Outcome6M == "Responder")
nResponderGroup_Neutral <- data_CohortROC_responseRate %>% dplyr::filter(imbalanceCategory == "Neutral" &
                                                                            Outcome6M == "Responder")

challenge.df = matrix(c(nResponderGroup_Positive$n, nResponderGroup_Negative$n,
                        nResponderGroup_Positive$sizePerGroup - nResponderGroup_Positive$n, nResponderGroup_Negative$sizePerGroup - nResponderGroup_Negative$n), nrow = 2)

pFisherExacTest <- fisher.test(challenge.df)$p.value

data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(Outcome6M = factor(Outcome6M, levels = c("Non-responder", "Responder")))

# change TSE status to include number of samples in each group
data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(imbalanceCategory = paste0(imbalanceCategory, ' (n=', sizePerGroup, ')')) %>%
  dplyr::mutate(imbalanceCategory = factor(imbalanceCategory, levels = c('Positive (n=15)', 'Neutral (n=14)', 'Negative (n=12)')))


# plot histograms of response rates with 95% CI from binomial distribution
plot_ResponseRate_imbalanceSig <- ggplot(data_CohortROC_responseRate, aes(x=imbalanceCategory, y=responderRate, fill = Outcome6M)) +
  geom_bar(stat='identity', color = 'black', width = 0.66) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%'))  +
  labs(title = paste0("p = ", format(round(pFisherExacTest, 5))),
       y = "Response (relative proportion)", x = 'TSE category') +
  annotate(geom="text", x=1, y=nResponderGroup_Positive$responderRate+4, size=2,
           label = paste0(nResponderGroup_Positive$responderRate, "%"), color="black") +
  annotate(geom="text", x=2, y=nResponderGroup_Neutral$responderRate+4, size=2,
           label = paste0(nResponderGroup_Neutral$responderRate, "%"), color="black") +
  annotate(geom="text", x=3, y=nResponderGroup_Negative$responderRate+4, size=2,
           label = paste0(nResponderGroup_Negative$responderRate, "%"), color="black") +
  scale_fill_manual(guide = guide_legend(title = 'Response to treatment', nrow = 1,
                                         title.position = 'top', title.hjust = 0.5, keywidth = 0.5, keyheight = 0.5),
                    name = NULL, values = c('Responder' = '#06A506', 'Non-responder' = '#E2C054')) +
  theme(
    legend.position = 'bottom',
    text=element_text(size=9, family='Helvetica'),
    plot.title = element_text(size=9, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )



# save ROC and Response rate Figure 5a
pdf(paste0(odir,"ROC_ResponseRate_SigImbalance.pdf"),width = 5, height = 3)#, width = 14, height = 21)
cowplot::plot_grid(plot_ResponseRate_imbalanceSig, plot_ROC_tmb_APOBECmut_stroma,
                   nrow = 1,
                   align = 'vh', axis = 'tblr')
dev.off()





# Plot OS curve for Imbalance signatures
data_CohortROC_survival <- data_CohortROC %>% dplyr::select(PFS_days = `PFS days`, PFS_status = `PFS status`, OS_status = `Survival status`, OS_days = `OS days`, imbalanceCategory = sigImbalanceCategory, response = OutcomeResponse, Outcome6M, TMB = Genome.TMB, APOBECmut = foldEnrichment, signatureImbalance_TcellsStroma, globalSigScore_Tcells, globalSigScore_stroma) %>%
  dplyr::mutate(response = 1 - response,
                OS_days = 12*OS_days/365,
                PFS_days = 12*PFS_days/365) %>%
  dplyr::mutate(imbalanceCategory = factor(imbalanceCategory, levels = c("Positive", "Neutral", "Negative")))

# get survival model
fit_OS_TSE <- survminer::surv_fit(survival::Surv(OS_days,  OS_status) ~ imbalanceCategory, data = data_CohortROC_survival)
fit_PFS_TSE <- survminer::surv_fit(survival::Surv(PFS_days,  PFS_status) ~ imbalanceCategory, data = data_CohortROC_survival)

# Multivariate cox regression analysis
coxRegressionAnalysis <- survival::coxph(survival::Surv(OS_days, OS_status) ~ signatureImbalance_TcellsStroma + TMB + APOBECmut, data_CohortROC_survival)
summary(coxRegressionAnalysis)

# Linear multivariate regression
multivariateReg <- lm(response ~ globalSigScore_Tcells + globalSigScore_stroma, data = data_CohortROC_survival)
summary(multivariateReg)

# change factors to calculate HR
data_CohortROC_survival$imbalanceCategory <- factor(data_CohortROC_survival$imbalanceCategory, levels = c("Negative", "Neutral", "Positive"))


# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(OS_days, OS_status) ~ imbalanceCategory, data =  dplyr::filter(data_CohortROC_survival, imbalanceCategory != "Neutral"))
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[2,1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[2,3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[2,4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))


# survival plot
plot_OSCurve <- survminer::ggsurvplot(fit_OS_TSE, data = data_CohortROC_survival,
                                       size = 1,                 # change line size
                                      palette =  c('#006D2C', '#FBE577', '#FF0000'),
                                       xlab = "Months",
                                      ylab = "OS probability (n = 41)",
                                      xlim = c(0,37),
                                       title = res.cox_HazardRatio,
                                       legend.title = "TSE",
                                      legend = c(0.7,0.6),
                                       pval.method.size = 5,
                                       break.time.by = 5,
                                       pval = TRUE,              # Add p-value (default log-rank)
                                       pval.size = 3,
                                       pval.coord = c(10, 0.9),
                                       risk.table = TRUE,        # Add risk table
                                      tables.y.text.col = TRUE, 
                                      risk.table.fontsize = 2,
                                      surv.median.line = "hv",
                                       risk.table.col = "strata",# Risk table color by groups
                                       legend.labs = c( "Positive", "Neutral", "Negative"),    # Change legend labels
                                       #ggtheme = theme_bw()
                                       ggtheme = theme(
                                         legend.position = 'bottom',
                                         text=element_text(size=10, family='Helvetica'),
                                         plot.title = element_text(size=8, hjust = 0),
                                         legend.title=element_text(size=8, hjust = 0.5), 
                                         legend.text=element_text(size=8),
                                         legend.key = element_rect(color = NA, fill = NA),
                                         legend.key.size = unit(0.25, "cm"),
                                         legend.spacing.x = unit(0.1, "line"),
                                         legend.spacing.y = unit(0, "line"),
                                         legend.background	= element_rect(fill = NA, colour = NA),
                                         panel.grid.major.x = element_blank(),
                                         panel.grid.minor.x = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         panel.background = element_rect(fill = 'white', colour = NA),
                                         panel.border = element_rect(fill = NA, colour = 'grey20'))
)      # Change ggplot2 theme)


# change theme of tables
plot_OSCurve$table <- plot_OSCurve$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )




# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(PFS_days, PFS_status) ~ imbalanceCategory, data =  dplyr::filter(data_CohortROC_survival, imbalanceCategory != "Neutral"))
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[2,1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[2,3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[2,4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))


plot_PFSCurve <- survminer::ggsurvplot(fit_PFS_TSE, data = data_CohortROC_survival,
                      size = 1,                 # change line size
                      palette =  c('#006D2C', '#FBE577', '#FF0000'),
                      xlab = "Months",
                      ylab = "PFS probability (n = 41)",
                      title = res.cox_HazardRatio,
                      legend.title = "TSE",
                      legend = c(0.7,0.6),
                      pval.method.size = 5,
                      break.time.by = 5,
                      pval = TRUE,              # Add p-value (default log-rank)
                      pval.size = 3,
                      pval.coord = c(10, 0.9),
                      risk.table = TRUE,        # Add risk table
                      tables.y.text.col = TRUE, 
                      risk.table.fontsize = 2,
                      surv.median.line = "hv",
                      risk.table.col = "strata",# Risk table color by groups
                      legend.labs = c("Positive", "Neutral", "Negative"),    # Change legend labels
                      risk.table.height = 0.25, # Useful to change when you have multiple groups
                      
                      #ggtheme = theme_bw()
                      ggtheme = theme(
                        legend.position = 'bottom',
                        #axis.title.y = element_text(size = 8),
                        text=element_text(size=10, family='Helvetica'),
                        plot.title = element_text(size=8, hjust = 0),
                        legend.title=element_text(size=8, hjust = 0.5), 
                        legend.text=element_text(size=8),
                        legend.key = element_rect(color = NA, fill = NA),
                        legend.key.size = unit(0.25, "cm"),
                        legend.spacing.x = unit(0.1, "line"),
                        legend.spacing.y = unit(0, "line"),
                        legend.background	= element_rect(fill = NA, colour = NA),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = 'white', colour = NA),
                        panel.border = element_rect(fill = NA, colour = 'grey20'))
                      )      # Change ggplot2 theme)


# change theme of tables
plot_PFSCurve$table <- plot_PFSCurve$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )



# save OS and PFS Figure 5b
pdf(paste0(odir,"Survival_SigImbalance.pdf"),width = 9, height = 3)#, width = 14, height = 21)

cowplot::plot_grid(
  
  cowplot::plot_grid(
    plot_OSCurve$plot,
    plot_OSCurve$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
  
  cowplot::plot_grid(
    plot_PFSCurve$plot,
    plot_PFSCurve$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
    
  align = 'v', axis = 'tblr', nrow = 1)


dev.off()





