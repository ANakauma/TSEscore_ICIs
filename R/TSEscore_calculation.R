# Author:  J. Alberto Nakauma Gonzalez
# Date:   09-05-2022
# Function: This script calculates the TSE score for the DR-176 cohort. It can be adjusted for other cohorts
# e-mail: j.nakaumagonzalez@erasmusmc.nl
# 

# Clean everything before starting and set working directory
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load normalized RNA counts (DESeq2) with the vst method
load("data/normalizedCountMatrix.RData")




#---------------------------------------------------------------------------------------------------
#----------------------------------Calculate TSE score----------------------------------------------
#---------------------------------------------------------------------------------------------------

# Load list signatures and their genes
gene_signatures <- readxl::read_xlsx("GeneSignaturesList.xlsx", sheet = 1)

# Indicate signatures for global T-cell and stromal signatures (see paper for details)  
selectedSignaturesForGlobalSig <- c("CAF", "EMT/stroma core genes", "Fibroblasts", "Stromal signature", "TBRS", 
                                 "IFN gamma", "tGE8", "T cell signature", "Immune gene signature", "T cell inflamed GEP", "Chemoattractants", "Cytotoxic CD8 T cell")


# Keep only the selected signatures and genes that are in the RNA-seq data set
gene_signatures <- gene_signatures %>%
  dplyr::filter(signatureID %in% selectedSignaturesForGlobalSig) %>%
  dplyr::filter(Gene %in% rownames(normalizedCountMatrix))


# Get expression of all genes
signatureScores <- normalizedCountMatrix[rownames(normalizedCountMatrix) %in% unique(gene_signatures$Gene), ]

# Center gene expression around the median across samples
signatureScores  <- signatureScores - rowMedians(signatureScores)

# Add column with gene names
signatureScores <- as.data.frame(signatureScores) %>% dplyr::mutate(geneName = rownames(signatureScores))

# Melt results
signatureScores <- reshape2::melt(signatureScores, id.vars = "geneName")

# Add marker type information
signatureScores <- gene_signatures %>% dplyr::right_join(signatureScores, by = c("Gene" = "geneName"))

# get signature score per gene signature
signatureScores <- signatureScores %>%
  dplyr::group_by(signatureID, variable) %>%
  dplyr::mutate(sigScore = mean(value)) %>% dplyr::ungroup() %>%
  dplyr::distinct(signatureID, variable, .keep_all = TRUE) %>%
  dplyr::mutate(Gene = NULL)


# Calculate global signature scores for T cells and stromal resident cells
signatureScores <- signatureScores %>%
  dplyr::group_by(MainCategory, variable) %>%
  dplyr::mutate(globalSigScore_mainCategories = mean(sigScore)) %>% dplyr::ungroup()

# Re-arrange data to get in columns the global signature scores for T cells and stromal resident cells
signatureScores_globalSigScore_Tcells <-  signatureScores %>% dplyr::filter(MainCategory == 'T cells') %>%
  dplyr::select(variable, globalSigScore_Tcells = globalSigScore_mainCategories) %>%
  dplyr::distinct(variable, globalSigScore_Tcells)
signatureScores_globalSigScore_stroma <-  signatureScores %>% dplyr::filter(MainCategory == 'Stromal resident cells') %>%
  dplyr::select(variable, globalSigScore_stroma = globalSigScore_mainCategories) %>%
  dplyr::distinct(variable, globalSigScore_stroma)


# Calculate TSE score and define categories
cutOff_TSEscore <- 0.5
TSEscore_results <- signatureScores_globalSigScore_Tcells %>%
  dplyr::full_join(signatureScores_globalSigScore_stroma, by = "variable") %>%
  dplyr::mutate(TSE_score = globalSigScore_Tcells - globalSigScore_stroma) %>%
  dplyr::mutate(TSE_score_category = ifelse(TSE_score > cutOff_TSEscore, 'Positive',
                                              ifelse(TSE_score < -cutOff_TSEscore, 'Negative', 'Neutral')))
  





#---------------------------------------------------------------------------------------------------
#---------------------------Identify  mUC transcriptomic subtypes-----------------------------------
#---------------------------------------------------------------------------------------------------

# Load genes associated with each transcriptomic subtype (see Nakauma, et. al: https://doi.org/10.1016/j.eururo.2022.01.026)
mUC_subtypesGenes <- readxl::read_xlsx("/Genes_mUC_RNAsubtypes.xlsx")

# Get expression of genes of interest
normMatrixRNA_mUC_subtype <- normalizedCountMatrix[rownames(normalizedCountMatrix) %in% unique(mUC_subtypesGenes$geneName), ]

# Normalize data by centering around the median each gene value across samples
normMatrixRNA_mUC_subtype  <- normMatrixRNA_mUC_subtype - rowMedians(normMatrixRNA_mUC_subtype)

# Add column with gene names
normMatrixRNA_mUC_subtype <- as.data.frame(normMatrixRNA_mUC_subtype) %>% dplyr::mutate(geneName = rownames(normMatrixRNA_mUC_subtype))

# Melt results
normMatrixRNA_mUC_subtype <- reshape2::melt(normMatrixRNA_mUC_subtype, id.vars = "geneName")

# Add subtype label of genes
normMatrixRNA_mUC_subtype <- normMatrixRNA_mUC_subtype %>%
  dplyr::left_join(dplyr::select(mUC_subtypesGenes, geneName, geneRNAsubtype_mUC = RNAcluster), by = "geneName")

# Calculate score per transcriptomic subtype
normMatrixRNA_mUC_subtype <- normMatrixRNA_mUC_subtype %>% dplyr::group_by(geneRNAsubtype_mUC, variable) %>%
  dplyr::mutate(meanExp = mean(value)) %>% dplyr::ungroup() %>%
  dplyr::distinct(variable, geneRNAsubtype_mUC, meanExp) %>%
  dplyr::select(sample = variable, geneRNAsubtype_mUC, subtypeScore = meanExp)

# Get score per transcriptomic mUC subtype
normMatrixRNA_mUC_subtype_LumA <- normMatrixRNA_mUC_subtype %>% 
  dplyr::filter(geneRNAsubtype_mUC == "Luminal-a") %>% dplyr::select(sample, LumA = subtypeScore)
normMatrixRNA_mUC_subtype_LumB <- normMatrixRNA_mUC_subtype %>%
  dplyr::filter(geneRNAsubtype_mUC == "Luminal-b") %>% dplyr::select(sample, LumB = subtypeScore)
normMatrixRNA_mUC_subtype_BaSq <-  normMatrixRNA_mUC_subtype %>%
  dplyr::filter(geneRNAsubtype_mUC == "Basal/Squamous") %>% dplyr::select(sample, BaSq = subtypeScore)
normMatrixRNA_mUC_subtype_Stroma <-  normMatrixRNA_mUC_subtype %>%
  dplyr::filter(geneRNAsubtype_mUC == "Stroma-rich") %>% dplyr::select(sample, Stroma = subtypeScore)
normMatrixRNA_mUC_subtype_NonSpe <- normMatrixRNA_mUC_subtype %>%
  dplyr::filter(geneRNAsubtype_mUC == "Non-specified") %>% dplyr::select(sample, NonSpe = subtypeScore)

# get all scores per sample
normMatrixRNA_mUC_subtypeSummary <- normMatrixRNA_mUC_subtype_LumA %>% dplyr::full_join(normMatrixRNA_mUC_subtype_LumB, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_mUC_subtype_BaSq, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_mUC_subtype_Stroma, by = "sample") %>%
  dplyr::full_join(normMatrixRNA_mUC_subtype_NonSpe, by = "sample")

# Based on the signature scores for each mUC subtype, we can now define a transcriptomic subtype per sample 
results_RNAmUCsubtype <- normMatrixRNA_mUC_subtype %>% dplyr::arrange(-subtypeScore) %>%
  dplyr::distinct(sample, geneRNAsubtype_mUC) %>% dplyr::distinct(sample, .keep_all = TRUE) %>%
  dplyr::full_join(normMatrixRNA_mUC_subtypeSummary, by = "sample")





