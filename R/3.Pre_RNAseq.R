# Author:  JA Nakauma Gonzalez
# Date:   28-06-2022
# Function: This script performs DGE and immune cell fraction from RNAseq for the DR-176 cohort
# e-mail: j.nakaumagonzalez@erasmusmc.nl

# Clean everything before starting a new project and set working directory (default R script path)---------------------------
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(DESeq2)
library(ConsensusClusterPlus)
library(immunedeconv)


# path to data and output dir
pathHMF <- "/DR176/postHMF/"
odir <- "/output/figures/"

# Import sample information
load(paste0(pathHMF, "RData/DR176.MetaData.RData"))
load(paste0(pathHMF, "RData/data.Cohort.RData"))
load(paste0(pathHMF, "RData/results.Cohort.RData"))



# Perform DESeq2 Responder vs Non-responder----------------------------------------------------------

# Import annotations
genesDB <- rtracklayer::import.gff('/annotation/hg19/GENCODE/gencode.v35lift37.annotation.gtf')

genesDB <- tibble::as_tibble(mcols(genesDB))
geneInfo <- genesDB %>% dplyr::distinct(gene_id, gene_name)

# Import counts.
counts <- readr::read_delim(paste0(pathHMF, 'RNAseq/counts/DR176_RNA.counts'), delim= '\t', comment = '#')
colnames(counts) <- gsub('_.*', '', gsub('.*/', '', colnames(counts)))

#Add geneNames
counts <- merge(counts, dplyr::mutate(geneInfo, Geneid = gene_id, gene_id = NULL), by = "Geneid")
counts$Geneid <- counts$gene_name
counts$gene_name <- NULL

# Convert to matrix.
countMatrix <- as.matrix(counts[7:ncol(counts)], rownames = 1)
rownames(countMatrix) <- counts$Geneid

# Only keep included samples.
countMatrix <- countMatrix[,colnames(countMatrix) %in% DR176.MetaData$sampleId]
DR176.MetaData <- DR176.MetaData[DR176.MetaData$sampleId %in% colnames(countMatrix), ]


####### Keep only Protein coding (mRNA)
RNA_typesToKeep <- genesDB %>% dplyr::filter(gene_type %in% c("protein_coding")) %>%
  dplyr::distinct(gene_id, gene_name)
countMatrix <- countMatrix[rownames(countMatrix) %in% RNA_typesToKeep$gene_name ,]

# prepare data for normalization
colData = DR176.MetaData[match(colnames(countMatrix), DR176.MetaData$sampleId),]
colData$Outcome6M <- factor(colData$Outcome6M, levels = c("non-responder", "responder"))

# DESeq2 with Wald Test between Outcome6M
DESeq2.DR176 <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix, colData = colData, design = ~Outcome6M)
DESeq2.DR176.WT <- DESeq2::DESeq(DESeq2.DR176, parallel = T)

# Normalize with DESeq2 with Variance Stabilizing Transformation
DESeq2.DR176.WT.vst <- DESeq2::vst(DESeq2.DR176.WT[rowSums(DESeq2::counts(DESeq2.DR176.WT, normalized = T)) > 10,], blind = F)


# Diff gene expression
DESeq2.DR176.WT.Diff <- tibble::as_tibble(DESeq2::results(DESeq2.DR176.WT, pAdjustMethod = 'BH', tidy = T))
DESeq2.DR176.WT.Diff <- DESeq2.DR176.WT.Diff %>% dplyr::filter(padj <= 0.05, baseMean >= 100, abs(log2FoldChange) > 1)

# Save RData frame generated
save(DESeq2.DR176.WT, file = paste0(pathHMF, "RData/DESeq2.DR176.WT.RespNonResp.RData"))





#============================== Immune cell fraction calculation =======================================================
##======================================================================================================================
#Load files
files.RSEM <- list.files(paste0(pathHMF, 'RNAseq/counts/RSEM_TPM/'), pattern = 'genes.results', full.names = T)

data_TPM <- tibble::as_tibble(do.call(cbind, lapply(files.RSEM, function(x){ y <- readr::read_delim(x, delim = '\t'); return(y$TPM)})))
colnames(data_TPM) <- gsub('_Al.*', '', basename(files.RSEM))

# Convert to matrix.
data_TPM <- as.matrix(data_TPM)

# Add annotation.
data.RSEM <- readr::read_delim(files.RSEM[[1]], delim = '\t')[,1]
data.RSEM <- data.RSEM %>% dplyr::inner_join(geneInfo, by = c('gene_id' = 'gene_id'))

# Combine.
rownames(data_TPM) <- data.RSEM$gene_name

# Only keep included samples.
data_TPM <- data_TPM[, colnames(data_TPM) %in% DR176.MetaData$sampleId]

# Only keep mRNA
data_TPM <- data_TPM[rownames(data_TPM) %in% RNA_typesToKeep$gene_name, ]

# Perform Quantiseq
immuneCFResults <- immunedeconv::deconvolute(data_TPM, method = 'quantiseq', tumor = T, scale_mrna = T)

# save results
save(immuneCFResults, file = paste0(pathHMF, "RData/immuneCFResults.RData"))






