# Author:  J. Alberto Nakauma-Gonzalez
# Date:   28-06-2022
# Function: This script pre-processes data for the DR176/CPCT02 mUC cohort
# using data pre-processed from HMF scripts and using the R2CPCT R package (https://github.com/J0bbie/R2CPCT).

# Libraries ---------------------------------------------------------------
pacman::p_load('plyr', 'dplyr', 'tidyr', 'ggplot2', 'BiocParallel', 'GenomicRanges', 'IRanges')
library(R2CPCT)
library(ShatterSeek)
library(graph)
library(CHORD)
library(ShatterSeek)
library(mutSigExtractor)

# defines pathways
path.output <- '/HMF/DR176/postHMF/RData/'
path.metaData <- '/HMF/DR176/metadata/'
pathHMFcombined <- '/HMF/DR176/postHMF/combinedData'

# load meta data
load(paste0(path.output, 'DR176.MetaData.RData'))

# Select same patients as previous release.
cpctIds <- as.character(DR176.MetaData$sampleId)


# Retrieve all relevant HMF files -----------------------------------------
R2CPCT::retrieveFilesWGS(
  pathHMF = '/DR176/somatics/',
  pathOutput = pathHMFcombined,
  cpctIds = cpctIds
)


# Perform VEP (Annotation) ------------------------------------------------
R2CPCT::performVEP(
  pathVCF = pathHMFcombined,
  cpctIds = cpctIds
)

# Run LINX (v1.11) to detect fusions and drivers ----------------------------------

R2CPCT::performLINX(
  pathHMFcombined,
  cpctIds = cpctIds,
  nThreads = 40,
  dryRun = T
)

# Import / convert WGS samples --------------------------------------------

# Import the WGS data of all samples in the cohort
data.Cohort <- R2CPCT::importWGSOfCohort(cpctIds = cpctIds, inputFolder = pathHMFcombined, nThreads = 20, performAggregation = T)
save(data.Cohort, file = file.path(path.output, 'data.Cohort.RData'))


# Perform GISTIC2 (v2.0.23) analysis ------------------------------------------------

# Generate the command to perform GISTIC2 (perform this in your Bash terminal).
R2CPCT::performGISTIC2(data.Cohort$copyNumbers, outputFolder = '/HMF/DR176/postHMF/GISTIC/')


# Cohort-wide analysis ----------------------------------------------------

# Initialize a list to contain the results of the WGS analysis.
results.Cohort <- list()

# Generate overview of mutational burden per sample and cohort.
results.Cohort$mutationalBurden <- R2CPCT::determineMutationalBurden(data.Cohort, minTAF.Muts = 0, minTAF.SV = 0.1)

# Incorporate the sample information.
results.Cohort$mutationalBurden <- results.Cohort$mutationalBurden %>% dplyr::full_join(DR176.MetaData %>% dplyr::distinct(sample = sampleId, hmfPatientId, hmfSampleId, tumorPurity, RESP_subject, gender, Outcome6M, cohortGroup), by = "sample")


# Perform dN/dS on cohort.
results.Cohort$dNdS <- R2CPCT::rundNdS(data.Cohort$somaticVariants)

# Import the GISTIC2-determined recurrent CNA peaks.
results.Cohort$GISTIC2 <- R2CPCT::importGISTIC2('/HMF/DR176/postHMF/GISTIC/')

# Generate the gene-level mutational overview.
results.Cohort$combinedReport <- R2CPCT::generateCombinedReport(data.Cohort, dNdS = results.Cohort$dNdS, GISTIC2 = results.Cohort$GISTIC2, nThreads = 12, mutantsOnly = T)

# Incorporate sample information.
results.Cohort$combinedReport <- results.Cohort$combinedReport %>% dplyr::inner_join(DR176.MetaData %>% dplyr::distinct(sample = sampleId, hmfSampleId, RESP_subject, Outcome6M), by = "sample")

# Determine mutational signatures.
results.Cohort$mutSigs <- R2CPCT::fitMutSigs(data.Cohort$somaticVariants, restrictiveFit = F)



#========================= Determine APOBEC enrichment ===================================

########### Add APOBEC mut_load Enrichment calculations ###########
########### if using this method, please also consider citing Nakauma-Gonzalez, et al. 2022: https://doi.org/10.1016/j.eururo.2022.01.026 ######
require("BSgenome.Hsapiens.UCSC.hg19")

#get A, T, C, G content
params <- new("BSParams", X = Hsapiens, FUN = alphabetFrequency, exclude = c("M","Y","random","hap", "chrUn"))
ATCG_content <- bsapply(params)

# transform to data frame
ATCG_content <-  as.data.frame(do.call(rbind, ATCG_content))

# Calculate TCW motifs
params <- new("BSParams", X = Hsapiens, FUN = countPDict, exclude = c("M","Y","random","hap", "chrUn"))

# Get TCW (WGA) content
pdict_TCA <- PDict(DNAStringSet("TCA"))
TCW_content_TCA <- bsapply(params, pdict = pdict_TCA)

pdict_TCT <- PDict(DNAStringSet("TCT"))
TCW_content_TCT <- bsapply(params, pdict = pdict_TCT)

pdict_AGA <- PDict(DNAStringSet("AGA"))
TCW_content_AGA <- bsapply(params, pdict = pdict_AGA)

pdict_TGA <- PDict(DNAStringSet("TGA"))
TCW_content_TGA <- bsapply(params, pdict = pdict_TGA)

# Create data frame with all TCW (WGA) DNA context
TCW_content <-  as.data.frame(do.call(rbind, TCW_content_TCA))
TCW_content <-  cbind(TCW_content, as.data.frame(do.call(rbind, TCW_content_TCT)))
TCW_content <-  cbind(TCW_content, as.data.frame(do.call(rbind, TCW_content_AGA)))
TCW_content <-  cbind(TCW_content, as.data.frame(do.call(rbind, TCW_content_TGA)))
colnames(TCW_content) <- c("V1", "V2", "V3", "V4")

# Get total number of TCW motifs and C (G) bp
total_Context_CorG <- sum(ATCG_content$C + ATCG_content$G)
total_Context_TCW <- sum(TCW_content$V1 + TCW_content$V2 + TCW_content$V3 + TCW_content$V4)

# get total mutations per patient and per tri-nucleotide context and define APOBEC TCW context mutations
SNV_data_totalMuts <- as.data.frame(results.Cohort[["mutSigs"]][["SNV"]][["mutMatrix"]])
SNV_data_totalMuts$mutType96 <- row.names(SNV_data_totalMuts)
SNV_data_totalMuts <- reshape2::melt(SNV_data_totalMuts, id.vars = "mutType96", value.name = "nMuts96", variable.name = "sample")
SNV_data_APOBEC <- SNV_data_totalMuts %>% dplyr::mutate(mutFromTo = gsub("\\].*", "", mutType96)) %>%
  dplyr::mutate(mutFromTo = gsub("*.\\[", "", mutFromTo),
                ref = gsub(">.*", "", mutFromTo),
                APOBECmuts = ifelse(mutType96 %in% c("T[C>T]A", "T[C>T]T", "T[C>G]A", "T[C>G]T"), nMuts96, 0),
                CtoT = ifelse(mutFromTo == "C>T", nMuts96, 0),
                CtoG = ifelse(mutFromTo == "C>G", nMuts96, 0)) %>%
  dplyr::group_by(sample) %>% dplyr::mutate(totalSNVs = sum(nMuts96 ),
                                            totalAPOBEC = sum(APOBECmuts),
                                            total_CtoT = sum(CtoT),
                                            total_CtoG = sum(CtoG),
                                            total_CtoT_CtoG = sum(CtoT) + sum(CtoG)) %>% dplyr::ungroup() %>%
  dplyr::distinct(sample, totalSNVs, totalAPOBEC, total_CtoT_CtoG)


# Calculate fold enrichment for APOBEC mutations
SNV_data_APOBEC <- SNV_data_APOBEC %>%
  dplyr::mutate(foldEnrichment = (totalAPOBEC * total_Context_CorG) / (total_CtoT_CtoG * total_Context_TCW))

# Initialize value for pValues calculation
SNV_data_APOBEC$E_pValue <- 1

# Calculate enrichment pValues (one-sided Fisher's exact test)
for (iPatient in c(1:nrow(SNV_data_APOBEC))) {
  challenge.df = matrix(c(SNV_data_APOBEC$totalAPOBEC[iPatient], total_Context_TCW,
                          SNV_data_APOBEC$total_CtoT_CtoG[iPatient]-SNV_data_APOBEC$totalAPOBEC[iPatient], total_Context_CorG - total_Context_TCW), nrow = 2)
  
  SNV_data_APOBEC$E_pValue[iPatient] <- fisher.test(challenge.df, alternative = "greater")$p.value
}

#Adjust p-values with Benjamini-Hochberg method
SNV_data_APOBEC$E_pAdj <- p.adjust(SNV_data_APOBEC$E_pValue, method = "BH")

# Tumor is APOBEC enriched if p_Adj < 0.01
APOBEC_enrich <- SNV_data_APOBEC %>% dplyr::mutate(APOBEC_enrich = ifelse(E_pAdj < 0.01, "Yes", "No"))
APOBEC_enrich$APOBEC_enrich <- factor(APOBEC_enrich$APOBEC_enrich, levels = c("No", "Yes"))

# classify tumors in low, medium, high or no APOBEC mutagenesis based on last result and fold enrichment
APOBEC_enrich <- APOBEC_enrich %>%
  dplyr::mutate(APOBEC_mutagenesis = ifelse(APOBEC_enrich == "Yes",
                                            ifelse(foldEnrichment > 3, "High",
                                                   ifelse(foldEnrichment > 2, "Medium", "Low")), "No"))
APOBEC_enrich$APOBEC_mutagenesis <- factor(APOBEC_enrich$APOBEC_mutagenesis, levels = c("No", "Low", "Medium", "High"))


# save data
data.Cohort$APOBEC_enrich <- APOBEC_enrich



# --------------------------------# --------------------------------# --------------------------------
# --------------------------------# --------------------------------# --------------------------------

# Save results.
save(results.Cohort, file = file.path(path.output, 'results.Cohort.RData'))

# --------------------------------# --------------------------------# --------------------------------
# --------------------------------# --------------------------------# --------------------------------






