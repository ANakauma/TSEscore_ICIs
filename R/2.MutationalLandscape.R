# Author:  J. Alberto Nakauma-Gonzalez
# Date:   28-06-2022
# Function: This script generates the mutational landscape figure for the DR176/HMF mUC cohort manuscript.

# Clean everything before starting a new project and set working directory (default R script path)---------------------------
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load libraries ----------------------------------------------------------
pacman::p_load('plyr', 'dplyr', 'tidyr', 'ggplot2', 'BiocParallel', 'GenomicRanges', 'IRanges', 'NMF', 'ggpubr',
               'ConsensusClusterPlus', 'R2CPCT', 'RColorBrewer', 'xlsx')


# path to data and output dir
path.hmf <- "/HMF/DR176/postHMF/RData/"
odir <- "/output/figures/"

# load pre-processed data
load(paste0(path.hmf,"DR176.MetaData.RData"))
load(paste0(path.hmf,"data.Cohort.RData"))
load(paste0(path.hmf,"results.Cohort.RData"))

# define genome size
genomeSize <- 2858.67 #(Mb)

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



##################################################################################################
#---------------------- Find genomic subtypes ----------------------

dataMutSigsSBS <- tibble::as_tibble(reshape2::melt(results.Cohort$mutSigs$SNV$relativeContribution, varnames = c('proposedAetiologyGrouped', 'sampleId')))
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::filter(sampleId %in% DR176.MetaData$sampleId)

# Get  relative contributions and combine
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::group_by(sampleId, proposedAetiologyGrouped) %>% dplyr::summarise(value = sum(value)) %>% dplyr::ungroup()

# Find clusters based on Mut Signatures 
mutSigsMatrix <- dataMutSigsSBS %>% reshape2::dcast(proposedAetiologyGrouped~sampleId, value.var = 'value')
rownames(mutSigsMatrix) <- mutSigsMatrix$proposedAetiologyGrouped
mutSigsMatrix$proposedAetiologyGrouped <- NULL
mutSigsMatrix <- as.matrix(mutSigsMatrix)
mutSigsMatrix[is.na(mutSigsMatrix)] <- 0


#+++++++++ Calculate Clusters
clusterDir <- paste0(odir, "GenS_SBScluster")
ConsensusMatrix <- ConsensusClusterPlus::ConsensusClusterPlus(mutSigsMatrix, pItem=0.8,
                                                              maxK = 10,reps = 1000, title = clusterDir, plot = "pdf", distance = "pearson")
icl = ConsensusClusterPlus::calcICL(ConsensusMatrix, title=clusterDir, plot="pdf")

#Choose num of clusters and Get Consensus Matrix
resultsConsensusClass <- as.data.frame(ConsensusMatrix[[6]]$consensusClass)
colnames(resultsConsensusClass)[1] <- "Cluster"
resultsConsensusClass$sample <- rownames(resultsConsensusClass)

#Order samples based on clusters size
resultsConsensusClass <- resultsConsensusClass %>%  dplyr::add_count(Cluster, name = "nSize") %>%
  dplyr::arrange(-nSize)

# name subtypes
# in the validation cohort cluster 1-3 are the same but 4 and 5 are different
resultsConsensusClass <- resultsConsensusClass %>%
  dplyr::mutate(genomicSubtype = ifelse(Cluster == 1, "GenS1",
                                     ifelse(Cluster == 2, "GenS2",
                                            ifelse(Cluster == 3, "GenS3",
                                                   ifelse(Cluster == 5, "GenS4",
                                                          ifelse(Cluster == 4, "GenS5","GenS6"))))))
                                                   #ifelse(Cluster == 4, "GenS6", "GenS7")))))
resultsConsensusClass$genomicSubtype <- factor(resultsConsensusClass$genomicSubtype,
                                               levels = c("GenS1", "GenS2", "GenS3", "GenS4", "GenS5", "GenS6"))
                                               #levels = c("GenS1", "GenS2", "GenS3", "GenS6"))

resultsConsensusClass$Cluster <- NULL

# Save MSig DBS Clusters
resultsGenomicSubtype <- resultsConsensusClass
save(resultsGenomicSubtype, file = paste0(path.hmf, "resultsGenomicSubtype.RData"))

#++++++++++++++++++++++++++++++ Finished Calculate Clusters ++++++++++++++++++++++++++++++++++++++
##################################################################################################


# get TMB
results.Cohort_mutationalBurden <- results.Cohort$mutationalBurden

# Sort by Responders and non-responders and decreasing sample TMB.
orderSamples <- (results.Cohort_mutationalBurden %>% dplyr::arrange(-Genome.TMB) %>% 
  dplyr::arrange(factor(Outcome6M, levels = c("responder", "non-responder"))) %>%
    dplyr::distinct(sample))$sample


# order samples for plot
results.Cohort_mutationalBurden <- results.Cohort_mutationalBurden %>% dplyr::mutate(sample = factor(sample, levels = orderSamples),
                                                                                     Outcome6M = factor(Outcome6M, levels = c("responder", "non-responder")))
DR176.MetaData$sampleId <- factor(DR176.MetaData$sampleId, levels = orderSamples) 

# plot Responder labels
plot.Responders <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('Response to treatment',
                    values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.Responders <- cowplot::get_legend(
  plot.Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.Responders <- plot.Responders + theme(legend.position = 'none')


# Melt the TMB
mutsPerMb <- reshape2::melt(results.Cohort_mutationalBurden %>% dplyr::select(sample, SNV = Muts.SNV, Indel = Muts.InDel, MNV = Muts.MNV, SV = totalSV))
mutsPerMb <- mutsPerMb %>% dplyr::mutate(value = value/genomeSize,
                                        variable = factor(variable, levels = c('SNV', 'Indel', 'MNV', 'SV')))

# order samples
mutsPerMb$sample <- factor(mutsPerMb$sample, levels = orderSamples)

# Barplot of total Mutations per Sample
plot_mutsPerMb <- ggplot(mutsPerMb, aes(x = sample, y = value, fill = variable)) +
  coord_trans(y = 'sqrt') +
  scale_y_continuous(breaks = c(0, 5, 10, 25, 50, 75, 100), expand = expand_scale(0,0)) +
  scale_x_discrete(expand = expand_scale(0,0)) +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2) +
  scale_fill_manual('Mutational category', values = c('SNV' = '#375D96', 'Indel' = '#FFBA00', 'MNV' = 'salmon', 'SV' = 'mediumseagreen'),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  labs(x = NULL, y = 'Mutations per\nmega base pair\n(Genome-wide)') +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 8),
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
legend_plot_mutsPerMb <- cowplot::get_legend(
  plot_mutsPerMb + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_mutsPerMb <- plot_mutsPerMb + theme(legend.position = 'none')


# order data and define TMB high (>10) TMB low (<10)
results.Cohort_mutationalBurden <- results.Cohort_mutationalBurden %>%
  dplyr::mutate(tmbStatus = ifelse(Genome.TMB >= 10, "High TMB (≥10)", "Low TMB (<10)")) %>%
  dplyr::mutate(tmbStatus = factor (tmbStatus, levels = c("High TMB (≥10)", "Low TMB (<10)")))

# plot TMB > 10
plot.TMB <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'TMB', fill = tmbStatus)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('#8E1B37', '#7A9BCF'),
                    guide = guide_legend(title = "Tumor mutational burden (TMB)", title.position = 'top', title.hjust = 0,
                                         nrow = 1, keywidth = 0.5, keyheight = 0.5)) + 
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.TMB <- cowplot::get_legend(
  plot.TMB + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.TMB <- plot.TMB + theme(legend.position = 'none')


#---------- Plot mutation load on protein missense mutation ------------------------------------
################################################################################################
# keep only missense mutations
somaticVariants <- data.frame(sample = data.Cohort$somaticVariants$sample,
                               overlappingGenes = data.Cohort$somaticVariants$ANN.SYMBOL,
                               chr = GenomeInfoDb::seqnames(data.Cohort$somaticVariants),
                               start = IRanges::start(data.Cohort$somaticVariants),
                               end = IRanges::end(data.Cohort$somaticVariants),
                               ref = data.Cohort$somaticVariants@ref,
                               alt = data.Cohort$somaticVariants@alt,
                               VAF = data.Cohort$somaticVariants$PURPLE_AF,
                               VCN = data.Cohort$somaticVariants$PURPLE_VCN,
                              IMPACT = data.Cohort$somaticVariants$ANN.IMPACT,
                              HGVSp = data.Cohort$somaticVariants$ANN.HGVSp,
                              mutType = data.Cohort$somaticVariants$ANN.VARIANT_CLASS,
                              aaChange = data.Cohort$somaticVariants$ANN.Amino_acids) %>%
  dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::filter(!(chr %in% c("chrM", "chrY"))) %>%
  dplyr::filter(!is.na(HGVSp)) %>% dplyr::filter(IMPACT != "LOW") %>%
  dplyr::filter(mutType == "SNV")


# count number of missense mutations
somaticVariants <- somaticVariants %>% dplyr::add_count(sample, name = "totalMissenseMuts")

# order samples
somaticVariants$sample <- factor(somaticVariants$sample, levels = orderSamples)

# plot total missense mutations
plot_missenseMuts <- ggplot(dplyr::distinct(somaticVariants, sample, .keep_all = TRUE),
                          aes(x = sample, y = totalMissenseMuts, fill = totalMissenseMuts)) +
  scale_y_continuous(expand = expand_scale(0,0)) +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2) +
  scale_fill_gradient(low = 'white', high = 'black') +
  labs(x = NULL, y = '# missense\nmutations') +
  theme(
    legend.position = 'none', axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'black', linetype = 'longdash'),
    panel.grid.minor.y = element_line(colour = 'black', linetype = 'dashed'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# define colors for Genomic subtypes
colDNAcluster <- c("GenS1" = '#ff6767', "GenS2" = '#b3e6ff', "GenS3" = 'darkseagreen4', "GenS4" = '#ffe596', "GenS5" = 'grey70')

resultsGenomicSubtype$sample <- factor(resultsGenomicSubtype$sample, levels = orderSamples)
color_GenS <- c("GenS1" = '#ff6767',
                "GenS2" = '#b3e6ff',
                "GenS3" = 'darkseagreen4',
                'GenS4' = '#ffe596',
                "GenS5" = 'grey70',
                "GenS6" = '#f4b0ff',
                'GenS7' = 'steelblue')

# Plot Matrix clusters
plot.GenomicSubtype <- ggplot(resultsGenomicSubtype, aes(sample, y = "Genomic subtype", fill = genomicSubtype)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = color_GenS,
                    guide = guide_legend(title = "Genomic subtype", title.position = 'top', title.hjust = 0,
                                         nrow = 2, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.GenomicSubtype <- cowplot::get_legend(
  plot.GenomicSubtype + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# delete legend from main plot
plot.GenomicSubtype <- plot.GenomicSubtype + theme(legend.position = 'none')



# define clonal when variantPloidy > 0.75
somaticVariants <- somaticVariants %>% dplyr::mutate(isclonal = ifelse(VCN >= 0.75, "clonal", "subclonal")) %>%
  dplyr::add_count(sample, name = "totalMutsSample") %>%
  dplyr::add_count(sample, isclonal, name = "numClonalMuts") %>%
  dplyr::mutate(relativeClonal = numClonalMuts/totalMutsSample) %>%
  dplyr::distinct(sample, isclonal, .keep_all = TRUE)


# order samples
somaticVariants <- somaticVariants %>% dplyr::inner_join(DR176.MetaData, by = c("sample" = "sampleId"))
somaticVariants$sample <- factor(somaticVariants$sample, levels = orderSamples)



# plot clonality results
plot_clonalFraction <- ggplot(dplyr::filter(somaticVariants, isclonal == "clonal"), aes(x = sample, y = "Clonal fraction", fill = relativeClonal)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = 'white', mid = "#f7d5e5", high = "#910449", midpoint = 0.6, limits = c(0.2,1),
                       guide = guide_colorbar(title = 'Clonal fraction',
                                              title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 3.5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))



# extract legends
legend_plot_clonalFraction <- cowplot::get_legend(
  plot_clonalFraction + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_clonalFraction <- plot_clonalFraction + theme(legend.position = 'none')


#========================================= plot mut signatures ==============

# get mut sig data
dataMutSigsSBS <- tibble::as_tibble(reshape2::melt(results.Cohort$mutSigs$SNV$relativeContribution, varnames = c('proposedAetiologyGrouped', 'sampleId')))
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::filter(sampleId %in% DR176.MetaData$sampleId)
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::mutate(sampleId  = factor(as.character(sampleId), levels = orderSamples))

# Get  relative contributions by Etiology group
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::group_by(sampleId, proposedAetiologyGrouped) %>% dplyr::summarise(value = sum(value)) %>% dplyr::ungroup()

# Define shorter names for Etiology groups
proposedAetiology_SBS <- R2CPCT::proposedAetiologyCOSMICv3.1 %>% dplyr::filter(grepl("SBS",Signature)) %>%
  dplyr::distinct(proposedAetiologyGrouped) %>%
  dplyr::mutate(proposedAetiologyGrouped_shorterName = c("Deamination of 5-methylcytosine (SBS1)",
                                                         "APOBEC activity (SBS2, SBS13)",
                                                         "Defective HR DNA damage repair (SBS3)",
                                                         "Tobacco smoking (SBS4)",
                                                         "Unknown (clock-like signature) (SBS5)",
                                                         "Defective DNA MMR (SBS6, SBS15, SBS21, SBS26, SBS44)",
                                                         "Ultraviolet light exposure (SBS7a-d)",
                                                         "Unknown (15 signatures)",
                                                         "Polimerase eta somatic hypermutation activity (SBS9)",
                                                         "POLE mutations (SBS10a-b)",
                                                         "Temozolomide treatment (SBS11)",
                                                         "Concurrent POLE mutations and defective DNA MMR (SBS14)",
                                                         "Damage by reactive oxygen species (SBS18)",
                                                         "Concurrent POLD1 mutations and defective DNA MMR (SBS20)",
                                                         "Aristolochic acid exposure (SBS22)",
                                                         "Aflatoxin exposure (SBS24)",
                                                         "Chemotherapy treatment (SBS25)",
                                                         "Tobacco chewing (SBS29)",
                                                         "Defective DNA BER due to NTHL1 mutations (SBS30)",
                                                         "Platinum chemotherapy treatment (SBS31, SBS35)",
                                                         "Azathioprine treatment (SBS32)",
                                                         "Defective DNA BER due to MUTYH mutations (SBS36)",
                                                         "Indirect effect of ultraviolet light (SBS38)",
                                                         "Haloalkane exposure (SBS42)",
                                                         "AID activity (SBS84)",
                                                         "Indirect effects of AID activity (SBS85)",
                                                         "Unknown chemotherapy treatment (SBS86)",
                                                         "Thiopurine chemotherapy treatment (SBS87)",
                                                         "Colibactin exposure (SBS88)",
                                                         "Duocarmycin exposure (SBS90)",
                                                         "Possible Sequencing Artifact (SBS27, SBS43, SBS45-60)"))

# assign short names
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::left_join(proposedAetiology_SBS, by = "proposedAetiologyGrouped")

# Filter signatures per sample below min. contrib.
minContribution <- 5 #5%
nameFilteredMutSigs <- sprintf('Filtered (<%d%% in all samples)', minContribution)
relevantSignatures <- unique(dplyr::filter(dataMutSigsSBS, value > minContribution)$proposedAetiologyGrouped_shorterName)
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::mutate(proposedAetiologyGrouped_shorterName = ifelse(proposedAetiologyGrouped_shorterName %in% relevantSignatures, proposedAetiologyGrouped_shorterName, nameFilteredMutSigs))

# Aggregate filtered signatures per sample together.
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::group_by(sampleId, proposedAetiologyGrouped_shorterName) %>% dplyr::summarise(value = sum(value)) %>% dplyr::ungroup()

# Reorder signatures.
dataMutSigsSBS$proposedAetiologyGrouped_shorterName <- factor(dataMutSigsSBS$proposedAetiologyGrouped_shorterName,
                                                  levels = c("Deamination of 5-methylcytosine (SBS1)",
                                                               "APOBEC activity (SBS2, SBS13)",
                                                               "Defective HR DNA damage repair (SBS3)",
                                                               "Tobacco smoking (SBS4)",
                                                               "Unknown (clock-like signature) (SBS5)",
                                                               "Defective DNA MMR (SBS6, SBS15, SBS21, SBS26, SBS44)",
                                                               "Ultraviolet light exposure (SBS7a-d)",
                                                               "Polimerase eta somatic hypermutation activity (SBS9)",
                                                               "POLE mutations (SBS10a-b)",
                                                               "Temozolomide treatment (SBS11)",
                                                               "Concurrent POLE mutations and defective DNA MMR (SBS14)",
                                                               "Damage by reactive oxygen species (SBS18)",
                                                               "Concurrent POLD1 mutations and defective DNA MMR (SBS20)",
                                                               "Aristolochic acid exposure (SBS22)",
                                                               "Aflatoxin exposure (SBS24)",
                                                               "Chemotherapy treatment (SBS25)",
                                                               "Tobacco chewing (SBS29)",
                                                               "Defective DNA BER due to NTHL1 mutations (SBS30)",
                                                               "Platinum chemotherapy treatment (SBS31, SBS35)",
                                                               "Azathioprine treatment (SBS32)",
                                                               "Defective DNA BER due to MUTYH mutations (SBS36)",
                                                               "Indirect effect of ultraviolet light (SBS38)",
                                                               "Haloalkane exposure (SBS42)",
                                                               "AID activity (SBS84)",
                                                               "Indirect effects of AID activity (SBS85)",
                                                               "Unknown chemotherapy treatment (SBS86)",
                                                               "Thiopurine chemotherapy treatment (SBS87)",
                                                               "Colibactin exposure (SBS88)",
                                                               "Duocarmycin exposure (SBS90)",
                                                               "Unknown (15 signatures)",
                                                             nameFilteredMutSigs))

# Sort samples on given sorting order.
dataMutSigsSBS <- dataMutSigsSBS %>% dplyr::mutate(sampleId  = factor(as.character(sampleId), levels = orderSamples))

colorsMutSigs <- c(
  'Deamination of 5-methylcytosine (SBS1)' = '#85660D',
  'APOBEC activity (SBS2, SBS13)' = '#c11d00',
  'Defective HR DNA damage repair (SBS3)' = '#0079c2',
  'Tobacco smoking (SBS4)' = '#16FF32',
  'Unknown (clock-like signature) (SBS5)' = '#EBB584',
  'Defective DNA MMR (SBS6, SBS15, SBS21, SBS26, SBS44)' = '#f9b320',
  'Ultraviolet light exposure (SBS7a-d)' = 'yellow',
  'Polimerase eta somatic hypermutation activity (SBS9)' = '#A088C3',
  'POLE mutations (SBS10a-b)' = '#5d052e',
  'Temozolomide treatment (SBS11)'  = '#91E4A6',
  'Concurrent POLE mutations and defective DNA MMR (SBS14)' = '#785EF0',
  'Damage by reactive oxygen species (SBS18)' = '#44AA99',
  'Concurrent POLD1 mutations and defective DNA MMR (SBS20)' = '#CC6677',
  'Aristolochic acid exposure (SBS22)' = '#D869C0',
  'Aflatoxin exposure (SBS24)' = '#68564E',
  'Chemotherapy treatment (SBS25)' = '#FF0074',
  'Tobacco chewing (SBS29)' = '#7ED7D1',
  'Defective DNA BER due to NTHL1 mutations (SBS30)' = '#525975',
  'Platinum chemotherapy treatment (SBS31, SBS35)' = '#88CCEE',
  'Azathioprine treatment (SBS32)' = '#FFB6C1',
  'Defective DNA BER due to MUTYH mutations (SBS36)' = '#B10DA1',
  'Indirect effect of ultraviolet light (SBS38)' = '#DEA0FD',
  'Haloalkane exposure (SBS42)' = '#90AD1C',
  'AID activity (SBS84)' = '#FFD599',
  'Indirect effects of AID activity (SBS85)' = '#FAE4C5',
  'Unknown chemotherapy treatment (SBS86)' = '#FA98C5',
  'Thiopurine chemotherapy treatment (SBS87)' = '#D4AEBF',
  'Colibactin exposure (SBS88)' = '#939e64',
  'Duocarmycin exposure (SBS90)' = '#99b8c7',
  'Unknown (15 signatures)'  = 'grey90')


# Add filtering color.
colorsMutSigs <- c(colorsMutSigs, '#D3D3D3')
names(colorsMutSigs)[length(names(colorsMutSigs))] <- nameFilteredMutSigs

# Order on colorsMutSigs
dataMutSigsSBS$proposedAetiologyGrouped_shorterName <- factor(dataMutSigsSBS$proposedAetiologyGrouped_shorterName, levels = names(colorsMutSigs))
dataMutSigsSBS$proposedAetiologyGrouped_shorterName <- base::droplevels(dataMutSigsSBS$proposedAetiologyGrouped_shorterName)

# only keep colours to be used
colorsMutSigs <- colorsMutSigs[names(colorsMutSigs) %in% unique(dataMutSigsSBS$proposedAetiologyGrouped_shorterName)]

# Generate plot.
plotMutSig_SBS <- ggplot(dataMutSigsSBS, aes(x = sampleId, y = value, fill = proposedAetiologyGrouped_shorterName)) +
  geom_bar(stat='identity', color = 'black', width = 1, size = 0.2) +
  labs(y = "SBS\nmut. signature\nrel. contribution", x = NULL) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%')) +
  scale_fill_manual(values = colorsMutSigs, guide = guide_legend(title = 'SBS mutational signatures (COSMIC v3.1)',
                                                                 title.position = 'top',
                                                                 title.hjust = 0,
                                                                 ncol = 1,
                                                                 keywidth = 0.5,
                                                                 keyheight = 0.5), name = NULL, drop = TRUE) +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )

# extract legends
legend_plotMutSig_SBS <- cowplot::get_legend(
  plotMutSig_SBS + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plotMutSig_SBS <- plotMutSig_SBS + theme(legend.position = 'none')



# plot APOBEC enrichment analysis
# order samples

APOBEC_enrichResults <- data.Cohort$APOBEC_enrich %>% dplyr::filter(sample %in% DR176.MetaData$sampleId) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples))

# plot APOBEC mutagenesis low, medium, high, and non
plot.APOBECenrich <- ggplot(APOBEC_enrichResults, aes(sample, y = 'APOBEC mutagenesis', fill = APOBEC_mutagenesis)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  scale_fill_manual(values = c('No' = 'white', 'Low' = '#F5BFAC', 'Medium' = '#EC7F59', 'High' = '#BB0303'),
                    guide = guide_legend(title = 'APOBEC mutagenesis',
                    title.position = 'top',
                    title.hjust = 0, nrow = 1,
                    keywidth = 0.5, keyheight = 0.5)) +
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.APOBECenrich <- cowplot::get_legend(
  plot.APOBECenrich + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.APOBECenrich <- plot.APOBECenrich + theme(legend.position = 'none')

# plot genome-wide ploidy
plot_meanPloidy <- ggplot(results.Cohort_mutationalBurden, aes(x = sample, y = genomePloidy, fill = genomePloidy)) +
  scale_y_continuous(expand = expand_scale(0,0)) +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2) +
  scale_fill_gradient(low = 'white', high = 'black',
                      guide = guide_colorbar(title = 'Mean ploidy', title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 3.5, barheight = 0.5))  +
  labs(x = NULL, y = 'Mean ploidy \ngenome-wide') +
  theme(
    legend.position = 'none', axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'black', linetype = 'longdash'),
    panel.grid.minor.y = element_line(colour = 'black', linetype = 'dashed'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )


# plot chromotripsis
plot.chromothripsis <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'Chromothripsis', fill = hasChromothripsis)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('white', 'black')) + 
  annotationTheme() + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# plot sex of patient
plot.sex <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'Female', fill = gender)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('black', 'white')) + 
  annotationTheme() + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# plot age at Biopsy
  plot.ageAtBiopsy <- ggplot(DR176.MetaData, aes(x = sampleId, y = "Age", fill = ageAtBiopsy)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient(low = 'white', high = "darkorchid4",
                      guide = guide_colorbar(title = 'Age at biopsy',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 3.5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.ageAtBiopsy <- cowplot::get_legend(
  plot.ageAtBiopsy + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.ageAtBiopsy <- plot.ageAtBiopsy + theme(legend.position = 'none')



# plot HR deficiency
plot.CHORD <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'HR status', fill = hr_status)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) + scale_fill_manual(values = c("cannot_be_determined" = 'grey60', "HR_deficient" = 'black', "HR_proficient" ='white'),
                                               labels = c("HR_deficient" = "Deficient",
                                                          "HR_proficient" = "Proficient",
                                                          "cannot_be_determined" = "Undefined"),
                                               guide = guide_legend(title = 'Homologous recombination (HR) status',
                                                                    title.position = 'top',
                                                                    title.hjust = 0, ncol = 3,
                                                                    keywidth = 0.5, keyheight = 0.5)) +
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.CHORD <- cowplot::get_legend(
  plot.CHORD + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.CHORD <- plot.CHORD + theme(legend.position = 'none')


# plot tumor purity
plot_tumorPurity <- ggplot(DR176.MetaData, aes(x = sampleId, y = "Tumor purity", fill = tumorPurity)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = 'white', mid = '#1d5c8f', high = '#1d5c8f', midpoint = 0.7,
                      guide = guide_colorbar(title = 'Tumor purity',
                                             title.position = 'top', title.hjust = 0.5, ncol = 1, barwidth = 3.5, barheight = 0.5))  +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_tumorPurity <- cowplot::get_legend(
  plot_tumorPurity + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_tumorPurity <- plot_tumorPurity + theme(legend.position = 'none')


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
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.preTreatment <- cowplot::get_legend(
  plot.preTreatment + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.preTreatment <- plot.preTreatment + theme(legend.position = 'none')



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
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))


# extract legends
legend_plot.biopsySite <- cowplot::get_legend(
  plot.biopsySite + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.biopsySite <- plot.biopsySite + theme(legend.position = 'none')



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
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.cancerSubtype <- cowplot::get_legend(
  plot.cancerSubtype + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.cancerSubtype <- plot.cancerSubtype + theme(legend.position = 'none')




# protein PD-L1 expression
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(PDL1exp = ifelse(PDL1 == "pos", "pos", 'non-pos'))
DR176.MetaData <- DR176.MetaData %>% dplyr::mutate(PDL1exp = factor(PDL1exp, levels = c('pos', 'non-pos')))

# plot
plot_PDL1exp <- ggplot(DR176.MetaData, aes(x = sampleId, y = "PD-L1 expression", fill = PDL1exp)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = F) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('pos' = 'black', 'non-pos' = 'white', na.value = "grey80"),
                    labels = c('pos' = 'Positive', 'non-pos' = 'Negative', 'na.value' = "NA"),
                    guide = guide_legend(title = 'PD-L1 protein expression',
                                         title.position = 'top', title.hjust = 0, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_PDL1exp <- cowplot::get_legend(
  plot_PDL1exp + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_PDL1exp <- plot_PDL1exp + theme(legend.position = 'none')






# export figures 

# Plot genomic Landscape Figure 2 ---------------------------------------------------

# Plot sample-level overview.
samplePlot <- cowplot::plot_grid(plot.Responders,
                                 plot.TMB,
                                 plot_mutsPerMb,
                                 plot_missenseMuts,
                                 
                                 
                                 plot_clonalFraction,
                                 plot.APOBECenrich,
                                 
                                 plot.GenomicSubtype,
                                 plotMutSig_SBS,
                                 plot_meanPloidy,
                                 plot_tumorPurity,
                                 
                                 plot.CHORD,
                                 plot.chromothripsis,
                                 plot_PDL1exp,
                                 plot.sex,
                                 plot.ageAtBiopsy,
                                 plot.biopsySite,
                                 plot.cancerSubtype,
                                 plot.preTreatment,
                                 ncol = 1, align = 'v', axis = 'tblr',
                                 rel_heights = c(0.1, 0.1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.5, 0.5, rep(0.1,9))
)

legend_samplePlot <- cowplot::plot_grid(legend_plot_mutsPerMb,
                                        legend_plot.Responders,
                                        legend_plot.TMB,

                                        legend_plot_PDL1exp,
                                        legend_plot_clonalFraction,
                                        legend_plot.GenomicSubtype,
                                        legend_plot.APOBECenrich,
                                        legend_plotMutSig_SBS,
                                        legend_plot_tumorPurity,
                                        legend_plot.CHORD,
                                        legend_plot.ageAtBiopsy,
                                        legend_plot.biopsySite,
                                        legend_plot.cancerSubtype,
                                        legend_plot.preTreatment,
                                        
                                        ncol = 1, align = 'v',
                                        rel_heights = c(0.1, 0.2, 0.15, 0.15, 0.15, 0.15, 0.5, 0.25, 0.25, 0.25, 0.3, 0.1, 0.1, 0.1, 0.1))


# export Figure 2
pdf(paste0(odir,"overviewMutLandscape.pdf"),width = 9, height = 6)
cowplot::plot_grid(samplePlot, legend_samplePlot, ncol = 2, rel_widths = c(0.6, 0.4), rel_heights = c(1, 0.8), align = 'h', axis = 'tblr')
dev.off()



# ------------------------------------------------------------------------------------------------------------------
# ---------------------------------- Calculating ROC curves for Genomic features (TMB and APOBEC) --------------
# ------------------------------------------------------------------------------------------------------------------

# plot ROC curve for responders on TMB-high 
data_CohortROC <- results.Cohort$mutationalBurden %>%
  dplyr::left_join(data.Cohort$APOBEC_enrich, by = "sample") %>%
  dplyr::left_join(dplyr::select(DR176.MetaData, sampleId, `PFS days`, `PFS status`, `OS days`, `Survival status`), by = c('sample' = 'sampleId')) %>%
  dplyr::mutate(OutcomeResponse = ifelse(Outcome6M == "responder", 1, 0),
                isAPOBEChigh = ifelse(APOBEC_mutagenesis == "High", 1, 0),
                isTMBhigh = ifelse(tmbStatus == "High TMB (≥10)", 1, 0))

# sort data from low TMB to high
data_CohortROC <- data_CohortROC %>% dplyr::arrange(Genome.TMB)

# calculate glm fit for TMB
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$Genome.TMB, family = binomial)
ROC_TMB <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# calculate glm fit for  APOBEC mutagenesis
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$Genome.TMB+data_CohortROC$foldEnrichment, family = binomial)
ROC_TMB_APOBECmut <- pROC::plot.roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# calculate glm fit for APOBEC mut
data_CohortROC <- data_CohortROC %>% dplyr::arrange(foldEnrichment)
gml_fit_ROC <- stats::glm(data_CohortROC$OutcomeResponse~data_CohortROC$foldEnrichment, family = binomial)
ROC_APOBECmut <- pROC::roc(data_CohortROC$OutcomeResponse, gml_fit_ROC$fitted.values, legacy.axes = TRUE, ci=TRUE)

# test sig between two ROC curves
pValueAUC_APOBEC <- round((pROC::roc.test(ROC_TMB, ROC_APOBECmut, method = 'delong'))$p.value, 2)
pValueAUC_TMBAPOBEC <- round((pROC::roc.test(ROC_TMB, ROC_TMB_APOBECmut, method = 'delong'))$p.value, 2)

# Multiple curves:
plot_ROC_tmb_APOBECmut <- pROC::ggroc(list(ROC_TMB, ROC_APOBECmut, ROC_TMB_APOBECmut), legacy.axes=TRUE) +
  labs(y = 'True positive rate', x = 'False positive rate') +
  scale_color_manual(values = c("#1762ad", "#f0a630", "#008837"),
                     labels = c('1' = 'TMB', '2' = 'APOBEC',
                                '3' = 'TMB + APOBEC'),
                    guide = guide_legend(title = 'ROC curve',
                                         title.position = 'top',
                                         title.hjust = 0, ncol = 1,
                                         keywidth = 0.5, keyheight = 0.5,
                                         override.aes = list(size=3))) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", size = 0.25) +
  annotate(geom="text", x=0.5, y=0.3, size=2,
           label = paste0("AUC = ", format(round(ROC_TMB$auc, 2))), color="#1762ad", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.225, size=2,
           label = paste0("AUC = ", format(round(ROC_APOBECmut$auc, 2)),
                          ", p = ", pValueAUC_APOBEC), color="#f0a630", hjust = 0) +
  annotate(geom="text", x=0.5, y=0.15, size=2,
           label = paste0("AUC = ", format(round(ROC_TMB_APOBECmut$auc, 2)),
                          ", p = ", pValueAUC_TMBAPOBEC), color="#008837", hjust = 0) +
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 9),
    text=element_text(size=9, family='Helvetica'),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )



# plot Response rate
data_CohortROC_responseRate <- dplyr::count(data_CohortROC, Outcome6M, isTMBhigh) %>%
  dplyr::group_by(isTMBhigh) %>% dplyr::mutate(totalTMBgroup = sum(n)) %>% dplyr::ungroup() %>%
  dplyr::mutate(responderRate = 100*n/totalTMBgroup,
                isTMBhigh = ifelse(isTMBhigh == 0, 'Low', 'High')) %>%
  dplyr::mutate(isTMBhigh = factor(isTMBhigh, levels = c('High', 'Low')))

# calculate fisher exact test
nResponderGroup1 <- data_CohortROC_responseRate %>% dplyr::filter(isTMBhigh == "High" &
                                                                    Outcome6M == "responder")
nResponderGroup2 <- data_CohortROC_responseRate %>% dplyr::filter(isTMBhigh == "Low" &
                                                                    Outcome6M == "responder")

challenge.df = matrix(c(nResponderGroup1$n, nResponderGroup2$n,
                        nResponderGroup1$totalTMBgroup - nResponderGroup1$n, nResponderGroup2$totalTMBgroup - nResponderGroup2$n), nrow = 2)
pFisherExacTest <- fisher.test(challenge.df)$p.value

# sort data
data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(Outcome6M = factor(Outcome6M, levels = c("non-responder", "responder")))

# change tmbStatus to include number of samples in each group
data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(isTMBhigh = paste0(isTMBhigh, ' (n=', totalTMBgroup, ')')) %>%
  dplyr::mutate(isTMBhigh = factor(isTMBhigh, levels = c('High (n=38)', 'Low (n=32)')))


# plot histograms of response rates with 95% CI from binomial distribution
plot_ResponseRate_TMB <- ggplot(data_CohortROC_responseRate, aes(x=isTMBhigh, y=responderRate, fill = Outcome6M)) +
  geom_bar(stat='identity', color = 'black', width = 0.6) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%'))  +
  labs(title = paste0("p = ", format(round(pFisherExacTest, 3))),
       y = "Response (cumulative cases)", x = 'TMB') +
  annotate(geom="text", x=1, y=nResponderGroup1$responderRate+4, size=2,
           label = paste0(round(nResponderGroup1$responderRate), "%"), color="black") +
  annotate(geom="text", x=2, y=nResponderGroup2$responderRate+4, size=2,
           label = paste0(round(nResponderGroup2$responderRate), "%"), color="black") +
  scale_fill_manual(guide = guide_legend(title = 'Response to treatment', nrow = 1,
                                         title.position = 'top', title.hjust = 0.5, keywidth = 0.5, keyheight = 0.5),
                    name = NULL, values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder') ) +
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


# plot Response rate using APOBEC high
data_CohortROC_responseRate <- dplyr::count(data_CohortROC, Outcome6M, isAPOBEChigh) %>%
  dplyr::group_by(isAPOBEChigh) %>% dplyr::mutate(totalAPOBECgroup = sum(n)) %>% dplyr::ungroup() %>%
  dplyr::mutate(responderRate = 100*n/totalAPOBECgroup,
                isAPOBEChigh = ifelse(isAPOBEChigh == 0, 'Non-high', 'High')) %>%
  dplyr::mutate(isAPOBEChigh = factor(isAPOBEChigh, levels = c('High', 'Non-high')))


# calculate fisher exact test
nResponderGroup1 <- data_CohortROC_responseRate %>% dplyr::filter(isAPOBEChigh == "High" &
                                                                    Outcome6M == "responder")
nResponderGroup2 <- data_CohortROC_responseRate %>% dplyr::filter(isAPOBEChigh == "Non-high" &
                                                                    Outcome6M == "responder")

challenge.df = matrix(c(nResponderGroup1$n, nResponderGroup2$n,
                        nResponderGroup1$totalAPOBECgroup - nResponderGroup1$n, nResponderGroup2$totalAPOBECgroup - nResponderGroup2$n), nrow = 2)

pFisherExacTest <- fisher.test(challenge.df)$p.value

# sort data
data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(Outcome6M = factor(Outcome6M, levels = c("non-responder", "responder")))

# change APOBEC mut status to include number of samples in each group
data_CohortROC_responseRate <- data_CohortROC_responseRate %>%
  dplyr::mutate(isAPOBEChigh = paste0(isAPOBEChigh, ' (n=', totalAPOBECgroup, ')')) %>%
  dplyr::mutate(isAPOBEChigh = factor(isAPOBEChigh, levels = c('High (n=29)', 'Non-high (n=41)')))


# plot histograms of response rates
plot_ResponseRate_APOBEC <- ggplot(data_CohortROC_responseRate, aes(x=isAPOBEChigh, y=responderRate, fill = Outcome6M)) +
  geom_bar(stat='identity', color = 'black', width = 0.6) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%'))  +
  labs(title = paste0("p = ", format(round(pFisherExacTest, 3))),
       y = "Response (cumulative cases)", x = 'APOBEC mutagenesis') +
  annotate(geom="text", x=1, y=nResponderGroup1$responderRate+4, size=2,
           label = paste0(round(nResponderGroup1$responderRate), "%"), color="black") +
  annotate(geom="text", x=2, y=nResponderGroup2$responderRate+4, size=2,
           label = paste0(round(nResponderGroup2$responderRate), "%"), color="black") +
  scale_fill_manual(guide = guide_legend(title = 'Response to treatment', nrow = 1,
                                         title.position = 'top', title.hjust = 0.5, keywidth = 0.5, keyheight = 0.5),
                    name = NULL, values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder') ) +
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


# save ROC and Response rate Figure S1
pdf(paste0(odir,"ROC_ResponseRate_TMB_APOBEC.pdf"),width = 7, height = 3)
cowplot::plot_grid(plot_ResponseRate_TMB, plot_ResponseRate_APOBEC, plot_ROC_tmb_APOBECmut,
                   nrow = 1, align = 'vh', axis = 'tblr')
dev.off()




# Plot OS and PFS curves (order factor according to numerator and denominator)
data_CohortROC_survival <- data_CohortROC %>% dplyr::select(PFS_days = `PFS days`, PFS_status = `PFS status`, OS_status = `Survival status`, OS_days = `OS days`, isAPOBEChigh, isTMBhigh, response = OutcomeResponse, Outcome6M) %>%
  dplyr::mutate(response = 1 - response,
                OS_days = 12*OS_days/365,
                PFS_days = 12*PFS_days/365,
                isAPOBEChigh = ifelse(isAPOBEChigh == 1, "High", "Non-high"),
                isTMBhigh = ifelse(isTMBhigh == 1, "High", "Low")) %>%
  dplyr::mutate(isAPOBEChigh = factor(isAPOBEChigh, levels = c("High", "Non-high")),
                isTMBhigh = factor(isTMBhigh, levels = c("High", "Low")))


# get OS survival functions
fit_OS_TMB <- survminer::surv_fit(survival::Surv(OS_days,  OS_status) ~ isTMBhigh, data = data_CohortROC_survival)
fit_OS_APOBEC <- survminer::surv_fit(survival::Surv(OS_days,  OS_status) ~ isAPOBEChigh, data = data_CohortROC_survival)
fit_PFS_TMB <- survminer::surv_fit(survival::Surv(PFS_days,  PFS_status) ~ isTMBhigh, data = data_CohortROC_survival)
fit_PFS_APOBEC <- survminer::surv_fit(survival::Surv(PFS_days,  PFS_status) ~ isAPOBEChigh, data = data_CohortROC_survival)



# change factors for TMB and APOBEC to calculate HR
data_CohortROC_survival$isTMBhigh <- factor(data_CohortROC_survival$isTMBhigh, levels = c("Low", "High"))
data_CohortROC_survival$isAPOBEChigh <- factor(data_CohortROC_survival$isAPOBEChigh, levels = c("Non-high", "High"))


# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(OS_days, OS_status) ~ isTMBhigh, data = data_CohortROC_survival)
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))

# plot survival curve
plot_OSCurve_TMB <- survminer::ggsurvplot(fit_OS_TMB, data = data_CohortROC_survival,
                                                size = 1,                 # change line size
                                                palette =  c('#3470D3', '#CD8B35'),
                                                xlab = "Months",
                                                ylab = "OS probability (n = 70)",
                                          title = res.cox_HazardRatio,
                                                legend.title = "TMB",
                                          legend = c(0.7,0.6),
                                          xlim = c(0,37),
                                                pval.method.size = 5,
                                                break.time.by = 5,
                                                pval = TRUE,              # Add p-value (default log-rank)
                                                pval.size = 3,
                                                pval.coord = c(10, 0.9),
                                          surv.median.line = "hv",
                                                risk.table = TRUE,        # Add risk table
                                                fontsize = 2,   # fontsize for table
                                          tables.y.text.col = TRUE,
                                          risk.table.col = "strata",
                                                legend.labs = c("High", "Low"),    # Change legend labels
                                          ggtheme = theme(
                                            legend.position = 'bottom',
                                            text=element_text(size=10, family='Helvetica'),
                                            legend.title=element_text(size=8, hjust = 0.5), 
                                            legend.text=element_text(size=8),
                                            plot.title = element_text(size=8, hjust = 0),
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
plot_OSCurve_TMB$table <- plot_OSCurve_TMB$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )



# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(OS_days, OS_status) ~ isAPOBEChigh, data = data_CohortROC_survival)
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))


# plot survival curve
plot_OSCurve_APOBEC <- survminer::ggsurvplot(fit_OS_APOBEC, data = data_CohortROC_survival,
                                                   size = 1,                 # change line size
                                                   palette =  c('#3470D3', '#CD8B35'),
                                                   xlab = "Months",
                                                   ylab = "OS probability (n = 70)",
                                             title = res.cox_HazardRatio,
                                                   legend.title = "APOBEC",
                                             legend = c(0.7,0.6),
                                             xlim = c(0,37),
                                                   pval.method.size = 5,
                                                   break.time.by = 5,
                                                   pval = TRUE,              # Add p-value (default log-rank)
                                                   pval.size = 3,
                                                   pval.coord = c(10, 0.9),
                                             surv.median.line = "hv",
                                                   risk.table = TRUE,        # Add risk table
                                                   fontsize = 2,   # fontsize for table
                                                   risk.table.col = "strata",# Risk table color by groups
                                                   legend.labs = c("High", "Non-high"),    # Change legend labels
                                                   risk.table.height = 0.25, # Useful to change when you have multiple groups
                                                   tables.y.text.col = FALSE, 
                                             ggtheme = theme(
                                               legend.position = 'bottom',
                                               text=element_text(size=10, family='Helvetica'),
                                               legend.title=element_text(size=8, hjust = 0.5), 
                                               legend.text=element_text(size=8),
                                               plot.title = element_text(size=8, hjust = 0),
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
plot_OSCurve_APOBEC$table <- plot_OSCurve_APOBEC$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )



# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(PFS_days, PFS_status) ~ isTMBhigh, data = data_CohortROC_survival)
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))

plot_PFSCurve_TMB <- survminer::ggsurvplot(fit_PFS_TMB, data = data_CohortROC_survival,
                                            size = 1,                 # change line size
                                            palette =  c('#3470D3', '#CD8B35'),
                                            xlab = "Months",
                                            ylab = "PFS probability (n = 70)",
                                           title = res.cox_HazardRatio,
                                            legend.title = "TMB",
                                           legend = c(0.7,0.6),
                                            pval.method.size = 5,
                                            break.time.by = 5,
                                            #conf.int = TRUE,          # Add confidence interval
                                            pval = TRUE,              # Add p-value (default log-rank)
                                            pval.size = 3,
                                            pval.coord = c(10, 0.9),
                                           surv.median.line = "hv",
                                            risk.table = TRUE,        # Add risk table
                                            fontsize = 2,   # fontsize for table
                                            risk.table.col = "strata",# Risk table color by groups
                                            legend.labs = c("High", "Low"),    # Change legend labels
                                            risk.table.height = 0.25, # Useful to change when you have multiple groups
                                            tables.y.text.col = FALSE, 
                                           ggtheme = theme(
                                             legend.position = 'bottom',
                                             text=element_text(size=10, family='Helvetica'),
                                             legend.title=element_text(size=8, hjust = 0.5), 
                                             legend.text=element_text(size=8),
                                             plot.title = element_text(size=8, hjust = 0),
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
plot_PFSCurve_TMB$table <- plot_PFSCurve_TMB$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )



# get cox hazard ratios
res.cox_HazardRatio <- survival::coxph(survival::Surv(PFS_days, PFS_status) ~ isAPOBEChigh, data = data_CohortROC_survival)
res.cox_HazardRatio <- summary(res.cox_HazardRatio)
res.cox_HazardRatio <- paste0("hr = ",
                              round(res.cox_HazardRatio$conf.int[1], 2), ", CI: ",
                              round(res.cox_HazardRatio$conf.int[3], 2), "-",
                              round(res.cox_HazardRatio$conf.int[4], 2),
                              ", ", ifelse(res.cox_HazardRatio$logtest["pvalue"] > 0.0001,
                                           paste0("p = ", format(round(res.cox_HazardRatio$logtest["pvalue"], 4), scientific=FALSE)),
                                           "p < 0.0001"))

plot_PFSCurve_APOBEC <- survminer::ggsurvplot(fit_PFS_APOBEC, data = data_CohortROC_survival,
                                                size = 1,                 # change line size
                                                palette =  c('#3470D3', '#CD8B35'),
                                                xlab = "Months",
                                                ylab = "PFS probability (n = 70)",
                                              title = res.cox_HazardRatio,
                                                legend.title = "APOBEC",
                                              legend = c(0.7,0.6),
                                              pval.method.size = 5,
                                              break.time.by = 5,
                                              pval = TRUE,              # Add p-value (default log-rank)
                                              pval.size = 3,
                                              pval.coord = c(10, 0.9),
                                              surv.median.line = "hv",
                                                risk.table = TRUE,        # Add risk table
                                                fontsize = 2,   # fontsize for table
                                                risk.table.col = "strata",# Risk table color by groups
                                                legend.labs = c("High", "Non-high"),    # Change legend labels
                                                risk.table.height = 0.25, # Useful to change when you have multiple groups
                                                tables.y.text.col = FALSE, 
                                              ggtheme = theme(
                                                legend.position = 'bottom',
                                                text=element_text(size=10, family='Helvetica'),
                                                legend.title=element_text(size=8, hjust = 0.5), 
                                                legend.text=element_text(size=8),
                                                plot.title = element_text(size=8, hjust = 0),
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
plot_PFSCurve_APOBEC$table <- plot_PFSCurve_APOBEC$table +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )


# save OS and PFS Figure 5c-d
pdf(paste0(odir,"Survival_TMB_APOBEC.pdf"),width = 9, height = 3)


cowplot::plot_grid(
  
  cowplot::plot_grid(
    plot_OSCurve_TMB$plot,
    plot_OSCurve_TMB$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
  
  cowplot::plot_grid(
    plot_PFSCurve_TMB$plot,
    plot_PFSCurve_TMB$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
  
  cowplot::plot_grid(
    plot_OSCurve_APOBEC$plot,
    plot_OSCurve_APOBEC$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
  
  cowplot::plot_grid(
    plot_PFSCurve_APOBEC$plot,
    plot_PFSCurve_APOBEC$table,
    align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(1, 0.35)),
  
  align = 'v', axis = 'tblr', nrow = 1)


dev.off()




#================================= Pathway analysis and DDR genes ============================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





#================================= Prepare data for Supplementary Fig 2 ============================

#++++++++++++++++++++++++++++++++++++ Plot oncoplot

# Melt GISTIC peaks
gisticPeaks <- reshape2::melt(data.frame(mcols(results.Cohort$GISTIC2$fullCohort$gisticNarrowPeaksWithAnno)), id.vars = c('Unique.Name', 'Descriptor', 'q.values', 'nGenes.GENCODE', 'overlapGenes.Final'))
gisticPeaks_aux <- data.frame(mcols(results.Cohort$GISTIC2$fullCohort$gisticNarrowPeaksWithAnno))
  
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

# rename 'dee' amplification/deletion' to just keep 'amplification' and 'deletion'
dataOncoplot <- dataOncoplot %>% dplyr::mutate(Consequence.CNA = gsub("Deep ", "", Consequence.CNA),
                                               Consequence.SV = gsub("Va", "va", Consequence.SV))

# Sometime results repeated for a sample (amplification and deletion in the same patient)
# discard these duplicate samples
dataOncoplot <- dataOncoplot %>% dplyr::distinct(sample, SYMBOL, .keep_all = TRUE)


# number of samples
nSamples <- length(unique(DR176.MetaData$sampleId))

# Add axis names showing number of unique mutated samples.
dataOncoplot <- dataOncoplot %>% dplyr::group_by(ENSEMBL) %>% dplyr::add_tally() %>%
  dplyr::mutate(axisName = sprintf('%s (%s%%)', SYMBOL, round(100*n/nSamples))) %>% dplyr::ungroup()
dataOncoplot$sample <- factor(dataOncoplot$sample, level = DR176.MetaData$sampleId)
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

# Get order of samples
orderSamples_drivers <- colnames(memoData)

# order samples by driver genes and by responders and non-responders
results.Cohort_mutationalBurden <- results.Cohort_mutationalBurden[match(orderSamples_drivers, results.Cohort_mutationalBurden$sample),]
results.Cohort_mutationalBurden <- results.Cohort_mutationalBurden %>%
  dplyr::arrange(factor(tmbStatus, levels = c("High TMB (≥10)", "Low TMB (<10)"))) %>%
  dplyr::arrange(factor(Outcome6M, levels = c("responder", "non-responder")))
  
# get order of samples
orderSamples_drivers <- as.character(results.Cohort_mutationalBurden$sample)

# Order samples
dataOncoplot$sample <- factor(dataOncoplot$sample, levels = orderSamples_drivers)



# Also order by p-value and Responder vs Non-resopnder

# prepare data for Fisher exact test
results.Cohort_mutationalBurden_test <- dplyr::select(results.Cohort_mutationalBurden, sample, Outcome6M) %>%
  dplyr::add_count(Outcome6M, name = "nSizeOutcome")
dataOncoplot_test <- dplyr::mutate(dataOncoplot, Outcome6M = NULL) %>%
  dplyr::full_join(results.Cohort_mutationalBurden_test, by = c("sample")) %>%
  dplyr::group_by(Outcome6M, axisName) %>%
  dplyr::mutate(nOutcome = sum(isMutant)) %>%
  dplyr::ungroup(Outcome6M, axisName) %>%
  dplyr::distinct(Outcome6M, axisName, .keep_all = TRUE)
dataOncoplot_test_Responder <- dplyr::filter(dataOncoplot_test, Outcome6M == "responder") %>% dplyr::select(axisName, groupSize1 = nSizeOutcome, nGenesCount1 = nOutcome)
dataOncoplot_test_NonResponder <- dplyr::filter(dataOncoplot_test, Outcome6M == "non-responder") %>% dplyr::select(axisName, groupSize2 = nSizeOutcome, nGenesCount2 = nOutcome)
dataOncoplot_test <- dplyr::full_join(dataOncoplot_test_Responder, dataOncoplot_test_NonResponder, by = c("axisName"))

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
  labs(x = NULL, y = "Driver genes") +
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
        axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5),
        axis.text.y = element_blank(),
        #plot.margin = unit(c(0, -1, -1, 0), "cm"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=10, family='Helvetica', colour = "black"))


# extract legends
legend_plot_pValues_oncoplot <- cowplot::get_legend(
  plot_pValues_oncoplot + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_pValues_oncoplot <- plot_pValues_oncoplot + theme(legend.position = 'none')




#++++++++++++++++++++++++++++++++++++ Plot Fusion genes
Fusions_data <- data.Cohort$fusionsLINX

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
Fusions_data_summary$sample <- factor(Fusions_data_summary$sample, levels = orderSamples_drivers)
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
    #axis.text.x=element_text(angle = 45, hjust = 1),
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
                snv_id = paste0(chr, ":", pos, ref))

# Find HotSpot SNVs 
SNV_data_summary_all <- SNV_data %>%
  #dplyr::filter(geneName %in% c('TERT', 'ADGRG6', 'PLEKHS1', 'TBC1D12', 'LEPROTL1')) %>%
  dplyr::group_by(snv_id) %>% dplyr::add_count(name = "totalHotspots", sort = TRUE) %>%
  dplyr::filter(totalHotspots>1) %>% dplyr::ungroup()


# keep only protein coding and intergenic variants
SNV_data_summary_all <- SNV_data_summary_all %>% dplyr::filter(biotype == 'protein_coding' | consequence == 'intergenic_variant')

# keep only most frequent ones (>= 4 of the cohort)
SNV_data_summary <- SNV_data_summary_all %>% dplyr::filter(totalHotspots>=7)

# Change ID of hotspot mutations to include gene name
SNV_data_summary <- SNV_data_summary %>% dplyr::mutate(snv_id_aux = ifelse(is.na(HGVSp), snv_id, gsub(".*:", "", HGVSp))) %>%
  dplyr::mutate(snv_id_aux = ifelse(is.na(geneName), snv_id_aux, paste0(geneName, ", ", snv_id_aux)))

# Get % of hotspot muts in the cohort
SNV_data_summary <- SNV_data_summary %>%
  dplyr::mutate(axisName = paste0(snv_id_aux, " (", round(100 * totalHotspots/length(unique(sample))), "%)"))

# order axis of hotspots
HotSpot_order <- rev(unique(SNV_data_summary$axisName))
SNV_data_summary$axisName <- factor(SNV_data_summary$axisName, levels = HotSpot_order)

# order samples
SNV_data_summary$sample <- factor(SNV_data_summary$sample, levels = orderSamples_drivers)


# Define consequence of hotspot
SNV_data_summary <- SNV_data_summary %>%
  dplyr::mutate(consequence = ifelse(consequence == "5_prime_UTR_variant", "5' UTR",
                                     ifelse(consequence == "intron_variant", "Intron variant",
                                            ifelse(consequence == "missense_variant", "Missense variant",
                                                   ifelse(consequence == "upstream_gene_variant", "Upstream gene variant",
                                                          ifelse(consequence == "downstream_gene_variant", "Downstream gene variant", "Intergenic variant"))))))
SNV_data_summary$consequence <- factor(SNV_data_summary$consequence, levels = c("Upstream gene variant", "Downstream gene variant", "5' UTR", "Intron variant", "Missense variant", "Intergenic variant"))




# Order by p-value Responder vs Non-resopnder
SNV_data_summary$isMutant <- 1
SNV_data_summary <- SNV_data_summary %>% tidyr::complete(sample, axisName)

# prepare data for Fisher exact test
results.Cohort_mutationalBurden_test <- dplyr::select(results.Cohort_mutationalBurden, sample, Outcome6M) %>%
  dplyr::add_count(Outcome6M, name = "nSizeOutcome")

SNV_data_summary_test <- dplyr::mutate(SNV_data_summary, Outcome6M = NULL) %>%
  dplyr::left_join(results.Cohort_mutationalBurden_test, by = c("sample")) %>%
  dplyr::group_by(Outcome6M, axisName) %>%
  dplyr::mutate(nOutcome = sum(isMutant, na.rm = TRUE)) %>%
  dplyr::ungroup(Outcome6M, axisName) %>%
  dplyr::distinct(Outcome6M, axisName, .keep_all = TRUE)
SNV_data_summary_test_Responder <- dplyr::filter(SNV_data_summary_test, Outcome6M == "responder") %>% dplyr::select(axisName, groupSize1 = nSizeOutcome, nGenesCount1 = nOutcome)
SNV_data_summary_test_NonResponder <- dplyr::filter(SNV_data_summary_test, Outcome6M == "non-responder") %>% dplyr::select(axisName, groupSize2 = nSizeOutcome, nGenesCount2 = nOutcome)
SNV_data_summary_test <- dplyr::full_join(SNV_data_summary_test_Responder, SNV_data_summary_test_NonResponder, by = c("axisName"))

# Calculate Pvalues
SNV_data_summary_test$pValue <- 1

for (iPatient in c(1:nrow(SNV_data_summary_test))) {
  challenge.df = matrix(c(SNV_data_summary_test$nGenesCount1[iPatient], SNV_data_summary_test$nGenesCount2[iPatient],
                          SNV_data_summary_test$groupSize1[iPatient]-SNV_data_summary_test$nGenesCount1[iPatient], SNV_data_summary_test$groupSize2[iPatient]-SNV_data_summary_test$nGenesCount2[iPatient]), nrow = 2)
  
  SNV_data_summary_test$pValue[iPatient] <- fisher.test(challenge.df)$p.value
}

# Adjust pValues for multiple testing
SNV_data_summary_test <- SNV_data_summary_test %>% dplyr::mutate(pAdj = p.adjust(pValue, method = "BH")) %>% dplyr::ungroup()

# order from low to high pValue
SNV_data_summary_test <- SNV_data_summary_test %>% dplyr::arrange(pValue)

# order data
SNV_data_summary$axisName <- factor(SNV_data_summary$axisName, levels = as.character(rev(SNV_data_summary_test$axisName)))





# Colors for Consequence of hotpot
# Color of the mutations.
colorHotSpot <- c('Missense variant' = '#000000', # Black
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

# extract legends
legend_plot_HotSpot_Cohort <- cowplot::get_legend(
  plot_HotSpot_Cohort + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_HotSpot_Cohort <- plot_HotSpot_Cohort + theme(legend.position = 'none')




# Plot P values for HotSpots
SNV_data_summary_test$axisName <- factor(SNV_data_summary_test$axisName, levels = rev(SNV_data_summary_test$axisName))
SNV_data_summary_test <- SNV_data_summary_test %>%
  dplyr::mutate(pValue_tmp = ifelse(pValue < 0.01, 0.01, pValue))

# plot p-value
plot_pValues_hotspots <- ggplot(SNV_data_summary_test, aes(x = "P-value", y = axisName, fill = -log10(pValue_tmp))) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_gradient2(low = "white", mid = "#FFFFBF", high = '#9E0142', midpoint = 1,
                       breaks=c(0, 1, 2),labels=c("1", "0.1", "<0.01"),
                       limits = c(0, 2),
                       guide = guide_colorbar(title = 'P-value', reverse = TRUE,
                                              title.position = 'top', title.hjust = 0.5, barwidth = 3.5, barheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        #axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5),
        axis.text.y = element_blank(),
        #plot.margin = unit(c(0, -1, -1, 0), "cm"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=10, family='Helvetica', colour = "black"))








# Order responder data 
results.Cohort_mutationalBurden$sample <- factor(results.Cohort_mutationalBurden$sample,
                                                 levels = orderSamples_drivers)

# plot Responder labels
plot_Responders <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = "Response to treatment", title.position = 'top', title.hjust = 0,
                                         ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=10, family='Helvetica', colour = "black"))

# extract legends
legend_plot_Responders <- cowplot::get_legend(
  plot_Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_Responders <- plot_Responders + theme(legend.position = 'none')


# plot TMB > 10
plot_TMB <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'Tumor mutational burden', fill = tmbStatus)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('#8E1B37', '#7A9BCF'),
                    guide = guide_legend(title = "Tumor mutational\nburden (TMB)", title.position = 'top', title.hjust = 0.5,
                                         ncol = 1, keywidth = 0.5, keyheight = 0.5)) + 
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 10),
        text=element_text(size=10, family='Helvetica', colour = "black"))

# extract legends
legend_plot_TMB <- cowplot::get_legend(
  plot_TMB + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_TMB <- plot_TMB + theme(legend.position = 'none')



# Plot oncoplot Figure S2 ---------------------------------------------------

# Plot sample-level overview.
oncoplotOverview <- cowplot::plot_grid(plot_Responders, 
                                       plot_TMB, 
                                       plot_oncoResponderCohort, 
                                       plot_HotSpot_Cohort, 
                                       plot_fusionGenes, 
                                       ncol = 1,
                                       align = 'v',
                                       axis = 'tblr',
                                       rel_heights = c(0.07, 0.07, 2.3, 0.4, 0.2)
)

oncoplotOverview_pValue <- cowplot::plot_grid(NULL,
                                       NULL,
                                       plot_pValues_oncoplot,
                                       plot_pValues_hotspots,
                                       NULL,
                                       ncol = 1,
                                       align = 'v',
                                       axis = 'tblr',
                                       rel_heights = c(0.09, 0.09, 2.3, 0.4, 0.2)
)

legend_oncoplotOverview <- cowplot::plot_grid(legend_plot_TMB,
                                              legend_plot_Responders,
                                              legend_plot_oncoResponderCohort,
                                              legend_plot_pValues_oncoplot,
                                              legend_plot_HotSpot_Cohort,
                                              ncol = 1, align = 'v')



pdf(paste0(odir,"oncoplotResponderCohort.pdf"),width = 10, height = 9)#, width = 14, height = 21)
cowplot::plot_grid(oncoplotOverview, oncoplotOverview_pValue, legend_oncoplotOverview,
                   ncol = 3, rel_widths = c(0.8, 0.02, 0.2), align = 'h', axis = 'tblr')
dev.off()






#================================= Plot Figure S3a-b - Pathway analysis + DDR genes ============================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# load ddr genes from Taber et. al, 2020 and pathway analysis from Nakauma-Gonzalez, et al., 2022
DDRgenes <- readxl::read_xlsx("/databases/DDRgenes.xlsx")
pathwayGenes <- xlsx::read.xlsx(file = "/databases/GeneAndPathways_TCGA.xlsx", sheetIndex = 2)

# One gene changed name from the list
pathwayGenes <- pathwayGenes %>% dplyr::mutate(Gene = ifelse(Gene == "NOV", "CCN3", as.character(Gene)))

# add DDR genes as extra pathway
DDRgenes <- DDRgenes %>% dplyr::mutate(Pathway = "DDR") %>% dplyr::select(Pathway, geneName)
pathwayGenes <- pathwayGenes %>% dplyr::select(Pathway, geneName = Gene) %>%
  rbind(DDRgenes) %>% dplyr::group_by(Pathway) %>% dplyr::add_count(name = 'sizePathway')

#Get all Genes for pathway analysis from combined report
dataMutPathway <- results.Cohort$combinedReport %>% dplyr::filter(SYMBOL %in% pathwayGenes$geneName)

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

# order samples
orderSamples <- (results.Cohort$mutationalBurden %>% dplyr::arrange(-Genome.TMB) %>%
                   dplyr::arrange(factor(Outcome6M, levels = c("responder", "non-responder"))) %>%
                   dplyr::distinct(sample))$sample


dataMutPathway <- dataMutPathway %>% dplyr::mutate(sample = factor(sample, levels = orderSamples))

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





# Melt the TMB
mutsPerMb <- reshape2::melt(results.Cohort_mutationalBurden %>% dplyr::select(sample, SNV = Muts.SNV, Indel = Muts.InDel, MNV = Muts.MNV, SV = totalSV))
mutsPerMb <- mutsPerMb %>% dplyr::mutate(value = value/genomeSize,
                                         variable = factor(variable, levels = c('SNV', 'Indel', 'MNV', 'SV')))

# order samples
mutsPerMb$sample <- factor(mutsPerMb$sample, levels = orderSamples)

# Barplot of total Mutations per Sample
plot_mutsPerMb <- ggplot(mutsPerMb, aes(x = sample, y = value, fill = variable)) +
  coord_trans(y = 'sqrt') +
  scale_y_continuous(breaks = c(0, 10, 50, 100), expand = expand_scale(0,0)) +
  scale_x_discrete(expand = expand_scale(0,0)) +
  geom_bar(stat = 'identity', col = 'grey20', size = 0.2) +
  scale_fill_manual('Mutational category', values = c('SNV' = '#375D96', 'Indel' = '#FFBA00', 'MNV' = 'salmon', 'SV' = 'mediumseagreen'),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, nrow = 2, keywidth = 0.5, keyheight = 0.5)) +
  labs(x = NULL, y = 'TMB') +
  theme(
    legend.position = 'right',
    axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey60', linetype = 'longdash'),
    #panel.grid.minor.y = element_line(colour = 'grey60', linetype = 'dashed'),
    panel.grid.minor.y = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )

# extract legends
legend_plot_mutsPerMb <- cowplot::get_legend(
  plot_mutsPerMb + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_mutsPerMb <- plot_mutsPerMb + theme(legend.position = 'none')



# plot TMB (order samples)
results.Cohort_mutationalBurden <- results.Cohort_mutationalBurden %>%
  dplyr::mutate(sample = factor (sample, levels = orderSamples),
                tmbStatus = factor (tmbStatus, levels = c('High TMB (≥10)', 'Low TMB (<10)')))

# plot TMB > 10
plot.TMB <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = 'TMB category', fill = tmbStatus)) + 
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) + 
  labs(y = NULL, x = NULL) + 
  scale_fill_manual(values = c('#8E1B37', '#7A9BCF'),
                    #labels = c('High TMB (≥10)' = "High TMB (>10)", 'Medium TMB (≥5-10)' = 'Medium TMB (5-10)', 'Low TMB (0-5)' = 'Low TMB (0-5)'),
                    guide = guide_legend(title = "TMB category", title.position = 'top', title.hjust = 0,
                                         ncol = 1, keywidth = 0.5, keyheight = 0.5)) + 
  annotationTheme() + 
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size = 8, family='Helvetica', colour = "black"))

# extract legends
legend_plot.TMB <- cowplot::get_legend(
  plot.TMB + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.TMB <- plot.TMB + theme(legend.position = 'none')



# plot Responder labels
plot_Responders <- ggplot(results.Cohort_mutationalBurden, aes(sample, y = "Response to treatment", fill = Outcome6M)) +
  geom_tile(colour = 'grey20', size = 0.2, na.rm = T) +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual(values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder'),
                    guide = guide_legend(title = "Response to treatment", title.position = 'top', title.hjust = 0,
                                         ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"))

# extract legends
legend_plot_Responders <- cowplot::get_legend(
  plot_Responders + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_Responders <- plot_Responders + theme(legend.position = 'none')



# plot proportion of responders with high and Low TMB per pathway
dataMutPathway <- dataMutPathway %>% dplyr::left_join(dplyr::select(DR176.MetaData, sampleId, Outcome6M), by = c('sample' = 'sampleId')) %>%
  dplyr::left_join(dplyr::select(results.Cohort_mutationalBurden, sample, Genome.TMB), by = 'sample') %>%
  dplyr::mutate(tmbCategory = ifelse(Genome.TMB > 10, 'TMB > 10', 'TMB < 10'))

# number of patients with TMB high,low and number of mutant pathways per pathway
dataMutPathway <- dataMutPathway %>% dplyr::group_by(Pathway, isMutant) %>%
  dplyr::add_count(name = 'categorySize') %>% dplyr::ungroup() %>%
  dplyr::group_by(Pathway, isMutant, Outcome6M) %>%
  dplyr::add_count(name = 'responseRate') %>% dplyr::ungroup() %>%
  dplyr::mutate(isMutant = ifelse(isMutant == 1, 'Mutant', 'Non-mutant'))

# prepare data for Fisher exact test
dataMutPathway_test_MutantPathway <- dplyr::filter(dataMutPathway, isMutant == 'Mutant') %>% dplyr::select(Pathway, Outcome6M, groupSize1 = categorySize, responseRate_1 = responseRate)
dataMutPathway_test_nonMutantPathway <- dplyr::filter(dataMutPathway, isMutant == 'Non-mutant') %>% dplyr::select(Pathway, Outcome6M, groupSize2 = categorySize, responseRate_2 = responseRate)
dataMutPathway_test <- dplyr::full_join(dataMutPathway_test_MutantPathway, dataMutPathway_test_nonMutantPathway, by = c("Pathway", "Outcome6M"))
dataMutPathway_test <- dataMutPathway_test %>% dplyr::distinct() %>% dplyr::filter(Outcome6M == "responder")


# Calculate Pvalues
dataMutPathway_test$pValue <- 1

for (iPathway in c(1:nrow(dataMutPathway_test))) {
  challenge.df = matrix(c(dataMutPathway_test$responseRate_1[iPathway], dataMutPathway_test$responseRate_2[iPathway],
                          dataMutPathway_test$groupSize1[iPathway]-dataMutPathway_test$responseRate_1[iPathway], dataMutPathway_test$groupSize2[iPathway]-dataMutPathway_test$responseRate_2[iPathway]), nrow = 2)
  
  dataMutPathway_test$pValue[iPathway] <- fisher.test(challenge.df)$p.value
}

# Adjust pValues for multiple testing
dataMutPathway_test <- dataMutPathway_test %>% dplyr::mutate(pAdj = p.adjust(pValue, method = "BH")) %>% dplyr::ungroup()


# get proportion of responder per pathway
dataMutPathway_barPlots <- dataMutPathway %>% dplyr::distinct(Pathway, isMutant, Outcome6M, responseRate, categorySize) %>%
  dplyr::mutate(responseRatePerc = 100*responseRate/categorySize,
                Pathway = factor(Pathway, levels = rev(levels(dataMutPathway$Pathway))),
                isMutant_aux = paste0(isMutant, ' (n=', categorySize, ')'))




# plot histograms of response rates with 95% CI from binomial distribution
plot_ResponseRate_Pathway <- ggplot(dataMutPathway_barPlots, aes(x=isMutant, y=responseRatePerc, fill = Outcome6M)) +
  geom_bar(stat='identity', color = 'black', width = 0.8) +
  facet_wrap(~Pathway, scales = "fixed", nrow = 1) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%')) +
  labs(y = 'Response rate', x = NULL) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c('0%', '25%', '50%','75%', '100%'))  +
  scale_fill_manual(guide = guide_legend(title = 'Response to treatment', nrow = 1,
                                         title.position = 'top', title.hjust = 0.5, keywidth = 0.5, keyheight = 0.5),
                    name = NULL, values = c('responder' = '#06A506', 'non-responder' = '#E2C054'),
                    labels = c('responder' = 'Responder', 'non-responder' = 'Non-responder')) +
  theme(
    legend.position = 'bottom',
    text=element_text(size=8, family='Helvetica'),
    plot.title = element_text(size=8, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )




# all plots in one Figure
samplePlot <- cowplot::plot_grid(
  plot_Responders,
  plot_mutsPerMb,
  plot.TMB,
  plot_pathwayAlterations,
  
  ncol = 1, align = 'v', axis = 'tblr',
  rel_heights =  c(0.1,0.4, 0.1, 1)
)

# Plot sample-level overview.
legend_samplePlot <- cowplot::plot_grid(
  legend_plot_Responders,
  legend_plot_mutsPerMb,
  legend_plot.TMB,
  legend_plot_pathwayAlterations,
  
  ncol = 1, align = 'vh',  axis = 'tblr',
  rel_heights =  c(2, 1.5, 1.5, 1)
)

# export
pdf(paste0(odir,"pathwayMutations.pdf"), width = 6.5, height = 4.8)
cowplot::plot_grid(
  cowplot::plot_grid(samplePlot, legend_samplePlot, nrow = 1,
                     rel_widths = c(1, 0.4), align = 'hv', axis = 'tblr'),
  plot_ResponseRate_Pathway, ncol = 1, rel_heights = c(1, 0.9), align = 'v', axis = 'tblr')
dev.off()







