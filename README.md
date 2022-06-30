# T Cell-to-Stroma Enrichment (TSE) score: a gene expression metric that predicts response to immune checkpoint inhibitors in patients with urothelial cancer

The TSE score is a transcriptomic marker that correlates with response to ICIs. The study reported by Rijnders, Nakauma-Gonzalez, et al., entitled [T Cell-to-Stroma Enrichment (TSE) score: a gene expression metric that predicts response to immune checkpoint inhibitors in patients with urothelial cancer]((https://doi.org/10.1101/2022.05.30.493997)) is available as preprint at bioRxiv. It has been submitted for publication and it is under review.  This workflow is dependent on pre-processed Hartwig Medical Foundation (HMF) data which are used to generate figures and process data for the manuscript. Access to HMF data is controlled and must be requested under the request number DR-176. WGS data, RNA-seq data and corresponding clinical data are freely available for academic use from the HMF through standardized procedures. Request forms can be found at https://www.hartwigmedicalfoundation.nl.

The `5.TSEscoreCalculation_example.R` script calculates the TSE score from normalized counts. A test data set is provided. For new data sets, they have to be a Matrix with samples as columns and gene names as rows (see test data set).


