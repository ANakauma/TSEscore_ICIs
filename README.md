# T Cell-to-Stroma Enrichment (TSE) score: a gene expression metric that predicts response to immune checkpoint inhibitors (ICI) in patients with urothelial cancer

------------------------------------------------------------------------

## Description

The **TSE score** is a transcriptomic metric that predicts response to ICIs. The study reported by Rijnders, Nakauma-Gonzalez, et al., entitled [T Cell-to-Stroma Enrichment (TSE) score: a gene expression metric that predicts response to immune checkpoint inhibitors in patients with urothelial cancer]((https://doi.org/10.1101/2022.05.30.493997)) is available as preprint at bioRxiv. It has been submitted for publication and it is under review.

We identified the abundance of T cells as a predictor of response, while the abundance of stromal cells/products is a predictor of non-response to ICI. Both predictors are affected by tumor purity, but when combined into the TSE score, the effect of purity disappears making the TSE score a robust biomarker even for immune cell-riched organs such as lymph nodes (Figure). Additionally, the TSE score is able to predict responders and non-responders because it combines the positive and negative predictive values of T and stromal cells into one single metric.

![The TSE score is a robust biomarker that outperforms the tumor mutational burden (TMB). Patients with a positive TSE score live twice as long as other patients.](images/Figure_RESPONDER.jpg)

------------------------------------------------------------------------

## Calculate the TSE score

To identify the TSE score category, you will have to clone this repository. A TSE classifier is available in the `classifier` folder. Load the classifier:

```         
source("classifier/TSE_classify.R")
load("data/centroids_TSE.RData")
```

Apply the classifier on a RNA count matrix (matrixCounts; rows = gene symbols, columns = sampleId)

```         
TSEclass <- TSE_classify(x = matrixCounts, centroids_TSE = centroids_TSE)
```

The result is a table with the TSE category per sample.

If you are interested in calculating the exact value of the TSE score, which is independent of the classifier. You may follow the example provided in `5.TSEscoreCalculation_example.R`. A test data set is provided. For new data sets, they have to be a Matrix with samples as columns and gene names as rows (see test data set).

The other scripts were used to process the data provided by the Hartwig Medical Foundation (HMF) and to generate figures and results for the manuscript. Access to HMF data is controlled and must be requested under the request number DR-176. WGS data, RNA-seq data and corresponding clinical data are freely available for academic use from the HMF through standardized procedures. Request forms can be found at <https://www.hartwigmedicalfoundation.nl>.

------------------------------------------------------------------------

## Reference

Coming soon...
