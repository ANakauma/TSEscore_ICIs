# function to call TSE groups
TSE_classify <- function (x, minCor = 0.2, centroids_TSE) 
{
  TSEclasses <- c("TSE_positive", "TSE_neutral", "TSE_negative")
  genesToKeep <- intersect(rownames(centroids_TSE), rownames(x))
  if (length(genesToKeep) == 0) {
    stop("Genes provided are not in the list of genes used for classification.\n Make sure that gene names are symbols") }
  if (length(gkeep) < 0.6 * nrow(centroids_TSE))  {
    warning("Less than 60% of the genes used for classification are present in the data provided. Results may not be relevant") }
  cor.dat <- as.data.frame(cor(x[genesToKeep, ], centroids_TSE[match(genesToKeep, rownames(centroids_TSE)), TSEclasses], use = "complete.obs"), row.names = colnames(x))
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y) {
    TSEclasses[which.max(y)]
  })
  cor.dat$corToNearest <- apply(cor.dat[, TSEclasses], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp) {
    cor.test(x[genesToKeep, smp], centroids_TSE[match(genesToKeep, rownames(centroids_TSE)), cor.dat[smp, "nearestCentroid"]])$p.value
  })
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - 
                                        cor.dat[, TSEclasses], 1, function(x) {
                                          sort(x)[2]
                                        })
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, TSEclasses], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  cor.dat$TSE_category <- cor.dat$nearestCentroid
  try(cor.dat[which(cor.dat$corToNearest < minCor), "TSE_category"] <- NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <- NA)
  cor.dat <- cor.dat[, c("TSE_category", "cor_pval", "separationLevel", TSEclasses)]
  cor.dat$sampleId <- rownames(cor.dat)
  return(cor.dat)
}