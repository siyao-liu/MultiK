#' Perform SigClust tests pairwisely on the terminals of the dendrogram
#'
#' Generate a dendrogram of cluster centroids on the left, with nodes indicating the significance of the P-values from the pairwise SigClust tests. On the right is a heatmap of P-values from the pairese SigClust tests.
#' @param seu A Seurat object
#' @param clusters A vector of clustering labels
#' @return A matrix of P values
#' @export
CalcSigClust <- function(seu, clusters) {
  suppressPackageStartupMessages(library(sigclust))
  seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
  seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000,
                              loess.span = 0.3, clip.max = "auto",
                              num.bin = 20, binning.method = "equal_width", verbose = F)
  hvg <- VariableFeatures(object=seu)
  norm.hvg <- seu@assays$RNA@data[hvg, ]
  ClustAssign <- as.character(clusters)
  n <- length(unique(ClustAssign))
  pval <- matrix(NA, ncol = n, nrow = n)
  rownames(pval) <- c(paste("cluster", c(0: (n-1)), sep = ""))
  colnames(pval) <- c(paste("cluster", c(0: (n-1)), sep = ""))

  # set up data to run in SigClust
  x.list <- list()
  l.list <- list()
  for (i in 1: n) {
    x.list[[paste0("cluster", i-1)]] <- as.matrix(norm.hvg[, ClustAssign == i-1])
    l.list[[i]] <- ClustAssign[ClustAssign == i-1]
  }

  # run SigClust
  for (i in 1: (n-1)) {
    for (j in (i+1): n) {
      pval[i, j] <- RunSigClust(x1 = x.list[[i]],
                                x2 = x.list[[j]],
                                l1 = l.list[[i]],
                                l2 = l.list[[j]])@pval
    }
  }
  return(pval)
}
