#' Plot cluster centroids dendrogram and overlay P-values from pairwise SigClust tests
#'
#' Generate a dendrogram of cluster centroids on the left, with nodes indicating the significance of the P-values from the pairwise SigClust tests. On the right is a heatmap of P-values from the pairese SigClust tests. The assignment of classes and subclasses is determined by the shape of the nodes. Classes are indicated by closed circles, and subclasses are indicated by open circles.
#' @param seu A Seurat object
#' @param clusters A vector of clustering labels
#' @param pval A matrix of P values from pairwise SigClust tests
#' @return A diagnostic plot that shows the SigClust results
#' @export
PlotSigClust <- function(seu, clusters, pval) {

  seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000,
                              loess.span = 0.3, clip.max = "auto",
                              num.bin = 20, binning.method = "equal_width", verbose = FALSE)
  hvg <- VariableFeatures(object=seu)
  norm.hvg <- seu@assays$RNA@data[hvg, ]
  ClustAssign <- as.character(clusters)
  tmp.list <- list()
  for (i in names(table(ClustAssign))) {
    tmp.list[[i]] <- rowMeans(as.matrix(norm.hvg[, ClustAssign == i]))
  }
  clustMean.mat <- do.call("cbind", tmp.list)

  # Run hierarchical clustering on the cluster means
  hc <- hclust(as.dist(1 - cor(clustMean.mat)))

  suppressPackageStartupMessages(library(scales))
  # plot pairwise p value heatmap
  hp <- PWSig_Heatmap(pval=pval, order=hc$order)

  # plot cluster mean dendrogram
  nodes_shapes <- assign_nodes_shapes(hc, pval)
  dend <- as.dendrogram(hc)
  plt_dendr <- PlotDendro(dend, nodes_shapes)

  lm=rbind(c(1,2),c(1,2))

  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  grid.arrange(grobs = list(plt_dendr, hp[[4]]), layout_matrix = lm)
}
