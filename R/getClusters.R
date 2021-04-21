#' Get clustering labels at optimal K
#'
#' Find the clustering labels at the selected optimal K levels.
#' @param seu A gene expression matrix
#' @param optK A vector of selected optimal Ks
#' @return A matrix of clustering labels for each optimal K and the resolution parameters used to identifky the optimal K clusters.
#' @export
getClusters <- function(seu, optK) {

  suppressPackageStartupMessages(library(Seurat))

  # normalizing the data
  seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  # Find HVG genes ~ 2000 genes
  seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000,
                                   loess.span = 0.3, clip.max = "auto",
                                   num.bin = 20, binning.method = "equal_width", verbose = FALSE)

  # Scaling unwanted variation
  all.genes <- rownames(x = seu)
  seu <- ScaleData(object = seu, features = all.genes, verbose = FALSE)

  # Run PCA to reduce dimensions
  seu <- RunPCA(object = seu, features = VariableFeatures(object = seu), npcs = 50, verbose = FALSE)

  # Run Clustering
  k.param <- 20
  nPC <- 30
  seu <- FindNeighbors(object = seu, k.param = k.param, reduction = "pca", dims = 1: nPC, verbose = FALSE)
  resolution <- seq(0.05, 2.0, by = 0.05)
  seu <- FindClusters(seu, resolution = resolution, verbose = F)
  meta.data <- seu@meta.data[, grep("RNA_snn_res.", colnames(seu@meta.data)) ]

  ks <- apply(meta.data, 2, function(x) length(table(x)))

  clusters <- matrix(NA, nrow=ncol(seu), ncol=length(optK))
  colnames(clusters) <- paste("K", optK, sep="_")
  rownames(clusters) <- rownames(meta.data)
  res <- c()

  for (i in 1: length(optK)) {
    #print(optK[i])
    # first find out if optK exist in the range of tested resolution, if yes, use the first occurence
    if ( any(optK[i]==ks) ) {
      clusters[, i] <- as.numeric(as.character(meta.data[, which(ks==optK[i])[1]]))
      res[i] <- gsub("RNA_snn_res.", "", names(which(ks==optK[i])[1]))
    }
    else {
      # optK not exist in the range of test resolution, find the small window
      res[i] <- findResol(seu, ks, optK[i])
      seu <- FindClusters(seu, resolution = res[i], verbose = F)
      clusters[, i] <- as.numeric(as.character(seu@meta.data$seurat_clusters))
    }
  }
  names(res) <- paste("K", optK, sep="_")
  return(list("clusters"=clusters, "resolution"=res))
}
