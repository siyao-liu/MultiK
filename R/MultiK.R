#' MultiK main algorithm
#'
#' MultiK main algorithm: takes a preprocessed gene expression matrix as input. Then subsamples 80\% of the cells and applies standard Seurat pipeline on the subsampled data matrix 100 times over multiple resolution parameters.
#' @param seu A Seurat object with normalized count
#' @param resolution A vector Seurat resolution parameters. Default is from 0.05 to 2 with step size of 0.05
#' @param nPC Number of principal components to use in clustering
#' @param reps Number of subsampling runs. Integer value. Default is 100
#' @param pSample Proportion of cells to sample. Numerical value. Default is 0.8
#' @param seed Optional numerical value. This sets a random seed for generating reproducible results
#' @return A list with components: k is a vector of number of runs for each K. clusters is a list containing the clustering labels for each subsampling run at each resolution parameter. consensus is a list containing a consensus matrix for each K.
#' @export
MultiK <- function(seu, resolution = seq(0.05, 2, 0.05), nPC = 30, reps = 100, pSample = 0.8, seed = NULL) {
  # setting seed for reproducibility
  if (is.null(seed) == TRUE) {
    seed <- timeSeed <- as.numeric(Sys.time())
  }
  set.seed(seed)

  # step 1: subsampling
  subcol <- list()
  for (i in 1: reps) {
    subcol[[i]] <- sample(round(ncol(seu) * pSample))
  }

  # step 2: loop over subsampling runs, with each run subsampling 80% of cells, reselect genes for clustering
  clusters <- list()
  messages <- c()
  ks <- c()
  count <- 1

  suppressPackageStartupMessages(library(Seurat))

  for(i in 1: reps) {

    print(paste("Rep: ", i, sep=""))
    # subsample the columns (the cells) from the full matrix
    subX <- seu[, subcol[[i]] ]

    # normalizing the data
    #subX <- NormalizeData(object = subX, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

    # Find HVG genes ~ 2000 genes
    subX <- FindVariableFeatures(object = subX, selection.method = "vst", nfeatures = 2000,
                                 loess.span = 0.3, clip.max = "auto",
                                 num.bin = 20, binning.method = "equal_width", verbose = F)

    # Scaling unwanted variation
    all.genes <- rownames(x = subX)
    subX <- ScaleData(object = subX, features = all.genes, verbose=F)
    # Run PCA to reduce dimensions
    subX <- RunPCA(object = subX, features = VariableFeatures(object = subX), npcs = 50, verbose=F)
    # Run Clustering
    subX <- FindNeighbors(object = subX,
                          k.param = 20, # default is 20-nearest neighbors
                          reduction = "pca", dims = 1: nPC, verbose=F)

    for (res in resolution) {
      print(paste("Rep", i, "Res", res, sep=" "))
      subX <- FindClusters(subX, resolution = res, verbose = F)
      subX.clusters <- Idents(subX)
      clusters[[count]] <- subX.clusters
      messages <- c(messages, paste("Rep_", i, "_res_", res, sep = ""))
      count <- count + 1
      ks <- c(ks, length(unique(subX.clusters)))
    }
    names(clusters) <- messages

  }

  # step 3: calculate consensus matrix across subsampling runs for each unique K
  mInit <- matrix(0, ncol = ncol(seu), nrow = ncol(seu))

  ml <- list()
  res <- list()
  all.clusters.by.K <- list()
  m.count <- list()
  unique.ks <- unique(ks)[order(unique(ks))]

  count.k <- 1
  for(k in unique.ks) {
    print(paste("k =", k, sep=" "))
    idx <- which(ks == k)
    cluster.k <- clusters[idx]
    all.clusters.by.K[[count.k]] <- cluster.k

    for (s in 1: length(cluster.k) ) {
      print(paste("run", s, sep = ""))
      sampleKey <- as.numeric(sapply(names(cluster.k[[s]]), function(x){which(colnames(seu) == x)}))
      if (s == 1){
        ml[[count.k]] <- connectivityMatrix(cluster.k[[s]], mInit, sampleKey)
        m.count[[count.k]] <- connectivityMatrix(rep(1, length(sampleKey)), mInit, sampleKey)
      }else{
        ml[[count.k]] <- connectivityMatrix(cluster.k[[s]], ml[[count.k]], sampleKey)
        m.count[[count.k]] <- connectivityMatrix(rep(1, length(sampleKey)), m.count[[count.k]], sampleKey)
      }
    }

    res[[count.k]] <- triangle(ml[[count.k]], mode = 3)/triangle(m.count[[count.k]], mode = 3)
    res[[count.k]][which(triangle(m.count[[count.k]], mode = 3) == 0)] = 0
    print(paste(k, " finished", sep = ""))
    count.k <- count.k + 1
  }

  return(list("consensus" = res, "k" = ks))
}
