#########################################
#Helper Functions
###########################################

###########################
# connectivityMatrix() #
###########################
connectivityMatrix = function( clusterAssignments, m, sampleKey) {
  ## input: named vector of cluster assignments, matrix to add connectivities
  ## output: connectivity matrix
  names( clusterAssignments ) <- sampleKey
  cls <- lapply( unique( clusterAssignments ),
                 function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId

  for ( i in 1: length( cls ) ) {
    nelts <- 1: ncol( m )
    ## produces a binary vector
    cl <- as.numeric( nelts %in% cls[[i]] )
    ## product of arrays with * function;
    ## with above indicator (1/0) statement updates all cells to indicate the sample pair was observed in the same cluster;
    updt <- outer( cl, cl, FUN = "*" )
    # same as %*%
    #updt <- cl %*% t(cl)
    m <- m + updt
  }

  return(m)
}

#######################
# triangle() #
#######################
triangle = function(m, mode = 1) {
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m


  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half

  fm = t(nm)+nm
  diag(fm) = diag(m)

  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]

  if(mode==1){
    return(vm) #vector
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }

}


##########################
# calculate rPAC score #
##########################
CalcPAC <- function(x1 = 0.1, x2 = 0.9, # threshold defining the intermediate sub-interval
                    xvec, ml) {

  PAC <- rep(NA, length(xvec))
  names(PAC) <- xvec #

  prop_zeroes <- rep(NA, length(xvec))
  names(prop_zeroes) <- xvec

  rPAC <- rep(NA, length(xvec))
  names(rPAC) <- xvec #

  for(i in 1: length(xvec)) {
    M <- ml[[i]]
    Fn <- ecdf(M[lower.tri(M)])
    # calculate PAC
    PAC[i] <- Fn(x2) - Fn(x1)
    # calculate proportion of zeroes in M
    prop_zeroes[i] <- sum(M==0)/length(M)
    # calculate relative PAC
    rPAC[i] <- PAC[i]/(1-prop_zeroes[i])
  }
  return(list("min PAC" = which.min(PAC),
              "min rPAC" = which.min(rPAC),
              "PAC" = PAC,
              "rPAC" = rPAC))
}


############################
# MultiK main algorithm #
############################
MultiK = function(seu, resolution = seq(0.05, 2, 0.05), nPC = 30, reps = 100, pSample = 0.8, seed = NULL) {

  if (is.null(seed) == TRUE) {
    seed <- timeSeed <- as.numeric(Sys.time())
  }
  set.seed(seed)

  # step 1: subsampling
  subcol <- list()
  for (i in 1: reps) {
    subcol[[i]] <- sample(x=ncol(seu), size=round(ncol(seu) * pSample), replace=FALSE)
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
    subX <- NormalizeData(object = subX, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

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


#########################
# Find optimal K #
#########################
findOptK = function(tog) {

  tog.f <- tog[tog$Freq > 100 | tog$Freq ==100, ]
  hpts <- chull(tog.f[, c("one_minus_rpac", "Freq")]) # in clockwise order
  hpts <- c(hpts, hpts[1])
  ch.df <- tog.f[hpts, ]

  df <- ch.df[ , c("ks", "one_minus_rpac", "Freq")]
  colnames(df) <- c("k", "x", "y")
  b <- c()
  end_points <- c()
  for (i in 1: (nrow(df)-1)) {
    end_points[i] <- paste(as.character(df[i, ]$k), as.character(df[(i+1),]$k), sep="-")
    b[i] <- (df[(i+1), ]$y - df[i, ]$y)/(df[(i+1), ]$x - df[i, ]$x)
  }

  # put in data frame for each line segment
  lineseg.df <- data.frame("end_points"=end_points, "slope"=b)
  lineseg.df$p1 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 1]
  lineseg.df$p2 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 2]

  # step 1: find k with largest # of runs
  which.k <- as.character(ch.df[which.max(ch.df$Freq), ]$ks)

  # step 2: see if a line segment with negative slope coming out
  if ( all(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope > 0) ) {
    optK <- which.k
  }
  else {

    # follow the line segment with the negative slope
    tmp <- which(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope < 0)
    tmp <- lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ][tmp, ]
    which.k2 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, tmp$p2)!=which.k)])

    # check if the slope becomes more negative
    lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k, ]

    if ( #any(lineseg.df[lineseg.df$p1==which.k2 | lineseg.df$p2==which.k2, ]$slope > 0)
      lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope > tmp$slope ) {
      optK <- c(which.k, which.k2)
    }

    else {
      tmp <- which(lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ]$slope < 0)
      tmp <- lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ][tmp, ]
      which.k3 <- as.character(c(tmp$p1, tmp$p2)[ which(c(tmp$p1, tmp$p2)!=which.k & c(tmp$p1, tmp$p2)!=which.k2)])

      lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k
                                   & lineseg.df$p1!=which.k2 & lineseg.df$p2 !=which.k2, ]

      if ( lineseg.df.sub[lineseg.df.sub$p1==which.k3 | lineseg.df.sub$p2 == which.k3, ]$slope > tmp$slope ) {
        optK <- c(which.k, which.k2, which.k3)
      }
      else {
        optK <- c(which.k, which.k2, which.k3)
      }
    }
  }
  return(optK)
}


###########################
# MultiK diagnositc plot #
##########################
DiagMultiKPlot = function(ks, res) {

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(cowplot))

  # set a data frame for the data to plot
  tog <- as.data.frame(table(ks)[table(ks) > 1])

  # calculate rpac
  pacobj <- CalcPAC(x1=0.1, x2=0.9, xvec = tog$ks, ml = res)
  tog$rpac <- pacobj$rPAC
  tog$one_minus_rpac  <- 1-tog$rpac

  # Plot Freq of # of runs
  freqPlot <- ggplot(data=tog, aes(x=ks, y=Freq)) +
    geom_bar(stat="identity") +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    geom_text(aes(label=Freq),  vjust=-0.3, size=3.5) +
    scale_x_discrete("K") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black")

  # Plot rPAC for each K
  rpacPlot <- ggplot(data=tog, aes(x=ks, y=rpac,group=1)) +
    geom_point(shape=21, color="black", fill="black", size=2) +
    geom_line() +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_discrete("K") +
    scale_y_continuous("rPAC")

  # Plot (1-rPAC) Vs freq for each K
  # first find optimal K using the convex hull
  optK <- findOptK(tog)
  cat("Optimal K: ", optK)

  scatPlot <- ggplot(data=tog, aes(x=one_minus_rpac, y=Freq)) +
    geom_point(shape=21, color="black", fill="black", size=1.5) +
    geom_path(color="grey", alpha=0.75, linetype=2) +
    theme_bw(base_size=14)+
    theme(axis.text=element_text(color="black"),
          panel.background=element_rect(color="black"),
          strip.text = element_text(size=12),
          strip.background = element_rect(fill="white")) +
    scale_x_continuous("1 - rPAC") +
    scale_y_continuous("Number of clustering runs") +
    geom_hline(yintercept=100, linetype="dashed", color = "black") +
    geom_label_repel(aes(label = ks), segment.color = 'grey50', size=3)  +
    geom_path(data=tog[match(findOptK(tog), tog$ks), ])

  plot_grid(freqPlot, rpacPlot, scatPlot, ncol=3)


}


#######################################################################
# Find the resolution parameter that gives optimal K in the full data #
#######################################################################
findResol <- function(seu, ks, optK) {
  L <- as.numeric(gsub("RNA_snn_res.", "", names(which(ks < optK))))
  R <- as.numeric(gsub("RNA_snn_res.", "", names(which(ks > optK)[1])))
  k_L <- ks[which(ks < optK)]
  k_L
  k_R <- ks[which(ks > optK)[1]]
  k_R
  repeat {
    mid = L+(R-L)/2
    seu <- FindClusters(seu, resolution = mid, verbose = F)
    mid_k <- length(table(seu@meta.data$seurat_clusters))
    mid_k
    if (mid_k < optK) {
      # use mid as L
      L <- mid
    }
    else {
      # use mid as R
      R <- mid
    }

    if (mid_k==optK)
      break
  }
  return(mid[length(mid)])
}

#####################################
# get cluster labels at Optimal K #
#####################################
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
  return(list("clusters"=clusters, "resolution"=res))
}


#######################
# Run SigClust #
#######################
RunSigClust <- function(x1, x2, l1, l2) {
  suppressPackageStartupMessages(library(sigclust))
  sig.dat <- as.matrix(t(cbind(x1, x2)))
  l1 <- rep(1, length(l1))
  l2 <- rep(2, length(l2))
  sig.label <- c(l1, l2)
  sig <- sigclust(x = sig.dat,
                  nsim = 100,
                  nrep = 1, labflag = 1, # uses known label
                  label = sig.label,
                  icovest = 2)

  return(sig)
}

#########################################################
# Calculate SigClust pvalues following the dendrogram #
#########################################################
CalcSigClust = function(seu, clusters) {

  suppressPackageStartupMessages(library(Seurat))
  seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
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

#####################################
# Plot PW SigClust p value heatmap #
#####################################
# only plot the lower triangle, set upper triangle NA and plot in grey
PWSig_Heatmap <- function(pval, order) {
  suppressPackageStartupMessages(library(pheatmap))
  # make pvalue matrix symmetric (let lower tri = upper tri)
  m <- as.matrix(pval)
  tm <- as.matrix(t(pval))

  mm <- matrix(NA, nrow = dim(m)[[1]], ncol = dim(m)[[1]])
  mm[upper.tri(mm)] <- m[upper.tri(m)]
  mm[lower.tri(mm)] <- tm[lower.tri(tm)]
  rownames(mm) <- rownames(m)
  colnames(mm) <- colnames(m)
  mm[is.na(mm)] = 1

  mat <- mm[order, order]
  mat[upper.tri(mat)] <- NA
  # put diagonals as NA as well (no need to show a cluster has pval=1 with itself)
  diag(mat) <- NA

  # create cell note matrix
  mlabel <- mat
  mlabel[is.na(mlabel)] <- ""

  # change color scale - any value above 0.05 is
  # create some vectors for breaks and color in the heatmap
  colors <- c(seq(0, 0.05, length=100), seq(0.05001, 1, length=100))
  my_palette <- colorRampPalette(c("white",  "blue"))(n = 199)

  pheatmap(mat,
           cluster_rows = F, cluster_cols = F,
           fontsize = 12, clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_rownames = TRUE,
           show_colnames = TRUE,
           display_numbers = mlabel,
           number_color = "black",
           col = rev(my_palette), main="Pairwise SigClust p values",
           breaks = colors,
           na_col = "grey")

}

# create my own function to plot dendrogram
PlotDendro <- function(dend, nodes_shapes) {
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(dendextend))
  suppressPackageStartupMessages(library(ggdendro))
  suppressPackageStartupMessages(library(ggplot2))

  dend_data <- dendro_data(dend)
  # Setup the data, so that the layout is inverted (this is more
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data),
    data.frame(x = y, y = x, xend = yend, yend = xend))

  # get the nodes position using the horizontal lines in ggdendro plot
  nodes_sub <- segment_data[which(segment_data$y==segment_data$yend), ]
  nodes_sub <- nodes_sub[nodes_sub$xend!=0, ]
  # manually add the top node
  nodes_pos <- nodes_sub[, c("xend", "yend")]

  topn_x <- max(nodes_sub$x)
  if ( length(nodes_sub[nodes_sub$x==max(nodes_sub$x), ]$y) > 1) {
    topn_y <- mean((nodes_sub[nodes_sub$x==max(nodes_sub$x), ]$y -1))+1
  }
  else {
    topn_y <- (nodes_sub[nodes_sub$x==max(nodes_sub$x), ]$y -1)/2 + 1
  }

  nodes_pos <- rbind(nodes_pos, c(topn_x, topn_y))
  nodes_pos <- nodes_pos[order(nodes_pos$xend, decreasing=TRUE), ]

  plt_dendr <- ggplot() +
    geom_segment(data=segment_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() +
    labs(x = "", y = "", colour = "", size = "") +
    geom_text(aes(x = y, y = x, label = label),
              angle = -90, hjust = 0,
              data= label(dend_data)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    geom_point(data=nodes_pos, aes(x=xend, y=yend),
               shape=nodes_shapes,
               size=3) + ggtitle("Cluster Centroids Dendrogram")
  print(plt_dendr)
}

##########################
# assign nodes shapes #
##########################
assign_nodes_shapes = function(hc, pval) {
  A = hc$merge
  labelList = c()
  for (i in seq(1, nrow(A))){
    if((A[i,1]<0) & (A[i,2]<0)){
      labelList[[i]] = hc$labels[c(-A[i,])]
    }
    else if((A[i,1]<0) & (A[i,2]>0)){
      labelList[[i]] = c( hc$labels[-A[i,1]], labelList[[A[i,2]]])
    }
    else if((A[i,1]>0) & (A[i,2]<0)){
      labelList[[i]] = c(hc$labels[-A[i,2]], labelList[[A[i,1]]])
    }
    else if((A[i,1]>0) & (A[i,2]>0)){
      labelList[[i]] =  c(labelList[[A[i,1]]], labelList[[A[i,2]]])
    }
  }

  pch <- list()
  for (g in 1: length(labelList)) {
    cluster_names <- paste0("cluster", labelList[[g]])
    cluster_names
    tmp <- as.vector(pval[ which(colnames(pval) %in% cluster_names), which(colnames(pval) %in% cluster_names)])
    tmp
    if ( any(tmp[!is.na(tmp)] < 0.05) ) {
      pch[[g]] = 19
    } else {
      pch[[g]] = 1
    }
  }
  return(rev(pch))
}

########################################################################
# Plot Cluster Mean dendrograms and overlay pairwise SigClust pvalues #
########################################################################
PlotSigClust = function(seu, clusters, pval) {

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

  # plot pairwise p value heatmap
  hp <- PWSig_Heatmap(pval=pval, order=hc$order)

  # plot cluster mean dendrogram
  nodes_shapes <- assign_nodes_shapes(hc, pval)
  dend <- as.dendrogram(hc)
  plt_dendr <- PlotDendro(dend, nodes_shapes)

  lm=rbind(c(1,2),
           c(1,2))

  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  grid.arrange(grobs = list(plt_dendr, hp[[4]]), layout_matrix = lm)
}




