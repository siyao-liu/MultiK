#' Make MultiK diagnostic plots
#'
#' Takes the output from MultiK() function, and generate diagnostic plots for determining optimal Ks
#' @param ks A vector of the number of clustering runs for each K
#' @param res A list of consensus matrix for each K
#' @return A bar plot of the frequency of K, a plot of rPAC score for each K, and a scatterplot of (1-rPAC) vs. the frequency of K
#' @export
DiagMultiKPlot <- function(ks, res) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))

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
