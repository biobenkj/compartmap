#' Denoise a correlation matrix using a Random Matrix Theory approach
#'
#' @name getDenoisedMatrix
#'
#' @param obj SummarizedExperiment object with rowRanges for each feature and colnames
#' @param res The resolution desired (default is a megabase 1e6)
#' @param chr Which chromosome to perform the denoising
#' @param genome Which genome (default is hg19)
#' @param iter How many iterations to perform denoising
#' @param targets Samples/cells to shrink towards
#' @param prior.means The means of the bin-level prior distribution (default will compute them for you)
#' @param assay What assay type this is ("rna", "atac", "bisulfite", "array")
#'
#' @return A denoised correlation matrix object for plotting with plotCorMatrix
#' 
#' @import covmat
#' @import doParallel
#' @import SummarizedExperiment
#' @import GenomicRanges
#' 
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' denoised_cor_mat <- getDenoisedCorMatrix(k562_scrna_chr14, genome = "hg19", assay = "rna")

getDenoisedCorMatrix <- function(obj, res = 1e6, chr = "chr14",
                                 genome = c("hg19", "hg38", "mm9", "mm10"),
                                 iter = 2, targets = NULL, prior.means = NULL,
                                 assay = c("rna", "atac",
                                           "bisulfite", "array")) {
  ## this is a wrapper to give back a denoised correlation matrix to plot
  #match the assay args
  assay <- match.arg(assay)
  
  #match the genome if given
  genome <- match.arg(genome)
  
  #double check the obj class is compatible
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")
  
  #subset to selected chromosome(s)
  obj <- keepSeqlevels(obj, chr, pruning.mode = "coarse")
  
  #get the prior means
  if (is.null(prior.means)) {
    prior.means <- getGlobalMeans(obj=obj, targets=targets,
                                  assay=assay)
  } else {
    message("Assuming the prior means passed were derived from the full sample set.")
  }
  
  #shrink bins
  message("Shrinking bins with the JSE.")
  bin.mat <- shrinkBins(x = obj, original.x = obj, prior.means = prior.means,
                        chr = chr, res = res, targets = targets,
                        jse = TRUE, assay = assay, genome = genome)
  #get the raw correlation matrix
  cor.mat <- getCorMatrix(binmat = bin.mat, squeeze = FALSE)
  #denoise with RMT
  message("Denoising the correlation matrix using RMT.")
  cor.mat.denoise <- covmat::estRMT(cor.mat$binmat.cor,
                                    parallel = FALSE)$cov
  #iterate?
  if (iter >= 2) {
    for (i in 2:iter) {
      message("Iterative denoising. Iteration: ", i)
      cor.mat.denoise <- covmat::estRMT(cor.mat.denoise,
                                        parallel = FALSE)$cov
    }
  }
  
  #rescale
  cor.mat.denoise <- scales::rescale(cor.mat.denoise,
                                     to = c(0,1))
  #tidy up
  colnames(cor.mat.denoise) <- rownames(cor.mat.denoise) <- as.character(granges(cor.mat$gr.cor))
  
  return(cor.mat.denoise)
}

#' Plot a denoised correlation matrix
#'
#' @name plotCorMatrix
#'
#' @param denoised.cor.mat The denoised correlation matrix object from getDenoisedMatrix
#' @param midpoint The midpoint for the coloring (default is 0.3)
#' @param return.plot.obj Whether to return the ggplot object
#' @param uppertri Whether to keep the upper triangle of the matrix
#' @param lowertri Whether to keep the lower triangle of the matrix
#'
#' @return Either a ggplot object or plot
#' 
#' @import ggplot2
#' @import reshape2
#' 
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(10000), ncol=25)
#' set.seed(1000)
#' my_plot <- plotCorMatrix(dummy, return.plot.obj = TRUE)

plotCorMatrix <- function(denoised.cor.mat,
                          midpoint = 0.3,
                          return.plot.obj = FALSE,
                          uppertri = FALSE,
                          lowertri = FALSE) {
  ## upper tri
  if (uppertri) {
    denoised.cor.mat[upper.tri(denoised.cor.mat)] <- 0
  }
  if (lowertri) {
    denoised.cor.mat[lower.tri(denoised.cor.mat)] <- 0
  }
  diag(denoised.cor.mat) <- 1
  ## melt for plotting
  cor.mat.melt <- reshape2::melt(denoised.cor.mat)
  ## plot
  p <- ggplot(cor.mat.melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "white", mid = "white", high = "red3", midpoint = midpoint) +
    theme_minimal() +
    theme(
      legend.position = "None",
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.margin= grid::unit(c(0, 0, 0, 0), "in")
      )
  
  if (return.plot.obj) {
    return(p)
  } else {
    p
  }
}
