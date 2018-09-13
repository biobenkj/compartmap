#' @title Estimate A/B compartments from WGBS data
#'
#' @description 
#' \code{getWGBSABsignal} returns estimated A/B compartments from whole genome bisulfite sequencing data
#'
#' @details 
#'
#' @param mat 
#' @param res 
#' @param globalMeanSet 
#' @param noMean 
#' @param targets 
#' @param parallel 
#' @param allchrs 
#' @param chr 
#' @param regions 
#' @param genome 
#' @param preprocess 
#' @param ... 
#'
#' @return compartment estimates
#' @export
#'
#' @examples

getWGBSABsignal <- function(mat, res=1e6, globalMeanSet = NULL, noMean = FALSE, targets = NULL, parallel=FALSE, allchrs=FALSE, chr = NULL, regions = NULL, genome = "hg19", preprocess = TRUE, ...) {
  
  # Run preprocessing if specified
  # Note: Dirichlet smoothing uses Jeffrey's prior as a value for k
  if (preprocess == TRUE) mat <- .preprocess_wgbs(mat, k = 0.5, minCov = 3, minSamp = 2)
  
  columns <- colnames(mat)
  names(columns) <- columns
  
  if (any(is.na(mat))) stop("The input matrix cannot have any NA values in it. Either impute the missing values or remove loci/rows that contain NAs.")
  if (allchrs == TRUE & !is.null(chr)) stop("You've requested genome wide analysis in addition to specific chromosomes. Set chrs = 'all' for all chromosomes.")
  if (allchrs == FALSE & is.null(chr)) stop("You need to specify which chromosomes to use for analysis...")
  if (is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeans(mat, targets)
  
  #Calculate things once... 
  if (is.null(regions)) regions <- .buildRanges(mat)
  
  #This assumes human... 
  #TODO: fix this to support more organisms...
  if (allchrs == TRUE) chr <- paste0("chr", c(seq(1,22)))
  
  #If we are squeezing to a global mean
  if (parallel & noMean == FALSE) {
    do.call(cbind, 
            mclapply(columns,.getWGBSShrunkenCompartments,mat=mat,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr))
  } else if (noMean == FALSE) { 
    do.call(cbind, 
            lapply(columns,.getWGBSShrunkenCompartments,mat=mat,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr))
  } else {
    #Not squeezing to global mean
    .getWGBSShrunkenCompartments(mat=mat,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr)
  } 
}

.getGlobalMeans <- function(mat, targets = NULL) {
  if (!is.null(targets)){
    stargets <- .getShrinkageTargets(mat, targets)
    message(paste0("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets..."))
    meanBeta <- matrix(rowMeans(stargets, na.rm=TRUE), ncol=1)
  }
  else (meanBeta <- matrix(rowMeans(mat, na.rm=TRUE), ncol=1))
  colnames(meanBeta) <- "globalMean"
  return(meanBeta) 
}

#Subset a group of samples for binning
.getShrinkageTargets <- function(mat, group) {
  if (all(group %in% colnames(mat))) stargets.mat <- mat[,group]
  else (stop(paste0("Could not find ", group, " in the colnames of the input matrix...")))
  return(stargets.mat)
}


#A helper function to reconstruct the GRanges object for binning if not supplied from the input matrix
#These *should* be the rownames of the input matrix post biscuiteer
.buildRanges <- function(mat) {
  if (!is.null(rownames(mat))) {
    message("Constructing rowRanges from the input matrix...")
    snames <- sapply(strsplit(rownames(mat), ":"), `[`, 1)
    loc <- as.numeric(sapply(strsplit(rownames(mat), ":"), `[`, 2))
    myranges <- GRanges(seqnames = snames,
                        strand = "*",
                        ranges = IRanges(start = loc, width = 1))
  }
  else (stop("Could not reconstruct the rowRanges from the input matrix..."))
  return(myranges)
}

#Realize and coerce to a matrix from bsseq object
#The ... are additional arguments to pass to getMvals or getLogitFracMeth
.coerceBSseq <- function(bsseq, k = 0.5, minCov = 3, minSamp = 2) {
  if (is(bsseq, "BSseq")) {
    message("Assuming this hasn't been smoothed yet since it is still a BSseq object...")
    message("Smoothing...")
    bsseq.smooth <- getLogitFracMeth(bsseq, k = k, minCov = minCov, minSamp = minSamp)
    message("Realizing DelayedArray object...")
    message("This could take a little while...")
    bsseq.smooth <- realize(bsseq.smooth)
    message("Coercing to a matrix for compartment inference...")
    bsseq.mat <- as.matrix(bsseq.smooth)
  }
  else if (is(bsseq, "matrix")) {
    message("Realizing DelayedArray object...")
    message("This could take a little while...")
    #Add a try statement here in case it's just a standard matrix at this point
    bsseq <- try(realize(bsseq))
    message("Trying to coerce to a matrix for compartment inference...")
    bsseq.mat <- as.matrix(bsseq)
  }
  else (stop(paste0("Not sure what do with ", class(bsseq))))
  return(bsseq.mat)
}

#Filter and impute missing values prior to binning
#Binning WILL fail with NAs
#The ... are other arguments you can pass to impute.knn
.filterAndImpute <- function(mat, fracNA = 0.5, k = 0.1, minCov = 3, minSamp = 2) {
  if (is(mat, "BSseq")) mat <- .coerceBSseq(mat, k = k, minCov = minCov, minSamp = minSamp)
  message(paste0("Removing loci with greater than ", fracNA*100, "% NA values prior to knn imputation..."))
  mat.clean <- mat[rowMeans(is.na(mat)) < fracNA,]
  message("Starting knn imputation...")
  #This will spew all kinds of output but some users like seeing things are working I suppose...
  mat.clean.imputed <- impute.knn(mat.clean)$data
  return(mat.clean.imputed)
}



.getWGBSShrunkenCompartments <- function(column, mat, regions = NULL, globalMeanSet=NULL, noMean = FALSE, res=1e6, targets = NULL, chr = NULL, genome = "hg19", ...) {
  if (is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeans(mat, targets)
  chrs <- chr
  names(chrs) <- chrs
  if (!is.null(globalMeanSet) & isFALSE(noMean)) {
    getMeanPairedChr <- function(chr) { 
      message("Computing shrunken eigenscores for ", column, " on ", chr, "...") 
      binmat <- getBinMatrix(cbind(mat[,column], globalMeanSet), regions, chr = chr, res = res, FUN = mean, genome = genome)
      cormat <- getCorMatrix(binmat)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cormat$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cormat$binmat.cor)))
      }
      else {
        absig <- getABSignal(cormat)$pc
      }
      
      return(absig)
    }
    unlist(lapply(chrs, getMeanPairedChr))
  }
  else if (isTRUE(noMean)) {
    getSmoothedPairedChr <- function(chr) { 
      message("Computing eigenscores on ", chr, "...") 
      binmat <- getBinMatrix(mat, regions, chr = chr, res = res, FUN = mean, genome = genome)
      cormat <- getCorMatrix(binmat)
      #Stupid check for NAs
      if (any(is.na(cormat$binmat.cor))) {
        warning("NAs were found in the correlation matrix. Eigenvector may not handle NAs.")
        absig <- try(getABSignal(cormat)$pc)
      }
      else {
        absig <- getABSignal(cormat)$pc
      }
      
      return(absig)
    }
    unlist(lapply(chrs, getSmoothedPairedChr))
  }
}

.preprocess_wgbs <- function(mat, k = 0.5, minCov = 3, minSamp = 2, ...) {
  mat <- .filterAndImpute(mat, k = k, minCov = minCov, minSamp = minSamp)
  return(mat)
}