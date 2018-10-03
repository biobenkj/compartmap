#' @title Estimate A/B compartments
#'
#' @description 
#' \code{getCompartments} returns estimated A/B compartments from ATAC-seq, whole genome bisulfite sequencing, and methylation array data
#'
#' @details 
#' This is a wrapper function to perform A/B compartment inference. Compartmentalizer implements a Stein estimator to shrink per-sample compartment estimates towards a global mean. The expected input for this function can be generated using packages like minfi, biscuiteer, and ATACseeker.
#'
#' @param obj The object with which to perform compartment inference
#' @param type The type of data that obj represents (e.g. atac, wgbs, or array)
#' @param res Resolution of compartments in base pairs (default is 1e6)
#' @param parallel Should the estimates be done in parallel (default is FALSE)
#' @param chrs Chromosomes to operate on (can be individual chromosomes, a list of chromosomes, or all)
#' @param shrink.targets Target samples to shrink towards (e.g. normal/control samples - default is the global mean)
#' @param regions GRanges object that contains corresponding genomic locations of the loci
#' @param genome Genome to use (default is hg19)
#' @param preprocess Should the data be preprocessed (currently only supports WGBS data)
#' @param gmean Squeeze towards a global mean? (default is TRUE)
#' @param ... Other parameters to pass to internal functions
#'
#' @return A p x n matrix (samples as columns and compartments as rows) to pass to embed_compartments
#' 
#' @import biscuiteer 
#' @import impute 
#' @import gtools 
#' @import parallel
#' @import Homo.sapiens
#' @import Mus.musculus
#' 
#' @export
#'
#' @examples
#' 

getCompartments <- function(obj, type = c("atac", "wgbs", "array"), res = 1e6, parallel = FALSE,
                             chrs = "chr1", shrink.targets = NULL, regions = NULL, genome = "hg19",
                             preprocess = FALSE, gmean = TRUE, ...) {
  
  # Perform initial check the input data type
  if (type %in% c("atac", "wgbs", "array")) {
    type <- match.arg(type)
  } else {
    stop("Don't know what to do with the supplied type argument. Can be 'atac', 'wgbs', or 'array'.")
  }
  
  #Pre-check the chromosomes to be analyzed
  if (chrs == "all") {
    allchrs <- TRUE
    chrs <- NULL
    message("Proceeding with all chromosomes...")
  } else {
    allchrs <- FALSE
    message(paste0("Proceeding with chromosome(s): ", paste(shQuote(chrs), collapse = ", ")))
  }
  
  # WGBS
  # Check whether the input object is a BSseq object for WGBS data
  if (is(obj, "BSseq") & type == "wgbs") {
    
    #Check for parallel and run A/B inference
    if (parallel == TRUE & gmean == TRUE) {
      warning("Running inference in parallel is memory hungry. No guarantee this won't error due to lack of memory.")
      compartments <- getWGBSABsignal(obj = obj, res = res, globalMeanSet = NULL, noMean = FALSE,
                                      targets = shrink.targets, parallel = TRUE, allchrs = allchrs, chr = chrs,
                                      regions = regions, genome = genome, preprocess = preprocess, ...)
    }
    # Not parallel and shrink towards global mean
    else if (parallel == FALSE & gmean == TRUE) {
      compartments <- getWGBSABsignal(obj = obj, res = res, globalMeanSet = NULL, noMean = FALSE,
                                      targets = shrink.targets, parallel = FALSE, allchrs = allchrs, chr = chrs,
                                      regions = regions, genome = genome, preprocess = preprocess, ...)
    }
    # Not shrinking towards global mean
    else {
      compartments <- getWGBSABsignal(obj = obj, res = res, globalMeanSet = NULL, noMean = TRUE,
                                      targets = shrink.targets, parallel = FALSE, allchrs = allchrs, chr = chrs,
                                      regions = regions, genome = genome, preprocess = preprocess, ...)
    }
    return(compartments)
  } else {
    stop("obj needs to be a BSseq object. Can be generated using the biscuiteer package.")
  }
  
  # ATACseq
  # Check whether the input object is a RSE object for ATACseq data
  if (is(obj, "RangedSummarizedExperiment") & type == "atac") {
    stop("Not implemented yet")
  } else {
    stop("obj needs to be a RangedSummarizedExperiment object. Can be generated using the ATACseeker package.")
  }
  
  # Methylation array (e.g. 450k or EPIC)
  # Check whether the input object is an GRset object for array data
  if (is(obj, "GenomicRatioSet") & type == "array") {
    if (allchrs == TRUE) {
      compartments <- getArrayABsignal(obj = obj, res = res, parallel = parallel, allchrs = allchrs, ...)
      } else {
      compartments <- getArrayABsignal(obj = obj, res = res, parallel = parallel, allchrs = FALSE, ...)
      }
    return(compartments)
  } else {
    stop("obj needs to be an GenomicRatioSet object. Can be generated using the minfi package.")
  }
}