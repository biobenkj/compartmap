#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description 
#' \code{getArrayABsignal} returns estimated A/B compartments from methylation array data.
#'
#' @details This function estimates A/B compartments shrinking towards a global mean of targets or across samples
#' 
#'
#' @param obj Input GenomicRatioSet object 
#' @param res Compartment resolution (in bp)
#' @param parallel Should the inference be done in parallel?
#' @param allchrs Whether all autosomes should be used for A/B inference
#' @param ... Additional arguments
#'
#' @return A p x n matrix (samples as columns and compartments as rows) of compartments
#' @export
#'
#' @examples
#' 

getATACABsignal <- function(obj, res=1e6, parallel=FALSE, allchrs=F, ...) {
  globalMeanSet <- .getGlobalMeans(obj)
  columns <- colnames(obj)
  names(columns) <- columns 
  
  getComp <- .getPaired
  if (allchrs == TRUE) getComp <- .getPairedAllChrs
  
  if (parallel) {
    options(mc.cores=detectCores()/2) # RAM blows up otherwise 
    do.call(cbind, 
            mclapply(columns,getComp,obj=obj,globalMeanSet=globalMeanSet))
  } else { 
    do.call(cbind, 
            lapply(columns,getComp,obj=obj,globalMeanSet=globalMeanSet))
  } 
}

.getGlobalMeans <- function(obj) { 
  meanBeta <- matrix(rowMeans(assay(obj), na.rm=TRUE), ncol=1) 
  sub_rse <- obj[,1]
  assays(sub_rse)$counts <- meanBeta
  colnames(sub_rse) <- "globalMean"
  return(sub_rse) 
}

.getPaired <- function(column, obj, globalMeanSet=NULL, res=1e6, ...) {
  message("Computing shrunken compartment eigenscores for ", column, "...") 
  if(is.null(globalMeanSet)) globalMeanSet <- getGlobalMeans(obj)
  binmat <- getBinMatrix(as.matrix(cbind(assay(obj)[,column], assay(globalMeanSet))), rowRanges(obj), chr = "chr1", res = res)
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

.getPairedAllChrs <- function(column, obj, globalMeanSet=NULL, res=1e6, ...) {
  if (is.null(globalMeanSet)) globalMeanSet <- getGlobalMeans(obj)
  chrs <- paste0("chr", c(seq(1,22)))
  names(chrs) <- chrs
  getPairedChr <- function(chr) { 
    message("Computing shrunken eigenscores for ", column, " on ", chr, "...") 
    binmat <- getBinMatrix(cbind(assay(obj)[,column], assay(globalMeanSet)), rowRanges(obj), chr = chr, res = res)
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
  unlist(lapply(chrs, getPairedChr))
}