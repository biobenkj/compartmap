#' Impute missing values/NAs with KNN
#'
#' @name imputeKNN
#'
#' @param obj Input SummarizedExperiment object
#' @param rowmax Maximum fraction of NAs that can exist in a row
#' @param colmax Maximum fraction of NAs that can exist in a column/sample
#' @param k Number of neighbors to be used in the imputation
#' @param maxp Largest block of regions/loci imputed using KNN
#' @param in.place Whether to modify the Beta/counts in place (default: TRUE)
#' @param drop.sparse.samps Whether to drop samples that are too sparse (default: TRUE)
#' @param assay The type of assay ("array", "bisulfite")
#'
#' @return Imputed data matrix that is added to the assays slot
#' @import impute
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' #impute
#' imputed <- imputeKNN(array.data.chr14, assay = "array")

imputeKNN <- function(obj, rowmax = 0.5, colmax = 0.8, k = 10,
                      maxp = 1500, in.place = TRUE, drop.sparse.samps = TRUE,
                      assay = c("array", "atac", "bisulfite")) {
  #match the assay args
  assay <- match.arg(assay)

  #double check the obj class is compatible
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")
  
  #check the names of the assays
  if (!any(getAssayNames(obj) %in% c("Beta", "counts"))) {
    stop("The assay slot should contain either 'Beta' for arrays or 'counts' for atac/bisulfite.")
  }
  
  #stop early if there aren't any NAs to impute
  if (!any(is.na(assay(obj)))) {
    message("No NAs found. Nothing to impute.")
    return(obj)
    }
  
  #filter out missing data based on chosen rowmax
  #otherwise imputation will blow up
  obj.clean <- cleanAssayRows(obj, rowmax = rowmax, assay = assay)
  
  #drop samples that have too sparse of data to use
  #this is the way to filter to samples with sufficient signal
  #before getting single cell imputation up
  if (drop.sparse.samps) {
    message("Dropping samples with >", colmax*100, "% NAs.")
    obj.clean <- cleanAssayCols(obj.clean, colmax = colmax, assay = assay)
    #stop if all the samples are now gone...
    if (ncol(obj.clean) == 0) {
      message("No samples left after sparisty filtering.")
      stop("Consider increasing the value of colmax closer to 1 and increasing maxp.")
    }
  } else {
    warning("Imputation may not work with samples that are too sparse!")
  }
  
  ## check if we are in beta land
  is.beta <- ifelse(min(assays(obj.clean)$Beta, na.rm = TRUE) < 0, FALSE, TRUE)
  message("Imputing missing data with kNN.")
  if (!is.beta) {
    if (assay == "array") {
      # assumes these are M-values
      imputed.data <- impute.knn(assays(obj.clean)$Beta, k = k,
                                 rowmax = rowmax, colmax = colmax,
                                 maxp = maxp)$data
    } else {
      # assumes the assay is bisulfite-seq
      # calculated as M-values
      imputed.data <- impute.knn(assays(obj.clean)$counts, k = k,
                                 rowmax = rowmax, colmax = colmax,
                                 maxp = maxp)$data
    }
  } else {
    if (assay == "array") {
      # assumes beta values and use squeezed M-values
      impute.knn(flogit(assays(obj.clean)$Beta), k = k,
                 rowmax = rowmax, colmax = colmax,
                 maxp = maxp)$data
    } else {
      # assumes that bisulfite-seq was given as betas
      impute.knn(flogit(assays(obj.clean)$counts), k = k,
                 rowmax = rowmax, colmax = colmax,
                 maxp = maxp)$data
    }
  }
  
  #add on another counts matrix to the assays slot if not in.place
  if (!in.place) {
    assays(obj.clean)$imputed.data <- imputed.data
  } else {
    switch(assay,
           #send M-values back to beta
           array = assays(obj.clean)$Beta <- fexpit(imputed.data),
           bisulfite = assays(obj.clean)$counts <- fexpit(imputed.data))
  }
  return(obj.clean)
}
