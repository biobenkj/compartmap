### verify errors {{{
err.verifySE <- "Input needs to be a SummarizedExperiment"
err.verifyCoords <- paste(
  "The SummarizedExperiment you have provided has no coordinates.\n",
  "Compartment extraction will fail.\n",
  "Please provide rowRanges with genomic coordinates for the object."
)
err.verifyAssayNames.beta <- "The 'assays' slot should contain 'Beta' for array data"
err.verifyAssayNames.atac_counts <- "The 'assays' slot should contain 'counts' for atac data"
err.verifyAssayNames.rna_counts <- "The 'assays' slot should contain 'counts' for rna data"
# }}}
