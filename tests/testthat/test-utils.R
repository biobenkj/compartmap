library(SummarizedExperiment)
library(RaggedExperiment)
library(GenomicRanges)

se <- SummarizedExperiment()
re <- RaggedExperiment()
df <- data.frame()
chrs <- c("chr2", "chr2", "chr1", "chr3")
gr <- GRanges(
  Rle(chrs, c(1, 3, 2, 4)),
  IRanges(1:10, width=10:1)
)

# extractOpenClosed {{{
gr.pc <- gr
pc.list <- list('pc' =  c(1,        -1,       1,         -1,       -2,       2,       -5,        5,      -10,     10))
pc.state.cutoff.zero <- c("open",   "closed", "open",   "closed", "closed", "open",   "closed", "open", "closed", "open")
pc.state.cutoff.two <-  c("closed", "closed", "closed", "closed", "closed", "closed", "closed", "open", "closed", "open")
pc.state.array <-  lapply(pc.state.cutoff.zero, function(i) ifelse(i == "open", "closed", "open")) |> unlist()
mcols(gr.pc) <- pc.list

test_that("extractOpenClosed", {
  expect_error(extractOpenClosed(se))
  expect_error(extractOpenClosed(re))
  expect_error(extractOpenClosed(df))
  expect_error(extractOpenClosed(gr))
  expect_no_error(extractOpenClosed(gr.pc))
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "rna", cutoff = 2),
    pc.state.cutoff.two
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "atac"),
    pc.state.cutoff.zero
  )
  expect_equal(
    extractOpenClosed(gr.pc, assay = "array"),
    pc.state.array
  )
})
# }}}

