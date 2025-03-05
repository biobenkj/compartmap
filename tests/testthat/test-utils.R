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

# checkAssayType{{{
test_that("checkAssayType", {
  expect_true(checkAssayType(se))
  expect_false(checkAssayType(re))
  expect_false(checkAssayType(gr))
  expect_false(checkAssayType(df))
})
# }}}

# getAssayNames{{{

nrows <- 200
ncols <- 6
counts1 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts2 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
counts3 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(
  Treatment = rep(c("a", "b"), 3),
  row.names = LETTERS[1:6]
)
se <- SummarizedExperiment(
  assays = SimpleList(one = counts1, two = counts2, three = counts3),
  colData = colData
)

test_that("getAssayNames", {
  expect_equal(
    getAssayNames(se),
    c("one", "two", "three")
  )
})
# }}}

# getChrs {{{
test_that("getChrs", {
  expect_equal(
    getChrs(gr),
    unique(chrs)
  )
})
# }}}

# getGenome {{{
test_that("getGenome", {
  bundled_genomes <- c("hg19", "hg38", "mm9", "mm10")
  lapply(bundled_genomes, function(g) {
    expect_equal(
      getGenome(g),
      get(paste0(g, ".gr"))
    )
    expect_equal(
      getGenome(g, type = "openseas"),
      get(paste0("openSeas.", g))
    )
  })
  expect_error(getGenome("NA50"))
})
# }}}

# getSeqLengths{{{
forseqlengths.gr <- GRanges(
  Rle(c("chr1", "chr2", "chr3"), c(1, 3, 6)),
  IRanges(1:10, width=10:1, names=head(letters, 10)),
  Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score=1:10,
  GC=seq(1, 0, length=10)
)
seqlengths.list <- c(chr1 = 100, chr2 = 200, chr3 = 300)
seqinfo(forseqlengths.gr) <- GenomeInfoDb::Seqinfo(paste0("chr", 1:3), seqlengths.list, NA, "mock1")
# to get seqlengths(gr) set to
# chr1 chr2 chr3
#  100  200  300

seqlengths.list <- c(chr1 = 100, chr2 = 200, chr3 = 300)
seqinfo(forseqlengths.gr) <- GenomeInfoDb::Seqinfo(paste0("chr", 1:3), seqlengths.list, NA, "mock1")
test_that("getSeqLengths", {
  lapply(seq(1:3), function(i) {
    expect_equal(
      getSeqLengths(forseqlengths.gr, paste0("chr", i)),
      seqlengths.list[i]
    )
  })
  # missing chromosome
  expect_error(getSeqLengths(forseqlengths.gr, "chr0"))
})
# }}}

# cleanAssay {{{

checkCleanAssay <- function(se.array, se.bisulfite, cleanFun, dimFun) {
  expect_equal(
    dimFun(assays(cleanFun(se.array, assay = "array"))$Beta),
    0
  )
  expect_equal(
    dimFun(assays(cleanFun(se.bisulfite, assay = "bisulfite"))$counts),
    0
  )
  expect_equal(
    dimFun(assays(cleanFun(se.array, assay = "array", na.max = 0.9))$Beta),
    10
  )
  expect_equal(
    dimFun(assays(cleanFun(se.bisulfite, assay = "bisulfite", na.max = 0.9))$counts),
    10
  )
}

test_that("cleanAssayRows", {
  cleanFun <- compartmap:::cleanAssayRows
  dimFun <- nrow
  m <- matrix(1:10, 10, 10)
  m.allNA <- apply(m, 1, function(i) ifelse(i <= 2, i, NA))
  se.array <- SummarizedExperiment(assays = SimpleList(Beta = m.allNA))
  se.bisulfite <- SummarizedExperiment(assays = SimpleList(counts = m.allNA))
  checkCleanAssay(se.array, se.bisulfite, cleanFun, dimFun)
})

test_that("cleanAssayCols", {
  cleanFun <- compartmap:::cleanAssayCols
  dimFun <- ncol
  m <- matrix(1:10, 10, 10)
  m.allNA <- apply(m, 2, function(i) ifelse(i <= 2, i, NA))
  se.array <- SummarizedExperiment(assays = SimpleList(Beta = m.allNA))
  se.bisulfite <- SummarizedExperiment(assays = SimpleList(counts = m.allNA))
  checkCleanAssay(se.array, se.bisulfite, cleanFun, dimFun)
})
### }}}

