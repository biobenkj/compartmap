chrs <- c("chr2", "chr2", "chr1", "chr3")
gr <- GRanges(
  Rle(chrs, c(1, 3, 2, 4)),
  IRanges(1:10, width=10:1)
)

gr.flip <- gr
mcols(gr) <- data.frame(conf.est = rep(0.7, 10), score = rep(0.02, 10))
mcols(gr.flip) <- data.frame(flip.conf.est = rep(0.7, 10), flip.score = rep(0.02, 10))

# filterer {{{
test_that("filterer", {
  expect_equal(
    gr[compartmap:::filterer(gr, min.conf = 0.7, min.eigen = 0.02),],
    gr
  )
  expect_equal(
    gr[compartmap:::filterer(gr, min.conf = 1, min.eigen = 1),],
    gr[0]
  )
  expect_equal(
    gr.flip[compartmap:::filterer(gr, min.conf = 0.7, min.eigen = 0.02),],
    gr.flip
  )
  expect_equal(
    gr.flip[compartmap:::filterer(gr, min.conf = 1, min.eigen = 1),],
    gr.flip[0]
  )
})
# }}}

# filterCompartments messages {{{
test_that("filterCompartments messages", {
  expect_message(
    filterCompartments(gr),
    "Filtering compartments based on a minimum confidence of 70%"
  )
  expect_message(
    filterCompartments(gr),
    "Filtering compartments based on a minimum absolute eigen value of 0.02"
  )

  expect_message(
    filterCompartments(gr, min.conf = 1),
    "Filtering compartments based on a minimum confidence of 100%"
  )
  expect_message(
    filterCompartments(gr, min.eigen = 5),
    "Filtering compartments based on a minimum absolute eigen value of 5"
  )
})
# }}}

# filterCompartments {{{
test_that("filterCompartments", {
  grlist <- list(gr, gr)
  expect_equal(
    filterCompartments(grlist),
    grlist
  )
  expect_equal(
    filterCompartments(grlist, min.conf = 1, 1),
    list(gr[0], gr[0])
  )
})
# }}}
