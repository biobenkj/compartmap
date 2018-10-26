# turn a matrix with appropriate rownames back into a GRanges
#
# note: the granges() method is already documented and this brings no change
#
setAs("matrix", "GRanges",
      function(from) { 
        gr <- as(rownames(from), "GRanges")
        mcols(gr) <- as(from, "DataFrame")
        return(gr) 
      })

# wrap
setMethod("granges", "matrix", function(x) as(x, "GRanges"))

# turn it back 
setAs("GRanges", "matrix", 
      function(from) { 
        mat <- as(mcols(from), "matrix")
        rownames(mat) <- as.character(from) 
        return(mat)  
      }) 

# wrap
setMethod("as.matrix", "GRanges", function(x, ...) as(x, "matrix"))

# a test:
if (FALSE) { 

  testMat <- matrix(rnorm(4), nrow=2, ncol=2, 
                    dimnames=list(c("chr1:1-1000","chr2:2-2000"), c("A","B")))
  show(testMat)
  #                        A            B
  # chr1:1-1000 -0.974000464 -0.705903477
  # chr2:2-2000 -0.251981621  1.421097044
  #
  identical(as(as(x, "GRanges"), "matrix"), x)
  # [1] TRUE 

  testGr <- as(testMat, "GRanges") 
  show(testGr) 
  # 
  # GRanges object with 2 ranges and 2 metadata columns:
  #     seqnames    ranges strand |                  A                  B
  #        <Rle> <IRanges>  <Rle> |          <numeric>          <numeric>
  # [1]     chr1    1-1000      * | -0.974000463966715 -0.705903476640502
  # [2]     chr2    2-2000      * | -0.251981620673546   1.42109704420397
  # -------
  # seqinfo: 2 sequences from an unspecified genome; no seqlengths
  # 
  identical(as.matrix(testGr), testMat)
  # [1] TRUE 
}

