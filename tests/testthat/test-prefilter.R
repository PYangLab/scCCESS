test_that("prefilter removes zero-count genes when requested", {
  counts <- matrix(
    c(
      0, 0, 0,
      5, 2, 1,
      1, 0, 4
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("zero_gene", "gene_a", "gene_b")

  filtered <- prefilter(
    counts,
    minReads = 1,
    minGene = 1,
    minCountsperGene = 1,
    removeZeroGene = TRUE
  )

  expect_false("zero_gene" %in% rownames(filtered))
  expect_equal(nrow(filtered), 2)
})

test_that("prefilter can retain zero-count genes when requested", {
  counts <- matrix(
    c(
      0, 0,
      2, 3
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(counts) <- c("zero_gene", "gene_a")

  filtered <- prefilter(
    counts,
    minReads = 1,
    minGene = 1,
    minCountsperGene = 1,
    removeZeroGene = FALSE
  )

  expect_true("zero_gene" %in% rownames(filtered))
})
