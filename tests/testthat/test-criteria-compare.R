test_that("criteria_compare returns perfect agreement for identical labels", {
  skip_if_not_installed("aricode")

  labels <- c(1, 1, 2, 2, 3)

  expect_equal(scCCESS:::criteria_compare(labels, labels, "ARI"), 1)
  expect_equal(scCCESS:::criteria_compare(labels, labels, "NMI"), 1)
})
