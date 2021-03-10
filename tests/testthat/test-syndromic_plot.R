
context("syndromic_plot")

library(Gifi)
# library(syndRomics)

pca<-prcomp(mtcars, center = TRUE, scale. = TRUE)
nlpca<-princals(mtcars, ndim=ncol(mtcars))

load_df<-extract_loadings(pca, mtcars)

testthat::test_that("warning", {
  testthat::expect_warning(
    syndromic_plot(pca, mtcars, cutoff = 0.1, ndim=13),
    "The specified ndim is bigger than number of dimensions, using ndim = 11"
  )
})

testthat::test_that("errors", {
  testthat::expect_error(
    extract_syndromic_plot(pca, mtcars, cutoff = 0.1,pc = "PC11"),
    "There is no loading above cutoff for PC11"
  )
})
