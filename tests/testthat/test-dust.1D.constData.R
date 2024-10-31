library(testthat)
library(dust)

############################################
#### constant data = no change
############################################

test_that("constant data = no change",
          {
            data <- rep(0, 100)
            res <- dust.1D(data, penalty = 1, model = "gauss")
            expect_equal(res$changepoints, 100)
          })


test_that("constant data = no change",
          {
            data <- rep(0, 100)
            res <- dust.1D(data, penalty = 1, model = "poisson")
            expect_equal(res$changepoints, 100)
          })

test_that("constant data = no change",
          {
            data <- rep(1, 100)
            res <- dust.1D(data, penalty = 1, model = "exp")
            expect_equal(res$changepoints, 100)
          })



test_that("constant data = no change",
          {
            data <- rep(1,100)
            res <- dust.1D(data, penalty = 1, model = "variance")
            expect_equal(res$changepoints, 100)
          })



test_that("constant data = no change",
          {
            data <- rep(1,100)
            res <- dust.1D(data, penalty = 1, model = "geom")
            expect_equal(res$changepoints, 100)
          })


test_that("constant data = no change",
          {
            data <- rep(0,100)
            res <- dust.1D(data, penalty = 1, model = "bern")
            expect_equal(res$changepoints, 100)
          })


test_that("constant data = no change",
          {
            data <- rep(0,100)
            res <- dust.1D(data, penalty = 1, model = "binom")
            expect_equal(res$changepoints, 100)
          })


test_that("constant data = no change",
          {
            data <- rep(0,100)
            res <- dust.1D(data, penalty = 1, model = "negbin")
            expect_equal(res$changepoints, 100)
          })

