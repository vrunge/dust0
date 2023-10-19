library(testthat)
library(dust)

test_that("test gauss model data = 0 constant size 100",
          {
            res <- dataGenerator_1D(chpts = 100, sdNoise = 0)
            op <- OP_R_1D(data = res)
            expect_equal(op$changepoints, 100)
            expect_equal(op$costQ, rep(0,100))
          })


test_that("test poisson model data = 0 constant size 100",
          {
            res <- rep(1,100)
            op <- OP_R_1D(data = res, type = "poisson")
            expect_equal(op$changepoints, 100)
            expect_equal(op$costQ, 1:100)
          })


test_that("test exp model data = 0 constant size 100",
          {
            res <- rep(0,100)
            op <- OP_R_1D(data = res, type = "exp")
            expect_equal(op$changepoints, 100)
            expect_equal(op$costQ, rep(-Inf, 100))
          })

