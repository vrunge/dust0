library(testthat)
library(dust0)


########## OP_R_1D with constant limit values ##########
########## OP_R_1D with constant limit values ##########
########## OP_R_1D with constant limit values ##########

test_that("test gauss model data = 0 constant size 100",
          {
            res <- dataGenerator_1D(chpts = 100, parameters = 0, sdNoise = 0)
            op <- OP_R_1D(data = res, type = "gauss")
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


########## OP_R_1D versus OP_R_MD ##########
########## OP_R_1D versus OP_R_MD ##########
########## OP_R_1D versus OP_R_MD ##########

test_that("OP_R_1D and OP_R_MD same result if p copy of univariate data and penalty multiply by p",
          {
            p <- 10
            data <- dataGenerator_1D(100, sdNoise = 2) # noise != 1 to generate changes at random positions
            pen <- 2*log(100)
            op1D <- OP_R_1D(data = data, penalty = pen, type = "gauss")
            dataM <- matrix(data, p, 100, byrow = T)
            opMD <- OP_R_MD(data = dataM, penalty = p*pen, type = "gauss")

            expect_equal(op1D$changepoints, opMD$changepoints)
            expect_equal(op1D$costQ, opMD$costQ/p)
          })







