library(testthat)
library(dust0)

############################################
#### constant data = no change
############################################


test_that("MD = 1D gauss",
          {
          d <- 2
          data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(1,0.01,1), type = "gauss")
          dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

          r1 <- dust.1D(data, penalty = 0.1*log(100), model = "gauss")
          rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "gauss", method = "detIndex_Eval0")
          expect_equal(rM$costQ, d*r1$costQ)
          expect_equal(r1$changepoints, rM$changepoints)
          })


test_that("MD = 1D poisson",
          {
            d <- 30
            type <- "poisson"
            data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(1,0.01,1), type = type)
            data <- data_normalization_1D(data, type = type)
            dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

            r1 <- dust.1D(data, penalty = 0.1*log(100), model = type)
            rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = type, method = "detIndex_Eval4")
            expect_equal(rM$costQ, d*r1$costQ)
            expect_equal(r1$changepoints, rM$changepoints)
          })


test_that("MD = 1D exp",
          {
            d <- 5
            type <- "exp"
            data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(1,0.01,1), type = type)
            data <- data_normalization_1D(data, type = type)
            dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

            r1 <- dust.1D(data, penalty = 0.1*log(100), model = type)
            rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = type, method = "detIndex_Eval4")
            expect_equal(rM$costQ, d*r1$costQ)
            expect_equal(r1$changepoints, rM$changepoints)
          })


test_that("MD = 1D bern",
          {
            d <- 5
            type <- "bern"
            data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(0.5,0.4,0.8), type = type)
            dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

            r1 <- dust.1D(data, penalty = log(100), model = type)
            rM <- dust.MD(dataMD, penalty =  d*log(100), model = type, method = "detIndex_Eval4")
            expect_equal(rM$costQ, d*r1$costQ)
            expect_equal(r1$changepoints, rM$changepoints)
          })


test_that("MD = 1D binom",
          {
            d <- 5
            type <- "binom"
            data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(0.5,0.4,0.8), type = type)
            data <- data_normalization_1D(data, type = type)
            dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

            r1 <- dust.1D(data, penalty = log(100), model = type)
            rM <- dust.MD(dataMD, penalty =  d*log(100), model = type, method = "detIndex_Eval1")
            expect_equal(rM$costQ, d*r1$costQ)
            expect_equal(r1$changepoints, rM$changepoints)
          })


test_that("MD = 1D variance",
          {
            d <- 10
            type <- "variance"
            data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(5,1,4), type = type)
            data <- data_normalization_1D(data, type = type)
            dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

            r1 <- dust.1D(data, penalty = log(100), model = type)
            rM <- dust.MD(dataMD, penalty =  d*log(100), model = type, method = "detIndex_Eval1")
            expect_equal(rM$costQ, d*r1$costQ)
            expect_equal(r1$changepoints, rM$changepoints)
          })
