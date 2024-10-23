library(testthat)
library(dust)



############################################
############## outcome size ################
############################################

test_that("costQ size = nb size = data size",
          {
            data <- dataGenerator_1D(parameters = 5, type = "poisson")
            res <- dust.1D(data)
            expect_equal(length(res$costQ), length(data))
            expect_equal(length(res$nb), length(data))
          })


########################################################
############## constant dat => no change ###############
########################################################

test_that("constant dat => no change",
          {
            data <- dataGenerator_1D(sdNoise = 0)
            res <- dust.1D(data)
            expect_equal(res$changepoints == 100, TRUE)
          })

test_that("constant dat => no change",
          {
            data <- dataGenerator_1D(parameters = 0, type = "bern")
            res <- dust.1D(data)
            expect_equal(res$changepoints == 100, TRUE)
            expect_equal(all(res$costQ == rep(0,100)), TRUE)
          })








