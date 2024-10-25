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

test_that("costQ size = nb size = data size",
          {
            data <- dataGenerator_1D(parameters = 5, type = "variance")
            res <- dust.1D(data)
            expect_equal(length(res$costQ), length(data))
            expect_equal(length(res$nb), length(data))
          })




