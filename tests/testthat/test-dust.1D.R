library(testthat)
library(dust)


############################################
#### zero penalty => each index is a change
### default data size = 100
############################################

test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "gauss")
            res <- dust.1D(data, penalty = 0, model = "gauss")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })


test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "poisson")
            res <- dust.1D(data, penalty = 0, model = "poisson")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })

test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "exp")
            res <- dust.1D(data, penalty = 0, model = "exp")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })

#test_that("zero penalty => each index is a change",
#          {
#            data <- dataGenerator_1D(type = "geom")
#            res <- dust.1D(data, penalty = 0, model = "geom")
#            expect_equal(all(res$changepoints == 1:100), TRUE)
#          })

test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "bern")
            res <- dust.1D(data, penalty = 0, model = "bern")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })

#test_that("zero penalty => each index is a change",
#          {
#            data <- dataGenerator_1D(type = "binom")
#            res <- dust.1D(data, penalty = 0, model = "binom")
#            expect_equal(all(res$changepoints == 1:100), TRUE)
#          })

#test_that("zero penalty => each index is a change",
#          {
#            data <- dataGenerator_1D(type = "negbin")
#            res <- dust.1D(data, penalty = 0, model = "negbin")
#            expect_equal(all(res$changepoints == 1:100), TRUE)
#          })

#test_that("zero penalty => each index is a change",
#          {
#            data <- dataGenerator_1D(type = "variance")
#            res <- dust.1D(data, penalty = 0, model = "variance")
#            expect_equal(all(res$changepoints == 1:100), TRUE)
#          })


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








