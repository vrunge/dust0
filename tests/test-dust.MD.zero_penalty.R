library(testthat)
library(dust)

############################################
#### CONTINUOUS DISTRIBUTION ONLY
#### zero penalty => each index is a change
### default data size = 100
############################################


test_that("zero penalty => each index is a change",
          {
            d <- 10
            data <- dataGenerator_MD(chpts = c(40,60,100),
                                     parameters =data.frame(matrix(abs(rnorm(3*d)), ncol = d)),
                                     type = "gauss")
            data <- data_normalization_MD(data, type = "gauss")
            res <- dust.MD(data, penalty = 0, model = "gauss")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })


test_that("zero penalty => each index is a change",
          {
            d <- 10
            data <- dataGenerator_MD(chpts = c(40,60,100),
                                     parameters =data.frame(matrix(abs(rnorm(3*d)), ncol = d)),
                                     type = "exp")
            data <- data_normalization_MD(data, type = "exp")
            res <- dust.MD(data, penalty = 0, model = "exp")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })


test_that("zero penalty => each index is a change",
          {
            d <- 100
            data <- dataGenerator_MD(chpts = c(40,60,100),
                                     parameters =data.frame(matrix(abs(rnorm(3*d)), ncol = d)),
                                     type = "variance")
            data <- data_normalization_MD(data, type = "variance")
            res <- dust.MD(data, penalty = 0, model = "variance")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })

