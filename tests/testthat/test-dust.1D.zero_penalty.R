library(testthat)
library(dust0)

############################################
#### CONTINUOUS DISTRIBUTION ONLY
#### zero penalty => each index is a change
### default data size = 100
############################################

test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "gauss")
            data <- data_normalization_1D(data, type = "gauss")
            res <- dust.1D(data, penalty = 0, model = "gauss")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })


test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "exp")
            data <- data_normalization_1D(data, type = "exp")
            res <- dust.1D(data, penalty = 0, model = "exp")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })


test_that("zero penalty => each index is a change",
          {
            data <- dataGenerator_1D(type = "variance")
            data <- data_normalization_1D(data, type = "variance")
            res <- dust.1D(data, penalty = 0, model = "variance", method = "detIndex_Eval1")
            expect_equal(all(res$changepoints == 1:100), TRUE)
          })

