library(testthat)
library(dust)


### RETURNS ###

test_that("return is a vector",
          {expect_equal(is.vector(dataGenerator_1D()), TRUE)
          })
test_that("return is a matrix",
          {expect_equal(is.matrix(dataGenerator_MD()), TRUE)
          })
test_that("returns is a vector",
          {expect_equal(is.vector(dataGenerator_meanVar()), TRUE)
          })
test_that("returns is a dataframe",
          {expect_equal(is.data.frame(dataGenerator_Reg()), TRUE)
          })
test_that("returns is a vector",
          {expect_equal(is.vector(dataGenerator_AR1()), TRUE)
          })


### ERRORS CHECK ###
### 1D GAUSS ###

test_that("chpts non strictly increasing",
  {expect_error(dataGenerator_1D(chpts = c(10,20,20,30)))
  })

test_that("chpts non increasing",
  {expect_error(dataGenerator_1D(chpts = c(20,10,30)))
  })

test_that("non positive index in chpts",
  {expect_error(dataGenerator_1D(chpts = c(0,10)))
  })

test_that("negative index in chpts",
  {expect_error(dataGenerator_1D(chpts = -30))
  })

test_that("negative index in chpts",
  {expect_error(dataGenerator_1D(chpts = -c(1,2)))
          })

test_that("more parameters than changes",
  {expect_error(dataGenerator_1D(chpts = c(5,10), parameters = c(1,2,3)))
          })

test_that("less parameters than changes",
    {expect_error(dataGenerator_1D(chpts = c(5,10), parameters = 2))
          })

test_that("negative sdNoise",
    {expect_error(dataGenerator_1D(sdNoise = -1))
          })

test_that("gamma is zero",
          {expect_error(dataGenerator_1D(gamma = 0))
          })



