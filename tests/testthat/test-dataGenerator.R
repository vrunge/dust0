library(testthat)
library(dust)

test_that("returns is a vector",
  {expect_equal(is.vector(dataGenerator_1D()), TRUE)
  })

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

