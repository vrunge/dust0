library(testthat)
library(dust0)



##################################################
############## object dust = dust ################
##################################################

test_that("object dust = dust GAUSS",
          {
            obj_dust <- new(DUST_1D, "gauss", "randIndex_Eval3", 5)
            penalty <- 2*log(5000)
            data_all <- NULL
            for(i in 1:5)
            {
              data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(0,1), type = "gauss")
              data_all <- c(data_all, data)
              obj_dust$append_c(data, penalty)
              obj_dust$update_partition()
            }
            resObject <- obj_dust$get_partition()
            res <- dust.1D(data = data_all,
                           penalty = penalty,
                           model = "gauss",
                           method = "randIndex_Eval3", nbLoops = 5)
            expect_equal(all(res$changepoints == resObject$changepoints), TRUE)
            expect_equal(all(res$costQ == resObject$costQ), TRUE)
          })



test_that("object dust = dust VARIANCE many changes detected",
          {
            obj_dust <- dust.object.1D("variance", "randIndex_Eval4", 5)
            penalty <- 2*log(5)
            data_all <- NULL
            for(i in 1:5)
            {
              data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(1,1.5), type = "variance")
              data_all <- c(data_all, data)
              obj_dust$append_c(data, penalty)
              obj_dust$update_partition()
            }
            resObject <- obj_dust$get_partition()
            res <- dust.1D(data = data_all,
                           penalty = penalty,
                           model = "variance",
                           method = "randIndex_Eval4", nbLoops = 5)
            expect_equal(all(res$changepoints == resObject$changepoints), TRUE)
            expect_equal(all(res$costQ == resObject$costQ), TRUE)
          })




test_that("object dust = dust POISSON one by one",
          {
            obj_dust <- dust.object.1D("poisson", "randIndex_Eval4", 9)
            penalty <- 2
            data_all <- NULL
            for(i in 1:500)
            {
              data <- rpois(1, 100)
              data_all <- c(data_all, data)
              obj_dust$append_c(data, penalty)
              obj_dust$update_partition()
            }
            resObject <- obj_dust$get_partition()
            res <- dust.1D(data = data_all,
                           penalty = penalty,
                           model = "poisson")
            expect_equal(all(res$changepoints == resObject$changepoints), TRUE)
            expect_equal(all(res$costQ == resObject$costQ), TRUE)
          })

