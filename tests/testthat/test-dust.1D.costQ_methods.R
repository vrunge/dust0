
library(testthat)
library(dust0)


##############################################################
############## all methods  - all costs as OP ################
##############################################################

test_that("costQ all the same, for all types and methods",
          {
            types <- c( "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")
            for(i in types)
            {
            data <- dataGenerator_1D(parameters = 0.5, type = i)
            data <- data_normalization_1D(data, type = i)
            res0 <- dust.1D(data, model = i, method = "randIndex_Eval0")
            res1 <- dust.1D(data, model = i, method = "randIndex_Eval1")
            res2 <- dust.1D(data, model = i, method = "randIndex_Eval2")
            res3 <- dust.1D(data, model = i, method = "randIndex_Eval3")
            res4 <- dust.1D(data, model = i, method = "randIndex_Eval4")
            res5 <- dust.1D(data, model = i, method = "randIndex_Eval5")
            expect_equal(res0$costQ, res1$costQ)
            expect_equal(res0$costQ, res2$costQ)
            expect_equal(res0$costQ, res3$costQ)
            expect_equal(res0$costQ, res4$costQ)
            expect_equal(res0$costQ, res5$costQ)
            res0 <- dust.1D(data, model = i, method = "detIndex_Eval0")
            res1 <- dust.1D(data, model = i, method = "detIndex_Eval1")
            res2 <- dust.1D(data, model = i, method = "detIndex_Eval2")
            res3 <- dust.1D(data, model = i, method = "detIndex_Eval3")
            res4 <- dust.1D(data, model = i, method = "detIndex_Eval4")
            res5<- dust.1D(data, model = i, method = "detIndex_Eval5")
            expect_equal(res0$costQ, res1$costQ)
            expect_equal(res0$costQ, res2$costQ)
            expect_equal(res0$costQ, res3$costQ)
            expect_equal(res0$costQ, res4$costQ)
            expect_equal(res0$costQ, res5$costQ)
            }
          })

