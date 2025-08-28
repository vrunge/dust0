library(testthat)
library(dust0)


test_that("Infinity A",
          {
            expect_equal(A(type = "exp")(0), Inf)
            expect_equal(A(type = "exp")(10^(-15)), Inf)
            expect_equal(A(type = "geom")(0), Inf)
            expect_equal(A(type = "geom")(10^(-15)), Inf)
            expect_equal(A(type = "negbin")(0), Inf)
            expect_equal(A(type = "negbin")(10^(-15)), Inf)
            })

test_that("minus Infinity B",
          {
            expect_equal(B(type = "exp")(0), -Inf)
            expect_equal(B(type = "exp")(-10^(-15)), -Inf)
            expect_equal(B(type = "poisson")(0), -Inf)
            expect_equal(B(type = "poisson")(-10^(-15)), -Inf)
            expect_equal(B(type = "geom")(1), -Inf)
            expect_equal(B(type = "geom")(1-10^(-15)), -Inf)
            expect_equal(B(type = "bern")(1), Inf)          # + INF !!!!
            expect_equal(B(type = "bern")(1+10^(-15)), Inf) # + INF !!!!
            expect_equal(B(type = "bern")(0), -Inf)
            expect_equal(B(type = "bern")(-10^(-15)), -Inf)
            expect_equal(B(type = "binom")(1), Inf)           # + INF !!!!
            expect_equal(B(type = "binom")(1+10^(-15)), Inf)  # + INF !!!!
            expect_equal(B(type = "binom")(0), -Inf)
            expect_equal(B(type = "binom")(-10^(-15)), -Inf)
            expect_equal(B(type = "negbin")(0), -Inf)
            expect_equal(B(type = "negbin")(-10^(-15)), -Inf)
          })


########## mu_max ##########
########## mu_max ##########  cumsum(S) instead cumsum(statistics(S)) => all statistics here = identity
########## mu_max ##########

test_that("case Mt = Ms. Ratio R(mu) = const, mu_max = 1",
          {
            S <- rep(1,3)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)  # Ratio R(mu) = 1
            expect_equal(mu_max(S,s1,s2,t,type = "exp"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "poisson"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "geom"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "negbin"), 1)

            S <- rep(0,3)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)  # Ratio R(mu) = 0
            expect_equal(mu_max(S,s1,s2,t,type = "exp"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "poisson"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "geom"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "negbin"), 1)
          })

test_that("case Mt = min possible or max possible, mu_max = 0",
          {
            #  Mt = 0
            S <- c(1,2,0)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)
            expect_equal(mu_max(S,s1,s2,t,type = "exp"), 0)
            expect_equal(mu_max(S,s1,s2,t,type = "poisson"), 0)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), 0)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), 0)
            expect_equal(mu_max(S,s1,s2,t,type = "negbin"), 0)

            # Mt = 1, for geom, values in S >=1, mu_max = 0
            S <- c(3,2,1)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)
            expect_equal(mu_max(S,s1,s2,t,type = "geom"), 0)

            # for bern and binom case Mt = 1, mu_max = 0
            S <- c(1,0,1) # Ms !=1 here (if = 1, mu_max = 1)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), 0)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), 0)

          })


test_that("case Ms = 0, mu_max = 1",
          {
            S <- c(1,0,1)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)
            expect_equal(mu_max(S,s1,s2,t,type = "exp"), 1)
            expect_equal(mu_max(S,s1,s2,t,type = "poisson"), 1)
            ### case not possible for geom, where all values >= 1
            expect_equal(mu_max(S,s1,s2,t,type = "negbin"), 1)
          })

test_that("case Ms = 0 or 1, mu_max = formula for bern + binom",
          {
            S <- c(1,1,1,0,0)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(5)
            Mt <- 1/3; Ms <- 1
            res <-  min(Mt / Ms, (1 - Mt) / (1 - Ms), na.rm = TRUE)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), res)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), res)

            S <- c(1,0,1,0,0)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(5)
            Mt <- 1/3; Ms <- 0
            res <-  min(Mt / Ms, (1 - Mt) / (1 - Ms), na.rm = TRUE)
            expect_equal(mu_max(S,s1,s2,t,type = "bern"), res)
            expect_equal(mu_max(S,s1,s2,t,type = "binom"), res)
          })


test_that("case Ms = 1, mu_max = 1, geom",
          {
            S <- c(1,1,3)
            S <- c(0, cumsum(S))
            s1 <- shift(2); s2 <- shift(1); t <- shift(3)
            expect_equal(mu_max(S,s1,s2,t,type = "geom"), 1)
          })




########## min_cost ##########
########## min_cost ##########
########## min_cost ##########


test_that("min_cost testing the limit cases",
          {
            s <- 0; t <- 1

            type <- "exp"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 0)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), -Inf)

            type <- "poisson"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 0)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)

            type <- "geom"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 1)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)

            type <- "bern"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 0)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)
            S <- c(0, 1)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)

            type <- "binom"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 0)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)
            S <- c(0, 1)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)

            type <- "negbin"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, 0)
            expect_equal(min_cost(A_, B_, S, shift(s), shift(t), 0), 0)
          })



########## evalDual ##########
########## evalDual ########## don't use it at mu_max value
########## evalDual ##########


test_that("evalDual at mu max point",
          {
            s1 <- 3; s2 <- 2; t <- 4
            type <- "gauss"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, cumsum(rnorm(3)))
            expect_equal(evalDual(mu_max(S,s1,s2,t,type = type), A_, B_, S, s1, s2, t, 0, 0), NaN)

            type <- "poisson"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, cumsum(c(3,2,5)))
            mu_max(S,s1,s2,t,type = type) # 1
            expect_equal(evalDual(mu_max(S,s1,s2,t,type = type), A_, B_, S, s1, s2, t, 0, 0), NaN)

            type <- "geom"
            A_ <- A(type = type); B_ <- B(type = type)
            S <- c(0, cumsum(c(1,1,2)))
            mu_max(S,s1,s2,t,type = type) # < 1
            expect_equal(evalDual(mu_max(S,s1,s2,t,type = type), A_, B_, S, s1, s2, t, 0, 0), NaN)
          })
