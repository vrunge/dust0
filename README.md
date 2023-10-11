<a id="top"></a>

# dust Vignette

### Vincent Runge
#### LaMME, Evry University
### October 11, 2023

___ 

> [Introduction](#intro)

> [Functions in R code](#Rcode)

> [DuST Algorithms](#dust)

> [Pruning Capacity](#pruning)


___ 

<a id="intro"></a>

## Introduction

The `dust` package contains methods for detecting multiple change-points within time-series based on the optimal partitioning algorithm. A few models from the exponential family are considered (Gauss, Poisson, Exponential...).

The proposed algorithm is a pruned dynamic programming algorithm optimizing a penalized likelihood **using a new pruning rule**, different from PELT or FPOP. We called this method the DuST pruning rule, standing for **Du**ality **S**imple **T**est.

Indeed, indices for potential last change-point are discarded by considering some constrained optimization problems. For each potential last change-point index, evaluating its associated dual function at a random testing point enables a fast and efficient test.

[Back to Top](#top)

___ 

<a id="Rcode"></a>

## Functions in R code

### Data Generators

**dataGenerator_1D** is used to generate data with a given vector of change-point (e.g. `chpts = c(50,100)`), parameter vector (e.g. `parameters = c(0,1)`), a shared variance for all data (usually `sdNoise = 1`) and a type of probability distribution in `type`. We have the following choices for type:

- `type = "gauss"`

- `type = "poisson"`

- `type = "exp"`

- `type = "bern"`

- ...

**dataGenerator_MultiD** concatenates `p` copies of `dataGenerator_1D` function.

**dataGenerator_MV** is used for change in mean and variance for the Gaussian problem


**dataGenerator_Reg** generates 2-dimensional data (x,y) following a simple linear regression link (y = Ax + B + noise) with A and B changing over time (after each change-point)



### Dual Functions

`Dual_1D` returns the value of the dual at a point `mu` when comparing index `s1` with the constraint from index `s2`. with option `OP = TRUE` the optimal partitioning algorithm is used to have the true constants in the cost functions.

`dual_1D <- function(mu, x, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(x))`

- `x` is raw data

- If `OP` is `true`, we run the OP algorithm to have the optimal cost vector to use in cost functions. See the function `OP_R`.

- at time `t`, we evaluate the dual function at point `mu` when trying to remove index `s1` using function linked to index `s2` (we have a unique constraint, which means that the dual is a one-parametric function)

- Depending on the `type`, different functions `A`, `B` and `mu_max` are used (see the code in file `functions_by_type.R`)


### dust_R

The function `dust_R` has the following parameters:

`dust_R <- function(data, penalty, type = "gauss")`

and returns a list of two elements

- the change-points found by the algo

- the number of indices to consider in the minimization at each time step

the function uses values obtained by `dual_1D` (with `OP = FALSE`) to discard indices with DuST rule. 


**... ADD OPTIONS to change the sampling?**

**TO BE DONE**


[Back to Top](#top)

___ 

<a id="dust"></a>

## DuST Algorithms



[Back to Top](#top)


<a id="pruning"></a>

___ 

## Pruning Capacity


[Back to Top](#top)

