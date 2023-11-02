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

The `dust` package contains methods for detecting multiple change-points within time-series based on the optimal partitioning algorithm with pruning. A few models from the exponential family are considered (Gauss, Poisson, Exponential...).

The proposed algorithm is a pruned dynamic programming algorithm optimizing a penalized likelihood **using a new pruning rule**, different from PELT or FPOP. We called this method the **DUST** pruning rule, standing for **Du**ality **S**imple **T**est.

Indeed, indices for potential last change-point are discarded by considering some constrained optimization problems. For each potential last change-point index, evaluating its associated dual function at a random testing point enables a fast and efficient test.

[Back to Top](#top)

___ 

<a id="Rcode"></a>

## Functions in R code

### Data Generators

**dataGenerator_1D** is used to generate data with a given vector of change-point (e.g. `chpts = c(50,100)` for one change at position `50` and data length `100`), parameter vector (e.g. `parameters = c(0,1)`) and a type of probability distribution in `type`. We have the following choices for type:
  
- `type = "gauss"` (additional parameters `sdNoise` and `gamma`)

- `type = "exp"` 

- `type = "poisson"`

- `type = "geom"` 

- `type = "bern"` 

- `type = "binom"` (additional parameter `nbTrials`)

- `type = "negbin"` (additional parameter `nbSuccess`)



**dataGenerator_MultiD** concatenates `p` copies of `dataGenerator_1D` function.

**dataGenerator_MV** is used for change in mean and variance for the Gaussian problem

**dataGenerator_Reg** generates 2-dimensional data frame `(x,y)` following a simple linear regression link (`y = Ax + B + noise`) with `A` and `B` changing over time (after each change-point)


### OP in R

The base function `OP_R` is used to compute the change-point vector with the simplest dynamic programming algorithm with no pruning. This method is of quadratic time complexity. We propose 3 such functions:`OP_R_1D`, `OP_R_MultiD`, `OP_R_2param`.

`OP_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss")`

Example: `OP_R_1D(dataGenerator_1D(chpts = c(200,400), c(0,1), type = "gauss"))`


`OP_R_2param` is used for:

- `type = "meanVar"` change in Gaussian data in mean and variance

- `type = regression` change in simple regression model



### Dual Functions
  
`dual_1D` returns the value of the dual at a point `mu` when comparing index `s1` with the constraint from index `s2` at time `t`. With option `OP = TRUE` the optimal partitioning algorithm is used to have the true constants in the cost functions with penalty `penalty` and a pruning option `pruningOpt`.

`dual_1D <- function(mu, data, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(data), pruningOpt = 3)`

- `data` is raw data

- If `OP` is `true`, we run the OP algorithm to have the optimal cost vector to use in cost functions. See the function `OP_R`.

- at time `t`, we evaluate the dual function at point `mu` when trying to remove index `s1` using function linked to index `s2` (we have a unique constraint, which means that the dual is a one-parametric function)

- Depending on the `type`, different functions `A`, `B`, `statistic`, `mu_max` and `evalDual` are used (see the code in file `functions_by_type.R`)

Function `dual_1D` allows us to study the shape of the dual. 

### dust_R 

We propose a few R functions computing the change-point location with dust method: `dust_R_1D`, `dust_R_MultiD`, `dust_R_2Dquad`.

The function `dust_R_1D` has the following parameters:

`dust_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss", pruningOpt = 2)`


We have different type of possible pruning:

- `pruningOpt == 0`: nothing

- `pruningOpt == 1`: PELT

- `pruningOpt == 2`: dust

- `pruningOpt == 3`: dust + PELT
  

and returns a list of two elements:

- `changepoints`: the change-points found by the algo

- `nb`: the number of indices to consider in the minimization at each time step

- `lastIndexSet`: the vecotr of indices saved by the algo in the dynamic programming algorithm at the last iteration

- `costQ`: the vector of optimal cost (of size `length(data)`)


### Plot functions 


`plot_dual_1D` is a function displaying the result of `dual_1D` for a vector of mu values.

`plot_dual_1D <- function(mu =  (1:99)/100,`
                         `data, s1, s2,`
                         `type = "gauss",`
                         `OP = FALSE,`
                         `penalty = 2*length(data))`
                         

We use the function `plot_dual_1D` with `OP = TRUE` to plot the true dual function seen by the dynamic programming algorithm. 

What we called the "pruning interval" is the interval of values between the vertical green lines for which the dual function takes a value higher than the pruning threshold (horizontal line in red), so that the index considered `s1` is pruned by `s2` at time `n`. 


Using function `barplot_dual_1D` we can repeat the generation of the pruning interval `nb` and count the number of time each value mu is in this interval.

We add the values in the bar plot only if at the final time step `n`, the index `s1` has not been removed by the algorithm (the pruning option is given by option `pruningOpt`).

`barplot_dual_1D <- function(nb = 1000, s1 = 18, s2 = 15,`
                            `n = 20,`
                            `oneParam = 0,`
                            `type = "gauss",`
                            `penalty = 2*log(n),`
                            `pruningOpt = 3)`
                            
                        
[Back to Top](#top)

___ 

<a id="dust"></a>

## DuST Algorithms



[Back to Top](#top)


<a id="pruning"></a>

___ 

## Pruning Capacity


[Back to Top](#top)

