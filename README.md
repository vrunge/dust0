<a id="top"></a>

# dust Vignette

### Vincent Runge and Simon Querné

#### LaMME, Evry University, August 31, 2024

<center>
<img src="man/figures/dust.png" alt="" style="width:30%;"/>
</center>

> [Quick start](#start)

> [Models And Data Generators](#Models)

> [dust 1D and MD Algorithms](#dust1D)

> [dust 2D Algorithms](#dust2D)

> [Rcpp Structure and new cost integration](#rcpp)

> [Pruning Capacity](#pruning)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="start"></a>

## Quick start

### Introduction

The `dust` package contains methods **for detecting multiple
change-points within time-series** based on the optimal partitioning
algorithm, which is a dynamic programming (DP) algorithm. Our algorithms
optimize a penalized likelihood and the DP algorithm is encoded with
pruning rules for reducing execution time. The novelty of the `dust`
package consists in its pruning step. We use a **new pruning rule**,
different from the two standard ones: [PELT rule
(2012)](\doi%7B10.1080/01621459.2012.737745%7D) and [FPOP rule
(2017)](\doi%7B10.1007/s11222-016-9636-3%7D).

We called this method the **DUST** pruning rule, standing for
**Du**ality **S**imple **T**est. This method is based on considering
some optimization problems under inequality constraints and its dual for
discarding indices in the search for the last change-point index.

Data are modeled by a cost derived from the **exponential family**
(Gauss, Poisson, Exponential…). We provide a polymorphic structure in
Rcpp **allowing any User to easily add a new cost** of their own choice.
The steps to follow are described in details [in this Section](#rcpp).
The User only have to provide few Rcpp functions: the minimal cost
function, the dual function and its derivative. In addition we need some
information on the domain of the dual.

Various tests and simulations are provided in this **readme file** and
in the **simulations folder** and show that the dust dual approach is
**highly efficient in all regimes** (many or few changes) with improved
time efficient comparing to PELT and FPOP. Furthermore, unlike these 2
methods, the DUST method is capable of reducing time for multivariate
cost functions (See Section [pruning](#pruning)).

### Installing the dust Rcpp package

**REQUIREMENTS:** - R &gt;= 3.4 - devtools :
`install.packages('devtools')`

The package can then be installed from the github repo with the
following command:

    devtools::install_github("vrunge/dust")

and imported with:

    library(dust)

### A simple example

We generate some 1D time series of length `400` from the Gaussian model
and one change in the middle of the sequence.

`data <- dataGenerator_1D(chpts = c(200,400), c(0,1), type = "gauss")`

We segment data using the dust 1D method coded in Rcpp. We give data,
the penalty value and the type of cost to be used. It can be one of the
following: `"gauss"` (additional parameters `sdNoise` and `gamma`),
`"exp"`, `"poisson"`, `"geom"`, `"bern"`, `"binom"` (additional
parameter `nbTrials`), `"negbin"` (additional parameter `nbSuccess`) and
`variance`. See next [Section](#Models).

`dust.1D(data, penalty = 2*log(length(data)), model = "gauss")`

The result is a list whose elements are:

-   `changepoints`A vector of change points (the index ending each of
    the segments)

-   `lastIndexSet`The list of non-pruned indices at the end of the
    analysis

-   `nb` The number of indices to consider at each time step (its length
    is equal to data length)

-   `costQ` The minimal cost of the optimization problem at each time
    step.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="Models"></a>

## Models And Data Generators

### Data Generators in 1D and MultiD

**dataGenerator\_1D** is used to generate data with a given vector of
change-point (e.g. `chpts = c(50,100)` for one change at position `50`
and data length `100`), parameter vector (e.g. `parameters = c(0,1)`)
and a type of probability distribution (from the exponential family) in
`type`. The following types are available in the current package
version:

-   `type = "gauss"` (additional parameters `sdNoise` and `gamma`)

-   `type = "exp"`

-   `type = "poisson"`

-   `type = "geom"`

-   `type = "bern"`

-   `type = "binom"` (additional parameter `nbTrials`)

-   `type = "negbin"` (additional parameter `nbSuccess`)

-   `type = "variance"`

We show two data examples with Gaussian and Exponential models
(`"gauss"` and `"exp"`)

<center>
<img src="man/figures/cost1.png" alt="" style="width:80%;"/>
</center>
and some other examples with integer-valued cost (`"poisson"`,`"geom"`,
`"binom"`, `"negbin"`):
<center>
<img src="man/figures/cost2.png" alt="" style="width:80%;"/>
</center>

**dataGenerator\_MD** concatenates `p` copies of `dataGenerator_1D`
function.

Additional information and examples are easily accessible in the help of
these functions (e.g. run `?dataGenerator_MD`).

### Data Generators in 2D

**dataGenerator\_meanVar** is used for change in mean and variance for
the Gaussian problem

**dataGenerator\_Reg** generates 2-dimensional data frame `(x,y)`
following a simple linear regression link (`y = Ax + B + noise`) with
`A` and `B` changing over time (after each change-point)

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="dust1D"></a>

## dust 1D and MD Algorithms

One simple example with the `variance` model:

    library(dust)
    set.seed(2)
    data <- dataGenerator_1D(chpts = c(50,100), parameters = c(2,1), type = "variance")
    dust.1D(data = data, model = "variance")

    ## $changepoints
    ## [1]  46 100
    ## 
    ## $lastIndexSet
    ## [1] 97 62 46  0
    ## 
    ## $nb
    ##   [1] 1 2 2 1 2 2 2 3 3 3 2 3 3 2 2 3 3 4 3 4 4 3 2 3 2 3 3 2 2 3 3 4 2 3 3 3 3
    ##  [38] 4 5 5 5 4 4 3 3 3 3 3 3 3 3 2 3 3 3 3 4 4 4 4 4 2 3 4 4 4 3 4 5 4 4 5 6 6
    ##  [75] 7 5 4 4 4 4 4 4 5 4 4 5 5 5 5 5 5 6 5 4 4 3 4 4 4 4
    ## 
    ## $costQ
    ##   [1]   1.084352   1.517155   3.749380   5.068406   5.781117   6.401682
    ##   [7]   7.286523   7.837598  10.816833  11.502406  12.231993  13.410460
    ##  [13]  14.104168  15.358092  17.759648  20.902828  22.014772  22.795887
    ##  [19]  24.005131  24.840711  27.469206  28.849124  30.647503  32.890257
    ##  [25]  33.750418  36.694660  37.672722  38.681489  39.775446  40.668828
    ##  [31]  41.718728  42.595164  43.865596  44.719010  45.768521  46.713655
    ##  [37]  48.706706  49.840570  50.771791  51.592334  52.438003  54.804666
    ##  [43]  55.891739  58.134740  59.103193  61.451178  62.317484  63.141338
    ##  [49]  63.965161  65.344485  66.217675  67.448376  68.273378  69.227719
    ##  [55]  70.121835  71.306815  72.087016  72.943729  73.841793  74.903528
    ##  [61]  76.009501  77.219493  78.017466  78.754792  79.511626  80.311312
    ##  [67]  81.504454  82.245029  82.951726  83.750042  84.549661  85.250900
    ##  [73]  85.934074  86.632574  87.300303  88.071311  88.956656  89.689978
    ##  [79]  90.489295  91.407495  92.185645  93.224911  93.899025  94.796064
    ##  [85]  96.117586  97.218507  97.907215  98.541455  99.078332  99.288808
    ##  [91] 100.382570 101.562323 102.230898 103.056630 104.036017 104.763639
    ##  [97] 106.266099 106.459378 106.890089 107.196194

The dust.1D function has the following main arguments:

-   `data`: the 1D time series to be segmented

-   `penalty`: the linear penalty value (often set to
    `2*sd(data)*log(n)`)

-   `model`: the data model to use: “gauss”, “poisson”, “exp”, “geom”,
    “bern”, “binom”, “negbin”, “variance”.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="dust2D"></a>

## dust 2D Algorithms

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="rcpp"></a>

## Rcpp Structure and new cost integration

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="pruning"></a>

## Pruning Capacity

Analysis of the pruning capacity (return field `nb`) for some of our
algorithms. We explore in particular the impact of the different choices
for the dual max evaluation.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

## Hidden Functions and Parameters (Package Development)

### OP in R

The base function `OP_R` is used to compute the change-point vector with
the simplest dynamic programming algorithm with no pruning. This method
is of quadratic time complexity. We propose 3 such functions:`OP_R_1D`,
`OP_R_MD`, `OP_R_2param`.

`OP_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss")`

Example:
`OP_R_1D(dataGenerator_1D(chpts = c(200,400), c(0,1), type = "gauss"))`

`OP_R_2param` is used for:

-   `type = "meanVar"` change in Gaussian data in mean and variance

-   `type = regression` change in simple regression model

### Dual Functions

`dual_1D` returns the value of the dual at a point `mu` when comparing
index `s1` with the constraint from index `s2` at time `t`. With option
`OP = TRUE` the optimal partitioning algorithm is used to have the true
constants in the cost functions with penalty `penalty` and a pruning
option `pruningOpt`.

`dual_1D <- function(mu, data, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(data), pruningOpt = 3)`

-   `data` is raw data

-   If `OP` is `true`, we run the OP algorithm to have the optimal cost
    vector to use in cost functions. See the function `OP_R`.

-   at time `t`, we evaluate the dual function at point `mu` when trying
    to remove index `s1` using function linked to index `s2` (we have a
    unique constraint, which means that the dual is a one-parametric
    function)

-   Depending on the `type`, different functions `A`, `B`, `statistic`,
    `mu_max` and `evalDual` are used (see the code in file
    `functions_by_type.R`)

Function `dual_1D` allows us to study the shape of the dual.

### dust\_R

We propose a few R functions computing the change-point location with
dust method: `dust_R_1D`, `dust_R_MD`, `dust_R_2Dquad`.

The function `dust_R_1D` has the following parameters:

`dust_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss", pruningOpt = 2)`

We have different type of possible pruning:

-   `pruningOpt == 0`: nothing

-   `pruningOpt == 1`: PELT

-   `pruningOpt == 2`: dust

-   `pruningOpt == 3`: dust + PELT

and returns a list of two elements:

-   `changepoints`: the change-points found by the algo

-   `nb`: the number of indices to consider in the minimization at each
    time step

-   `lastIndexSet`: the vecotr of indices saved by the algo in the
    dynamic programming algorithm at the last iteration

-   `costQ`: the vector of optimal cost (of size `length(data)`)

### Plot functions

`plot_dual_1D` is a function displaying the result of `dual_1D` for a
vector of mu values.

`plot_dual_1D <- function(mu =  (1:99)/100,` `data, s1, s2,`
`type = "gauss",` `OP = FALSE,` `penalty = 2*length(data))`

We use the function `plot_dual_1D` with `OP = TRUE` to plot the true
dual function seen by the dynamic programming algorithm.

What we called the “pruning interval” is the interval of values between
the vertical green lines for which the dual function takes a value
higher than the pruning threshold (horizontal line in red), so that the
index considered `s1` is pruned by `s2` at time `n`.

Using function `barplot_dual_1D` we can repeat the generation of the
pruning interval `nb` and count the number of time each value mu is in
this interval.

We add the values in the bar plot only if at the final time step `n`,
the index `s1` has not been removed by the algorithm (the pruning option
is given by option `pruningOpt`).

`barplot_dual_1D <- function(nb = 1000, s1 = 18, s2 = 15,` `n = 20,`
`oneParam = 0,` `type = "gauss",` `penalty = 2*log(n),`
`pruningOpt = 3)`

[(Back to Top)](#top)
