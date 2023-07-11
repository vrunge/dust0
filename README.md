
# dust Vignette
### Vincent Runge
#### LaMME, Evry University
### June 19, 2023

> [Introduction](#intro)

The `dust` package contains methods for detecting change-points within univariate time-series based on the optimal partitioning algorithm. A few models from the exponential family are considered (Gauss, Poisson, Exponential...).
The proposed algorithm is a pruned dynamic programming algorithm with an original pruning rule, different from PELT or FPOP. 

Indices for potential last change-point are discarded by considering some constrained optimization problems. Evaluating the dual function at a random testing point (for each of these problems) enables a quick and efficient test.

[Back to Top](#top)

