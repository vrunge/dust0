
# dust Vignette

### Vincent Runge
#### LaMME, Evry University
### October 11, 2023


> [Introduction](#intro)

> [Functions in R code](#Rcode)

> [DuST Algorithms](#dust)


## Introduction

The `dust` package contains methods for detecting multiple change-points within time-series based on the optimal partitioning algorithm. A few models from the exponential family are considered (Gauss, Poisson, Exponential...).

The proposed algorithm is a pruned dynamic programming algorithm optimizing a penalized likelihood *using an original pruning rule*, different from PELT or FPOP. We called this method, the DuST pruning rule, standing for Duality Sample Test.

Indeed, indices for potential last change-point are discarded by considering some constrained optimization problems. For each potential last change-point index, evaluating its associated dual function at a random testing point enables a fast and efficient test.

A

A

A
A


A

A

A

A

A

A

A

A
A
A
A
A
A


[Back to Top](#top)

<a id="Rcode"></a>

## Functions in R code

A
A
A
A
A

A

A

A

A

A

a

A

A

[Back to Top](#top)

<a id="dust"></a>

## DuST Algorithms


A

A

A
A


A

A

A

A

A

A

A

A
A
A
A
A
A


A

A

A
A


A

A

A

A

A

A

A

A
A
A
A
A
A


A

A

A
A


A

A

A

A

A

A

A

A
A
A
A
A
A


A

A

A
A


A

A

A

A

A

A

A

A
A
A
A
A
A

[Back to Top](#top)

