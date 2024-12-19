# MOODE


Multi-objective Optimal Design of experiments (MOODE) for targeting the
experimental objectives directly, ensuring as such that the full set of
research questions is answered as economically as possible.

## Installation

You can install the development version of `MOODE` from
[GitHub](https://github.com/vkstats/MOODE) with:

``` r
# install.packages("devtools")
devtools::install_github("vkstats/MOODE")
```

## Example

As a basic example, consider an experiment with `K=2` factors, each
having `Levels = 3` levels. The primary (assumed) model contains
first-order terms, and the potential model also contains squared terms.
The experiment will have `Nruns = 24` runs. An optimal compound design
will be sought combining $DP_S$-, $LoF$-$D$- and $MSE(D)$-optimality;
see Egorova (2017, PhD Thesis). We define the parameters for this
expriment using the `mood` function.

``` r
library("MOODE")
ex.mood <- mood(K = 2, Levels = 3, Nruns = 24, 
                model_terms = list(primary.terms = c("x1", "x2"), 
                                   potential.terms = c("x12", "x22")), 
                criterion.choice = "MSE.D", 
                kappa = list(kappa.DP = 1 / 3, kappa.LoF = 1 / 3, 
                             kappa.mse = 1 / 3))
```

The `kappa` list defines weights for each criterion, with
$\kappa_i\ge 0$ and $\sum \kappa_i = 1$.

Optimal designs are found using a point exchange algorithm, via the
`Search` function.

``` r
search.ex <- Search(ex.mood)
```

    ✔ Design search complete. Final compound objective function value = 0.19728

The best design found is available as element `X.design`, ordered here
by treatment number.

``` r
fd <- search.ex$X.design[order(search.ex$X1[, 1]),]
cbind(fd[1:12, ], fd[13:24, ])
```

          x1 x2 x1 x2
     [1,] -1 -1  0  1
     [2,] -1 -1  0  1
     [3,] -1 -1  1 -1
     [4,] -1 -1  1 -1
     [5,] -1  0  1 -1
     [6,] -1  0  1 -1
     [7,] -1  1  1  0
     [8,] -1  1  1  0
     [9,] -1  1  1  1
    [10,]  0 -1  1  1
    [11,]  0 -1  1  1
    [12,]  0  0  1  1

The `path` element records the compound objective function value from
each of the (by default) 10 attempts of the algorithm from different
random starting designs.

``` r
search.ex$path
```

     [1] 0.1979960 0.1971942 0.1979960 0.1971696 0.1971662 0.1990148 0.1971457
     [8] 0.1979960 0.1971299 0.1971564
