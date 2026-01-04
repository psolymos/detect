# Simulated example for occupancy model

Simulated example for occupancy model, see code below.

## Usage

``` r
data(datocc)
```

## Format

A data frame with 1000 observations on the following 6 variables.

- `Y`:

  true occupancy

- `W`:

  observations

- `x1`:

  random variables used as covariates

- `x2`:

  random variables used as covariates

- `x3`:

  random variables used as covariates

- `x4`:

  random variables used as covariates

- `p.occ`:

  probability of occurrence

- `p.det`:

  probability of detection

## Details

This simulated example corresponds to the ZI Binomial model implemented
in the function
[`svocc`](https://peter.solymos.org/detec/reference/svocc.md).

## Source

Simulated example.

## References

Lele, S.R., Moreno, M. and Bayne, E. (2012) Dealing with detection error
in site occupancy surveys: What can we do with a single survey? *Journal
of Plant Ecology*, **5(1)**, 22–31. \<doi:10.1093/jpe/rtr042\>

## Examples

``` r
data(datocc)
str(datocc)
#> 'data.frame':    1000 obs. of  8 variables:
#>  $ Y    : num  1 1 1 1 1 1 0 0 0 1 ...
#>  $ W    : num  1 0 0 0 1 1 0 0 0 1 ...
#>  $ x1   : num  -0.773 0.245 0.219 0.247 0.722 ...
#>  $ x2   : Factor w/ 2 levels "0","1": 2 1 1 1 2 1 1 1 2 2 ...
#>  $ x3   : num  -1.205 0.301 -1.539 0.635 0.703 ...
#>  $ x4   : num  -0.9738 -0.0996 -0.1107 1.1922 -1.6559 ...
#>  $ p.occ: num  0.71 0.872 0.869 0.873 0.927 ...
#>  $ p.det: num  0.605 0.591 0.457 0.615 0.562 ...
if (FALSE) { # \dontrun{
## simulation
n <- 1000
set.seed(1234)
x1 <- runif(n, -1, 1)
x2 <- as.factor(rbinom(n, 1, 0.5))
x3 <- rnorm(n)
x4 <- rnorm(n)
beta <- c(0.6, 0.5)
theta <- c(0.4, -0.5, 0.3)
X <- model.matrix(~ x1)
Z <- model.matrix(~ x1 + x3)
mu <- drop(X %*% beta)
nu <- drop(Z %*% theta)
p.occ <- binomial("cloglog")$linkinv(mu)
p.det <- binomial("logit")$linkinv(nu)
Y <- rbinom(n, 1, p.occ)
W <- rbinom(n, 1, Y * p.det)
datocc <- data.frame(Y, W, x1, x2, x3, x4, p.occ, p.det)
} # }
```
