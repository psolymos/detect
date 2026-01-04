# Simulated example for abundance model

Simulated example for abundance model, see code below.

## Usage

``` r
data(databu)
```

## Format

A data frame with 1000 observations on the following 11 variables.

- `N`:

  true counts

- `Y`:

  observed counts

- `x1`:

  random variables used as covariates

- `x2`:

  random variables used as covariates

- `x3`:

  random variables used as covariates

- `x4`:

  random variables used as covariates

- `x5`:

  random variables used as covariates

- `x6`:

  random variables used as covariates

- `p`:

  probability of detection

- `lambda`:

  mean of the linear predictor

- `A`:

  occupancy

- `phi`:

  zero inflation probabilities

## Details

This simulated example corresponds to the Binomial - ZIP model
implemented in the function
[`svabu`](https://peter.solymos.org/detec/reference/svabu.md).

## Source

Simulated example.

## References

Solymos, P., Lele, S. R. and Bayne, E. 2012. Conditional likelihood
approach for analyzing single visit abundance survey data in the
presence of zero inflation and detection error. *Environmetrics*,
**23**, 197–205. \<doi:10.1002/env.1149\>

## Examples

``` r
data(databu)
str(databu)
#> 'data.frame':    1000 obs. of  12 variables:
#>  $ N     : num  8 12 3 3 0 0 11 8 9 5 ...
#>  $ Y     : num  7 2 3 2 0 0 0 8 6 4 ...
#>  $ x1    : num  0.114 0.622 0.609 0.623 0.861 ...
#>  $ x2    : num  0.985 -1.225 0.71 -0.109 1.783 ...
#>  $ x3    : num  0.927 -0.586 -0.828 -0.568 -0.521 ...
#>  $ x4    : num  -0.6699 0.4302 -0.0794 0.2769 -0.0882 ...
#>  $ x5    : num  1 1 0 0 1 0 1 1 1 0 ...
#>  $ x6    : num  0 0 1 0 1 0 1 1 1 1 ...
#>  $ p     : num  0.922 0.125 0.918 0.686 0.983 ...
#>  $ lambda: num  11.12 7.41 4.54 4.49 6.12 ...
#>  $ A     : num  1 1 1 1 0 0 1 1 1 1 ...
#>  $ phi   : num  0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
if (FALSE) { # \dontrun{
## simulation
n <- 1000
set.seed(1234)
x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)
x3 <- runif(n,-1,1)
x4 <- runif(n,-1,1)
x5 <- rbinom(n,1,0.6)
x6 <- rbinom(n,1,0.4)
x7 <- rnorm(n,0,1)
X <- model.matrix(~ x1 + x5)
Z <- model.matrix(~ x2 + x5)
Q <- model.matrix(~ x7)
beta <- c(2,-0.8,0.5)
theta <- c(1, 2, -0.5)
phi <- 0.3
p <- drop(binomial("logit")$linkinv(Z %*% theta))
lambda <- drop(exp(X %*% beta))
A <- rbinom(n, 1, 1-phi)
N <- rpois(n, lambda * A)
Y <- rbinom(n, N, p)
databu <- data.frame(N=N, Y=Y, x1, x2, x3, x4, x5, x6, p=p, lambda=lambda, A, phi)
} # }
```
