# Single-bin QPAD (SQPAD)

Single bin QPAD (SQPAD) approach by Lele and Solymos (2025).

## Usage

``` r
sqpad(formula, data, dis, dur,
    type = c("full", "approx"), det = c("joint", "pq"),
    Nmax = NULL, K = NULL, A = NULL,
    montecarlo = FALSE, np = 1000, distcorr = 2/3, ...)

sqpad.fit(Y, dis, dur, X = NULL, Z = NULL, 
    type = c("full", "approx", "conv"), det = c("joint", "pq"),
    init = NULL, method = "Nelder-Mead", hessian = TRUE,
    tt1 = NULL, A = NULL, Nmax = NULL, K = NULL, 
    montecarlo = FALSE, np = 1000, distcorr = 2/3, dislist = NULL, ...)

# S3 method for class 'sqpad'
print(x, digits, ...)
# S3 method for class 'sqpad'
vcov(object, ...)
# S3 method for class 'sqpad'
fitted(object, ...)
# S3 method for class 'sqpad'
logLik(object, ...)
# S3 method for class 'sqpad'
summary(object, ...)
# S3 method for class 'summary.sqpad'
print(x, digits, ...)
# S3 method for class 'sqpad'
predict(object, newdata = NULL, 
    type = c("link", "response"), ...)

rtriang(n, r = 1)
```

## Arguments

- formula:

  formula, LHS takes a vector, RHS is either `1` or some covariates, see
  Examples.

- data:

  data.

- dis, dur:

  a vector with distance and duration values.

- A:

  a vector with area values. Area is calculated from `dis`, use `A` if
  distance is not directly related to area (e.g. for ARU data).

- type:

  character, one of `"full"` (full likelihood), `"approx"` (approximate
  method), and `"conv"` (convolution likelihood).

- det:

  character, the detection model: `"joint"` or independence (`"pq"`).

- Y:

  response vector, non-negative integers.

- X, Z:

  design matrices, `X` is the matrix with covariates for the density
  parameters. `Z` is the matrix with covariates for the cue rate model.
  `NULL` means constant density and detection.

- K:

  truncation value, `K=1` implements the occupancy version of SQPAD.
  `K>1` implements the truncated count models.

- Nmax:

  maximum abundance value for numerical integration under the full and
  convolution likelihoods.

- montecarlo, np, distcorr:

  should mean probability (`TRUE`) using `np` points or approximation be
  used with distance correction `distcorr` in the detection model.

- dislist:

  distance list for the convolution likelihood.

- tt1:

  vector with time to 1st detection values, same units as for `dur`.

- init:

  optional initial values for
  [`optim`](https://rdrr.io/r/stats/optim.html).

- method:

  method for [`optim`](https://rdrr.io/r/stats/optim.html).

- hessian:

  logical, should the Hessian matrix be computed by
  [`optim`](https://rdrr.io/r/stats/optim.html).

- object, x:

  fitted model object.

- newdata:

  optionally, a data frame in which to look for variables with which to
  predict. If omitted, the fitted linear predictors are used.

- digits:

  digits to use when providing summaries.

- ...:

  additional options for [`optim`](https://rdrr.io/r/stats/optim.html)
  or the methods.

- n:

  number of random variates to generate from the triangular
  distribution.

- r:

  maximum value for the triangular distribution.

## Details

Single bin QPAD (SQPAD) approach for robust analysis of point count data
with detection error by Lele and Solymos (2025).

## Value

An object of class 'sqpad'.

## References

Solymos, P., Lele, S. R., 2025. Single bin QPAD (SQPAD) approach for
robust analysis of point count data with detection error.
*Ornithological Applications*, **xx**, xx–xx.
\<doi:10.1093/ornithapp/duaf078\>

Supporting info for the SQPAD method:
<https://github.com/psolymos/sqpad-paper>,
\<doi:10.5281/zenodo.16172209\>

## Author

Peter Solymos and Subhash Lele

## Examples

``` r

set.seed(0)
n <- 100
x <- rnorm(n)
D <- exp(-2 + 0.5 * x)
phi <- 0.25
tau <- 1
dur <- sample(1:10, n, replace=TRUE)
dis <- sample(seq(0.5, 2, 0.25), n, replace=TRUE)
A <- dis^2 * pi
dcorr <- 2/3
p <- 1 - exp(-dur * phi * exp(-(dis*dcorr)^2/tau^2))
N <- rpois(n, D*A)
Y <- rbinom(n, N, p)

df <- data.frame(x = x, y = Y)

m <- sqpad(y ~ x | 1, data = df, dis = dis, dur = dur, type = "full", det = "joint", K = NULL)

print(m)
#> 
#> Call:
#> sqpad(formula = y ~ x | 1, data = df, dis = dis, dur = dur, type = "full", 
#>     det = "joint", K = NULL)
#> 
#> Single-bin QPAD model with dependence
#> Full Maximum Likelihood estimates
#> 
#> Coefficients:
#>     log.D_(Intercept)                log.D_x  log.delta_(Intercept)  
#>               -2.6490                 0.5793                 0.6531  
#>               log.phi  
#>               -0.2066  
#> 
summary(m)
#> 
#> Call:
#> sqpad(formula = y ~ x | 1, data = df, dis = dis, dur = dur, type = "full", 
#>     det = "joint", K = NULL)
#> 
#> Single-bin QPAD model with dependence
#> Full Maximum Likelihood estimates
#> 
#> Coefficients:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> log.D_(Intercept)      -2.6490     0.2559 -10.354  < 2e-16 ***
#> log.D_x                 0.5793     0.2028   2.857  0.00428 ** 
#> log.delta_(Intercept)   0.6531     1.2882   0.507  0.61218    
#> log.phi                -0.2066     0.2579  -0.801  0.42313    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Log-likelihood: -61.23 
#> BIC = 140.9 
#> 

coef(m)
#>     log.D_(Intercept)               log.D_x log.delta_(Intercept) 
#>            -2.6490003             0.5792880             0.6530735 
#>               log.phi 
#>            -0.2065773 
nobs(m)
#> [1] 100
vcov(m)
#>                       log.D_(Intercept)      log.D_x log.delta_(Intercept)
#> log.D_(Intercept)           0.065460969 -0.016823675           -0.11176285
#> log.D_x                    -0.016823675  0.041119454           -0.01041559
#> log.delta_(Intercept)      -0.111762852 -0.010415591            1.65946167
#> log.phi                     0.005060442  0.002226058           -0.29981180
#>                            log.phi
#> log.D_(Intercept)      0.005060442
#> log.D_x                0.002226058
#> log.delta_(Intercept) -0.299811804
#> log.phi                0.066510636
confint(m)
#>                            2.5 %     97.5 %
#> log.D_(Intercept)     -3.1504638 -2.1475368
#> log.D_x                0.1818478  0.9767282
#> log.delta_(Intercept) -1.8717539  3.1779010
#> log.phi               -0.7120453  0.2988906

logLik(m)
#> 'log Lik.' -61.23112 (df=4)
AIC(m)
#> [1] 130.4622
BIC(m)
#> [1] 140.8829

fitted(m)
#>          1          2          3          4          5          6          7 
#> 0.14699069 0.05854361 0.15279417 0.14779970 0.08992282 0.02898216 0.04129937 
#>          8          9         10         11         12         13         14 
#> 0.05962214 0.07048600 0.28478709 0.11006769 0.04451822 0.03637679 0.05980405 
#>         15         16         17         18         19         20         21 
#> 0.05946710 0.05572181 0.08184809 0.04218547 0.09102563 0.03453122 0.06210579 
#>         22         23         24         25         26         27         28 
#> 0.08800342 0.07640095 0.11268680 0.06842058 0.09467870 0.13265186 0.04739391 
#>         29         30         31         32         33         34         35 
#> 0.03360255 0.07266231 0.06169562 0.05163844 0.05502257 0.04854658 0.10774345 
#>         36         37         38         39         40         41         42 
#> 0.13783311 0.12565016 0.05514374 0.14490664 0.06015551 0.19579871 0.09786495 
#>         43         44         45         46         47         48         49 
#> 0.05440536 0.04367441 0.03598041 0.03814792 0.02858479 0.13820291 0.11452005 
#>         50         51         52         53         54         55         56 
#> 0.06199577 0.08251047 0.05685679 0.29090833 0.04461297 0.06850900 0.08174943 
#>         57         58         59         60         61         62         63 
#> 0.10117947 0.06399188 0.01950120 0.03401353 0.08705693 0.07027081 0.04101132 
#>         64         65         66         67         68         69         70 
#> 0.06613240 0.04410854 0.08137721 0.03097599 0.08742141 0.08166761 0.07344785 
#>         71         72         73         74         75         76         77 
#> 0.07151105 0.08209097 0.04855957 0.06600443 0.10390539 0.13382502 0.07686419 
#>         78         79         80         81         82         83         84 
#> 0.06605857 0.04169598 0.03075272 0.04456775 0.14623724 0.11061411 0.06227700 
#>         85         86         87         88         89         90         91 
#> 0.05529417 0.05548123 0.12600196 0.06027998 0.14640132 0.10285967 0.15011940 
#>         92         93         94         95         96         97         98 
#> 0.04264392 0.07106565 0.04245635 0.09989910 0.07580058 0.06005705 0.16438177 
#>         99        100 
#> 0.08075527 0.12596964 
predict(m, type = "link")
#>         1         2         3         4         5         6         7         8 
#> -1.917386 -2.837983 -1.878664 -1.911897 -2.408804 -3.541075 -3.186908 -2.819728 
#>         9        10        11        12        13        14        15        16 
#> -2.652341 -1.256013 -2.206660 -3.111857 -3.313824 -2.816682 -2.822332 -2.887384 
#>        17        18        19        20        21        22        23        24 
#> -2.502890 -3.165680 -2.396614 -3.365891 -2.778916 -2.430380 -2.571760 -2.183143 
#>        25        26        27        28        29        30        31        32 
#> -2.682082 -2.357266 -2.020027 -3.049262 -3.393153 -2.621932 -2.785542 -2.963489 
#>        33        34        35        36        37        38        39        40 
#> -2.900012 -3.025231 -2.228002 -1.981712 -2.074254 -2.897812 -1.931666 -2.810822 
#>        41        42        43        44        45        46        47        48 
#> -1.630668 -2.324167 -2.911293 -3.130993 -3.324781 -3.266284 -3.554881 -1.979032 
#>        49        50        51        52        53        54        55        56 
#> -2.167005 -2.780689 -2.494830 -2.867220 -1.234747 -3.109731 -2.680790 -2.504096 
#>        57        58        59        60        61        62        63        64 
#> -2.290859 -2.748999 -3.937279 -3.380997 -2.441193 -2.655399 -3.193907 -2.716097 
#>        65        66        67        68        69        70        71        72 
#> -3.121102 -2.508660 -3.474543 -2.437015 -2.505098 -2.611180 -2.637903 -2.499927 
#>        73        74        75        76        77        78        79        80 
#> -3.024964 -2.718033 -2.264274 -2.011222 -2.565715 -2.717214 -3.177351 -3.481777 
#>        81        82        83        84        85        86        87        88 
#> -3.110745 -1.922525 -2.201708 -2.776163 -2.895088 -2.891710 -2.071458 -2.808755 
#>        89        90        91        92        93        94        95        96 
#> -1.921404 -2.274390 -1.896324 -3.154871 -2.644151 -3.159279 -2.303595 -2.579649 
#>        97        98        99       100 
#> -2.812460 -1.805564 -2.516332 -2.071714 
predict(m, newdata = df[1:10,], type = "response")
#>          1          2          3          4          5          6          7 
#> 0.14699069 0.05854361 0.15279417 0.14779970 0.08992282 0.02898216 0.04129937 
#>          8          9         10 
#> 0.05962214 0.07048600 0.28478709 
predict(m, newdata = df[1:10,], type = "link")
#>         1         2         3         4         5         6         7         8 
#> -1.917386 -2.837983 -1.878664 -1.911897 -2.408804 -3.541075 -3.186908 -2.819728 
#>         9        10 
#> -2.652341 -1.256013 

if (FALSE) { # \dontrun{
m0 <- sqpad(y ~ 1 | 1, data = df, dis = dis, dur = dur, type = "full", det = "joint", K = NULL)
m1 <- sqpad(y ~ 1 | x, data = df, dis = dis, dur = dur, type = "full", det = "joint", K = NULL)
m2 <- sqpad(y ~ x | x, data = df, dis = dis, dur = dur, type = "full", det = "joint", K = NULL)

AIC(m, m0, m1, m2) # m2 is best
BIC(m, m0, m1, m2) # this is needed!

# average probability

set.seed(5)
n <- 1000
x <- rnorm(n)
D <- exp(-2 + 0.5 * x)
phi <- 0.25
tau <- 1
dur <- sample(1:10, n, replace=TRUE)
dis <- sample(seq(0.5, 2, 0.25), n, replace=TRUE)
A <- dis^2 * pi

dcorr <- 2/3
p <- 1 - exp(-dur * phi * exp(-(dis*dcorr)^2/tau^2))
N <- rpois(n, D*A)
Y <- rbinom(n, N, p)
df <- data.frame(x = x, y = Y)

np <- 1000
dij <- rtriang(as.integer(np), 1)
dij <- matrix(dij, nrow = n, ncol = np, byrow = TRUE)
p2 <- rowMeans(1 - exp(-dur * phi * exp(-(dis*dij)^2/tau^2)))
Y2 <- rbinom(n, N, p2)
df2 <- data.frame(x = x, y = Y2)

plot(p, p2)
abline(0, 1)

mAP <- sqpad(y ~ x | 1, data = df, dis = dis, dur = dur, type = "full", 
    det = "joint", K = NULL, montecarlo = FALSE, hessian = FALSE)
mMC <- sqpad(y ~ x | 1, data = df, dis = dis, dur = dur, type = "full", 
    det = "joint", K = NULL, montecarlo = TRUE, hessian = FALSE)
m2AP <- sqpad(y ~ x | 1, data = df2, dis = dis, dur = dur, type = "full", 
    det = "joint", K = NULL, montecarlo = FALSE, hessian = FALSE)
m2MC <- sqpad(y ~ x | 1, data = df2, dis = dis, dur = dur, type = "full", 
    det = "joint", K = NULL, montecarlo = TRUE, hessian = FALSE)

cbind(mAP=coef(mAP), mMC=coef(mMC), m2AP=coef(m2AP), m2MC=coef(m2MC))

# convolution likelihood

D <- 1
phi <- 0.25
tau <- 1.2
n <- 200

dur <- sample(1:10, n, replace=TRUE)
dis <- sample(seq(0.25, 2, 0.25), n, replace=TRUE)
A <- dis^2 * pi
N <- rpois(n, D*A)
dislist <- lapply(seq_len(n), function(i) {
    dij <- rtriang(N[i], dis[i])
    # pij is based on time of 1st detection and distance at 1st detection
    pij <- 1 - exp(-dur[i] * phi * exp(-(dij)^2/tau^2))
    Yij <- rbinom(length(pij), 1, pij)
    dij[Yij > 0]
})
Y <- sapply(dis, length)

m4 <- sqpad.fit(Y = Y, dis = dis, dur = dur, dislist = dislist, 
    type = "conv", det = "joint", Nmax = 20)
cbind(parameters = c(D = D, phi = phi, tau = tau), estimates = exp(m4$coef))
} # }
```
