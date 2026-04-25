# ZI Binomial model with single visit

ZI Binomial model with single visit

## Usage

``` r
svocc(formula, data, link.sta = "cloglog", link.det = "logit",
    penalized = FALSE, method = c("optim", "dc"), inits,
    model = TRUE, x = FALSE, ...)
svocc.fit(Y, X, Z, link.sta = "cloglog", link.det = "logit",
    penalized = FALSE, auc = FALSE, method = c("optim", "dc"),
    inits, ...)

extractMLE(object, ...)
svocc.step(object, model, trace = 1, steps = 1000,
    criter = c("AIC", "BIC", "cAUC"), test = FALSE, k = 2,
    control, ...)
```

## Arguments

- formula:

  formula of the form `y ~ x | z`, where `y` is a vector of
  observations, `x` is the set of covariates for the occurrence model,
  `z` is the set of covariates for the detection model

- Y, X, Z:

  vector of observation, design matrix for occurrence model, and design
  matrix for detection model

- data:

  data

- link.sta, link.det:

  link function for the occurrence (true state) and detection model

- penalized:

  logical, if penalized likelihood estimate should be computed

- method:

  optimization or data cloning to be used as optimization

- inits:

  initial values

- model:

  a logical value indicating whether model frame should be included as a
  component of the returned value, or true state or detection model

- x:

  logical values indicating whether the response vector and model matrix
  used in the fitting process should be returned as components of the
  returned value

- auc:

  logical, if AUC should be calculated

- object:

  a fitted model object

- trace:

  info returned during the procedure

- steps:

  max number of steps

- criter:

  criterion to be minimized (cAUC=1-AUC)

- test:

  logical, if decrease in deviance should be tested

- k:

  penalty to be used with AIC

- control:

  controls for optimization, if missing taken from object

- ...:

  other arguments passed to the functions

## Details

See Examples.

The right hand side of the formula must contain at least one continuous
(i.e. non discrete/categorical) covariate. This is the necessary
condition for the single-visit method to be valid and parameters to be
identifiable. See References for more detailed description.

## Value

An object of class 'svocc'.

## References

Lele, S.R., Moreno, M. and Bayne, E. 2012. Dealing with detection error
in site occupancy surveys: What can we do with a single survey? *Journal
of Plant Ecology*, **5(1)**, 22–31. \<doi:10.1093/jpe/rtr042\>

Moreno, M. and Lele, S. R. 2010. Improved estimation of site occupancy
using penalized likelihood. *Ecology*, **91**, 341–346.
\<doi:10.1890/09-1073.1\>

Solymos, P., Lele, S. R. 2016. Revisiting resource selection probability
functions and single-visit methods: clarification and extensions.
*Methods in Ecology and Evolution*, **7**, 196–205.
\<doi:10.1111/2041-210X.12432\>

## Author

Peter Solymos and Monica Moreno

## Examples

``` r
data(datocc)
## MLE
m00 <- svocc(W ~ x1 | x1 + x3, datocc)
## PMLE
m01 <- svocc(W ~ x1 | x1 + x3, datocc, penalized=TRUE)

## print
m00
#> 
#> Call:
#> svocc(formula = W ~ x1 | x1 + x3, data = datocc)
#> 
#> Single visit site-occupancy model
#> Maximum Likelihood estimates (optim method)
#> 
#> Coefficients for occurrence (cloglog link):
#> (Intercept)           x1  
#>     -0.1323       0.1952  
#> Coefficients for detection (logit link):
#> (Intercept)           x1           x3  
#>      1.7159      -0.8027       0.6823  
#> 
## summary
summary(m00)
#> 
#> Call:
#> svocc(formula = W ~ x1 | x1 + x3, data = datocc)
#> 
#> 
#> Single visit site-occupancy model
#> Maximum Likelihood estimates (optim method)
#> 
#> Occupancy model coefficients with cloglog link:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  -0.1323     0.2473  -0.535    0.593
#> x1            0.1952     0.2743   0.711    0.477
#> Detection model coefficients with logit link:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)   1.7159     1.0181   1.685   0.0919 .
#> x1           -0.8027     0.6527  -1.230   0.2187  
#> x3            0.6823     0.4615   1.479   0.1392  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Log-likelihood: -685.3 on 5 Df
#> AIC =  1381 
#> 
## coefficients
coef(m00)
#> sta_(Intercept)          sta_x1 det_(Intercept)          det_x1          det_x3 
#>      -0.1322958       0.1951644       1.7158595      -0.8026997       0.6823315 
## state (occupancy) model estimates
coef(m00, "sta")
#> (Intercept)          x1 
#>  -0.1322958   0.1951644 
## detection model estimates
coef(m00, "det")
#> (Intercept)          x1          x3 
#>   1.7158595  -0.8026997   0.6823315 
## compare estimates
cbind(truth=c(0.6, 0.5, 0.4, -0.5, 0.3),
mle=coef(m00), pmle=coef(m01))
#>                 truth        mle        pmle
#> sta_(Intercept)   0.6 -0.1322958  0.02430915
#> sta_x1            0.5  0.1951644  0.25655519
#> det_(Intercept)   0.4  1.7158595  1.17728804
#> det_x1           -0.5 -0.8026997 -0.65455395
#> det_x3            0.3  0.6823315  0.48569843

## AIC, BIC
AIC(m00)
#> [1] 1380.536
BIC(m00)
#> [1] 1405.075
## log-likelihood
logLik(m00)
#> 'log Lik.' -685.2682 (df=5)
## variance-covariance matrix
vcov(m00)
#>                 sta_(Intercept)      sta_x1 det_(Intercept)      det_x1
#> sta_(Intercept)      0.06115389  0.05059639     -0.24356451 -0.02790893
#> sta_x1               0.05059639  0.07525105     -0.18873669 -0.12264430
#> det_(Intercept)     -0.24356451 -0.18873669      1.03655409  0.02910525
#> det_x1              -0.02790893 -0.12264430      0.02910525  0.42595958
#> det_x3              -0.10172898 -0.09312257      0.42897959  0.06620209
#>                      det_x3
#> sta_(Intercept) -0.10172898
#> sta_x1          -0.09312257
#> det_(Intercept)  0.42897959
#> det_x1           0.06620209
#> det_x3           0.21295780
vcov(m00, model="sta")
#>             (Intercept)         x1
#> (Intercept)  0.06115389 0.05059639
#> x1           0.05059639 0.07525105
vcov(m00, model="det")
#>             (Intercept)         x1         x3
#> (Intercept)  1.03655409 0.02910525 0.42897959
#> x1           0.02910525 0.42595958 0.06620209
#> x3           0.42897959 0.06620209 0.21295780
## confidence intervals
confint(m00)
#>                       2.5%     97.5%
#> sta_(Intercept) -0.6169815 0.3523898
#> sta_x1          -0.3424914 0.7328203
#> det_(Intercept) -0.2796054 3.7113243
#> det_x1          -2.0818816 0.4764821
#> det_x3          -0.2221399 1.5868030
confint(m00, model="sta")
#>                   2.5%     97.5%
#> (Intercept) -0.6169815 0.3523898
#> x1          -0.3424914 0.7328203
confint(m00, model="det")
#>                   2.5%     97.5%
#> (Intercept) -0.2796054 3.7113243
#> x1          -2.0818816 0.4764821
#> x3          -0.2221399 1.5868030

## fitted values
## (conditional probability of occurrence given detection history:
## if W=1, fitted=1,
## if W=0, fitted=(phi*(1-delta)) / ((1-delta) + phi * (1-delta))
summary(fitted(m00))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.01269 0.17830 0.43880 0.58491 1.00000 1.00000 
## estimated probabilities: (phi*(1-delta)) / ((1-delta) + phi * (1-delta))
summary(m00$estimated.probabilities)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.01269 0.10232 0.17443 0.19693 0.27211 0.57908 
## probability of occurrence (phi)
summary(m00$occurrence.probabilities)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5137  0.5494  0.5850  0.5849  0.6206  0.6551 
## probability of detection (delta)
summary(m00$detection.probabilities)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2468  0.7581  0.8493  0.8200  0.9084  0.9887 

if (FALSE) { # \dontrun{
## model selection
m02 <- svocc(W ~ x1 | x3 + x4, datocc)
m03 <- drop1(m02, model="det")
## dropping one term at a time, resulting change in AIC
m03
## updating the model
m04 <- update(m02, . ~ . | . - x4)
m04
## automatic model selection
## part of the model (sta/det) must be specified
m05 <- svocc.step(m02, model="det")
summary(m05)

## nonparametric bootstrap
m06 <- bootstrap(m01, B=25)
attr(m06, "bootstrap")
extractBOOT(m06)
summary(m06, type="mle")
summary(m06, type="pmle") ## no SEs! PMLE!!!
summary(m06, type="boot")
## vcov
#vcov(m06, type="mle") ## this does not work with PMLE
vcov(m06, type="boot") ## this works
## confint
confint(m06, type="boot") ## quantile based

## parametric bootstrap
## sthis is how observations are simulated
head(simulate(m01, 5))
m07 <- bootstrap(m01, B=25, type="param")
extractBOOT(m07)
summary(m07)

data(oven)
ovenc <- oven
ovenc[, c(4:8,10:11)][] <- lapply(ovenc[, c(4:8,10:11)], scale)
ovenc$count01 <- ifelse(ovenc$count > 0, 1, 0)
moven <- svocc(count01 ~ pforest | julian + timeday, ovenc)
summary(moven)
drop1(moven, model="det")
moven2 <- update(moven, . ~ . | . - timeday)
summary(moven)

BIC(moven, moven2)
AUC(moven, moven2)
rocplot(moven)
rocplot(moven2, col=2, add=TRUE)
} # }
```
