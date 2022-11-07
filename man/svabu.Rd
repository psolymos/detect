\name{svabu}
\alias{svabu}
\alias{svabu.fit}
\alias{svabu_nb.fit}
\alias{zif}
\alias{is.present}
\alias{predictMCMC}
\alias{svabu.step}
\title{
Single visit N-mixture abundance models
}
\description{
Binomial-Poisson, Binomial-NegBin, Binomial-ZIP, and Binomial-ZINB models with single visit.
}
\usage{
svabu(formula, data, zeroinfl = TRUE, area = 1, N.max = NULL,
    inits, link.det = "logit", link.zif = "logit",
    model = TRUE, x = FALSE, distr = c("P", "NB"), ...)

svabu.fit(Y, X, Z, Q = NULL, zeroinfl = TRUE, area = 1, N.max = NULL,
    inits, link.det = "logit", link.zif = "logit", ...)
svabu_nb.fit(Y, X, Z, Q = NULL, zeroinfl = TRUE, area = 1, N.max = NULL,
    inits, link.det = "logit", link.zif = "logit", ...)

zif(x)
is.present(object, ...)
predictMCMC(object, ...)
svabu.step(object, model, trace = 1, steps = 1000,
    criter = c("AIC", "BIC"), test = FALSE, k = 2, control, ...)
}
\arguments{
  \item{formula}{
formula of the form \code{y ~ x | z}, where \code{y} is a vector of observations,
\code{x} is the set of covariates for the occurrence model,
\code{z} is the set of covariates for the detection model.
\code{x} can further expanded as \code{x1 + zif(x2)} into terms for
the nonzero count data part (\code{x1}) and the zero inflation
component (\code{zif(x2)}) using the \code{zif} special.
}
  \item{Y, X, Z, Q}{
vector of observation, design matrix for abundance model, design matrix for detection
and design matrix for zero inflation model
}
  \item{data}{
data
}
  \item{area}{
area
}
  \item{N.max}{
maximum of true count values (for calculating the integral)
}
  \item{zeroinfl}{
logical, if the Binomial-ZIP model should be fitted
}
  \item{inits}{
initial values used by \code{link{optim}}
}
  \item{link.det, link.zif}{
link function for the detection and zero inflation parts of the model
}
  \item{model}{
a logical value indicating whether model frame should be included as a component of the returned value,
or true state or detection model
}
  \item{x}{
logical values indicating whether the response vector and model matrix used in the
fitting process should be returned as components of the returned value.
For the function \code{zif} it is any object to be returned.
}
  \item{object}{
a fitted object.
}
  \item{trace}{
info returned during the procedure
}
  \item{steps}{
max number of steps
}
  \item{criter}{
criterion to be minimized (cAUC=1-AUC)
}
  \item{test}{
logical, if decrease in deviance should be tested
}
  \item{k}{
penalty to be used with AIC
}
  \item{control}{
controls for optimization, if missing taken from object
}
  \item{distr}{
character, abundance distribution: \code{"P"} for Poisson,
\code{"NB"} for Negative Binomial.
}
  \item{\dots}{
other arguments passed to the functions
}
}
\details{
See Examples.

The right hand side of the formula must contain at least one
continuous (i.e. non discrete/categorical) covariate.
This is the necessary condition for the single-visit
method to be valid and parameters to be identifiable.
See References for more detailed description.

The Binomial-Poisson model is the single visit special case
of the \emph{N}-mixture model proposed by Royle (2004)
and explained in Solymos et a. (2012) and Solymos and Lele (2016).
}
\value{
An object of class 'svabu'.
}
\references{
Royle, J. A. 2004.
\emph{N}-Mixture Models for Estimating Population Size from Spatially Replicated Counts.
\emph{Biometrics}, \bold{60(1)}, 108--115.
<doi:10.1111/j.0006-341X.2004.00142.x>

Solymos, P., Lele, S. R. and Bayne, E. 2012.
Conditional likelihood approach for analyzing single visit
abundance survey data in the presence of zero inflation and
detection error.
\emph{Environmetrics}, \bold{23}, 197--205.
<doi:10.1002/env.1149>

Solymos, P., Lele, S. R. 2016.
Revisiting resource selection probability functions and single-visit methods:
clarification and extensions.
\emph{Methods in Ecology and Evolution}, \bold{7}, 196--205.
<doi:10.1111/2041-210X.12432>

Denes, F., Solymos, P., Lele, S. R., Silveira, L. & Beissinger, S. 2017.
Biome scale signatures of land use change on raptor abundance:
insights from single-visit detection-based models.
\emph{Journal of Applied Ecology}, \bold{54}, 1268--1278.
<doi:10.1111/1365-2664.12818>
}
\author{
Peter Solymos and Subhash Lele
}
\examples{
data(databu)

## fit BZIP and BP models
m00 <- svabu(Y ~ x1 + x5 | x2 + x5, databu[1:200,])

## print method
m00
## summary: CMLE
summary(m00)
## coef
coef(m00)
coef(m00, model="sta") ## state (abundance)
coef(m00, model="det") ## detection
coef(m00, model="zif") ## zero inflation (this is part of the 'true state'!)

\dontrun{
## Diagnostics and model comparison

m01 <- svabu(Y ~ x1 + x5 | x2 + x5, databu[1:200,], zeroinfl=FALSE)
## compare estimates (note, zero inflation is on the logit scale!)
cbind(truth=c(2,-0.8,0.5, 1,2,-0.5, plogis(0.3)),
"B-ZIP"=coef(m00), "B-P"=c(coef(m01), NA))

## fitted
plot(fitted(m00), fitted(m01))
abline(0,1)

## compare models
AIC(m00, m01)
BIC(m00, m01)
logLik(m00)
logLik(m01)
## diagnostic plot
plot(m00)
plot(m01)

## Bootstrap

## non parametric bootstrap
## - initial values are the estimates
m02 <- bootstrap(m00, B=25)
attr(m02, "bootstrap")
extractBOOT(m02)
summary(m02)
summary(m02, type="cmle")
summary(m02, type="boot")
## vcov
vcov(m02, type="cmle")
vcov(m02, type="boot")
vcov(m02, model="sta")
vcov(m02, model="det")
## confint
confint(m02, type="cmle") ## Wald-type
confint(m02, type="boot") ## quantile based
## parametric bootstrap
simulate(m00, 5)
m03 <- bootstrap(m00, B=5, type="param")
extractBOOT(m03)
summary(m03)

## Model selection

m04 <- svabu(Y ~ x1 + x5 | x2 + x5 + x3, databu[1:200,], phi.boot=0)
m05 <- drop1(m04, model="det")
m05
m06 <- svabu.step(m04, model="det")
summary(m06)
m07 <- update(m04, . ~ . | . - x3)
m07

## Controls

m00$control
getOption("detect.optim.control")
getOption("detect.optim.method")
options("detect.optim.method"="BFGS")
m08 <- svabu(Y ~ x1 + x5 | x2 + x5, databu[1:100,])
m08$control ## but original optim method is retained during model selection and bootstrap
## fitted models can be used to provide initial values
options("detect.optim.method"="Nelder-Mead")
m09 <- svabu(Y ~ x1 + x5 | x2 + x5, databu[1:100,], inits=coef(m08))

## Ovenbirds dataset

data(oven)
ovenc <- oven
ovenc[, c(4:8,10:11)][] <- lapply(ovenc[, c(4:8,10:11)], scale)
moven <- svabu(count ~ pforest | observ + pforest + julian + timeday, ovenc)
summary(moven)
drop1(moven, model="det")
moven2 <- update(moven, . ~ . | . - timeday)
summary(moven2)
moven3 <- update(moven2, . ~ . | ., zeroinfl=FALSE)
summary(moven3)
BIC(moven, moven2, moven3)
}
}
\keyword{ htest }
\keyword{ models }
