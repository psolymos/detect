\name{sqpad}
\alias{sqpad}
\alias{sqpad.fit}
\alias{fitted.sqpad}
\alias{predict.sqpad}
\alias{vcov.sqpad}
\alias{logLik.sqpad}
\alias{print.sqpad}
\alias{summary.sqpad}
\alias{print.summary.sqpad}
\title{
Single-bin QPAD (SQPAD)
}
\description{
Single bin QPAD (SQPAD) approach by Lele and Solymos (2025).
}
\usage{
sqpad(formula, data, dis, dur,
    type = c("full", "approx"), det = c("joint", "pq"),
    Nmax = NULL, K = NULL, A = NULL,
    montecarlo = FALSE, np = 1000, distcorr = 2/3, ...)

sqpad.fit(Y, dis, dur, X = NULL, Z = NULL, 
    type = c("full", "approx", "conv"), det = c("joint", "pq"),
    init = NULL, method = "Nelder-Mead", hessian = TRUE,
    tt1 = NULL, A = NULL, Nmax = NULL, K = NULL, 
    montecarlo = FALSE, np = 1000, distcorr = 2/3, dislist = NULL, ...)

\method{print}{sqpad}(x, digits, ...)
\method{vcov}{sqpad}(object, ...)
\method{fitted}{sqpad}(object, ...)
\method{logLik}{sqpad}(object, ...)
\method{summary}{sqpad}(object, ...)
\method{print}{summary.sqpad}(x, digits, ...)
\method{predict}{sqpad}(object, newdata = NULL, 
    type = c("link", "response"), ...)
}
\arguments{
  \item{formula}{
formula, LHS takes a vector, RHS is either \code{1} or some covariates, 
see Examples.
}
  \item{data}{
data.
}
  \item{dis, dur}{
a vector with distance and duration values.
}
  \item{A}{
a vector with area values. Area is calculated from \code{dis}, use \code{A}
if distance is not directly related to area (e.g. for ARU data).
}
  \item{type}{
character, one of \code{"full"} (full likelihood), \code{"approx"} (approximate 
method), and \code{"conv"} (convolution likelihood).
}
  \item{det}{
character, the detection model: \code{"joint"} or independence (\code{"pq"}).
}
  \item{Y}{
response vector, non-negative integers.
}
  \item{X, Z}{
design matrices, \code{X} is the matrix with covariates for the
density parameters. \code{Z} is the matrix with covariates for the cue rate model.
\code{NULL} means constant density and detection.
}
  \item{K}{
truncation value, \code{K=1} implements the occupancy version of SQPAD.
\code{K>1} implements the truncated count models.
}
  \item{Nmax}{
maximum abundance value for numerical integration under the full 
and convolution likelihoods.
}
  \item{montecarlo, np, distcorr}{
should mean probability (\code{TRUE}) using \code{np} points or 
approximation be used with distance correction \code{distcorr} 
in the detection model.
}
  \item{dislist}{
distance list for the convolution likelihood.
}
  \item{tt1}{
vector with time to 1st detection values, same units as for \code{dur}.
}
  \item{init}{
optional initial values for \code{\link{optim}}.
}
  \item{method}{
method for \code{\link{optim}}.
}
  \item{hessian}{
logical, should the Hessian matrix be computed by \code{\link{optim}}.
}
  \item{object, x}{
fitted model object.
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to predict.
If omitted, the fitted linear predictors are used.
}
  \item{digits}{
digits to use when providing summaries.
}
  \item{\dots}{
additional options for \code{\link{optim}} or the methods.
}
}
\details{
Single bin QPAD (SQPAD) approach for robust analysis of point count data 
with detection error by Lele and Solymos (2025).
}
\value{
An object of class 'sqpad'.
}
\references{
Solymos, P., Lele, S. R., 2025. Single bin QPAD (SQPAD) approach for robust 
analysis of point count data with detection error.
\emph{Ornithological Applications}, \bold{xx}, xx--xx.

Supporting info for the SQPAD method:
\url{https://github.com/psolymos/sqpad-paper},
<doi:10.5281/zenodo.16172209>.
}
\author{
Peter Solymos and Subhash Lele
}
\examples{

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
summary(m)

coef(m)
nobs(m)
vcov(m)
confint(m)

logLik(m)
AIC(m)
BIC(m)

fitted(m)
predict(m, type = "link")
predict(m, newdata = df[1:10,], type = "response")
predict(m, newdata = df[1:10,], type = "link")

\dontrun{
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
dij <- extraDistr::rtriang(as.integer(np), 0, 1, 1)
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
    dij <- extraDistr::rtriang(N[i], 0, dis[i], dis[i])
    # pij is based on time of 1st detection and distance at 1st detection
    pij <- 1 - exp(-dur[i] * phi * exp(-(dij)^2/tau^2))
    Yij <- rbinom(length(pij), 1, pij)
    dij[Yij > 0]
})
Y <- sapply(dis, length)

m4 <- sqpad.fit(Y = Y, dis = dis, dur = dur, dislist = dislist, 
    type = "conv", det = "joint", Nmax = 20)
cbind(parameters = c(D = D, phi = phi, tau = tau), estimates = exp(m4$coef))
}
}
\keyword{ htest }
