source("~/repos/detect/R/zi.fit.R")

## Lognormal

n <- 2000
phi <- 0.7
beta <- c(1,0.5)
sigma <- 0.5
x <- rnorm(n)
X <- model.matrix(~x)
mu <- drop(X %*% beta)
o <- rnorm(n, mu, sigma)
z <- rbinom(n, 1, 1-phi)
Z <- matrix(1, n, 1)
Y <- z * exp(o)
hist(Y)

m2 <- zi.fit(Y, X, Z, distr="lognorm", type="ML")
m2$ML$coef
c(beta, log(sigma), qlogis(phi))
lapply(zi.fit(Y, X, Z, distr="lognorm", type=c("ML","CL","PL")), "[[", "coef")



## Poisson

n <- 2000
phi <- 0.7
beta <- c(1,0.5)
x <- rnorm(n)
X <- model.matrix(~x)
mu <- drop(X %*% beta)
z <- rbinom(n, 1, 1-phi)
Z <- matrix(1, n, 1)
Y <- rpois(n, exp(mu) * z)
plot(table(Y))

library(pscl)
m1 <- zeroinfl(Y ~ x | 1, dist="poisson")
m2 <- zi.fit(Y, X, Z, distr="pois", type="ML")
coef(m1)
m2$ML$coef
c(beta, qlogis(phi))

lapply(zi.fit(Y, X, Z, distr="pois", type=c("ML","CL","PL")), "[[", "coef")

## testing fit
clpl <- zi.fit(Y, X, Z, distr="pois", type=c("CL","PL"))
cl <- zi.fit(Y, X, Z, distr="pois", type="CL")
pl <- zi.fit(Y, X, Z, distr="pois", type="PL", fit=cl)
pl0 <- zi.fit(Y, X, Z, distr="pois", type="PL", fit=-exp(mu))
pl$PL$coef
pl0$PL$coef
plogis(pl0$PL$coef)
1-plogis(-pl0$PL$coef)

## NegBin

n <- 2000
phi <- 0.7
beta <- c(1,0.5)
x <- rnorm(n)
X <- model.matrix(~x)
mu <- drop(X %*% beta)
z <- rbinom(n, 1, 1-phi)
Z <- matrix(1, n, 1)
theta <- 0.1

lambda <- drop(exp(X %*% beta))
A <- rbinom(n, 1, 1-phi)
#Y <- stats:::rnbinom(n, size=1/gvar, prob=1/(1+gvar*lambda*A))
Y <- stats:::rnbinom(n, size=theta, mu=lambda*A)
table(Y)

library(pscl)
m1 <- zeroinfl(Y ~ x | 1, dist="negbin")
m2 <- zi.fit(Y, X, Z, distr="negbin", type="ML")
coef(m1)
log(m1$theta)
m2$ML$coef
c(beta, log(theta), qlogis(phi))

lapply(zi.fit(Y, X, Z, distr="negbin", type=c("ML","CL","PL")), "[[", "coef")

## Beta

n <- 2000
phi <- 0.7
beta <- c(1,0.5)
x <- rnorm(n)
X <- model.matrix(~x)
mu <- plogis(drop(X %*% beta))
z <- rbinom(n, 1, 1-phi)
Z <- matrix(1, n, 1)
gamma <- 2
a <- mu * gamma
b <- (1 - mu) * gamma
rb <- rbeta(n, a, b)
Y <- ifelse(z == 0, z, rb)
hist(Y)

m2 <- zi.fit(Y, X, Z, distr="beta", type="ML")
m2$ML$coef
c(beta, log(gamma), qlogis(phi))
lapply(zi.fit(Y, X, Z, distr="beta", type=c("ML","CL","PL")), "[[", "coef")

## Binomial

n <- 2000
phi <- 0.7
beta <- c(1,0.5)
x <- rnorm(n)
X <- model.matrix(~x)
mu <- plogis(drop(X %*% beta))
z <- rbinom(n, 1, 1-phi)
Z <- matrix(1, n, 1)
N <- rep(4, n)
W <- rbinom(n, N, mu)
Y <- ifelse(z == 0, z, W)
hist(Y)

m2 <- zi.fit(Y, X, Z, distr="binom", type="ML", N=N)
m2$ML$coef
c(beta, qlogis(phi))
lapply(zi.fit(Y, X, Z, distr="binom", type=c("ML","CL","PL"), N=N), "[[", "coef")

