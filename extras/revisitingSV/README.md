Revisiting resource selection probability functions and single-visit methods: Clarification and extensions
=====

## Introduction

Supporting Information for CITATION GOES HERE.

Here we give rationale and code for the simulation used in the paper.
The organization of this archive is as follows:

* [Examples for the RSPF condition](#examples-for-the-rspf-condition)
* [Quasi-Bayesian single-visit occupancy model](#quasi-bayesian-single-visit-occupancy-model)
* [B-B-ZIP CL identifiability](#b-b-zip-cl-identifiability)

## Examples for the RSPF condition

The following figure from the paper demosntrates that the RSPF condition
defines a class of functions that are fairly flexible.
Here we are generating regression coefficients (`cf`) randomly but users can supply their own. 

```R
m <- 5 # how many lines to draw
x <- seq(-1, 1, by=0.01) # predictor
col <- 1:m

## Regression coefficients. 
set.seed(12345)
cf <- sapply(1:m, function(i) 
    c(rnorm(1, 0, 1), rnorm(1, 0, 1+i/2),
    rnorm(1, 0, 2), rnorm(1, 0, 4)))

## plot the response curves: 
## Linear or Quadratic or Cubic terms included 
## in the link (logit or cloglog).

op <- par(mfrow=c(3,2))

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="logit, linear")
for (i in 1:m)
    lines(x, plogis(cf[1,i]+cf[2,i]*x), col=col[i])

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="cloglog, linear")
for (i in 1:m)
    lines(x, binomial("cloglog")$linkinv(cf[1,i]+cf[2,i]*x), col=col[i])

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="logit, quadratic")
for (i in 1:m)
    lines(x, plogis(cf[1,i]+cf[2,i]*x+cf[3,i]*x^2), col=col[i])

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="cloglog, quadratic")
for (i in 1:m)
    lines(x, binomial("cloglog")$linkinv(cf[1,i]+cf[2,i]*x+cf[3,i]*x^2), col=col[i])

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="logit, cubic")
for (i in 1:m)
    lines(x, plogis(cf[1,i]+cf[2,i]*x+cf[3,i]*x^2+cf[4,i]*x^3), col=col[i])

plot(0, type="n", ylim=c(0,1), xlim=c(-1,1),
    xlab="x", ylab="p", main="cloglog, cubic")
for (i in 1:m)
    lines(x, binomial("cloglog")$linkinv(cf[1,i]+cf[2,i]*x+cf[3,i]*x^2+cf[4,i]*x^3), col=col[i])

par(op)
```

In the following, we illustrate the estimation of the abundance 
and probability of detection when the logit link has cubic terms 
included in the model. The probability of detection 
does not reach 1 for any covariate value. 
This illustrates that the critical condition is that the 
model satisfy the RSPF condition and not that probability 
of detection is 1 for some combination of the covariates. 

We simulate data under BP N-mixture with logit 
link for detection with cubic terms. 

```R
n <- 1000
k <- 3 # order of polynomial
link <- "logit"
beta <- c(1,1)
theta <- cf[1:(1+k),3]

x1 <- rnorm(n, 0, 1)
X <- model.matrix(~ x1)
lambda <- exp(X %*% beta)

z1 <- sort(runif(n, -1, 0.75))
z2 <- z1^2
z3 <- z1^3
Z <- model.matrix(~ z1 + z2 + z3)[,1:(1+k)]
p <- binomial(link)$linkinv(Z %*% theta)
summary(p)    # This shows that the probability of detection does not reach 1 for any covariate value. 
N <- rpois(n, lambda)
Y <- rbinom(n, N, p)
table(Y)
```

Now we estimate the coefficients for BP N-mixture. Things to observe:

* AIC for the true model is the lowest, as it should be.
* The best fitting model approximates the true probabilities of detection quite well. 

```R
library(detect)
m1 <- svabu(Y ~ x1 | z1, zeroinfl=FALSE, link.det=link)
m2 <- svabu(Y ~ x1 | z1 + z2, zeroinfl=FALSE, link.det=link)
m3 <- svabu(Y ~ x1 | z1 + z2 + z3, zeroinfl=FALSE, link.det=link)
maic <- AIC(m1, m2, m3)
maic$dAIC <- maic$AIC - min(maic$AIC)
## AIC for the true model is the lowest, as it should be. 
maic

## How do the fitted probabilities of detection compare with the true ones?
plot(z1, p, type="l", ylim=c(0,1))
lines(z1, m1$detection.probabilities, col=2)
lines(z1, m2$detection.probabilities, col=3)
lines(z1, m3$detection.probabilities, col=4)
legend("bottomright", lty=1, col=1:4,
    legend=paste(c("truth", "linear", "quadratic", "cubic"),
        "AIC =", c("NA",round(maic$AIC, 2))))
```

Now let us see if we can estimate the absolute probability of selection if the 
RSPF condition is satisfied. We generate the use-available 
data under a cubic model (same as the probability of 
detection model above). User can use any other model 
that satisfies the RSPF condition. 

We simulate data under RSPF model.
The `simulateUsedAvail` function generates the used points  under probability sampling with replacement. 
The data frame `dat` also includes a sample of available points that is obtained using simple random sampling with replacement. (See Lele and Keim 2006: Animal as a sampler description in the discussion).

```R
library(ResourceSelection)

n <- 3000
k <- 3 # order of polynomial
link <- "logit"
beta <- c(1,1)
theta <- cf[1:(1+k),3]

x1 <- rnorm(n, 0, 1)
X <- model.matrix(~ x1)
lambda <- exp(X %*% beta)

z1 <- sort(runif(n, -1, 0.75))
z2 <- z1^2
z3 <- z1^3

Z <- model.matrix(~ z1 + z2 + z3)[,1:(1+k)]
p <- binomial(link)$linkinv(Z %*% theta)
summary(p)     # Notice that the probability of selection never reaches 1. 

z <- data.frame(z1=z1, z2=z2, z3=z3)[,1:k]   # This is the distribution of the covariates: The available distribution. 

dat <- simulateUsedAvail(z, theta, n, m=10, link=link)  
```

Estimate coefficients for RSPF. The best fitting model approximates the true probabilities quite well. 
As usual, larger the sample size, better will be the fit.

```R
r1 <- rspf(status ~ z1, dat, m=0, B=0, link=link)
r2 <- rspf(status ~ z1 + z2, dat, m=0, B=0, link=link)
r3 <- rspf(status ~ z1 + z2 + z3, dat, m=0, B=0, link=link)
raic <- AIC(r1, r2, r3)
raic$dAIC <- raic$AIC - min(raic$AIC)
raic   # The best fitting model is the cubic model as it should be. 

## How do the estimated probability of selection compare 
## with the true probabilities of selection? 
plot(z1, p, type="l", ylim=c(0,1))
points(dat$z1, fitted(r1), col=2, pch=".")
points(dat$z1, fitted(r2), col=3, pch=".")
points(dat$z1, fitted(r3), col=4, pch=".")
legend("bottomright", lty=1, col=1:4,
    legend=paste(c("truth", "linear", "quadratic", "cubic"),
        "AIC =", c("NA",round(raic$AIC, 2))))
```

## Quasi-Bayesian single-visit occupancy model

### Generate data 


```R
library(dclone)
library(boot)


OccData.fun = function(psi,p1){
	
	N = length(p1)
	Y = rbinom(N,1,psi)
	W = rbinom(N,1,Y*p1)
	out = list(Y=Y,W=W)
return(out)
}
```

### Penalized estimation

This is basically a Bayesian approach but prior is based on the 
currently observed data. The prior on beta is centered at the 
naive estimator. This forces the estimator to not veer off 
too far from the naive estimator. This penalty function is 
somewhat easier and simpler than the convoluted penalty 
function used in Lele et al. (2012). We write it using the 
conditional distribution form. These estimators are stable 
but biased for small samples (as are all other estimators used here).

`beta.naive` are from the naive estimator. 
We want to keep it close to this estimator unless the 
information in the data forces us to move. 

```R
OccPenal.est = function(){

# Data list includes: W, beta.naive, penalty,X1,Z1,nS
	
	for (i in 1:nS){
		W[i] ~ dbin(p1[i]*Y[i],1)
		p1[i] <- exp(inprod(Z1[i,],theta))/(1+exp(inprod(Z1[i,],theta)))	
		Y[i] ~ dbin(psi[i],1)
		psi[i] <- exp(inprod(X1[i,],beta))/(1+exp(inprod(X1[i,],beta)))
	}
for (i in 1:c1){
		beta[i] ~ dnorm(beta.naive[i],penalty)     
	}	
	for (i in 1:c2){
		theta[i] ~ dnorm(0,0.01)
	}	
}
```

### Data generation and analysis

There are two ways to increase the information in the data under the 
regression setting: 

* by increasing the sample size, 
* and by increasing the variability on the covariate space. 

The other is much more cost effective and is under the control of the researcher.

```R
N = 1000
beta = c(1,2)   # Occupancy parameters
theta = c(-1,0.6)  # Detection parameters

X = runif(N,-2,2)  # Occupancy covariate
Z = runif(N,-2,2)  # Detection covariate

X1 = model.matrix(~X)
Z1 = model.matrix(~Z)

psi = inv.logit(X1 %*% beta)
p1 = inv.logit(Z1 %*% theta)

mean(psi)   # Check the average occupancy probability
mean(p1)    # Check the average detection probability

occ.data = OccData.fun(psi,p1)

occ.data = list (Y = occ.data$Y, W = occ.data$W, X1 = X1, Z1= Z1)
```

### Using only the regularized likelihood function


```R
S = 10      # Number of simulations
out1 = matrix(0,S,5)   # Parameter estimates
out2 = matrix(0,S,3)

for (s in 1:S){

nS = 400    # Sample size from the population
sample.index = sample(seq(1:N),nS,replace=F)

beta.naive = glm(occ.data$W[sample.index] ~ X1[sample.index,2], family="binomial")
penalty = 0.5    # This is the precision for the beta parameters. This should be at least as small as the 1/SE for beta.naive but it could be substantially smaller than that (but not too small). 

dat = list(W = occ.data$W[sample.index], X1 = occ.data$X1[sample.index,], Z1 = occ.data$Z1[sample.index,], c1 = ncol(occ.data$X1), c2 = ncol(occ.data$Z1),nS=nS, beta.naive = beta.naive$coef, penalty=penalty)
inits = list(Y = rep(1,nS))

OccPenalty.fit = jags.fit(dat,c("beta","theta"),OccPenal.est, inits=inits, n.adapt=10000, n.update=10000, n.iter=1000)

diag = gelman.diag(OccPenalty.fit)  # Check for convergence.  

parms.est = summary(OccPenalty.fit)$quantiles[,c(1,3,5)]  # Median value and the 95% credible interval. 

OccPenalty.fit = jags.fit(dat,c("p1","psi"),OccPenal.est, inits=inits, n.adapt=10000, n.update=10000, n.iter=1000)

tmp = summary(OccPenalty.fit)
naive.medianOcc = summary(beta.naive$fitted)[3]
true.medianOcc = summary(psi[sample.index])[3]
Adj.medianOcc = summary(tmp$quantiles[(nS+1):(2*nS),3])[3]

out1[s,] = c(parms.est[,2],diag$mpsrf)   # Only spits out the median for simulation study.
out2[s,] = c(naive.medianOcc,true.medianOcc,Adj.medianOcc)

}

out1
out2
```

## B-B-ZIP CL identifiability

The goal is to establish identifiability for svabuRD models:
Binomial-Binomial-Poisson for p*p_max issue.

We simulate under the single visit distance sampling model 
to prove identifiability (constant and variable EDR): 
to demonstrate that the concern of 0.5*p can be addressed 
if there is a mechanism that causes that.

### B-B-ZIP simulation

```R
R <- 100
n <- 2000
K <- 250

setwd(set_wd)
source("svabuRD.R")

set.seed(4321)
link <- "logit"
linkinvfun.det <- binomial(link=make.link(link))$linkinv

x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)
x3 <- runif(n,-1,1)

X <- model.matrix(~ x1)
ZR <- model.matrix(~ x3)
ZD <- model.matrix(~ x2)

beta <- c(-1,1)
thetaR <- c(0.9, -1.5) # for singing rate
edr <- 0.65

## point count radius
r <- sample(c(0.5,1), n, replace=TRUE)

q <- (edr / r)^2 * (1 - exp(-(r / edr)^2))

D <- exp(drop(X %*% beta))
lambda <- D * pi*r^2

p <- linkinvfun.det(drop(ZR %*% thetaR))

summary(D)
summary(p)
summary(q)

res2 <- list()
for (i in 1:R) {
cat(i, "\n");flush.console()

Y <- rpois(n, lambda * p * q)

m <- svabuRD.fit(Y, X, ZR, NULL, Q=NULL, zeroinfl=FALSE, r=r, N.max=K,
    inits=c(beta, thetaR, log(edr))+rnorm(5, 0, 0.2))

res2[[i]] <- coef(m)

}
save.image("out-sim-2.Rdata")
```

### Results from B-B-ZIP simulation

```R
load("out-sim-2.Rdata")

True2 <- c(beta, thetaR, "logedr"=log(edr))
cf2 <- t(sapply(res2, unlist))
bias2 <- t(sapply(res2, unlist) - True2)

lamTrue <- mean(exp(X %*% beta))
pTrue <- mean(plogis(drop(X %*% thetaR)))
qTrue <- mean((edr / r)^2 * (1 - exp(-(r / edr)^2)))

lam_hat <- apply(cf2[,1:2], 1, function(z) {
    mean(exp(X %*% z))
})
p_hat <- apply(cf2[,3:4], 1, function(z) {
    mean(plogis(drop(X %*% z)))
})
q_hat <- sapply(cf2[,5], function(z) {
    mean((exp(z) / r)^2 * (1 - exp(-(r / exp(z))^2)))
})

rr <- 0:150 / 100
qq <- sapply(cf2[,5], function(z) 
    (exp(z) / rr)^2 * (1 - exp(-(rr / exp(z))^2)))
qq_true <- (edr / rr)^2 * (1 - exp(-(rr / edr)^2))

est <- cbind(lam=lam_hat, p=p_hat, q=q_hat)
bias <- cbind(lam=(lam_hat-lamTrue)/lamTrue, 
    p=(p_hat-pTrue)/pTrue, q=(q_hat-qTrue)/qTrue)

## figure
pdf("svabuRD.pdf", width=12, height=4)
par(mfrow=c(1,3), mar=c(5, 5, 4, 2) + 0.1)
boxplot(bias2,
    col="grey",
    ylab=expression(hat(theta)-theta),
    names=c(expression(alpha[0]),
      expression(alpha[1]), 
      expression(beta[0]),
      expression(beta[1]),
      expression(log(tau))))
abline(h=0, lty=2, col=1)
matplot(rr*100, qq, lty=1, col="grey", type="l",
    ylab="q", xlab="Distance from observer (m)")
lines(rr*100, qq_true, col=1, lwd=2)
boxplot(bias,
    col="grey",
    ylab="Relative bias",
    names=c(expression(bar(lambda)),
      expression(bar(p)), 
      expression(bar(q))))
abline(h=0, lty=2, col=1)
dev.off()
```

