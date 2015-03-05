Revisiting resource selection probability functions and single-visit methods: Clarification and extensions
=====

## Introduction

Supporting Information for CITATION GOES HERE.

Here we give rationale and code for the simulation used in the paper.
The organization of this archive is as follows:

* [Examples for the RSPF condition](#examples-for-the-rspf-condition)
* [B-B-ZIP CL identifiability](#b-b-zip-cl-identifiability)

## Examples for the RSPF condition


## B-B-ZIP CL identifiability

The goal is to establish identifiability for svabuRD models:
Binomial-Binomial-Poisson for p*p_max issue.

We simulate under the single visit distance sampling model 
to prove identifiability (constant and variable EDR): 
to demonstrate that the concern of 0.5*p can be addressed 
if there is a mechanism that causes that.

### B-B-ZIP simulation

```
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

```
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

