Revisiting the single-visit method for occupancy and abundance studies: A rebuttal of criticisms and some extensions
=====

## Introduction

Supporting Information for CITATION GOES HERE.

Here we give rationale and code for the simulation used in the rebuttal paper.
The organization of this archve is as follows:

| Simulation | Results |
|--------------|--------------|
| [Simulation 1](#simulation-1) | [log1 identifiability](#log1-identifiability) |
| [Simulation 1](#simulation-1) | [Model selection](#model-selection) |
| [Simulation 1](#simulation-1) | [Link function bias for log1](#link-function-bias-for-log1) |
| [Simulation 2](#simulation-2) | [B-B-ZIP CL identifiability](#b-b-zip-cl-identifiability) |
| [Simulation 3](#simulation-3) | [Half-logit bias](#half-logit-bias) |

## Simulation 1

The goal is to prove that N-mixture selects the correct model based on AIC 
and that the constrained log-link (log1) is identifiable under single-visit.

```
## number of locations
n <- 500

## number of runs for the simulations
R <- 100

## set working directory here
setwd("/set/your/working/directory/here/")
source("svabu_link.R")

## loading detect package
library(detect)

## ensuring reproducibility
set.seed(1234)

## predictors
x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)
x <- data.frame(x1=x1, x2=x2)

## model matrices for abundance (X) and detection (Z)
X <- model.matrix(~ x1)
Z <- model.matrix(~ x2)

## this is the constrained exponential function
## used as inverse log link
log1 <- function(x) pmin(binomial("log")$linkinv(x), 1-.Machine$double.eps)
loglog <- function(x) exp(-exp(-x))

## true coefficients
beta <- c(2,-0.8)
theta <- c(-1, 2)

## N.max for integral
K <- 100

## loop for link functions to generate the data starts here
Links <- c("logit","probit","cloglog", "loglog", "log1")

estfun <- function(i, link) {

    ## define link function for simulating the data
    if (link %in% c("logit","probit","cloglog"))
        linkinv <- binomial(link)$linkinv else {
        if (link == "log1")
            linkinv <- log1
        if (link == "loglog")
            linkinv <- loglog
    }

    ## expected values are same for all R runs
    ## inverse link function for p is defined by 'link'
    p <- linkinv(drop(Z %*% theta))
    ## inverse link function for abundance is kept the same
    lambda <- poisson("log")$linkinv(drop(X %*% beta))

    ## simulate the observations
    N <- rpois(n, lambda)
    y <- rbinom(n, N, p)

    cat("\n", link, i);flush.console()

    ## fit logit link to simulated data
    cat(" logit");flush.console()
    m_logit <- svabu(y ~ x1 | x2, x, N.max=K,
        zeroinfl=FALSE, link.det="logit")

    ## fit probit link to simulated data
    cat(" probit");flush.console()
    m_probit <- svabu(y ~ x1 | x2, x, N.max=K,
        zeroinfl=FALSE, link.det="probit")

    ## fit cloglog link to simulated data
    cat(" cloglog");flush.console()
    m_cloglog <- svabu(y ~ x1 | x2, x, N.max=K,
        zeroinfl=FALSE, link.det="cloglog")

    ## fit loglog link to simulated data
    cat(" loglog");flush.console()
    m_loglog <- svabu_link(y ~ x1 | x2, x, N.max=K,
        zeroinfl=FALSE, link.det=loglog)

    ## fit log1 link to simulated data
    cat(" log1");flush.console()
    m_log1 <- svabu_link(y ~ x1 | x2, x, N.max=K,
        zeroinfl=FALSE, link.det=log1)

    ## summarize AIC and estimated coefficients
    aic <- AIC(m_logit, m_probit, m_cloglog, m_loglog, m_log1)
    daic <- aic$AIC - min(aic$AIC)
    cf <- sapply(list(m_logit, m_probit, m_cloglog, m_loglog, m_log1), coef)
    names(daic) <- Links
    colnames(cf) <- Links
    list(dAIC=daic, coef=cf, N=N, y=y)
} # loop for R ends here

library(snow)
cl <- makeCluster(10)

clusterSetupRNG (cl, type = "RNGstream")
#tmp <- clusterEvalQ(cl, setwd("c:/Dropbox/pkg/detect2/mee-rebuttal/gh/"))
tmp <- clusterEvalQ(cl, setwd("c:/Users/Peter/Dropbox/pkg/detect2/mee-rebuttal/gh/"))
tmp <- clusterEvalQ(cl, source("svabu_link.R"))
tmp <- clusterEvalQ(cl, library(detect))
clusterExport(cl, c("n","x1","x2","x","X","Z","log1","loglog","beta","theta","K","Links"))

r_1 <- parLapply(cl, 1:R, estfun, link="logit")
r_2 <- parLapply(cl, 1:R, estfun, link="probit")
r_3 <- parLapply(cl, 1:R, estfun, link="cloglog")
r_4 <- parLapply(cl, 1:R, estfun, link="loglog")
r_5 <- parLapply(cl, 1:R, estfun, link="log1")

stopCluster(cl)

save.image("out-sim-1.Rdata")
```

## Results from Simulation 1

```
load("out-sim-1.Rdata")

AIC1 <- list()
AIC1$logit <- lapply(r_1, "[[", "dAIC")
AIC1$probit  <- lapply(r_2, "[[", "dAIC")
AIC1$cloglog <- lapply(r_3, "[[", "dAIC")
AIC1$loglog <- lapply(r_4, "[[", "dAIC")
AIC1$log1 <- lapply(r_5, "[[", "dAIC")

CF1 <- list()
CF1$logit <- lapply(r_1, "[[", "coef")
CF1$probit  <- lapply(r_2, "[[", "coef")
CF1$cloglog <- lapply(r_3, "[[", "coef")
CF1$loglog <- lapply(r_4, "[[", "coef")
CF1$log1 <- lapply(r_5, "[[", "coef")

trueLink <- names(AIC1)
nLink <- length(trueLink)

AIC0 <- matrix(NA, nLink, nLink-0)
rownames(AIC0) <- trueLink
colnames(AIC0) <- paste("fit", trueLink[1:(nLink-0)])
LAM0 <- P0 <- AIC0

lamTrue <- mean(exp(X %*% beta))
ip <- drop(X %*% theta)
pTrue <- numeric(4)
pTrue[1] <- mean(binomial("logit")$linkinv(ip))
pTrue[2] <- mean(binomial("probit")$linkinv(ip))
pTrue[3] <- mean(binomial("cloglog")$linkinv(ip))
pTrue[4] <- mean(loglog(ip))
pTrue[5] <- mean(log1(ip))
names(pTrue) <- trueLink
pTrue

## summary stat can be changed here
f <- mean

cfun <- function(z) {
    lam <- colMeans(exp(apply(z[1:2,], 2, function(x) X %*% x)))
    m_ip <- apply(z[3:4,], 2, function(x) Z %*% x)
    m_ip[,1] <- binomial("logit")$linkinv(m_ip[,1])
    m_ip[,2] <- binomial("probit")$linkinv(m_ip[,2])
    m_ip[,3] <- binomial("cloglog")$linkinv(m_ip[,3])
    m_ip[,4] <- loglog(m_ip[,4])
    p <- colMeans(m_ip)
    rbind(lam, p)
}

BEST <- list()
for (Link in trueLink) {
    dAIC1 <- do.call(rbind, AIC1[[Link]])
    tmp1 <- lapply(CF1[[Link]], cfun)
    
    best <- sapply(1:R, function(i) tmp1[[i]][,dAIC1[i,]==0])
    lam1 <- apply(sapply(tmp1, function(z) z["lam",]) - lamTrue, 1, f)
    p1 <- apply(sapply(tmp1, function(z) z["p",]) - pTrue[Link], 1, f)

    AIC0[Link,] <- colMeans(dAIC1 == 0)
    LAM0[Link,] <- lam1
    P0[Link,] <- p1
    BEST[[Link]] <- best - c(lamTrue, pTrue[Link])
}

## bias needs to use the best fit per run
for (Link in trueLink) {
    dAIC1 <- do.call(rbind, AIC1[[Link]])
    best <- apply(dAIC1, 1, which.min)
}


```

### log1 identifiability

```
log1bias <- t(sapply(CF1$log1, function(z) z[,"log1"])-c(beta,theta))
log1lam <- (sapply(CF1$log1, function(z) mean(exp(drop(X %*% z[1:2,"log1"])))) - lamTrue)/lamTrue
log1p <- (sapply(CF1$log1, function(z) mean(log1(drop(X %*% z[3:4,"log1"])))) - pTrue["log1"])/pTrue["log1"]

pdf("log1bias.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(5, 5, 4, 2) + 0.1)
boxplot(log1bias,
    col="grey",
    ylab=expression(hat(theta)-theta),
    names=c(expression(alpha[0]),
      expression(alpha[1]), 
      expression(beta[0]),
      expression(beta[1])))
abline(h=0, lty=2, col=1)
boxplot(cbind(log1lam, log1p),
    col="grey",
    ylab="Relative bias",
    names=c(expression(bar(lambda)),
      expression(bar(p))))
abline(h=0, lty=2, col=1)
dev.off()
```

### Model selection

```
tab <- data.frame("True model" = rep(trueLink, each=nLink-0),
    "Fitted model" = rep(trueLink[1:(nLink-0)], nLink),
    "Freq"=as.numeric(t(AIC0)),
    "lam bias"=as.numeric(t(LAM0 / lamTrue)),
    "p bias"=as.numeric(t(P0 / pTrue)))
write.csv(AIC0, file="model-selection-SVaic.csv")
```

### Link function bias for log1

```
cfun2 <- function(i, j=NULL) {
    if (is.null(j))
        j <- which(AIC1[["log1"]][[i]]==0)
    cf <- CF1[["log1"]][[i]][3:4,j]
    m_ip <- drop(Z %*% cf)
    if (j == 1)
        m_p <- binomial("logit")$linkinv(m_ip)
    if (j == 2)
        m_p <- binomial("probit")$linkinv(m_ip)
    if (j == 3)
        m_p <- binomial("cloglog")$linkinv(m_ip)
    if (j == 4)
        m_p <- loglog(m_ip)
    if (j == 5)
        m_p <- log1(m_ip)
    m_p
}

mu <- drop(Z %*% theta)
ii <- order(mu)
p_log1 <- log1(mu)
best <- sapply(1:R, cfun2)
best1 <- sapply(1:R, cfun2, j=1)
best2 <- sapply(1:R, cfun2, j=2)
best3 <- sapply(1:R, cfun2, j=3)
best4 <- sapply(1:R, cfun2, j=4)
best5 <- sapply(1:R, cfun2, j=5)

cfun3 <- function(i, j=NULL) {
    if (is.null(j))
        j <- which(AIC1[["log1"]][[i]]==0)
    cf <- CF1[["log1"]][[i]][1:2,j]
    mean(exp(drop(X %*% cf)))
}
xbest <- (sapply(1:R, cfun3) - lamTrue)/lamTrue
xbest1 <- (sapply(1:R, cfun3, j=1) - lamTrue)/lamTrue
xbest2 <- (sapply(1:R, cfun3, j=2) - lamTrue)/lamTrue
xbest3 <- (sapply(1:R, cfun3, j=3) - lamTrue)/lamTrue
xbest4 <- (sapply(1:R, cfun3, j=4) - lamTrue)/lamTrue
xbest5 <- (sapply(1:R, cfun3, j=5) - lamTrue)/lamTrue

xtext <- -5
ytext <- 0.9

pdf("log1-fit.pdf", width=9, height=6)
par(mfrow=c(2,3), mar=c(5, 5, 4, 2) + 0.1)
matplot(mu[ii], cbind(p_log1[ii], best1[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="logit(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest1), 2), "%"))
matplot(mu[ii], cbind(p_log1[ii], best2[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="probit(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest2), 2), "%"))
matplot(mu[ii], cbind(p_log1[ii], best3[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="cloglog(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest3), 2), "%"))
matplot(mu[ii], cbind(p_log1[ii], best4[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="loglog(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest4), 2), "%"))
matplot(mu[ii], cbind(p_log1[ii], best5[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="log1(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest5), 2), "%"))
## this is needed only if cloglog is not 100%
matplot(mu[ii], cbind(p_log1[ii], best[ii,]), lty=1, col="grey", type="l",
    ylab="p", xlab="best fit h(p)")
lines(mu[ii], p_log1[ii], lwd=2)
text(xtext, ytext, paste0("Rel. Bias = ", round(mean(100*xbest), 2), "%"))
dev.off()
```

## Simulation 2

Here text

## Results from Simulation 2

### B-B-ZIP CL identifiability

Here text

## Simulation 3

Here text

## Results from Simulation 3

### Half-logit bias

Here text


