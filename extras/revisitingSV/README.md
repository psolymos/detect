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
## set wd only once: download files from the repo here
set_wd <- "/set/your/working/directory/here/"

## number of locations
n <- 500

## number of runs for the simulations
R <- 100

## set working directory here
setwd(set_wd)
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
tmp <- clusterEvalQ(cl, library(detect))
clusterExport(cl, c("n","x1","x2","x","X","Z","log1","loglog","beta","theta","K","Links",
    "svabu_link", svabu.fit_link"))

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

The goal is to establish identifiability for svabuRD models:
Binomial-Binomial-Poisson for p*p_max issue.

We simulate under the single visit distance sampling model 
to prove identifiability (constant and variable EDR): 
to demonstrate that the concern of 0.5*p can be addressed 
if there is a mechanism that causes that.

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

## Results from Simulation 2

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
```

### B-B-ZIP CL identifiability

```
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

## Simulation 3

We show the bias in multiple visit results when closure is violated, 
contrast this bias to the single visit-bias when using the half logistic link.

```
n <- 200
R <- 100

library(detect)
library(unmarked)
setwd(set_wd)
source("svabu_link.R")

## link functions defined here
p_max <- 0.5
linkinv_p1 <- function(eta) 1 * binomial("logit")$linkinv(eta)
linkinv_phalf <- function(eta) 0.5 * binomial("logit")$linkinv(eta)
linkinv_p <- if (p_max < 1)
    linkinv_phalf else linkinv_p1
linkinv_N <- poisson("log")$linkinv

beta <- c(2,-0.8)
theta <- c(-1, 2)
T <- 3 # from mallard example
K <- 100

set.seed(1234)
x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)

X <- model.matrix(~ x1)
Z <- model.matrix(~ x2)
x <- data.frame(x1, x2)
xstack <- rbind(x, x, x)

p <- linkinv_p(drop(Z %*% theta))
lambda <- linkinv_N(drop(X %*% beta))

res4 <- list()

for (i in 1:R) {

cat(i, " halflogit\n");flush.console()

N <- y0 <- y1 <- matrix(0, n, T)
colnames(N) <- paste0("N_", 1:T)
colnames(y0) <- paste0("y0_", 1:T)
colnames(y1) <- paste0("y1_", 1:T)
for (t in 1:T) {
    N[,t] <- rpois(n, lambda)
    y0[,t] <- rbinom(n, N[,1], p) # closure true
    y1[,t] <- rbinom(n, N[,t], p) # closure false
}
y1[,1] <- y0[,1] # keep things as consistent as possible

y1stack <- as.numeric(y1)

## multiple visit
visitMat <- matrix(as.character(1:T), n, T, byrow=TRUE)
umf0 <- unmarkedFramePCount(y=y0, siteCovs=x, obsCovs=list(visit=visitMat))
umf1 <- unmarkedFramePCount(y=y1, siteCovs=x, obsCovs=list(visit=visitMat))

m_mv0 <- pcount(~ x2 ~ x1, umf0, K=K, mixture="P")
m_mv1 <- pcount(~ x2 ~ x1, umf1, K=K, mixture="P")

## generalized N-mixture
umf1g <- unmarkedFramePCO(y=y1, siteCovs=x, obsCovs=list(visit=visitMat),
    numPrimary=T)
m_gmv1 <- pcountOpen(~ x1, ~ 1, ~ 1, ~ x2, umf1g, K=K, mixture="P",
    dynamics="constant")

## single visit with n obs
m_sv_halflogit <- svabu_link(y0[,1] ~ x1 | x2, x, 
    zeroinfl=FALSE, link.det=linkinv_phalf, N.max=K)
m_sv_logit <- svabu_link(y0[,1] ~ x1 | x2, x, 
    zeroinfl=FALSE, link.det=linkinv_p1, N.max=K)

## single visit with n obs * T obs
m_svstack_halflogit <- svabu_link(y1stack ~ x1 | x2, xstack, 
    zeroinfl=FALSE, link.det=linkinv_phalf, N.max=K)
m_svstack_logit <- svabu_link(y1stack ~ x1 | x2, xstack, 
    zeroinfl=FALSE, link.det=linkinv_p1, N.max=K)

cflam_mv0 <- coef(m_mv0)[grepl("lam\\(", names(coef(m_mv0)))]
cfp_mv0 <- coef(m_mv0)[grepl("p\\(", names(coef(m_mv0)))]
cflam_mv1 <- coef(m_mv1)[grepl("lam\\(", names(coef(m_mv1)))]
cfp_mv1 <- coef(m_mv1)[grepl("p\\(", names(coef(m_mv1)))]

cflam_gmv1 <- coef(m_gmv1)[grepl("lam\\(", names(coef(m_gmv1)))]
cfp_gmv1 <- coef(m_gmv1)[grepl("p\\(", names(coef(m_gmv1)))]

cflam_sv1 <- coef(m_sv_logit, "sta")
cfp_sv1 <- coef(m_sv_logit, "det")
cflam_svhalf <- coef(m_sv_halflogit, "sta")
cfp_svhalf <- coef(m_sv_halflogit, "det")

cflam_sv1stack <- coef(m_svstack_logit, "sta")
cfp_sv1stack <- coef(m_svstack_logit, "det")
cflam_svhalfstack <- coef(m_svstack_halflogit, "sta")
cfp_svhalfstack <- coef(m_svstack_halflogit, "det")

lam_hat <- cbind(cflam_true=beta, cflam_mv0, cflam_mv1, cflam_gmv1,
    cflam_sv1, cflam_svhalf, cflam_sv1stack, cflam_svhalfstack)
p_hat <- cbind(cfp_true=theta, cfp_mv0, cfp_mv1, cfp_gmv1, 
    cfp_sv1, cfp_svhalf, cfp_sv1stack, cfp_svhalfstack)

tmp <- rbind(lam_hat,p_hat)
attr(tmp, "D-M") <- coef(m_gmv1)[c("gamConst(Int)", "omega(Int)")]
res4[[i]] <- tmp

}
save.image("out-sim-3.Rdata")
```

## Results from Simulation 3

```
load("out-sim-3.Rdata")
```

### Half-logit bias

```
#### producing output for publication

True <- c(beta, theta)
cf1 <- t(sapply(res4, function(z) z[,"cflam_mv0"]))
bias1 <- t(sapply(res4, function(z) z[,"cflam_mv0"]) - True)
cf2 <- t(sapply(res4, function(z) z[,"cflam_mv1"]))
bias2 <- t(sapply(res4, function(z) z[,"cflam_mv1"]) - True)
cf3 <- t(sapply(res4, function(z) z[,"cflam_sv1"])) # logit
bias3 <- t(sapply(res4, function(z) z[,"cflam_sv1"]) - True)
cf4 <- t(sapply(res4, function(z) z[,"cflam_svhalf"]))
bias4 <- t(sapply(res4, function(z) z[,"cflam_svhalf"]) - True)

cf5 <- t(sapply(res4, function(z) z[,"cflam_gmv1"]))
bias5 <- t(sapply(res4, function(z) z[,"cflam_gmv1"]) - True)
cf6 <- t(sapply(res4, function(z) z[,"cflam_sv1stack"])) # logit
bias6 <- t(sapply(res4, function(z) z[,"cflam_sv1stack"]) - True)
cf7 <- t(sapply(res4, function(z) z[,"cflam_svhalfstack"]))
bias7 <- t(sapply(res4, function(z) z[,"cflam_svhalfstack"]) - True)

## MV bias is related to wrong model specification
## cannot separate 0.5 and p intercept

pdf("halflogit-bias-coefs.pdf", width=12, height=6)
par(mfrow=c(2,4))
boxplot(bias1, ylim=c(-5,5), main="MV N-mix closed")
abline(h=0, col=2)
boxplot(bias2, ylim=c(-5,5), main="MV N-mix open")
abline(h=0, col=2)
boxplot(bias5, ylim=c(-5,5), main="GMV N-mix open")
abline(h=0, col=2)
plot.new()
boxplot(bias3, ylim=c(-5,5), main="SV N-mix logit")
abline(h=0, col=2)
boxplot(bias4, ylim=c(-5,5), main="SV N-mix half-logit")
abline(h=0, col=2)
boxplot(bias6, ylim=c(-5,5), main="SV N-mix logit stack")
abline(h=0, col=2)
boxplot(bias7, ylim=c(-5,5), main="SV N-mix half-logit stack")
abline(h=0, col=2)
dev.off()

## estimate lam and p
lam1 <- apply(cf1, 1, function(z) mean(exp(X %*% z[1:2])))
lam2 <- apply(cf2, 1, function(z) mean(exp(X %*% z[1:2])))
lam3 <- apply(cf3, 1, function(z) mean(exp(X %*% z[1:2])))
lam4 <- apply(cf4, 1, function(z) mean(exp(X %*% z[1:2])))
lam5 <- apply(cf5, 1, function(z) mean(exp(X %*% z[1:2])))
lam6 <- apply(cf6, 1, function(z) mean(exp(X %*% z[1:2])))
lam7 <- apply(cf7, 1, function(z) mean(exp(X %*% z[1:2])))
p1 <- apply(cf1, 1, function(z) mean(plogis(X %*% z[3:4])))
p2 <- apply(cf2, 1, function(z) mean(plogis(X %*% z[3:4])))
p3 <- apply(cf3, 1, function(z) mean(plogis(X %*% z[3:4])))
p4 <- apply(cf4, 1, function(z) mean(0.5 * plogis(X %*% z[3:4])))
p5 <- apply(cf5, 1, function(z) mean(plogis(X %*% z[3:4])))
p6 <- apply(cf6, 1, function(z) mean(plogis(X %*% z[3:4])))
p7 <- apply(cf7, 1, function(z) mean(0.5 * plogis(X %*% z[3:4])))
lamTrue <- mean(exp(X %*% beta))
pTrue <- mean(plogis(X %*% theta))

SHOW <- c(1,2,6,7) # results dropped from publication

NAM <- c("MV\nclosed","MV\nopen","GMV\nopen",
"SV\nlogit","SV\nhalf-logit", "SV\nlogit","SV\nhalf-logit")[SHOW]
pdf("halflogit-bias-preds.pdf", width=6, height=10)
par(mfrow=c(2,1), mar=c(2, 5, 1, 2))
boxplot(cbind(MV0=lam1,MV1=lam2,GMV1=lam5,
    SVlogit=lam3,SVhalf=lam4,SSVlogit=lam6,SSVhalf=lam7)[,SHOW], 
    ylab=expression(hat(lambda)), names=NAM,
    ylim=c(0, max(lam2)), col="grey", axes=FALSE)
axis(2)
box()
abline(h=lamTrue, col=1, lty=2)
text(2, lamTrue*0.9, expression(bar(lambda)))
mtext(NAM, side=1, at=1:length(NAM), line=1.5)
par(mar=c(2, 5, 1, 2))
boxplot(cbind(MV0=p1,MV1=p2,GMV1=p5,
    SVlogit=p3,SVhalf=p4,SSVlogit=p6,SSVhalf=p7)[,SHOW], 
    ylab=expression(hat(p)), names=NAM, col="grey", axes=FALSE)
axis(2)
box()
abline(h=c(1, 0.5)*pTrue, col=1, lty=2)
text(2, pTrue*1.05, expression(bar(p)))
text(2, pTrue*0.5*1.1, expression(bar(p)*p[max]))
dev.off()
```


