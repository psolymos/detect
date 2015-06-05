## single-visit estimation
library(detect)
library(snow)
library(rlecuyer)

## for the sake of reproducibility
set.seed(1234)
B <- 100 # number of replicates
n <- 1000 # number of sites
x1 <- runif(n,0,1) # random predictors
x2 <- rnorm(n,0,1)
x3 <- runif(n,0,1)
x4 <- rnorm(n,0,1)
x <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4) # nice data frame
X <- model.matrix(~ x1 + x2, x)
Z1 <- model.matrix(~ x3 + x4, x) # no overlap with X
Z2 <- model.matrix(~ x3 + x2, x) # overlap with X

K <- 100 # N_max for numerical integration
T <- 3 # max number of visits

## true coefficients
beta <- c(2,-0.8, 0.6)
theta <- c(-1, 2, -0.5)
gvar <- 0.8

source("~/repos/detect/extras/revisitingSV/svabuScaledNB.R")

## runid: iterator, not used but needed for parallel functions
## q: 1/c, scaled logit constant 'c' needs to be >= 1
## overlap: is there overlap between abundance & detectability covariates?
fun <-
function(runid=NA, q=1, overlap=FALSE)
{

    ## model matrices
    Z <- if (!overlap)
        Z1 else Z2

    ## links and transformations
    linkinv <- binomial("logit")$linkinv
    p <- linkinv(drop(Z %*% theta))
    delta <- p * q
    lambda <- poisson("log")$linkinv(drop(X %*% beta))

    N <- rnbinom(n, size=1/gvar, prob=1/(1+gvar*lambda))
    Y <- rbinom(n, N, delta)

    m1 <- svabu_nb2.fit(Y, X, Z, Q = NULL,
        zeroinfl = FALSE, area = 1, N.max = max(K, max(N)+50),
        link.det = "logit", link.zif = "logit")

    cf1 <- unlist(coef(m1))
    cf0 <- c(beta, theta, qlogis(q))

    lam0 <- mean(lambda)
    lam1 <- mean(exp(drop(X %*% cf1[1:ncol(X)])))
    p0 <- mean(p)
    p1 <- mean(linkinv(drop(Z %*% cf1[(ncol(X)+1):(ncol(X)+ncol(Z))])))

    out <- cbind(truth=c(lam=lam0, p=p0, q=q, coef=cf0, log.gvar=log(gvar)),
        SV=c(lam1, p1, plogis(cf1[length(cf1)]), cf1, m1$var$est))

    attr(out, "overlap") <- overlap
    attr(out, "runid") <- runid
    out
}

## see how long it takes for a single run
#system.time(aa <- fun(q=1, overlap=FALSE))
#library(pbapply)
#aa1 <- pblapply(c(1, 0.75, 0.5, 0.25), function(qv) fun(q=qv, overlap=FALSE))
#aa2 <- pblapply(c(1, 0.75, 0.5, 0.25), function(qv) fun(q=qv, overlap=TRUE))


## parallel stuff
qval <- c(1, 0.75, 0.5, 0.25)
ncl <- 10
cl <- makeSOCKcluster(ncl)

## load pkgs on workers
clusterEvalQ(cl, library(detect))
## push data to workers
clusterExport(cl, c("n","x","X","Z1","Z2","K","T","B","beta","theta",
    "gvar", "svabu_nb2.fit"))
## set RNGs
clusterSetupRNG (cl, type = "RNGstream")

## magic happens here
resF <- list()
for (qv in qval) {
    cat("q =", qv, "no overlap\n")
    flush.console()
    res1 <- parLapply(cl, seq_len(B), fun,
            q=qv,
            overlap=FALSE)
    resF[[paste0("no overlap, q=",qv)]] <- res1
}
resT <- list()
for (qv in qval) {
    cat("q =", qv, "with overlap\n")
    flush.console()
    res1 <- parLapply(cl, seq_len(B), fun,
        q=qv,
        overlap=TRUE)
    resT[[paste0("with overlap, q=",qv)]] <- res1
}

## save results
save.image(file=paste0("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/",
    "MEE-rev2-simul-nbs-n", n, ".Rdata"))

## shutting down safely
stopCluster(cl)

load("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/MEE-rev2-simul-nbs.Rdata")

f <- function(xx)
    t(sapply(xx, function(z) z[,2] - z[,1]))

aa <- lapply(resF, f)
bb <- lapply(resT, f)

pf <- function(tp) {
    boxplot(tp, col="grey", ylim=c(-5,5))
    abline(h=0)
}
pf(aa[[1]])

g <- function(xx)
    sapply(xx, function(z) (z["lam",2]-z["lam",1]) / z["lam",1])

op <- par(mfrow=c(1,2))
boxplot(sapply(resF, g), col="grey", main="no overlap",
        names=c(1,0.75,05,0.25), ylim=c(-1,3))
abline(h=0)
boxplot(sapply(resT, g), col="grey", main="with overlap",
        names=c(1,0.75,05,0.25), ylim=c(-1,3))
abline(h=0)
par(op)
