## preliminaries for MPI
.Last <- function() {
    if (getOption("CLUSTER_ACTIVE")) {
        stopCluster(cl)
        cat("active cluster stopped by .Last\n")
    } else {
        cat("no active cluster found\n")
    }
}
options("CLUSTER_ACTIVE" = FALSE)
library(snow)
library(Rmpi)
library(rlecuyer)

## multiple-visits estimation
library(unmarked)
## single-visit estimation
library(detect)

## for the sake of reproducibility
set.seed(1234)
B <- 120 # number of replicates
n <- 200 # number of sites
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


## runid: iterator, not used but needed for parallel functions
## cvalue: scaled logit constant 'c' needs to be >= 1
## omega: constant survival omega
## overlap: is there overlap between abundance & detectability covariates?
fun <-
function(runid=NA, cvalue=1, omega=1, overlap=FALSE)
{
    ## sanity check
    if (cvalue < 1)
        stop("cvalue must be >= 1")
    if (omega < 0 || omega > 1)
        stop("omega must be in (0,1)")
    if (length(cvalue) > 1L)
        stop("cvalue length must be 1")
    if (length(omega) > 1L)
        stop("omega length must be 1")

    dmtype <- "constant"

    ## model matrices
    Z <- if (!overlap)
        Z1 else Z2

    ## links and transformations
    linkinv <- binomial("logit")$linkinv
    p <- linkinv(drop(Z %*% theta)) / cvalue
    lambda <- cvalue * poisson("log")$linkinv(drop(X %*% beta))
    if (dmtype == "notrend")
        gamma <- lambda * (1-omega)
    if (dmtype == "constant")
        gamma <- 1

    ## matrices to store simulated values
    N <- Y <- matrix(NA, n, T)
    S <- G <- matrix(NA, n, T-1)
    N[,1] <- rpois(n, lambda) # 1st visit
    Y[,1] <- rbinom(n, N[,1], p) # 1st visit

    ## simulate values for subsequent visits
    for(i in 1:(T-1)) {
        S[,i] <- rbinom(n, N[,i], omega) # survival from t-1
        G[,i] <- rpois(n, gamma) # arrivals
        N[,i+1] <- S[,i] + G[,i] # new N is survived + arrived
        Y[,i+1] <- rbinom(n, N[,i+1], p) # new Y as observed
    }

    ## generalized N-mixture, T=3
    visitMat <- matrix(as.character(1:3), n, 3, byrow=TRUE)
    umfg <- unmarkedFramePCO(y=Y[,1:3], siteCovs=x, obsCovs=list(visit=visitMat),
        numPrimary=3)
    m0 <- if (!overlap) {
        pcountOpen(~ x1 + x2, ~ 1, ~ 1, ~ x3 + x4, umfg, K=K, mixture="P",
            dynamics=dmtype)
    } else {
        pcountOpen(~ x1 + x2, ~ 1, ~ 1, ~ x3 + x2, umfg, K=K, mixture="P",
            dynamics=dmtype)
    }

    ## single visit with n obs, T=1
    m1 <- if (!overlap) {
        svabu(Y[,1] ~ x1 + x2 | x3 + x4, x,
            zeroinfl=FALSE, link.det="logit", N.max=K)
    } else {
        svabu(Y[,1] ~ x1 + x2 | x3 + x2, x,
            zeroinfl=FALSE, link.det="logit", N.max=K)
    }
    m2 <- if (!overlap) {
        svabu(Y[,2] ~ x1 + x2 | x3 + x4, x,
              zeroinfl=FALSE, link.det="logit", N.max=K)
    } else {
        svabu(Y[,2] ~ x1 + x2 | x3 + x2, x,
              zeroinfl=FALSE, link.det="logit", N.max=K)
    }
    m3 <- if (!overlap) {
        svabu(Y[,3] ~ x1 + x2 | x3 + x4, x,
              zeroinfl=FALSE, link.det="logit", N.max=K)
    } else {
        svabu(Y[,3] ~ x1 + x2 | x3 + x2, x,
              zeroinfl=FALSE, link.det="logit", N.max=K)
    }

    beta0 <- coef(m0)[grepl("lam\\(", names(coef(m0)))]
    beta1 <- coef(m1, "sta")
    beta2 <- coef(m2, "sta")
    beta3 <- coef(m3, "sta")

    theta0 <- coef(m0)[grepl("p\\(", names(coef(m0)))]
    theta1 <- coef(m1, "det")
    theta2 <- coef(m2, "det")
    theta3 <- coef(m3, "det")


    lamt1 <- mean(lambda)
    lamt2 <- mean(lamt1 * omega + gamma)
    lamt3 <- mean(lamt2 * omega + gamma)
    lam0 <- exp(drop(X %*% beta0))
    lam01 <- mean(lam0)
    lam02 <- mean(lam01 * plogis(coef(m0)["omega(Int)"])) +
        exp(coef(m0)["gamConst(Int)"])
    lam03 <- mean(lam02 * plogis(coef(m0)["omega(Int)"])) +
        exp(coef(m0)["gamConst(Int)"])
    lam1 <- mean(exp(drop(X %*% beta1)))
    lam2 <- mean(exp(drop(X %*% beta2)))
    lam3 <- mean(exp(drop(X %*% beta3)))
    pt <- mean(p)
    p0 <- mean(linkinv(drop(Z %*% theta0)))
    p1 <- mean(linkinv(drop(Z %*% theta1)))
    p2 <- mean(linkinv(drop(Z %*% theta2)))
    p3 <- mean(linkinv(drop(Z %*% theta3)))

    out <- cbind(truth=c(lamt1, lamt2, lamt3, pt, pt, pt),
        DM=c(lam01, lam02, lam03, p0, p0, p0),
        SV=c(lam1, lam2, lam3, p1, p2, p3))
    rownames(out) <- c("lam1","lam2","lam3","p1","p2","p3")

    attr(out, "DM") <- coef(m0)
    attr(out, "SV1") <- coef(m1)
    attr(out, "SV2") <- coef(m2)
    attr(out, "SV3") <- coef(m3)
    attr(out, "overlap") <- overlap
    attr(out, "cvalue") <- cvalue
    attr(out, "omega") <- omega
    attr(out, "dmtype") <- dmtype
    attr(out, "runid") <- runid
    out
}

## see how long it takes for a single run
#system.time(aa <- fun(cvalue=1, omega=1, overlap=FALSE))

## settings
vals <- expand.grid(cvalue=1/c((1:10)/10),
    omega=seq(0, 1, by=0.1),
    overlap=c(FALSE,TRUE))

(args <- commandArgs(trailingOnly = TRUE))
nodes <- as.numeric(args[1])
ncl <- nodes * 12

TEST <- FALSE
if (TEST) {
    B <- 2
    T <- 5
    vals <- vals[1:2,]
    ncl <- 2
}

## parallel stuff
cl <- makeMPIcluster(ncl)
options("CLUSTER_ACTIVE" = TRUE)

## load pkgs on workers
clusterEvalQ(cl, library(unmarked))
clusterEvalQ(cl, library(detect))
## push data to workers
clusterExport(cl, c("n","x","X","Z1","Z2","K","T","B","beta","theta"))
## set RNGs
clusterSetupRNG (cl, type = "RNGstream")

## magic happens here
res <- list()
for (rowid in seq_len(nrow(vals))) {
    cat(rowid, "of", nrow(vals), "\n")
    flush.console()
    res1 <- parLapply(cl, seq_len(B), fun,
            cvalue=vals[rowid,"cvalue"],
            omega=vals[rowid,"omega"],
            overlap=vals[rowid,"overlap"])
    save(res1, file=paste0("out/simres4-constant-", rowid, ".Rdata"))
    res[[rowid]] <- res1
}

## save results
save.image(file="MEE-rev2-simul-4-constant.Rdata")

## shutting down safely
stopCluster(cl)
options("CLUSTER_ACTIVE" = FALSE)
mpi.quit("no")
