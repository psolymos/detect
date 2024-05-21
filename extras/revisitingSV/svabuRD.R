## single visit B-B-ZIP abundance model
## with distance sampling
## r is PC radius, we **DO NOT** allow r=Inf unlimited radius

svabuRD.fit <-
function(Y, X, ZR=NULL, ZD=NULL, Q=NULL, zeroinfl=TRUE,
r=1, N.max=NULL, inits,
link.det = "logit", link.zif = "logit", ...)
{
    nll.P <- function(parms) {
        ## EDR
        edr <- exp(drop(ZD %*% parms[(np.abu+np.detR+1):(np.abu+np.detR+np.detD)]))
        deltaD <- (edr / r)^2 * (1 - exp(-(r / edr)^2))
        ## lambda is D*A, but the notation is not fully followed
        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each=N.max+1) * area.rep
        ## singing rate
        #srate <- exp(drop(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        #deltaR <- 1 - exp(-srate * t)
        deltaR <- drop(linkinvfun.det(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        delta <- deltaR * deltaD
        delta.rep <- rep(delta, each=N.max+1)
        dp <- dpois(N.rep, lambda.rep)
        db <- dbinom(Y.rep, N.rep, delta.rep)
        intsum <- colSums(matrix(dp * db, nrow=N.max+1))
        out <- -sum(log(intsum))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.ZIP1 <- function(parms) {
        edr1 <- exp(drop(ZD1 %*% parms[(np.abu+np.detR+1):(np.abu+np.detR+np.detD)]))
        deltaD1 <- (edr1 / r1)^2 * (1 - exp(-(r1 / edr1)^2))
        lambda1.rep <- rep(drop(exp(X1 %*% parms[1:np.abu])), each=N.max+1) * area1.rep
        #srate1 <- exp(drop(ZR1 %*% parms[(np.abu+1):(np.abu+np.detR)]))
        #deltaR1 <- 1 - exp(-srate1 * t1)
        deltaR1 <- drop(linkinvfun.det(ZR1 %*% parms[(np.abu+1):(np.abu+np.detR)]))
        delta1 <- deltaR1 * deltaD1
        delta1.rep <- rep(delta1, each=N.max+1)
        dp1 <- dpois(N1.rep, lambda1.rep)
        db1 <- dbinom(Y1.rep, N1.rep, delta1.rep)
        part1 <- colSums(matrix(dp1 * db1, nrow=N.max+1))
        part2 <- 1 / (1 - exp(-lambda1.rep * delta1))
        out <- -sum(log(part1*part2))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.ZIP0 <- function(parms, emlp, id1) {
        phi <- drop(linkinvfun.zif(Q %*% parms))
        tmp0 <- exp(sum(log(phi[!id1] + (1 - phi[!id1]) * emlp[!id1])))
        tmp1 <- prod((1 - phi[id1]))
        out <- -log(tmp0 * tmp1)
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.ZIP <- function(parms) {
        edr <- exp(drop(ZD %*% parms[(np.abu+np.detR+1):(np.abu+np.detR+np.detD)]))
        deltaD <- (edr / r)^2 * (1 - exp(-(r / edr)^2))
        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each=N.max+1) * area.rep
        #srate <- exp(drop(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        #deltaR <- 1 - exp(-srate * t)
        deltaR <- drop(linkinvfun.det(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        delta <- deltaR * deltaD
        delta.rep <- rep(delta, each=N.max+1)
        phi.rep <- rep(drop(linkinvfun.zif(Q %*% parms[(np.abu+np.detR+np.detD+1):length(parms)])), each=N.max+1)
        loglik0 <- log(phi.rep + exp(log(1 - phi.rep) - lambda.rep))
        loglik1 <- log(1 - phi.rep) + dpois(N.rep, lambda=lambda.rep, log=TRUE)
        dp <- exp(ifelse(id1.rep, loglik1, loglik0))
#        dp <- exp(ifelse(N1full.rep, loglik1, loglik0))
        db <- dbinom(Y.rep, N.rep, delta.rep)
        intsum <- colSums(matrix(dp * db, nrow=N.max+1))
        out <- -sum(log(intsum))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    ## preliminaries
    control <- getOption("detect.optim.control")
    ## negLogLik must be minimized
    control$fnscale <- 1
    method <- getOption("detect.optim.method")
    if (is.null(method))
        method <- "Nelder-Mead"
    Control <- list(optim.control=control, optim.method=method)
    ## from post of Spencer Graves to avoid: optim() non-finite finite-difference value
    ## https://tolstoy.newcastle.edu.au/R/help/05/04/3272.html
    good.num.limit = c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    n <- length(Y)
    if (is.null(ZR)) {
        ZR <- matrix(1, n, 1)
        colnames(ZR) <- "(Intercept)"
    }
    if (is.null(ZD)) {
        ZD <- matrix(1, n, 1)
        colnames(ZD) <- "(Intercept)"
    }
    nam.abu <- colnames(X)
    nam.detR <- paste0("R_", colnames(ZR))
    nam.detD <- paste0("D_", colnames(ZD))
    np.abu <- NCOL(X)
    np.detR <- NCOL(ZR)
    np.detD <- NCOL(ZD)

    ## zeroinfl and site specific phi
    if (zeroinfl && is.null(Q)) {
        Q <- matrix(1, n, 1)
        colnames(Q) <- "(Intercept)"
    }
    if (zeroinfl) {
        np.zif <- NCOL(Q)
        nam.zif <- colnames(Q)
    } else {
        np.zif <- 0
    }

    np <- np.abu + np.detR + np.detD + np.zif

    linkinvfun.det <- binomial(link=make.link(link.det))$linkinv
    linkinvfun.zif <- binomial(link=make.link(link.zif))$linkinv

    if(missing(inits)) {
#        inits <- rep(0, np)
        inits <- if (zeroinfl) {
            c(glm.fit(X,Y,family=poisson())$coef,
                glm.fit(ZR,Y,family=poisson())$coef,
                glm.fit(ZD,Y,family=poisson())$coef, # glm.fit(Z,ifelse(Y>0,1,0),family=binomial())$coef
                glm.fit(Q,ifelse(Y>0,0,1),family=binomial())$coef)
        } else {
            c(glm.fit(X,Y,family=poisson())$coef,
                glm.fit(ZR,Y,family=poisson())$coef,
                glm.fit(ZD,Y,family=poisson())$coef) # glm.fit(Z,ifelse(Y>0,1,0),family=binomial())$coef
        }
        inits[is.na(inits)] <- 0
    }

#    if (inherits(inits, "svisit"))
#        inits <- coef(inits)
    inits <- rep(inits, np)[1:np]

    if (is.null(N.max))
        N.max <- max(max(Y, na.rm = TRUE)+20, max(Y, na.rm = TRUE)*3)
    if(N.max <= max(Y, na.rm = TRUE))
        stop("N.max too small")
    N.rep <- rep(0:N.max, n)
    Y.rep <- rep(Y, each=N.max+1)

    ## infinite r leads to infinite area - use edr.tr
    ## edr.tr cannot be used in Poisson
    if (any(!is.finite(r)))
        stop("r must be finite")
    r.rep <- if (length(r) == 1)
        r else rep(r, each=N.max+1)
#    t.rep <- if (length(t) == 1)
#        t else rep(t, each=N.max+1)
    area <- pi*r^2
    area.rep <- if (length(area) == 1)
        area else rep(area, each=N.max+1)

    ## estimate the count part
    if (zeroinfl) {
        id1 <- Y > 0
        id1.rep <- rep(id1, each=N.max+1)
        Y1 <- Y[id1]
        n1 <- length(Y1)
        N1.rep <- rep(0:N.max, n1) # needed in nll.ZIP1
#        N1full.rep <- N.rep > 0 # needed in nll.ZIP
        Y1.rep <- rep(Y1, each=N.max+1)
        X1 <- data.matrix(X[id1,])
        ZR1 <- data.matrix(ZR[id1,])
        ZD1 <- data.matrix(ZD[id1,])
        if (length(r) == 1) {
            area1 <- area
            area1.rep <- area
        } else {
            area1 <- area[id1]
            area1.rep <- rep(area1, each=N.max+1)
        }
        r1 <- if (length(r) == 1)
            r else r[id1]
#        t1 <- if (length(t) == 1)
#            t else t[id1]
        results <- optim(inits[1:(np.abu + np.detR + np.detD)],
            nll.ZIP1, method=method, hessian=TRUE, control=control)
    } else {
        id1 <- Y >= 0
        id1.rep <- rep(id1, each=N.max+1)
        results <- optim(inits[1:(np.abu + np.detR + np.detD)],
            nll.P, method=method, hessian=TRUE, control=control)
    }
    estimates <- results$par
    par.state <- estimates[1:np.abu]
    par.det <- estimates[(np.abu+1):(np.abu+np.detR+np.detD)]
    names(par.state) <- nam.abu
    names(par.det) <- c(nam.detR, nam.detD)
    if (rcond(results$hessian) <= 1e-06)
        std.error <- rep(NA, np)
    if (rcond(results$hessian) > 1e-06) {
        ## due to negLogLik, we take H^-1 and not -H^-1
        opvar <- diag(solve(results$hessian))
        if (any(opvar < 0)) {
            opvar[opvar < 0] <- NA
            warning("negative variance values in optim, NAs produced")
        }
        std.error <- sqrt(opvar)
    }
    se.state <- std.error[1:np.abu]
    se.det <- std.error[(np.abu+1):(np.abu+np.detR+np.detD)]
    names(se.state) <- nam.abu
    names(se.det) <- c(nam.detR, nam.detD)
    lambda <- drop(exp(X %*% par.state))

    #singrate <- exp(drop(ZR %*% par.det[1:np.detR]))
    #deltaR <- 1 - exp(-singrate * t)
    deltaR <- drop(linkinvfun.det(ZR %*% par.det[1:np.detR]))
#    singrate <- -log(1-deltaR) / t
    edr <- exp(drop(ZD %*% par.det[(np.detR+1):(np.detR+np.detD)]))
    deltaD <- (edr / r)^2 * (1 - exp(-(r / edr)^2))
    delta <- deltaR * deltaD

    ## estimate the zero part
    zif.results <- NULL
    par.zif <- if (zeroinfl)
        0 else NULL
    se.zif <- NULL
    phi <- NULL
    lLik <- -results$value
    if (sum(!id1) > 0 && zeroinfl) {
        emlp <- exp(-lambda * area * delta)
        zif.results <- suppressWarnings(optim(inits[(np.abu+np.detR+np.detD+1):np],
            nll.ZIP0, emlp=emlp, id1=id1, method=method, hessian=TRUE, control=control))
        par.zif <- zif.results$par
        names(par.zif) <- nam.zif
        phi <- drop(linkinvfun.zif(Q %*% par.zif))
        lLik <- -nll.ZIP(c(par.state, par.det, par.zif))

        if (rcond(zif.results$hessian) <= 1e-06)
            se.zif <- rep(NA, np.zif)
        if (rcond(zif.results$hessian) > 1e-06) {
            ## due to negLogLik, we take H^-1 and not -H^-1
            opvar2 <- diag(solve(data.matrix(zif.results$hessian)))
            if (any(opvar2 < 0)) {
                opvar2[opvar2 < 0] <- NA
                warning("negative variance values in optim, NAs produced")
            }
            se.zif <- sqrt(opvar2)
        }
        names(se.zif) <- nam.zif
    }
    ## assembling return object
    Converged <- if (sum(!id1) > 0 && zeroinfl) {
        results$convergence == 0 && zif.results$convergence == 0
    } else results$convergence == 0
    out <- list(coefficients = list(sta = par.state, det = par.det),
        std.error = list(sta = se.state, det = se.det),
        fitted.values = lambda,
        detection.probabilities = delta,
#        singing.rate = singrate,
        edr = edr,
        zif.probabilities = phi,
        zeroinfl = zeroinfl,
        nobs = n,
        N.max = N.max,
        link = list(sta="log", det=link.det, zif=link.zif),
        df.null = n - 2,
        df.residual = n - np,
        inits = inits,
        loglik = lLik,
        results = list(count=results, zero=zif.results),
        converged = Converged,
        control = Control,
        t = t,
        r = r)
    if (zeroinfl) {
        out$coefficients$zif <- par.zif
        out$std.error$zif <- se.zif
    }
    if (!out$converged)
        warning("model did not converge")
    return(out)
}

if (FALSE) {

#library(detect)
#source("c:/Dropbox/pkg/detect/R/svabusra.R")
#source("c:/Dropbox/pkg/detect/R/svabusra.fit.R")
#svisitFormula <- detect:::svisitFormula

set.seed(4321)
n <- 2000
T <- 3
K <- 250
link <- "logit"
linkinvfun.det <- binomial(link=make.link(link))$linkinv

x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)
x3 <- runif(n,-1,1)
x4 <- runif(n,-1,1)
x5 <- rbinom(n,1,0.6)
x6 <- rbinom(n,1,0.4)
x7 <- rnorm(n,0,1)

X <- model.matrix(~ x1)
ZR <- model.matrix(~ x3)
ZD <- model.matrix(~ x2)
Q <- model.matrix(~ x7)

beta <- c(-1,1)
thetaR <- c(0.2, -2) # for singing rate
thetaD <- c(-0.3, 0.3) # for EDR
phi <- 0

r <- sample(c(0.5,1), n, replace=TRUE)
#t <- sample(c(3,5,4,5,6,7,8,9,10), n, replace=TRUE)

edr <- exp(drop(ZD %*% thetaD))
#edr.tr <- edr*sqrt(1-exp(-r^2/edr^2))
#Area <- ifelse(is.finite(r), pi*edr.tr^2, pi*edr^2)
Area <- ifelse(is.finite(r), pi*r^2, pi*edr^2)
q <- ifelse(is.finite(r), (edr / r)^2 * (1 - exp(-(r / edr)^2)), 1)

D <- exp(drop(X %*% beta))
lambda <- D * Area

p <- linkinvfun.det(drop(ZR %*% thetaR))
#srate <- -log(1-deltaR) / t

## using the marginal distribution, so that we can simulate unlimited counts
## the interest is in D

Y <- rpois(n, lambda * p * q)

#x <- data.frame(x1, x2, x3, x4, x5, x6)
summary(Y)
summary(Y/(p*q))
summary(D)
summary(p)
summary(q)
#summary(pi*r^2)
table(Y)


m <- svabuRD.fit(Y, X, ZR, ZD, Q=NULL, zeroinfl=phi!=0, r=r, N.max=K)

cbind(est=unlist(coef(m)), true=c(beta, thetaR, thetaD))

#                         est true
#sta.(Intercept)   -0.9205812 -1.0
#sta.x1             0.7472176  1.0
#det.R_(Intercept)  0.3954751  0.2
#det.R_x3          -2.5397480 -2.0
#det.D_(Intercept) -0.5145289 -0.5
#det.D_x2           1.2768396  1.2


## this is constant EDR model

set.seed(4321)
n <- 2000
T <- 3
K <- 250
link <- "logit"
linkinvfun.det <- binomial(link=make.link(link))$linkinv

x1 <- runif(n,0,1)
x2 <- rnorm(n,0,1)
x3 <- runif(n,-1,1)
x4 <- runif(n,-1,1)
x5 <- rbinom(n,1,0.6)
x6 <- rbinom(n,1,0.4)
x7 <- rnorm(n,0,1)

X <- model.matrix(~ x1)
ZR <- model.matrix(~ x3)
ZD <- model.matrix(~ x2)
Q <- model.matrix(~ x7)

beta <- c(-1,1)
thetaR <- c(0.2, -2) # for singing rate
thetaD <- c(-0.5, 1.2) # for EDR
phi <- 0

r <- sample(c(0.5,1), n, replace=TRUE)
#t <- sample(c(3,5,4,5,6,7,8,9,10), n, replace=TRUE)

edr <- 0.65 # exp(drop(ZD %*% thetaD))
Area <- ifelse(is.finite(r), pi*r^2, pi*edr^2)
q <- ifelse(is.finite(r), (edr / r)^2 * (1 - exp(-(r / edr)^2)), 1)

D <- exp(drop(X %*% beta))
lambda <- D * Area

p <- linkinvfun.det(drop(ZR %*% thetaR))
#srate <- -log(1-deltaR) / t

## using the marginal distribution, so that we can simulate unlimited counts
## the interest is in D

Y <- rpois(n, lambda * p * q)

summary(Y)
summary(Y/(p*q))
summary(D)
summary(p)
summary(q)
table(Y)


m <- svabuRD.fit(Y, X, ZR, NULL, Q=NULL, zeroinfl=FALSE, t=t, r=r, N.max=K)
cbind(est=unlist(coef(m)), true=c(beta, thetaR, log(edr)))

m <- svabuRD.fit(Y, X, ZR, NULL, Q=NULL, zeroinfl=TRUE, t=t, r=r, N.max=K)
cbind(est=unlist(coef(m))[-length(unlist(coef(m)))], true=c(beta, thetaR, log(edr)))

}
