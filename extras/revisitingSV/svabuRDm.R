## single visit Muntinom-Binom-ZIP abundance model
## with distance sampling
## Y is a matrix matching D
## D is is matrix with PC radius endpoints,
## we **DO NOT** allow r=Inf unlimited radius in D

svabuRDm.fit <-
function(Y, X, ZR=NULL, ZD=NULL, Q=NULL, zeroinfl=TRUE,
D, N.max=NULL, inits,
link.det = "logit", link.zif = "logit", ...)
{
    nll.P <- function(parms) {
        edr <- exp(drop(ZD %*% parms[(np.abu+np.detR+1):(np.abu+np.detR+np.detD)]))
        deltaD <- (edr / r)^2 * (1 - exp(-(D / edr)^2))
        deltaD <- deltaD - cbind(0, deltaD[, -ncol(deltaD), drop=FALSE])

        deltaR <- drop(linkinvfun.det(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        ## lambda is D*A, but the notation is not fully followed
        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each=N.max+1) * area.rep

        ## this is p*pi cell probs
        delta <- deltaR * deltaD
        ## unobserved (pi_0)
        delta0 <- 1 - rowSums(delta, na.rm=TRUE)

        dp <- dpois(N.rep, lambda.rep)

        delta.rep <- delta[id1.repx,]
        delta0.rep <- rep(delta0, each=N.max+1)
        #db <- dbinom(Y.rep, N.rep, delta.rep)
        db <- sapply(1:length(N.rep), function(i) {
            yvec <- c(YY.rep[i,Yok.rep[i,]], YY0.rep[i])
            if (N.rep[i] < sum(yvec))
                return(NA)
            dmultinom(yvec, N.rep[i],
                c(delta.rep[i,Yok.rep[i,]], delta0.rep[i]))})

        intsum <- colSums(matrix(dp * db, nrow=N.max+1), na.rm=TRUE)
        out <- -sum(log(intsum))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.ZIP1 <- function(parms) {
        edr <- exp(drop(ZD %*% parms[(np.abu+np.detR+1):(np.abu+np.detR+np.detD)]))
        deltaD <- (edr / r)^2 * (1 - exp(-(D / edr)^2))
        deltaD <- deltaD - cbind(0, deltaD[, -ncol(deltaD), drop=FALSE])

        lambda1.rep <- rep(drop(exp(X1 %*% parms[1:np.abu])), each=N.max+1) * area1.rep
        #srate1 <- exp(drop(ZR1 %*% parms[(np.abu+1):(np.abu+np.detR)]))
        #deltaR1 <- 1 - exp(-srate1 * t1)
        deltaR <- drop(linkinvfun.det(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        ## this is p*pi cell probs
        delta <- deltaR * deltaD
        ## unobserved
        delta0 <- 1 - rowSums(delta, na.rm=TRUE)
        delta1 <- rowSums(delta, na.rm=TRUE)[id1]
        delta1.rep <- delta[id1.repx,]
        delta01.rep <- delta0[id1.repx]

        dp1 <- dpois(N1.rep, lambda1.rep)

        #db1 <- dbinom(Y1.rep, N1.rep, delta1.rep)
        db1 <- sapply(1:length(N1.rep), function(i) {
            yvec <- c(YY1.rep[i,Yok1.rep[i,]], YY01.rep[i])
            if (N1.rep[i] < sum(yvec))
                return(NA)
            dmultinom(yvec, N1.rep[i],
                c(delta1.rep[i,Yok1.rep[i,]], delta01.rep[i]))})

        part1 <- colSums(matrix(dp1 * db1, nrow=N.max+1), na.rm=TRUE)
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
        deltaD <- (edr / r)^2 * (1 - exp(-(D / edr)^2))
        deltaD <- deltaD - cbind(0, deltaD[, -ncol(deltaD), drop=FALSE])

        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each=N.max+1) * area.rep
        #srate <- exp(drop(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))
        #deltaR <- 1 - exp(-srate * t)
        deltaR <- drop(linkinvfun.det(ZR %*% parms[(np.abu+1):(np.abu+np.detR)]))

        ## this is p*pi cell probs
        delta <- deltaR * deltaD
        ## unobserved (pi_0)
        delta0 <- 1 - rowSums(delta, na.rm=TRUE)

        delta.rep <- delta[i.rep,]
        delta0.rep <- rep(delta0, each=N.max+1)
        phi.rep <- rep(drop(linkinvfun.zif(Q %*% parms[(np.abu+np.detR+np.detD+1):length(parms)])), each=N.max+1)
        loglik0 <- log(phi.rep + exp(log(1 - phi.rep) - lambda.rep))
        loglik1 <- log(1 - phi.rep) + dpois(N.rep, lambda=lambda.rep, log=TRUE)
        dp <- exp(ifelse(id1.rep, loglik1, loglik0))
#        dp <- exp(ifelse(N1full.rep, loglik1, loglik0))
#        db <- dbinom(Y.rep, N.rep, delta.rep)
        db <- sapply(1:length(N.rep), function(i) {
            yvec <- c(YY.rep[i,Yok.rep[i,]], YY0.rep[i])
            if (N.rep[i] < sum(yvec))
                return(NA)
            dmultinom(yvec, N.rep[i],
                c(delta.rep[i,Yok.rep[i,]], delta0.rep[i]))})
        intsum <- colSums(matrix(dp * db, nrow=N.max+1), na.rm=TRUE)
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
    ## http://tolstoy.newcastle.edu.au/R/help/05/04/3272.html
    good.num.limit = c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    if (is.null(dim(Y)))
        stop("Y must be matrix")
    if (is.null(dim(D)))
        stop("D must be matrix")
    YY <- Y
    Yok <- !is.na(Y)
    Y <- rowSums(Y, na.rm=TRUE)
    if (any(!is.finite(rowSums(D, na.rm=TRUE))))
        stop("unlimited distances not supported")

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

    ## check YY and D pattern
    cmx <- if (ncol(ZD) < 2)
        NULL else ZD
    ledr <- try(cmulti.fit(YY, D, cmx, "dis"))
    edr.hat <- if (!inherits(ledr, "try-error"))
        ledr$coefficients else rep(0, ncol(ZD))

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
                edr.hat, # glm.fit(Z,ifelse(Y>0,1,0),family=binomial())$coef
                glm.fit(Q,ifelse(Y>0,0,1),family=binomial())$coef)
        } else {
            c(glm.fit(X,Y,family=poisson())$coef,
                glm.fit(ZR,Y,family=poisson())$coef,
                edr.hat) # glm.fit(Z,ifelse(Y>0,1,0),family=binomial())$coef
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
    i.rep <- rep(seq_len(n), each=N.max+1)
    YY.rep <- YY[i.rep,]
    Yok.rep <- Yok[i.rep,]
    D.rep <- D[i.rep,]
    r <- apply(D, 1, max, na.rm=TRUE) # truncation distance
    r.rep <- rep(r, each=N.max+1)
    area <- pi*r^2
    area.rep <- rep(area, each=N.max+1)
    ## unobserved count
    YY0.rep <- pmax(0, N.rep - Y.rep)

    ## estimate the count part
    if (zeroinfl) {
        id1 <- Y > 0
        id1.rep <- rep(id1, each=N.max+1)
        id1.repx <- rep(which(Y > 0), each=N.max+1)
        Y1 <- Y[id1]
        r1 <- r[id1]
        YY1 <- YY[id1,]
        Yok1 <- YY[id1,]
        n1 <- length(Y1)
        N1.rep <- rep(0:N.max, n1) # needed in nll.ZIP1
#        N1full.rep <- N.rep > 0 # needed in nll.ZIP
        Y1.rep <- rep(Y1, each=N.max+1)
        YY1.rep <- YY[id1.repx,]
        YY01.rep <- pmax(0, N1.rep - Y1.rep)
        #YY01.rep <- YY0.rep[id1.repx]
        D1.rep <- D[id1.repx,]
        Yok1.rep <- Yok[id1.repx,]
        X1 <- data.matrix(X[id1,])
        ZR1 <- data.matrix(ZR[id1,])
        ZD1 <- data.matrix(ZD[id1,])
        area1 <- area[id1]
        area1.rep <- rep(area1, each=N.max+1)
        #r1 <- r[id1]
        D1 <- D[id1,]
        results <- optim(inits[1:(np.abu + np.detR + np.detD)],
            nll.ZIP1, method=method, hessian=TRUE, control=control)
    } else {
        id1 <- Y >= 0
        id1.rep <- rep(id1, each=N.max+1)
        id1.repx <- rep(seq_len(n), each=N.max+1)
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
    deltaD <- (edr / r)^2 * (1 - exp(-(D / edr)^2))
    delta <- deltaR * rowSums(deltaD, na.rm=TRUE)

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

## this is constant EDR model

set.seed(1234)
library(detect)
library(unmarked)
source("~/repos/detect/extras/revisitingSV/svabuRDm.R")
n <- 200
T <- 3
K <- 50
B <- 100
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

beta <- c(0,1)
thetaR <- c(1, -1.5) # for singing rate
#thetaD <- c(-0.5, 1.2) # for EDR

edr <- 0.8 # exp(drop(ZD %*% thetaD))
Dm <- matrix(c(0.5, 1), n, 2, byrow=TRUE)
## truncation distance must be finite
r <- apply(Dm, 1, max, na.rm=TRUE)
Area <- pi*r^2
q <- (edr / r)^2 * (1 - exp(-(Dm / edr)^2))
q <- q - cbind(0, q[,-ncol(Dm), drop=FALSE])

## test case 1
D <- exp(drop(X %*% beta))
lambda <- D * Area
p <- linkinvfun.det(drop(ZR %*% thetaR))
delta <- cbind(p * q, 1-rowSums(p * q))
summary(delta)
summary(lambda)
summary(rowSums(q))
summary(p)

res_mn0 <- list()
res_mnp <- list()
res_sv <- list()
res_mn <- list()
for (i in 1:B) {
    cat("variable p, run", i, "of", B, ":\t");flush.console()

    N <- rpois(n, lambda)
    Y10 <- t(sapply(1:n, function(i)
        rmultinom(1, N[i], delta[i,])))
    Y <- Y10[,-ncol(Y10)]

    zi <- FALSE

    cat("mn_p,  ");flush.console()
    m0 <- try(svabuRDm.fit(Y, X, NULL, NULL, Q=NULL, zeroinfl=zi, D=Dm, N.max=K))
    res_mn0[[i]] <- try(cbind(est=unlist(coef(m0)), true=c(beta, mean(qlogis(p)), log(edr))))

    cat("mn_pi,  ");flush.console()
    m1 <- try(svabuRDm.fit(Y, X, ZR, NULL, Q=NULL, zeroinfl=zi, D=Dm, N.max=K))
    res_mnp[[i]] <- try(cbind(est=unlist(coef(m1)), true=c(beta, thetaR, log(edr))))

    cat("mn,  ");flush.console()
    umf <- unmarkedFrameDS(y=Y,
        siteCovs=data.frame(x1=x1,x2=x2,x3=x3),
        dist.breaks=c(0,50,100), unitsIn="m", survey="point")
    m <- distsamp(~1 ~x1, umf, output="abund")
    sig <- exp(coef(m, type="det"))
    ea <- 2*pi * integrate(grhn, 0, 100, sigma=sig)$value # effective area
    logedr <- log(sqrt(ea / pi)/100) # effective radius
    res_mn[[i]] <- cbind(est=c(coef(m)[1:2], logedr), true=c(beta, log(edr)))

    cat("b_pi.\n")
    yy <- rowSums(Y)
    m2 <- svabu.fit(yy, X, ZR, Q=NULL, zeroinfl=zi, N.max=K)
    res_sv[[i]] <- cbind(est=unlist(coef(m2)), true=c(beta, thetaR))

}
save.image(paste0("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom_final.Rdata"))

## test case 2
res2_mn0 <- list()
res2_mn <- list()
p <- 1
delta <- cbind(p * q, 1-rowSums(p * q))
for (i in 1:B) {
    cat("variable p, 2, run", i, "of", B, ":\t");flush.console()

    N <- rpois(n, lambda)
    Y10 <- t(sapply(1:n, function(i)
        rmultinom(1, N[i], delta[i,])))
    Y <- Y10[,-ncol(Y10)]

    zi <- FALSE

    cat("mn_p,  ");flush.console()
    m0 <- svabuRDm.fit(Y, X, NULL, NULL, Q=NULL, zeroinfl=zi, D=Dm, N.max=K)
    res2_mn0[[i]] <- cbind(est=unlist(coef(m0)), true=c(beta, mean(qlogis(p)), log(edr)))

    cat("mn\n");flush.console()
    umf <- unmarkedFrameDS(y=Y,
        siteCovs=data.frame(x1=x1,x2=x2,x3=x3),
        dist.breaks=c(0,50,100), unitsIn="m", survey="point")
    m <- distsamp(~1 ~x1, umf, output="abund")
    sig <- exp(coef(m, type="det"))
    ea <- 2*pi * integrate(grhn, 0, 100, sigma=sig)$value # effective area
    logedr <- log(sqrt(ea / pi)/100) # effective radius
    res2_mn[[i]] <- cbind(est=c(coef(m)[1:2], logedr), true=c(beta, log(edr)))

}
save.image(paste0("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom_all4_",
    n, ".Rdata"))

## test case 3
res3_mn0 <- list()
res3_mn <- list()
p <- 0.5
delta <- cbind(p * q, 1-rowSums(p * q))
for (i in 1:B) {
    cat("variable p, 3, run", i, "of", B, ":\t");flush.console()

    N <- rpois(n, lambda)
    Y10 <- t(sapply(1:n, function(i)
        rmultinom(1, N[i], delta[i,])))
    Y <- Y10[,-ncol(Y10)]

    zi <- FALSE

    cat("mn_p,  ");flush.console()
    m0 <- svabuRDm.fit(Y, X, NULL, NULL, Q=NULL, zeroinfl=zi, D=Dm, N.max=K)
    res3_mn0[[i]] <- cbind(est=unlist(coef(m0)), true=c(beta, mean(qlogis(p)), log(edr)))

    cat("mn\n");flush.console()
    umf <- unmarkedFrameDS(y=Y,
        siteCovs=data.frame(x1=x1,x2=x2,x3=x3),
        dist.breaks=c(0,50,100), unitsIn="m", survey="point")
    m <- distsamp(~1 ~x1, umf, output="abund")
    sig <- exp(coef(m, type="det"))
    ea <- 2*pi * integrate(grhn, 0, 100, sigma=sig)$value # effective area
    logedr <- log(sqrt(ea / pi)/100) # effective radius
    res3_mn[[i]] <- cbind(est=c(coef(m)[1:2], logedr), true=c(beta, log(edr)))

}
save.image(paste0("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom_all4_",
    n, ".Rdata"))

## unmarked script -- scaling biases lambda estimates

set.seed(1234)
library(unmarked)
n <- 200
T <- 3
K <- 50
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

beta <- c(1,1)
p <- 0.2
edr <- 50

Dm <- matrix(c(50, 100), n, 2, byrow=TRUE)
## truncation distance must be finite
r <- apply(Dm, 1, max, na.rm=TRUE)
q <- (edr / r)^2 * (1 - exp(-(Dm / edr)^2))
q <- q - cbind(0, q[,-ncol(Dm), drop=FALSE])
## unobserved
delta <- cbind(p * q, 1-rowSums(p * q))
lambda <- exp(drop(X %*% beta))


N <- rpois(n, lambda)
Y10 <- t(sapply(1:n, function(i)
    rmultinom(1, N[i], delta[i,])))
Y <- Y10[,-ncol(Y10)]

summary(N)
summary(Y)
summary(delta)
summary(rowSums(Y)/N)
summary(rowSums(delta[,-3]))

umf <- unmarkedFrameDS(y=Y,
 siteCovs=data.frame(x1=x1,x2=x2,x3=x3),
 dist.breaks=c(0,50,100), unitsIn="m", survey="point")

m <- distsamp(~1 ~x1, umf, output="abund")
summary(m)

## effective radius
sig <- exp(coef(m, type="det"))
ea <- 2*pi * integrate(grhn, 0, 100, sigma=sig)$value # effective area
edr
sqrt(ea / pi) # effective radius
mean(lambda)
mean(lambda*p)
mean(exp(drop(X %*% coef(m)[1:2])))

## plot the results

load("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom-cval_all4_200.Rdata")
load("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom_c-2_n-200.Rdata")
#"res_mn"         "res_mn0"        "res_mnp"        "res_sv"
#"res2_mn"        "res2_mn0"
#"res3_mn"        "res3_mn0"

load("~/Dropbox/pkg/detect2/mee-rebuttal/rev2/multinom_final.Rdata")
f <- function(res) {
    true <- res[[1]][,"true"]
    true[!is.finite(true)] <- 0
    est <- t(sapply(res, function(z) 
        if (inherits(z, "try-error")) rep(NA, length(true)) else z[,"est"]))
    bias <- t(t(est) - true)
    list(true=true, est=est, bias=bias)
}

## bias in parameters
par(mfrow=c(2,2))
ylim <- c(-2,2)
boxplot(f(res_mn0)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res_mnp)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res_sv)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res_mn)$bias, ylim=ylim);abline(h=0, col=2)

par(mfrow=c(2,2))
ylim <- c(-2,2)
boxplot(f(res2_mn)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res2_mn0)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res3_mn)$bias, ylim=ylim);abline(h=0, col=2)
boxplot(f(res3_mn0)$bias, ylim=ylim);abline(h=0, col=2)

## bias in lambda
#lam <- mean(exp(X %*% z[1:2,"true"]))
lam <- mean(lambda)
g0 <- function(z, A=1) {
    lam.hat <- mean(A * exp(X %*% z[1:2,"est"]))
    (lam.hat - lam) / lam
    #lam.hat / lam
}
g <- function(res, ...)
    sapply(res, g0, ...)

## Royle & SV abundance, MN is density (D * 1^2*pi)
bLam <- cbind(Royle=g(res_mn),
    MN0=g(res_mn0, pi), MNi=g(res_mnp, pi), SV=g(res_sv))
bLam2 <- cbind(Royle=g(res2_mn), MN0=g(res2_mn0, pi))
bLam3 <- cbind(Royle=g(res3_mn), MN0=g(res3_mn0, pi))

par(mfrow=c(1,3))
ylim <- c(-1,1)
boxplot(bLam, ylim=ylim);abline(h=0)
boxplot(bLam2, ylim=ylim);abline(h=0)
boxplot(bLam3, ylim=ylim);abline(h=0)

ylim <- c(-1,0.5)
toPlot <- bLam[,c(2,4,1)]
colnames(toPlot) <- c("Multinomial", "SV", "Distsamp")
boxplot(toPlot, ylim=ylim, col="grey", ylab="Relative bias")
abline(h=0)

ylim <- c(-1,0.5)
toPlot <- bLam[,c(2,3,4,1)]
colnames(toPlot) <- c("Multinomial, p", "Multinom, p_i", "SV", "Distsamp")
boxplot(toPlot, ylim=ylim, col="grey", ylab="Relative bias")
abline(h=0)

ylim <- c(-2,0.5)
toPlot <- cbind(f(res_mn0)$est[,1],
    #f(res_mnp)$est[,1],
    f(res_sv)$est[,1]-log(pi))#,
    #f(res_mn)$est[,1]-log(pi))
colnames(toPlot) <- c("Multinomial", #"Multinomial_p", 
    "SV")#, "Distsamp")
boxplot(toPlot, col="grey", ylab="Bias")
abline(h=0)
abline(h=log(1/2), lty=2)
text(0.6, log(1/2)+0.08, "log(q)", cex=0.8)



}
