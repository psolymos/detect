## initial values need to be supplied from outside
## this uses 0 if not specified
## type can be several of c("ML", "CL", "PL")
## fit can be a zi.fit (list) or numeric vector of logd0
zi.fit <-
function(Y, X, Z, offsetx, offsetz, weights,
distr=c("pois","negbin","binom","lognorm","beta"), linkx, linkz="logit",
type="ML", fit, N,
method="Nelder-Mead", inits, control=list(), hessian=FALSE, ...)
{

    ## >>> ZI-LOGNORM
    ## linkinvx=identity (because mu in dlnorm is on log scale)
    nll_ZILN_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        sigma <- exp(parms[kx + 1])
        phi <- as.vector(linkinvz(Z %*% parms[(kx+1+1):(kx+kz+1)] + offsetz))
        ## 0 mass of Ln leads to 0 density
        loglik0 <- log(phi)
        loglik1 <- log(1 - phi) + dlnorm(Y, mu, sigma, log = TRUE)
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] *
            loglik1[id1])
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    nll_ZILN_CL <- function(parms) {
        mu1 <- as.vector(linkinvx(X1 %*% parms[1:kx] + offsetx1))
        sigma <- exp(parms[kx + 1])
        loglik <- sum(weights1 * dlnorm(Y1, mu1, sigma, log = TRUE))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    logd0_ZILN <- function(parms) -Inf
    nll_ZILN_PL <- function(parms, logd0) {
        phi <- as.vector(linkinvz(Z %*% parms + offsetz))
        loglik <- sum(weights * dbinom(W, 1, 1-phi, log=TRUE))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    ## >>> ZI-BETA
    ## linkinvx=ilogit
    ## linkinv for gamma=exp
    nll_ZIB_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        gamma <- exp(parms[kx + 1]) # precision
        phi <- as.vector(linkinvz(Z %*% parms[(kx+1+1):(kx+kz+1)] + offsetz))
        alpha <- mu * gamma
        beta <- (1 - mu) * gamma
        loglik0 <- log(phi)
        loglik1 <- log(1 - phi) + suppressWarnings(dbeta(Y,
            alpha, beta, log = TRUE))
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] *
            loglik1[id1])
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    nll_ZIB_CL <- function(parms) {
        mu1 <- as.vector(linkinvx(X1 %*% parms[1:kx] + offsetx1))
        gamma <- exp(parms[kx + 1]) # precision
        alpha1 <- mu1 * gamma
        beta1 <- (1 - mu1) * gamma
        loglik <- sum(weights1 * suppressWarnings(dbeta(Y1,
            alpha1, beta1, log = TRUE)))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    logd0_ZIB <- logd0_ZILN
    nll_ZIB_PL <- nll_ZILN_PL

    ## >>> ZI-POIS
    ## linkinvx=exp
    nll_ZIP_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
        loglik0 <- log(phi + exp(log(1 - phi) - mu))
        loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] *
            loglik1[id1])
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    nll_ZIP_CL <- function(parms) {
        mu1 <- as.vector(linkinvx(X1 %*% parms[1:kx] + offsetx1))
        ## P(Y=y|Y>0)=f(y;theta)/(1-0)=f(y;theta)
        num <- dpois(Y1, mu1, log = TRUE)
        den <- log(1 - exp(-mu1))
        loglik <- sum(weights1 * (num - den))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    logd0_ZIP <- function(parms) {
        -as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
    }
    ## this function applies to all distributions
    ## because logd0 (log density for 0 obs) is plug-in
    ## logd0 needs to incorporate offsets but not weights
    nll_ZIP_PL <- function(parms, logd0) {
        phi <- as.vector(linkinvz(Z %*% parms + offsetz))
        loglik0 <- log(phi + exp(log(1 - phi) + logd0))
        loglik1 <- log(1 - phi) + log(1 - exp(logd0))
        loglik <- sum(weights * ifelse(Y==0, loglik0, loglik1))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    ## >>> ZI-NEGBIN
    ## linkinvx=exp
    ## theta is precision
    nll_ZINB_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        theta <- exp(parms[kx + 1])
        phi <- as.vector(linkinvz(Z %*% parms[(kx+1+1):(kx+kz+1)] + offsetz))
        loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0,
            size = theta, mu = mu, log = TRUE))))
        loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y,
            size = theta, mu = mu, log = TRUE))
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] *
            loglik1[id1])
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    nll_ZINB_CL <- function(parms) {
        mu1 <- as.vector(linkinvx(X1 %*% parms[1:kx] + offsetx1))
        theta <- exp(parms[kx + 1])
        ## P(Y=y|Y>0)=f(y;theta)/(1-0)=f(y;theta)
        num <- suppressWarnings(dnbinom(Y1,
            size = theta, mu = mu1, log = TRUE))
        den <- log(1 - exp(suppressWarnings(dnbinom(0,
            size = theta, mu = mu1, log = TRUE))))
        loglik <- sum(weights1 * (num - den))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    logd0_ZINB <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        theta <- exp(parms[kx + 1])
        suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))
    }
    nll_ZINB_PL <- nll_ZIP_PL


    ## >>> ZI-BINOM
    ## linkinvx=exp
    nll_ZIBIN_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
        loglik0 <- log(phi + (1 - phi) * (1-mu)^N)
        loglik1 <- log(1 - phi) + dbinom(Y, N, mu, log = TRUE)
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] *
            loglik1[id1])
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    nll_ZIBIN_CL <- function(parms) {
        mu1 <- as.vector(linkinvx(X1 %*% parms[1:kx] + offsetx1))
        num <- dbinom(Y1, N1, mu1, log = TRUE)
        den <- log(1 - (1-mu1)^N1)
        loglik <- sum(weights1 * (num - den))
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }
    logd0_ZIBIN <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        log((1-mu)^N)
    }
    nll_ZIBIN_PL <- nll_ZIP_PL


    good.num.limit = c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    kx <- ncol(X)
    kz <- ncol(Z)
    np <- kx + kz
    n <- length(Y)
    distr <- match.arg(distr)
    if (distr %in% c("negbin","lognorm","beta"))
        np <- np + 1
    if (missing(offsetx))
        offsetx <- rep(0, n)
    if (missing(offsetz))
        offsetz <- rep(0, n)
    if (missing(weights))
        weights <- rep(1, n)
    if (missing(linkx))
        linkx <- switch(distr,
            "pois" = "log",
            "negbin" = "log",
            "binom" = "logit",
            "lognorm" = "identity", # (because mu in dlnorm is on log scale)
            "beta" = "logit")
    linkinvx <- poisson(linkx)$linkinv
    linkinvz <- binomial(linkz)$linkinv
    id1 <- Y > 0
    id0 <- !id1
    W <- ifelse(id1, 1L, 0L)
    Y1 <- Y[id1]
    X1 <- X[id1,,drop=FALSE]
    Z1 <- Z[id1,,drop=FALSE]
    offsetx1 <- offsetx[id1]
    offsetz1 <- offsetz[id1]
    weights1 <- weights[id1]

    if (distr == "binom") {
        if (missing(N))
            stop("Binomial size N missing")
        if (length(N) != n)
            N <- rep(N, n)[1:n]
        N1 <- N[id1]
    }

    type <- match.arg(type, c("ML", "CL", "PL"), several.ok=TRUE)
    if (missing(inits))
        inits <- rep(0, np)
    res_ML <- res_CLPL <- res_CL <- res_PL <- NULL
    control$fnscale <- 1 # minimize
    if ("ML" %in% type) {
        nll_ML <- switch(distr,
            "pois" = nll_ZIP_ML,
            "negbin" = nll_ZINB_ML,
            "binom" = nll_ZIBIN_ML,
            "lognorm" = nll_ZILN_ML,
            "beta" = nll_ZIB_ML)
        t_ML <- system.time(res_ML <- optim(inits, nll_ML,
            method=method, hessian=hessian, control=control, ...))
        res_ML$coef <- res_ML$par
        res_ML$loglik <- -res_ML$value
        res_ML$vcov <- NULL
        if (hessian)
            res_ML$vcov <- solve(res_ML$hessian)
        res_ML$time <- t_ML
        res_ML$par <- res_ML$value <- NULL
    }
    if ("CL" %in% type) {
        nll_CL <- switch(distr,
            "pois" = nll_ZIP_CL,
            "negbin" = nll_ZINB_CL,
            "binom" = nll_ZIBIN_CL,
            "lognorm" = nll_ZILN_CL,
            "beta" = nll_ZIB_CL)
        t_CL <- system.time(res_CL <- suppressWarnings(optim(inits[1:(np-kz)], nll_CL,
            method=method, hessian=hessian, control=control, ...)))
        res_CL$coef <- res_CL$par
        res_CL$loglik <- -res_CL$value
        res_CL$vcov <- NULL
        if (hessian)
            res_CL$vcov <- solve(res_CL$hessian)
        res_CL$time <- t_CL
        res_CL$par <- res_CL$value <- NULL
    }
    if ("PL"  %in% type) {
        ## logd0 (log density for 0 obs) is plug-in
        ## logd0 needs to incorporate offsets but not weights
        logd0_fun <- switch(distr,
            "pois" = logd0_ZIP,
            "negbin" = logd0_ZINB,
            "binom" = logd0_ZIBIN,
            "lognorm" = logd0_ZILN,
            "beta" = logd0_ZIB)
        if (is.null(res_CL)) {
            if (missing(fit))
                stop("type PL needs fit")
            if (is.list(fit)) {
                ## fit is a previous object from zi.fit
                res_CL <- fit$CL
                logd0 <- logd0_fun(fit$CL$coef)
            } else {
                ## fit is the logdensity f(y=0)
                #res_CL <- NULL # already NULL
                logd0 <- fit
            }
        } else {
            if (!missing(fit))
                warning("type includes CL, fit argument ignored")
            logd0 <- logd0_fun(res_CL$coef)
        }
        nll_PL <- switch(distr,
            "pois" = nll_ZIP_PL,
            "negbin" = nll_ZINB_PL,
            "binom" = nll_ZIBIN_PL,
            "lognorm" = nll_ZILN_PL,
            "beta" = nll_ZIB_PL)
        t_PL <- system.time(res_PL <- suppressWarnings(optim(inits[(np-kz+1):np], nll_PL,
            logd0=logd0,
            method=method, hessian=hessian, control=control, ...)))
        res_PL$coef <- res_PL$par
        res_PL$loglik <- -res_PL$value
        res_PL$vcov <- NULL
        if (hessian)
            res_PL$vcov <- solve(res_PL$hessian)
        res_PL$time <- t_PL
        res_PL$par <- res_PL$value <- NULL
    }
    if (!is.null(res_CL) && !is.null(res_PL)) {
        nll_ML <- switch(distr,
            "pois" = nll_ZIP_ML,
            "negbin" = nll_ZINB_ML,
            "binom" = nll_ZIBIN_ML,
            "lognorm" = nll_ZILN_ML,
            "beta" = nll_ZIB_ML)
        res_CLPL <- list(coef = c(res_CL$coef, res_PL$coef))
        res_CLPL$loglik <- -nll_ML(res_CLPL$coef)
        res_CLPL$time <- res_CL$time + res_PL$time
        res_CLPL$vcov <- NULL
        res_CLPL$convergence <- c(CL=res_CL$convergence, PL=res_PL$convergence)
        if (hessian) {
            res_CLPL$vcov <- matrix(NA, np, np)
            res_CLPL$vcov[1:(np-kz),1:(np-kz)] <- res_CL$vcov
            res_CLPL$vcov[(np-kz+1):np,(np-kz+1):np] <- res_PL$vcov
        }

    }
    list(ML=res_ML, CL=res_CL, PL=res_PL, CLPL=res_CLPL)
}
