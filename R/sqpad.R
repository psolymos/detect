# SQPAD estimation

# replacing extraDistr::rtriang(n, 0, r, r)
rtriang <- function(n, r = 1) {
    o <- numeric(0L)
    while (length(o) < n) {
        m <- n - length(o)
        d <- sqrt(runif(m * 1.3, 0, r)^2 + runif(m * 1.3, 0, r)^2)
        d <- d[d <= r]
        o <- c(o, d)
    }
    o[seq_len(n)]
}

.sqpad.fit <- function(
    Y,
    dis,
    dur,
    X = NULL,
    Z = NULL,
    type = c("full", "approx", "tt1", "conv"),
    det = c("joint", "pq"),
    init = NULL,
    method = "Nelder-Mead",
    hessian = TRUE,
    tt1 = NULL,
    A = NULL,
    Nmax = 50,
    K = NULL,
    montecarlo = FALSE,
    np = 1000,
    distcorr = 2 / 3,
    dislist = NULL,
    check_diff = TRUE,
    ...
) {
    .pfun <- function(dur, dis, delta, tau, det, distcorr, montecarlo, dij) {
        if (montecarlo) {
            DMAT <- dij * dis # dij is a n by np matrix
            switch(
                det,
                "pq" = rowMeans(
                    (1 - exp(-dur * delta)) *
                        (tau^2 * (1 - exp(-(DMAT)^2 / tau^2)) / (DMAT)^2)
                ),
                "joint" = rowMeans(
                    1 - exp(-dur * delta * exp(-(DMAT)^2 / tau^2))
                ),
                stop("Unrecognized detection function")
            )
        } else {
            switch(
                det,
                "pq" = (1 - exp(-dur * delta)) *
                    (tau^2 *
                        (1 - exp(-(distcorr * dis)^2 / tau^2)) /
                        (distcorr * dis)^2),
                "joint" = 1 -
                    exp(-dur * delta * exp(-(distcorr * dis)^2 / tau^2)),
                stop("Unrecognized detection function")
            )
        }
    }
    .solvenear <- function(x) {
        xinv <- try(solve(x), silent = TRUE)
        if (inherits(xinv, "try-error")) {
            xinv <- as.matrix(solve(Matrix::nearPD(x)$mat))
        }
        xinv
    }

    type <- match.arg(type)
    det <- match.arg(det)
    distcorr <- as.numeric(distcorr[1L])
    if (!is.null(K)) {
        K <- as.integer(K[1L])
    }

    if (!montecarlo && distcorr <= 0) {
        stop("distcorr must be positive")
    }
    if (montecarlo && np < 1) {
        stop("np must be positive")
    }
    if (length(Y) != length(dis) || length(Y) != length(dur)) {
        stop("Y, dis, and dur must have same length.")
    }
    if (any(is.infinite(dis))) {
        stop("dis cannot be infinite.")
    }
    if (any(is.infinite(dur))) {
        stop("dur cannot be infinite.")
    }
    if (any(is.na(dis))) {
        stop("dis cannot be NA.")
    }
    if (any(is.na(dur))) {
        stop("dur cannot be NA.")
    }
    if (any(dis <= 0)) {
        stop("dis must be positive.")
    }
    if (any(dur <= 0)) {
        stop("dur must be positive.")
    }

    if (type == "conv") {
        if (is.null(dislist) || !is.list(dislist)) {
            stop(
                "dislist must be a list with individual distances if type='conv'."
            )
        }
        if (!is.null(X)) {
            stop("X must be NULL if type='conv'.")
        }
        if (!is.null(Z)) {
            stop("Z must be NULL if type='conv'.")
        }
        if (det != "joint") {
            stop("det must be joint if type='conv'.")
        }
    }

    n <- length(Y)
    if (is.null(A)) {
        A <- dis^2 * pi
    }
    if (is.null(X)) {
        X <- matrix(1, n, 1)
        colnames(X) <- "(Intercept)"
    }
    if (is.null(Z)) {
        Z <- matrix(1, n, 1)
        colnames(Z) <- "(Intercept)"
    }
    mx <- NCOL(X)
    mz <- NCOL(Z)
    ix <- seq_len(mx)
    iz <- seq_len(mz) + mx
    ma <- mx + mz
    ma <- ma + 1L

    if (is.null(colnames(X))) {
        colnames(X) <- paste0("v", seq_len(ncol(X)))
    }
    if (is.null(colnames(Z))) {
        colnames(Z) <- paste0("v", seq_len(ncol(Z)))
    }

    if (is.null(init)) {
        cf_tmp <- glm(Y ~ X - 1, family = poisson, offset = log(A))
        init <- c(unname(coef(cf_tmp)), rep(0, mz + 1))
    }

    if (montecarlo) {
        # dij <- extraDistr::rtriang(as.integer(np), 0, 1, 1)
        dij <- rtriang(as.integer(np), 1)
        dij <- matrix(dij, nrow = n, ncol = np, byrow = TRUE)
    } else {
        dij <- NULL
    }

    # convolution likelihood
    # empirical distances are used (no Monte Carlo or approximation needed)
    if (type == "conv") {
        obs <- NULL
        for (i in seq_len(n)) {
            nint <- 100
            dcats <- seq_len(ceiling(dis[i] * nint))
            d <- data.frame(r1 = (dcats - 1) / nint, r2 = dcats / nint, Y = 0)
            if (Y[i] > 0) {
                dc <- as.integer(dislist[[i]] * nint) + 1
                for (dcc in dc) {
                    d$Y[dcc] <- d$Y[dcc] + 1
                }
            }
            d$a <- d$r2^2 * pi - d$r1^2 * pi # area of annulus
            d$dur <- dur[i]
            obs <- rbind(obs, d)
        }
        nll1 <- function(parms, obs, n, Nmax) {
            D.hat <- exp(parms[1])
            phi <- exp(parms[2])
            tau <- exp(parms[3])
            r <- (obs$r1 + obs$r2) / 2
            p <- 1 - exp(-obs$dur * phi * exp(-r^2 / tau^2))
            M <- matrix(0, nrow(obs), Nmax + 1)
            for (i in seq_len(Nmax + 1)) {
                Ni <- i - 1
                M[, i] <- dpois(Ni, D.hat * obs$a) * dbinom(obs$Y, Ni, p)
            }
            d <- rowSums(M)
            rv <- -sum(log(d))
            if (is.na(rv) || is.infinite(rv)) {
                rv <- 10^18
            }
            rv
        }
        opt <- optim(
            init,
            nll1,
            method = method,
            hessian = hessian,
            obs = obs,
            n = n,
            Nmax = Nmax,
            ...
        )
        if (check_diff && abs(max(opt$par - init)) < 10^-6) {
            init <- rnorm(length(init), init, 0.05)
            opt <- optim(
                init,
                nll1,
                method = method,
                hessian = hessian,
                obs = obs,
                n = n,
                Nmax = Nmax,
                ...
            )
        }
    }

    if (type == "approx") {
        # approximate abundance model
        if (is.null(K)) {
            nll2 <- function(
                parms,
                Y,
                A,
                dur,
                dis,
                X,
                Z,
                det,
                distcorr,
                montecarlo,
                dij
            ) {
                D.hat <- exp(drop(X %*% parms[ix]))
                delta <- exp(drop(Z %*% parms[iz]))
                tau <- exp(parms[ma])

                p.hat <- .pfun(
                    dur,
                    dis,
                    delta,
                    tau,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij
                )
                lam.hat <- D.hat * A
                rv <- -sum(dpois(Y, lam.hat * p.hat, log = TRUE))
                if (is.na(rv) || is.infinite(rv)) {
                    rv <- 10^18
                }
                rv
            }
            opt <- optim(
                init,
                nll2,
                method = method,
                hessian = hessian,
                Y = Y,
                A = A,
                X = X,
                Z = Z,
                dur = dur,
                dis = dis,
                det = det,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij,
                ...
            )
            if (check_diff && abs(max(opt$par - init)) < 10^-6) {
                init <- rnorm(length(init), init, 0.05)
                opt <- optim(
                    init,
                    nll2,
                    method = method,
                    hessian = hessian,
                    Y = Y,
                    A = A,
                    X = X,
                    Z = Z,
                    dur = dur,
                    dis = dis,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij,
                    ...
                )
            }
        } else {
            # approximate occupancy (0/1) model
            if (K > 1) {
                stop("Use K=1 with the approximation")
            }
            nll3 <- function(
                parms,
                Y01,
                A,
                dur,
                dis,
                X,
                Z,
                det,
                distcorr,
                montecarlo,
                dij
            ) {
                D.hat <- exp(drop(X %*% parms[ix]))
                delta <- exp(drop(Z %*% parms[iz]))
                tau <- exp(parms[ma])

                p.hat <- .pfun(
                    dur,
                    dis,
                    delta,
                    tau,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij
                )
                lam.hat <- D.hat * A
                rv <- -sum(dbinom(
                    Y01,
                    1,
                    (1 - exp(-lam.hat)) * p.hat,
                    log = TRUE
                ))
                if (is.na(rv) || is.infinite(rv)) {
                    rv <- 10^18
                }
                rv
            }
            opt <- optim(
                init,
                nll3,
                method = method,
                hessian = hessian,
                Y01 = ifelse(Y > 0, 1, 0),
                A = A,
                X = X,
                Z = Z,
                dur = dur,
                dis = dis,
                det = det,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij,
                ...
            )
            if (check_diff && abs(max(opt$par - init)) < 10^-6) {
                init <- rnorm(length(init), init, 0.05)
                opt <- optim(
                    init,
                    nll3,
                    method = method,
                    hessian = hessian,
                    Y01 = ifelse(Y > 0, 1, 0),
                    A = A,
                    X = X,
                    Z = Z,
                    dur = dur,
                    dis = dis,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij,
                    ...
                )
            }
        }
    }

    if (type == "full") {
        # full likelihood abundance model
        if (is.null(K)) {
            nll4 <- function(
                parms,
                Y,
                A,
                dur,
                dis,
                X,
                Z,
                Nmax,
                det,
                distcorr,
                montecarlo,
                dij
            ) {
                D.hat <- exp(drop(X %*% parms[ix]))
                delta <- exp(drop(Z %*% parms[iz]))
                tau <- exp(parms[ma])

                p.hat <- .pfun(
                    dur,
                    dis,
                    delta,
                    tau,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij
                )
                lam.hat <- D.hat * A
                d <- sapply(seq_along(Y), function(i) {
                    sum(sapply(Y[i]:Nmax, function(Ni) {
                        dpois(Ni, lam.hat[i]) * dbinom(Y[i], Ni, p.hat[i])
                    }))
                })
                rv <- -sum(log(d))
                if (is.na(rv) || is.infinite(rv)) {
                    rv <- 10^18
                }
                rv
            }
            opt <- optim(
                init,
                nll4,
                method = method,
                hessian = hessian,
                Y = Y,
                A = A,
                X = X,
                Z = Z,
                dur = dur,
                dis = dis,
                det = det,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij,
                Nmax = Nmax,
                ...
            )
            if (check_diff && abs(max(opt$par - init)) < 10^-6) {
                init <- rnorm(length(init), init, 0.05)
                opt <- optim(
                    init,
                    nll4,
                    method = method,
                    hessian = hessian,
                    Y = Y,
                    A = A,
                    X = X,
                    Z = Z,
                    dur = dur,
                    dis = dis,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij,
                    Nmax = Nmax,
                    ...
                )
            }
        } else {
            # full likelihood occupancy model
            if (K < 1 || K > max(Y, na.rm = TRUE) - 1) {
                stop("K must be NULL or between 1 and max(Y)-1")
            }
            if (K >= Nmax) {
                stop("Nmax must be higher than K")
            }
            YK <- Y
            YK[YK > K] <- K
            nll5 <- function(
                parms,
                YK,
                A,
                dur,
                dis,
                X,
                Z,
                K,
                Nmax,
                det,
                distcorr,
                montecarlo,
                dij
            ) {
                D.hat <- exp(drop(X %*% parms[ix]))
                delta <- exp(drop(Z %*% parms[iz]))
                tau <- exp(parms[ma])

                p.hat <- .pfun(
                    dur,
                    dis,
                    delta,
                    tau,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij
                )
                lam.hat <- D.hat * A

                # Y = 0 case
                d0 <- sapply(which(YK == 0), function(i) {
                    sum(sapply(0:Nmax, function(Ni) {
                        dpois(Ni, lam.hat[i]) * dbinom(0, Ni, p.hat[i])
                    }))
                })
                # Y >= K case
                dk <- sapply(which(YK == K), function(i) {
                    sum(sapply(K:Nmax, function(Ni) {
                        dpois(Ni, lam.hat[i]) *
                            sum(
                                sapply(K:Ni, function(j) {
                                    dbinom(j, Ni, p.hat[i])
                                })
                            )
                    }))
                })
                d <- c(d0, dk)
                # Y = 1,...,K-1 cases
                if (K > 1) {
                    for (k in 1:(K - 1)) {
                        d1 <- sapply(which(YK == k), function(i) {
                            sum(sapply(k:Nmax, function(Ni) {
                                dpois(Ni, lam.hat[i]) * dbinom(k, Ni, p.hat[i])
                            }))
                        })
                        d <- c(d, d1)
                    }
                }

                rv <- -sum(log(d))
                if (is.na(rv) || is.infinite(rv)) {
                    rv <- 10^18
                }
                rv
            }
            opt <- optim(
                init,
                nll5,
                method = method,
                hessian = hessian,
                YK = YK,
                A = A,
                X = X,
                Z = Z,
                dur = dur,
                dis = dis,
                det = det,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij,
                Nmax = Nmax,
                K = K,
                ...
            )
            if (check_diff && abs(max(opt$par - init)) < 10^-6) {
                init <- rnorm(length(init), init, 0.05)
                opt <- optim(
                    init,
                    nll5,
                    method = method,
                    hessian = hessian,
                    YK = YK,
                    A = A,
                    X = X,
                    Z = Z,
                    dur = dur,
                    dis = dis,
                    det = det,
                    distcorr = distcorr,
                    montecarlo = montecarlo,
                    dij = dij,
                    Nmax = Nmax,
                    K = K,
                    ...
                )
            }
        }
    }

    if (type == "tt1") {
        # time to 1st perception data (only 0 or >0 data)
        if (det != "joint") {
            stop("det argument must be joing with type='tt1'")
        }
        if (!is.null(K) && K > 1) {
            warning("K argument is ignored")
        }
        if (is.null(tt1)) {
            stop("provide tt1")
        }
        tt1 <- ifelse(Y > 0, tt1, dur)
        if (any(tt1 > dur)) {
            stop("tt1 cannot be higher than dur")
        }
        nll6 <- function(
            parms,
            Y,
            A,
            tt1,
            dur,
            dis,
            X,
            Z,
            Nmax,
            distcorr,
            montecarlo,
            dij
        ) {
            D.hat <- exp(drop(X %*% parms[ix]))
            delta <- exp(drop(Z %*% parms[iz]))
            tau <- exp(parms[ma])

            lam.hat <- D.hat * A
            p.hat1 <- delta *
                exp(-(dis * distcorr)^2 / tau^2) *
                exp(-delta * exp(-(dis * distcorr)^2 / tau^2) * tt1)
            p.hat0 <- .pfun(
                dur,
                dis,
                delta,
                tau,
                det = det,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij
            )
            Int <- sapply(1:Nmax, function(j) {
                dpois(j, lam.hat) * dbinom(0, j, p.hat0)
            })
            rv <- -sum(log(
                ifelse(
                    Y > 0,
                    (1 - exp(-lam.hat)) * p.hat1, # Y>0
                    exp(-lam.hat) + sum(Int)
                )
            )) # Y=0
            if (is.na(rv) || is.infinite(rv)) {
                rv <- 10^18
            }
            rv
        }
        opt <- optim(
            init,
            nll6,
            method = method,
            hessian = hessian,
            Y = Y,
            A = A,
            X = X,
            Z = Z,
            tt1 = tt1,
            dur = dur,
            dis = dis,
            Nmax = Nmax,
            distcorr = distcorr,
            montecarlo = montecarlo,
            dij = dij,
            ...
        )
        if (check_diff && abs(max(opt$par - init)) < 10^-6) {
            init <- rnorm(length(init), init, 0.05)
            opt <- optim(
                init,
                nll6,
                method = method,
                hessian = hessian,
                Y = Y,
                A = A,
                X = X,
                Z = Z,
                tt1 = tt1,
                dur = dur,
                dis = dis,
                Nmax = Nmax,
                distcorr = distcorr,
                montecarlo = montecarlo,
                dij = dij,
                ...
            )
        }
    }

    cf <- opt$par
    # names(cf) <- c(paste0("X_", colnames(X)), paste0("Z_", colnames(Z)), "tau")
    names(cf) <- c(
        paste0("log.D_", colnames(X)),
        paste0("log.delta_", colnames(Z)),
        "log.phi"
    )
    ll <- -opt$value
    if (hessian) {
        S <- .solvenear(opt$hessian)
        se <- sqrt(diag(S))
    } else {
        S <- matrix(NA_real_, length(cf), length(cf))
        se <- rep(NA_real_, length(cf))
    }
    dimnames(S) <- list(names(cf), names(cf))
    names(se) <- names(cf)
    tstat <- cf / se
    pval <- 2 * pnorm(-abs(tstat))
    coefs <- cbind(cf, se, tstat, pval)
    colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(coefs) <- names(cf)
    out <- list(
        coefficients = cf,
        loglik = ll,
        vcov = S,
        nobs = length(Y),
        df = length(cf),
        summary = coefs,
        Y = Y,
        X = X,
        Z = Z,
        A = A,
        dur = dur,
        dis = dis,
        tt1 = tt1,
        type = type,
        det = det,
        Nmax = Nmax,
        K = K,
        dislist = dislist,
        montecarlo = montecarlo,
        np = np,
        dij = dij[1L, ],
        distcorr = distcorr,
        init = init,
        hessian = hessian,
        method = method,
        results = opt
    )
    out
}

sqpad.fit <- function(
    Y,
    dis,
    dur,
    X = NULL,
    Z = NULL,
    type = c("full", "approx", "conv"),
    det = c("joint", "pq"),
    init = NULL,
    method = "Nelder-Mead",
    hessian = TRUE,
    tt1 = NULL,
    A = NULL,
    Nmax = NULL,
    K = NULL,
    montecarlo = FALSE,
    np = 1000,
    distcorr = 2 / 3,
    dislist = NULL,
    ...
) {
    Nmaxmax <- 100
    type <- match.arg(type)

    if (type == "full" && is.null(init)) {
        qout <- .sqpad.fit(
            Y = Y,
            dis = dis,
            dur = dur,
            X = X,
            Z = Z,
            type = "approx",
            det = det,
            tt1 = tt1,
            A = A,
            init = init,
            method = method,
            hessian = hessian,
            Nmax = if (is.null(Nmax)) 50 else Nmax,
            K = if (is.null(K)) K else 1L,
            distcorr = distcorr,
            ...
        )
        init <- qout$coefficients
        if (is.null(Nmax)) {
            Dhat <- exp(drop(qout$X %*% init[seq_len(ncol(qout$X))]))
            Nmax <- min(
                Nmaxmax,
                max(rpois(10^4, max(Dhat * dis^2 * pi, na.rm = TRUE)))
            )
        }
    }

    if (is.null(Nmax)) {
        Nmax <- 50
    }
    out <- .sqpad.fit(
        Y = Y,
        dis = dis,
        dur = dur,
        X = X,
        Z = Z,
        type = type,
        det = det,
        init = init,
        method = method,
        hessian = hessian,
        tt1 = tt1,
        A = A,
        Nmax = Nmax,
        K = K,
        montecarlo = montecarlo,
        np = np,
        distcorr = distcorr,
        dislist = dislist,
        ...
    )

    out$call <- match.call()
    class(out) <- "sqpad"
    out
}


sqpad <- function(
    formula,
    data,
    dis,
    dur,
    type = c("full", "approx"),
    det = c("joint", "pq"),
    Nmax = NULL,
    K = NULL,
    A = NULL,
    montecarlo = FALSE,
    np = 1000,
    distcorr = 2 / 3,
    ...
) {
    type <- match.arg(type)
    det <- match.arg(det)

    ## parsing the formula
    if (missing(data)) {
        data <- NULL
    }
    d <- svisitFormula(formula, data, n = 0) ## global environment
    Y <- d$y
    X <- d$x$sta
    Z <- d$x$det
    Xlevels <- d$levels

    ## check variables
    if (length(Y) < 1) {
        stop("empty model")
    }
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001))))) {
        stop("invalid dependent variable, non-integer values")
    }
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0)) {
        stop("invalid dependent variable, negative counts")
    }
    if (length(Y) != NROW(X)) {
        stop("invalid dependent variable, not a vector")
    }

    if (!is.numeric(dur)) {
        stop("dur must be numeric")
    }
    if (any(dur <= 0)) {
        stop("dur must be positive")
    }
    if (any(is.na(dur))) {
        stop("dur should not have missing (NA) values")
    }
    if (length(dur) < length(Y)) {
        stop("dur must be a vector with the same length as the response")
    }

    if (!is.numeric(dis)) {
        stop("dis must be numeric")
    }
    if (any(dis <= 0)) {
        stop("dis must be positive")
    }
    if (any(is.na(dis))) {
        stop("dis should not have missing (NA) values")
    }
    if (length(dis) < length(Y)) {
        stop("dis must be a vector with the same length as the response")
    }
    if (!is.null(A)) {
        if (length(A) != length(dis)) {
            stop("Length of A must equal length of dis")
        }
        if (any(A <= 0)) {
            stop("A must be positive")
        }
    }

    ## fit
    fit <- sqpad.fit(
        Y = Y,
        X = X,
        Z = Z,
        dis = dis,
        dur = dur,
        type = type,
        det = det,
        Nmax = Nmax,
        K = K,
        A = A,
        montecarlo = montecarlo,
        np = np,
        distcorr = distcorr,
        ...
    )
    fit$call = match.call()
    fit$formula <- d$formula
    fit$tt1 <- NULL
    class(fit) <- c("sqpad")
    return(fit)
}

print.sqpad <- function(x, digits, ...) {
    if (missing(digits)) {
        digits <- max(3, getOption("digits") - 3)
    }
    cat(
        "\nCall:",
        deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
        "",
        sep = "\n"
    )
    cat(
        "Single-bin QPAD model with ",
        if (x$det == "pq") "independence" else "dependence",
        "\n",
        if (x$type == "full") "Full" else "Approximate",
        " Maximum Likelihood estimates\n\n",
        "Coefficients:\n",
        sep = ""
    )
    print.default(
        format(x$coefficients, digits = digits),
        print.gap = 2,
        quote = FALSE
    )
    cat("\n")
    invisible(x)
}

vcov.sqpad <- function(object, ...) {
    object$vcov
}

fitted.sqpad <- function(object, ...) {
    X <- object$X
    poisson("log")$linkinv(drop(X %*% coef(object)[seq_len(ncol(X))]))
}

logLik.sqpad <- function(object, ...) {
    structure(
        object$loglik,
        df = object$df,
        nobs = object$nobs,
        class = "logLik"
    )
}

summary.sqpad <- function(object, ...) {
    coefs <- object$summary
    out <- list(
        call = object$call,
        coefficients = coefs,
        loglik = object$loglik,
        fitted.values = object$fitted.values,
        bic = BIC(object),
        type = object$type,
        det = object$det
    )
    class(out) <- "summary.sqpad"
    out
}

print.summary.sqpad <- function(x, digits, ...) {
    if (missing(digits)) {
        digits <- max(3, getOption("digits") - 3)
    }

    cat(
        "\nCall:",
        deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
        "",
        sep = "\n"
    )
    cat(
        "Single-bin QPAD model with ",
        if (x$det == "pq") "independence" else "dependence",
        "\n",
        if (x$type == "full") "Full" else "Approximate",
        " Maximum Likelihood estimates\n\n",
        "Coefficients:\n",
        sep = ""
    )
    printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
    if (!any(is.na(array(x$coefficients)))) {
        if (getOption("show.signif.stars") & any(x$coefficients[, 4] < 0.1)) {
            cat(
                "---\nSignif. codes: ",
                "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                "\n"
            )
        }
    }
    cat(
        "\nLog-likelihood:",
        formatC(x$loglik, digits = digits),
        "\nBIC =",
        formatC(x$bic, digits = digits),
        "\n"
    )
    cat("\n")
    invisible(x)
}

predict.sqpad <- function(
    object,
    newdata = NULL,
    type = c("link", "response"),
    ...
) {
    type <- match.arg(type)
    if (is.null(newdata)) {
        rval <- fitted(object)
        if (type == "link") {
            rval <- poisson("log")$linkfun(rval)
        }
    } else {
        X <- model.matrix(formula(object$formula$sta, lhs = FALSE), newdata)
        rval <- drop(X %*% coef(object)[seq_len(ncol(X))])
        if (type == "response") {
            rval <- poisson("log")$linkinv(rval)
        }
    }
    rval
}
