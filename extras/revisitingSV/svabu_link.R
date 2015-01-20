svabu_link <- 
function (formula, data, zeroinfl = TRUE, area = 1, N.max = NULL, 
    inits, link.det = "logit", link.zif = "logit", model = TRUE, 
    x = FALSE, cm_ini=FALSE, ...) 
{
    if (missing(data)) 
        data <- NULL
    out <- detect:::svisitFormula(formula, data, n = 0)
    out$call = match.call()
    Y <- out$y
    X <- out$x$sta
    Z <- out$x$det
    Q <- out$x$zif
    Xlevels <- out$levels
    if (!zeroinfl && !is.null(Q)) {
        warning("'zeroinfl = FALSE': zero inflation part in formula ignored")
        Q <- NULL
    }
    if (length(Y) < 1) 
        stop("empty model")
    if (all(Y > 0) && zeroinfl) 
        warning("dependent variable has no 0 value, zero-inflated model might not be adequate")
    if (!any(Y > 1)) 
        warning("dependent variable does not contain any count > 1")
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 
        0.001))))) 
        stop("invalid dependent variable, non-integer values")
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0)) 
        stop("invalid dependent variable, negative counts")
    if (length(Y) != NROW(X)) 
        stop("invalid dependent variable, not a vector")
    if (setequal(colnames(Z), colnames(X))) 
        stop("at least one covariate should be separate for occupancy and detection parts of the formula")
    if (all(union(colnames(X), colnames(Z))[-1] %in% names(unlist(Xlevels)))) 
        stop("model must include at least one continuous covariate")
    links <- c("logit", "probit", "cloglog")

    ## fair warning
    if (!is.function(link.det)) {
        link.det <- match.arg(link.det, links)
    } else {
        warning("!!! unknown link function used, beware of the consequences when using methods !!!")
    }
    if (!is.function(link.zif)) {
        link.zif <- match.arg(link.zif, links)
    } else {
        warning("!!! unknown link function used, beware of the consequences when using methods !!!")
    }

    if (!(length(area) %in% c(1, length(Y)))) 
        stop("invalid length for 'area' argument")
    if (any(area == 0)) 
        stop("'area' cannot be 0")
    fit <- svabu.fit_link(Y, X, Z, Q, zeroinfl = zeroinfl, area, N.max, 
        inits, link.det, link.zif, cm_ini, ...)
    sclass <- "svabu_p"
    out <- c(fit, out)
    if (!model) 
        out$model <- NULL
    if (!x) 
        out$x <- NULL
    class(out) <- c(sclass, "svabu", "svisit")
    return(out)
}


svabu.fit_link <-
function (Y, X, Z, Q = NULL, zeroinfl = TRUE, area = 1, N.max = NULL, 
    inits, link.det = "logit", link.zif = "logit", cm_ini=FALSE, ...) 
{
    nll.P <- function(parms) {
        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each = N.max + 
            1) * area.rep
        delta.rep <- rep(drop(linkinvfun.det(Z %*% parms[(np.abu + 
            1):(np.abu + np.det)])), each = N.max + 1)
        dp <- dpois(N.rep, lambda.rep)
        db <- dbinom(Y.rep, N.rep, delta.rep)
        intsum <- colSums(matrix(dp * db, nrow = N.max + 1))
        out <- -sum(log(intsum))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.ZIP1 <- function(parms) {
        lambda1 <- drop(exp(X1 %*% parms[1:np.abu]))
        delta1 <- drop(linkinvfun.det(Z1 %*% parms[(np.abu + 
            1):(np.abu + np.det)]))
        lambda1.rep <- rep(lambda1, each = N.max + 1) * area1.rep
        delta1.rep <- rep(delta1, each = N.max + 1)
        dp1 <- dpois(N1.rep, lambda1.rep)
        db1 <- dbinom(Y1.rep, N1.rep, delta1.rep)
        part1 <- colSums(matrix(dp1 * db1, nrow = N.max + 1))
        part2 <- 1/(1 - exp(-lambda1 * delta1))
        out <- -sum(log(part1 * part2))
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
        lambda.rep <- rep(drop(exp(X %*% parms[1:np.abu])), each = N.max + 
            1) * area.rep
        delta.rep <- rep(drop(linkinvfun.det(Z %*% parms[(np.abu + 
            1):(np.abu + np.det)])), each = N.max + 1)
        phi.rep <- rep(drop(linkinvfun.zif(Q %*% parms[(np.abu + 
            np.det + 1):length(parms)])), each = N.max + 1)
        loglik0 <- log(phi.rep + exp(log(1 - phi.rep) - lambda.rep))
        loglik1 <- log(1 - phi.rep) + dpois(N.rep, lambda = lambda.rep, 
            log = TRUE)
        dp <- exp(ifelse(id1.rep, loglik1, loglik0))
        db <- dbinom(Y.rep, N.rep, delta.rep)
        intsum <- colSums(matrix(dp * db, nrow = N.max + 1))
        out <- -sum(log(intsum))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }

    ## Marginal Conditional for inits
    nll.mZIP1 <- function(parms) {
        lambda1 <- drop(exp(X1 %*% parms[1:np.abu]))
        delta1 <- drop(linkinvfun.det(Z1 %*% parms[(np.abu + 
            1):(np.abu + np.det)]))
        lamp <- lambda1 * area1 * delta1
        out <- -sum(dpois(Y1, lamp, log=TRUE) - log(1 - dpois(0, lamp)))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }
    nll.mP <- function(parms) {
        lambda <- drop(exp(X %*% parms[1:np.abu]))
        delta <- drop(linkinvfun.det(Z %*% parms[(np.abu + 
            1):(np.abu + np.det)]))
        lamp <- lambda * area * delta
        out <- -sum(dpois(Y, lamp, log=TRUE))
        out <- ifelse(abs(out) == Inf, good.num.limit[2], out)
        out <- ifelse(is.na(out), good.num.limit[2], out)
        out
    }

    control <- getOption("detect.optim.control")
    control$fnscale <- 1
    method <- getOption("detect.optim.method")
    if (is.null(method)) 
        method <- "Nelder-Mead"
    Control <- list(optim.control = control, optim.method = method)
    good.num.limit = c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
    n <- length(Y)
    nam.abu <- colnames(X)
    nam.det <- colnames(Z)
    np.abu <- NCOL(X)
    np.det <- NCOL(Z)
    if (zeroinfl && is.null(Q)) {
        Q <- matrix(1, n, 1)
        colnames(Q) <- "(Intercept)"
    }
    if (zeroinfl) {
        np.zif <- NCOL(Q)
        nam.zif <- colnames(Q)
    }
    else {
        np.zif <- 0
    }
    np <- np.abu + np.det + np.zif

    ## simple as this
    if (is.function(link.det)) {
        linkinvfun.det <- link.det
    } else {
        linkinvfun.det <- binomial(link = make.link(link.det))$linkinv
    }
    if (is.function(link.det)) {
        linkinvfun.zif <- link.zif
    } else {
        linkinvfun.zif <- binomial(link = make.link(link.zif))$linkinv
    }

    user_ini <- TRUE
    if (missing(inits)) {
        user_ini <- FALSE
        inits <- if (zeroinfl) {
            c(glm.fit(X, Y, family = poisson())$coef, glm.fit(Z, 
                Y, family = poisson())$coef, glm.fit(Q, ifelse(Y > 
                0, 0, 1), family = binomial())$coef)
        }
        else {
            c(glm.fit(X, Y, family = poisson())$coef, glm.fit(Z, 
                Y, family = poisson())$coef)
        }
        inits[is.na(inits)] <- 0
    }
    inits <- rep(inits, np)[1:np]
    if (is.null(N.max)) 
        N.max <- max(max(Y, na.rm = TRUE) + 20, max(Y, na.rm = TRUE) * 
            3)
    if (N.max <= max(Y, na.rm = TRUE)) 
        stop("N.max too small")
    N.rep <- rep(0:N.max, n)
    Y.rep <- rep(Y, each = N.max + 1)
    area.rep <- if (identical(area, 1)) 
        area
    else rep(area, each = N.max + 1)
    if (zeroinfl) {
        id1 <- Y > 0
        id1.rep <- rep(id1, each = N.max + 1)
        Y1 <- Y[id1]
        n1 <- length(Y1)
        N1.rep <- rep(0:N.max, n1)
        Y1.rep <- rep(Y1, each = N.max + 1)
        X1 <- data.matrix(X[id1, ])
        Z1 <- data.matrix(Z[id1, ])
        if (identical(area, 1)) {
            area1 <- area
            area1.rep <- area
        }
        else {
            area1 <- area[id1]
            area1.rep <- rep(area1, each = N.max + 1)
        }
        if (!user_ini && cm_ini) {
            ini00 <- try(optim(inits[1:(np.abu + np.det)], nll.mZIP1, 
                method = method, hessian = FALSE, control = control))
            if (!inherits(ini00, "try-error"))
                inits[1:(np.abu + np.det)] <- ini00$par
        }
        results <- optim(inits[1:(np.abu + np.det)], nll.ZIP1, 
            method = method, hessian = TRUE, control = control)
    } else {
        id1 <- Y >= 0
        id1.rep <- rep(id1, each = N.max + 1)
        if (!user_ini && cm_ini) {
            ini00 <- try(optim(inits[1:(np.abu + np.det)], nll.mP, 
                method = method, hessian = FALSE, control = control))
            if (!inherits(ini00, "try-error"))
                inits[1:(np.abu + np.det)] <- ini00$par
        }
        results <- optim(inits[1:(np.abu + np.det)], nll.P, method = method, 
            hessian = TRUE, control = control)
    }
    estimates <- results$par
    par.state <- estimates[1:np.abu]
    par.det <- estimates[(np.abu + 1):(np.abu + np.det)]
    names(par.state) <- nam.abu
    names(par.det) <- nam.det
    if (rcond(results$hessian) <= 1e-06) 
        std.error <- rep(NA, np)
    if (rcond(results$hessian) > 1e-06) {
        opvar <- diag(solve(results$hessian))
        if (any(opvar < 0)) {
            opvar[opvar < 0] <- NA
            warning("negative variance values in optim, NAs produced")
        }
        std.error <- sqrt(opvar)
    }
    se.state <- std.error[1:np.abu]
    se.det <- std.error[(np.abu + 1):(np.abu + np.det)]
    names(se.state) <- nam.abu
    names(se.det) <- nam.det
    lambda <- drop(exp(X %*% par.state))
    delta <- drop(linkinvfun.det(Z %*% par.det))
    zif.results <- NULL
    par.zif <- if (zeroinfl) 
        0
    else NULL
    se.zif <- NULL
    phi <- NULL
    lLik <- -results$value
    if (sum(!id1) > 0 && zeroinfl) {
        emlp <- exp(-lambda * delta)
        zif.results <- suppressWarnings(optim(inits[(np.abu + 
            np.det + 1):np], nll.ZIP0, emlp = emlp, id1 = id1, 
            method = method, hessian = TRUE, control = control))
        par.zif <- zif.results$par
        names(par.zif) <- nam.zif
        phi <- drop(linkinvfun.zif(Q %*% par.zif))
        lLik <- -nll.ZIP(c(par.state, par.det, par.zif))
        if (rcond(zif.results$hessian) <= 1e-06) 
            se.zif <- rep(NA, np.zif)
        if (rcond(zif.results$hessian) > 1e-06) {
            opvar2 <- diag(solve(data.matrix(zif.results$hessian)))
            if (any(opvar2 < 0)) {
                opvar2[opvar2 < 0] <- NA
                warning("negative variance values in optim, NAs produced")
            }
            se.zif <- sqrt(opvar2)
        }
        names(se.zif) <- nam.zif
    }
    Converged <- if (zeroinfl) {
        results$convergence == 0 && zif.results$convergence == 
            0
    }
    else results$convergence == 0
    out <- list(coefficients = list(sta = par.state, det = par.det), 
        std.error = list(sta = se.state, det = se.det), fitted.values = lambda, 
        detection.probabilities = delta, zif.probabilities = phi, 
        zeroinfl = zeroinfl, nobs = n, N.max = N.max, link = list(sta = "log", 
            det = link.det, zif = link.zif), df.null = n - 2, 
        df.residual = n - np, inits = inits, loglik = lLik, results = list(count = results, 
            zero = zif.results), converged = Converged, control = Control, 
        area = area)
    if (zeroinfl) {
        out$coefficients$zif <- par.zif
        out$std.error$zif <- se.zif
    }
    if (!out$converged) 
        warning("model did not converge")
    return(out)
}
