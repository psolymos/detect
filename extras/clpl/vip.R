## this model uses V as inflation point, therefore the name VIP

set.seed(1)
n <- 100
lam <- 2 # poisson mean, can be a vector of length n
phi <- 0.4 # V-inflation probability, can be a vector of length n
V <- 2 # V is the count value, can be 0, 2, etc
y <- y0 <- rpois(n, lam)
a <- rbinom(n, 1, phi)
y[a > 0] <- V
table(Poisson=y0, Vinflated=y)

Y <- y
id0 <- Y == V
id1 <- !id0
X <- Z <- matrix(1, n, 1)
kx <- ncol(X)
kz <- ncol(Z)
offsetx <- offsetz <- 0
weights <- rep(1, length(Y))
linkinvx <- poisson("log")$linkinv
linkinvz <- binomial("logit")$linkinv
good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

nll_VIP_ML <- function(parms) {
    mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    loglik0 <- log(phi + (1 - phi) * dpois(Y, lambda = mu, log = FALSE))
    loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
    loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
    if (!is.finite(loglik) || is.na(loglik))
        loglik <- -good.num.limit[2]
    -loglik
}

opt <- optim(c(0,0), nll_VIP_ML, hessian=TRUE, method="Nelder-Mead")

cbind(true=c(lam=lam, phi=phi), optim=c(lam=exp(opt$par[1]), phi=plogis(opt$par[2])))

## use covariates for the count part

set.seed(1)
n <- 100
x <- rnorm(n)
df <- data.frame(x=x)
X <- model.matrix(~x, df)
beta <- c(1.5,-1)
lam <- exp(X %*% beta) # poisson mean, can be a vector of length n
phi <- 0.4 # V-inflation probability, can be a vector of length n
V <- 2 # V is the count value, can be 0, 2, etc
y <- y0 <- rpois(n, lam)
a <- rbinom(n, 1, phi)
y[a > 0] <- V
table(Poisson=y0, Vinflated=y)

Y <- y
id0 <- Y == V
id1 <- !id0
Z <- matrix(1, n, 1)
kx <- ncol(X)
kz <- ncol(Z)
offsetx <- offsetz <- 0
weights <- rep(1, length(Y))
linkinvx <- poisson("log")$linkinv
linkinvz <- binomial("logit")$linkinv
good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

opt <- optim(rep(0, kx+kz), nll_VIP_ML, hessian=TRUE, method="Nelder-Mead")

cbind(true=c(beta=beta, phi=phi), optim=c(beta=opt$par[1:kx], phi=plogis(opt$par[kx+kz])))

## wrapping up

vip <-
function(Y, X, Z, V=0,
offsetx, offsetz, weights, linkz="logit",
conditional=FALSE, ...) {
    if (missing(Y))
        stop("C'mon, you must have some data?!")
    if (conditional && any(Y < 1))
        stop("Y must be >0 when conditional=TRUE")
    n <- length(Y)
    id0 <- Y == V
    id1 <- !id0
    if (missing(X)) {
        X <- matrix(1, n, 1)
        colnames(X) <- "(Intercept)"
    }
    if (missing(Z)) {
        Z <- matrix(1, n, 1)
        colnames(Z) <- "(Intercept)"
    }
    kx <- ncol(X)
    kz <- ncol(Z)
    if (missing(offsetx))
        offsetx <- 0
    if (missing(offsetz))
        offsetz <- 0
    if (missing(weights))
        weights <- rep(1, n)
    linkinvx <- poisson("log")$linkinv
    linkinvz <- binomial(linkz)$linkinv
    good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

    nll_VIP_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
        loglik0 <- log(phi + (1 - phi) * dpois(Y, lambda = mu, log = FALSE))
        loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
        if (conditional) {
            num <- ifelse(id0, loglik0, loglik1)
            den <- log(1 - (1-phi) * exp(-mu))
            loglik <- sum(weights * (num - den))
        } else {
            loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
        }
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    opt <- optim(rep(0, kx+kz), nll_VIP_ML, hessian=TRUE, method="Nelder-Mead")
    par <- opt$par
    names(par) <- c(paste0("P_", colnames(X)), paste0("V_", colnames(Z)))
    vc <- solve(opt$hessian)
    dimnames(vc) <- list(names(par), names(par))
    out <- list(call=match.call(),
        coefficients=par, loglik=-opt$value, vcov=vc, nobs=n,
        conditional=conditional)
    class(out) <- "vip"
    out
}
vcov.vip <- function(object, ...) object$vcov
logLik.vip <- function (object, ...)
    structure(object$loglik, df = object$nobs - length(object$coef),
        nobs = object$nobs, class = "logLik")
summary.vip <- function (object, ...) {
    k <- length(object$coefficients)
    coefs <- coef(object)
    se <- sqrt(diag(vcov(object)))
    tstat <- coefs/se
#    pval <- 2 * pt(abs(tstat), object$df.residual, lower.tail = FALSE)
    ## z test because no overdspersion
    pval <- 2 * pnorm(-abs(tstat))
    coefs <- cbind(coefs, se, tstat, pval)
#    colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    coefs <- coefs[1:k, , drop = FALSE]
    rownames(coefs) <- names(coef(object))
    out <- list(call = object$call, coefficients=coefs, loglik = object$loglik,
        bic=BIC(object), conditional=object$conditional)
    class(out) <- "summary.vip"
    return(out)
}
print.summary.vip <- function (x, digits, ...)
{
    if (missing(digits))
        digits <- max(3, getOption("digits") - 3)
    cat("\nCall:", deparse(x$call,
        width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    cat("V-Inflated", if (x$conditional) "(Zero Conditioned)" else "", "Poisson Model\n\n")
    cat(paste("Coefficients:\n", sep = ""))
    printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
    if (!any(is.na(array(x$coefficients)))) {
        if (getOption("show.signif.stars") & any(x$coefficients[,4] < 0.1))
            cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    }
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "\nBIC =", formatC(x$bic, digits = digits), "\n")
    cat("\n")
    invisible(x)
}
confint.vip <-
function (object, parm, level = 0.95, ...)
{
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) {
        parm <- pnames
    } else {
        if (is.numeric(parm))
            parm <- pnames[parm]
    }
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%", sep="")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    fac <- qnorm(a)
    ses <- sqrt(diag(vcov(object, model, type)))
    ci[] <- cf[parm] + ses[parm] %o% fac
    ci
}



mod <- vip(Y, X, V=2)
coef(mod)
vcov(mod)
summary(mod)
confint(mod)

nll_VIP_CML <- function(parms) {
    mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    loglik0 <- log(phi + (1 - phi) * dpois(Y, lambda = mu, log = FALSE))
    loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
    num <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
    den <- log(1 - (1-phi) * exp(-mu))
    loglik <- sum(weights1 * (num - den))
    if (!is.finite(loglik) || is.na(loglik))
        loglik <- -good.num.limit[2]
    -loglik
}


set.seed(1)
n <- 1000
x <- rnorm(n)
z <- runif(n, -1, 1)
df <- data.frame(x=x, z=z)
X <- model.matrix(~x, df)
Z <- model.matrix(~z, df)
beta <- c(-0.5, -0.5)
alpha <- c(0, 0.5)
lam <- exp(X %*% beta)
phi <- plogis(Z %in% alpha)
V <- 2 # V is the count value, cannot be 0
y <- y0 <- rpois(n, lam)
a <- rbinom(n, 1, phi)
y[a > 0] <- V
keep <- y0>0
y <- y[keep] # conditioning
y0 <- y0[keep]
X <- X[keep,]
Z <- Z[keep,]
table(Poisson=y0, Vinflated=y)

Y=y
kx <- ncol(X)
kz <- ncol(Z)
offsetx <- offsetz <- 0
weights <- rep(1, length(Y))
linkinvx <- poisson("log")$linkinv
linkinvz <- binomial("logit")$linkinv
good.num.limit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
    id0 <- Y == V
    id1 <- !id0

    nll_VIP_ML <- function(parms) {
        mu <- as.vector(linkinvx(X %*% parms[1:kx] + offsetx))
        phi <- as.vector(linkinvz(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
        loglik0 <- log(phi + (1 - phi) * dpois(Y, lambda = mu, log = FALSE))
        loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
        loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
        if (conditional) {
            num <- ifelse(id0, loglik0, loglik1)
            den <- log(1 - (1-phi) * exp(-mu))
            loglik <- sum(weights * (num - den))
        } else {
            loglik <- sum(weights[id0] * loglik0[id0]) + sum(weights[id1] * loglik1[id1])
        }
        if (!is.finite(loglik) || is.na(loglik))
            loglik <- -good.num.limit[2]
        -loglik
    }

    opt <- optim(rep(0, kx+kz), nll_VIP_ML, hessian=TRUE, method="Nelder-Mead")

mod <- vip(Y=y, X=X, Z=Z, V=2, conditional=TRUE)
summary(mod)
cbind(True=c(beta=beta, logit_phi=qlogis(phi)),
      coef(mod))
