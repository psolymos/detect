## response to observer

## visualize the diffs
avoidDiff <-
function(r=seq(0,3,len=100), sigma=1, param=1,
type=c("exp","norm", "thres"), int=TRUE)
{
    if (is.function(type)) {
        hfun <- type
        type <- "custom"
    } else {
        type <- match.arg(type)
        hfun <- switch(type,
            "exp"   = function(r, param) 1-exp(-r/param),
            "norm"  = function(r, param) 1-exp(-(r/param)^2),
            "thres" = function(r, param) ifelse(r < param,0,1))
    }
#    require(MASS)
    ## this calculates the product c(r)*g(r)*h(r) for the integral
    pfun <- function(r, sigma, param, standard=FALSE) {
        aa <- exp(-(r/sigma)^2)
        if (standard)
            aa else aa * hfun(r, param)
    }
    tmpfun <- function(r, sigma, param, standard=FALSE) {
        aa <- pi*2*r * exp(-(r/sigma)^2)
        if (standard)
            aa else aa * hfun(r, param)
    }
#    require(MASS)
    ## function for the integral
    intfun1 <- function(r, ...) {
        suppressWarnings(MASS:::area(f=tmpfun, a=0, b=r, ...))
    }
    intfun2 <- function(r, ...) {
        integrate(f=tmpfun, lower=0, upper=r, ...)$value
    }
    intfun <- if (int) intfun2 else intfun1
    a0  <-   pfun(r, sigma, param, standard=TRUE)
    p0  <- tmpfun(r, sigma, param, standard=TRUE)
    cp0 <- sapply(r, function(z) intfun(z, sigma, param, standard=TRUE))
    a   <-   pfun(r, sigma, param, standard=FALSE)
    p   <- tmpfun(r, sigma, param, standard=FALSE)
    cp  <- sapply(r, function(z) intfun(z, sigma, param, standard=FALSE))
    op <- par(mfrow=c(1,3))
    plot(r,a0,ylim=c(0,max(a,a0)),
        xlab="r", ylab="circ prob", type="l", col=4, lty=2)
    lines(r,a,col=2, lty=1)
    plot(r,p0,ylim=c(0,max(p,p0)),
        xlab="r", ylab="circ prob", type="l", col=4, lty=2)
    lines(r,p,col=2, lty=1)
    plot(r,cp0,ylim=c(0,max(cp,cp0)),
        xlab="r", ylab="int dens", type="l", col=4, lty=2)
    lines(r,cp,col=2, lty=1)
    par(op)
    invisible(cbind(r=r,p0=p0,cp0=cp0,p=p,cp=cp))
}
if (FALSE) {
r <- seq(0,3,len=3000)
avoidDiff(r, sigma=1, param=1, type=function(r,param) 1)
avoidDiff(r, sigma=1, param=0.5, type="thres")
avoidDiff(r, sigma=1, param=0.5, type="exp")
avoidDiff(r, sigma=1, param=0.5, type="norm")
}

## type can be a custom function with formals (r, param)
## threshold model is unreasonable without another parameter (nu)
##    and should be restricted to 1st band that might vary
##    -- not generalizable !!!!!!!!!!!!!!
cmultiAvoid.fit <-
function(Y, D, X=NULL, type=c("exp","norm"),
inits=NULL, method="Nelder-Mead", ...)
{
    Ysum <- rowSums(Y, na.rm=TRUE)
    Y <- Y[Ysum > 0,,drop=FALSE]
    D <- D[Ysum > 0,,drop=FALSE]
    if (!is.null(X))
        X <- X[Ysum > 0,,drop=FALSE]
    Ysum <- Ysum[Ysum > 0]

    if (any(D[!is.na(D)] == Inf))
        stop("unlimited distance not supported in D")
    nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
        type <- match.arg(type)
    hfun <- switch(type,
        "exp"   = function(r, param) 1-exp(-r/param),
        "norm"  = function(r, param) 1-exp(-(r/param)^2))
    ## this calculates the product c(r)*g(r)*h(r) for the integral
    tmpfun <- function(r, sigma, param) {
        pi*2*r * exp(-(r/sigma)^2) * hfun(r, param)
    }
#    require(MASS)
    ## function for the integral
    intfun <- function(r, ...) {
        integrate(f=tmpfun, lower=0, upper=r, ...)$value
#        suppressWarnings(MASS:::area(f=tmpfun, a=0, b=r, ...))
    }
#    Ysum <- rowSums(Y, na.rm=TRUE)
    Yok <- !is.na(Y)
    n <- nrow(Y)
    k <- ncol(Y)
    if (is.null(inits)) {
        v0 <- cmulti.fit(Y, D, X, "dis")$coef
        ## param is on log scale, with no avoidance effect
        inits <- c(v0, -7)
    }
    R <- row(D)

    pifun0 <- function(r, sigma, param) {
        rr <- r
        dim(rr) <- NULL
        ii <- which(!duplicated(rr))
        vals <- sapply(ii, function(i)
            if (is.na(r[i]))
              NA else intfun(r=r[i], sigma=sigma, param=param))
        out <- vals[match(rr, rr[ii])]
        dim(out) <- dim(r)
        out
    }
    pifun1 <- function(r, sigma, param) {
        out <- sapply(1:length(r), function(i)
            if (is.na(r[i]))
                NA else intfun(r=r[i], sigma=sigma[R[i]], param=param))
        dim(out) <- dim(r)
        out
    }
    ## parameter is fixed, w is known
    if (is.null(X)) {
        pifun <- pifun0
        nll.fun <- function(param) {
            sig <- poisson("log")$linkinv(param[1])
            theta <- poisson("log")$linkinv(param[2])
            CP <- pifun(D, sig, theta)
            P <- CP - cbind(0, CP[, -k, drop=FALSE])
            Psum <- rowSums(P, na.rm=TRUE)
            PPsum <- P / Psum
            nll <- -sum(sapply(1:n, function(i) {
                logdmultinom(Y[i,Yok[i,]], Ysum[i], PPsum[i,Yok[i,]])
            }))
            if (nll %in% c(NA, NaN, Inf, -Inf))
                nlimit[2] else nll
        }
    } else {
        pifun <- pifun1
        nll.fun <- function(param) {
            sig <- poisson("log")$linkinv(drop(Xt %*% param[-length(param)]))
            theta <- poisson("log")$linkinv(param[length(param)])
            CP <- pifun(D, sig, theta)
            P <- CP - cbind(0, CP[, -k, drop=FALSE])
            Psum <- rowSums(P, na.rm=TRUE)
            PPsum <- P / Psum
            nll <- -sum(sapply(1:n, function(i) {
                logdmultinom(Y[i,Yok[i,]], Ysum[i], PPsum[i,Yok[i,]])
            }))
            if (nll %in% c(NA, NaN, Inf, -Inf))
                nlimit[2] else nll
        }
    }
    res <- optim(inits, nll.fun, method=method, hessian=TRUE, ...)
    rval <- list(coefficients=unname(res$par),
        vcov=try(solve(res$hessian)),
        loglik=-res$value)
    if (inherits(rval$vcov, "try-error"))
        rval$vcov <- matrix(NA, length(rval$coef), length(rval$coef))
    rval$vcov <- unname(rval$vcov)
    rval
}

simfunAvoid <-
function(n = 10, sigma=0.8, theta=NULL, type=c("exp","norm"),
D, lam=10, Nfixed=FALSE)
{
    type <- match.arg(type)
    hfun <- switch(type,
        "exp"   = function(r, param) 1-exp(-r/param),
        "norm"  = function(r, param) 1-exp(-(r/param)^2))
    tmpfun <- function(r, sigma, param) {
        pi*2*r * exp(-(r/sigma)^2) * hfun(r, param)
    }
    ## function for the integral
    intfun <- function(r, ...) {
        integrate(f=tmpfun, lower=0, upper=r, ...)$value
#        suppressWarnings(MASS:::area(f=tmpfun, a=0, b=r, ...))
    }
#    require(MASS)

    pifun0 <- function(r, sigma, param) {
        rr <- r
        dim(rr) <- NULL
        ii <- which(!duplicated(rr))
        vals <- sapply(ii, function(i)
            if (is.na(r[i]))
              NA else intfun(r=r[i], sigma=sigma, param=param))
        out <- vals[match(rr, rr[ii])]
        dim(out) <- dim(r)
        out
    }
    pifun1 <- function(r, sigma, param) {
        out <- sapply(1:length(r), function(i)
            if (is.na(r[i]))
                NA else intfun(r=r[i], sigma=sigma[R[i]], param=param))
        dim(out) <- dim(r)
        out
    }
    pifun <- if (length(sigma) == 1)
        pifun0 else pifun1
#    sigma <- rep(sigma, n)[1:n]
    if (missing(D)) {
        Dparts <- matrix(c(0.5, 1, NA,
                      0.5, 1, NA,
#                      0.5, 1, Inf,
                      0.25, 0.5, 1), 3, 3, byrow=TRUE)
        D <- Dparts[sample.int(3, n, replace=TRUE),]
    }
    R <- row(D)
    k <- ncol(D)

    CP <- pifun(D, sigma, theta)
    P <- CP - cbind(0, CP[, -k, drop=FALSE])
#    P[P<0] <- 0
    ## integral is not exact, small negative P might happen
    ## a rescaling seems more appropriate than hard thresholding at 0
#    P <- (P - min(P, na.rm=TRUE)) / (max(P, na.rm=TRUE) - min(P, na.rm=TRUE))
    Psum <- rowSums(P, na.rm=TRUE)
    PPsum <- P / Psum
    Pok <- !is.na(PPsum)

    lam <- rep(lam, n)[1:n] ## lam can be a vector, recycled
    N <- if (Nfixed)
        rep(round(lam), n) else rpois(n, lam)
    Y <- matrix(NA, ncol(PPsum), nrow(PPsum))

    Ypre <- sapply(1:n, function(i)
        rmultinom(1, N, PPsum[i,Pok[i,]]))
    Y[t(Pok)] <- unlist(Ypre)
    Y <- t(Y)
    list(Y=Y, D=D)
}

if (FALSE) {
n <- 1000

edr <- exp(-0.2)
D <- matrix(c(0.25, 0.5, 1, 1.5), n, 4, byrow=TRUE)

theta <- 0.001 # ----------- no effect

cmulti.fit <- detect::cmulti.fit

## CDS with exp
#avoidDiff(seq(0,1.5,len=200), sigma=edr, param=theta, type="exp")
vv <- simfunAvoid(n=n, sigma=edr, theta=theta, type="exp", D=D)
colMeans(vv$Y)
m0 <- cmulti.fit(vv$Y, vv$D, NULL, "dis")
m1 <- cmultiAvoid.fit(vv$Y, vv$D, NULL, type="exp")
cbind(true=c(log.edr=log(edr), log.theta=log(theta)), cds=c(m0$coef, NA), cdsAv=m1$coef)

## CDS with norm
#avoidDiff(seq(0,1.5,len=200), sigma=edr, param=theta, type="norm")
vv <- simfunAvoid(n=n, sigma=edr, theta=theta, type="norm", D=D)
colMeans(vv$Y)
m0 <- cmulti.fit(vv$Y, vv$D, NULL, "dis")
m1 <- cmultiAvoid.fit(vv$Y, vv$D, NULL, type="exp")
cbind(true=c(log.edr=log(edr), log.theta=log(theta)), cds=c(m0$coef, NA), cdsAv=m1$coef)

theta <- 0.25 # ----------- with effect

## CDS with exp
#avoidDiff(seq(0,1.5,len=200), sigma=edr, param=theta, type="exp")
vv <- simfunAvoid(n=n, sigma=edr, theta=theta, type="exp", D=D)
colMeans(vv$Y)
m0 <- cmulti.fit(vv$Y, vv$D, NULL, "dis")
m1 <- cmultiAvoid.fit(vv$Y, vv$D, NULL, type="exp")
cbind(true=c(log.edr=log(edr), log.theta=log(theta)), cds=c(m0$coef, NA), cdsAv=m1$coef)

## CDS with norm - might not be identifiable ??? (very close to CDS)
#avoidDiff(seq(0,1.5,len=200), sigma=edr, param=theta, type="norm")
vv <- simfunAvoid(n=n, sigma=edr, theta=theta, type="norm", D=D)
colMeans(vv$Y)
m0 <- cmulti.fit(vv$Y, vv$D, NULL, "dis")
m1 <- cmultiAvoid.fit(vv$Y, vv$D, NULL, type="norm")
cbind(true=c(log.edr=log(edr), log.theta=log(theta)), cds=c(m0$coef, NA), cdsAv=m1$coef)

}

