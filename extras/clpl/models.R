source("~/repos/clpl/R/clpl.R")
n <- 1000
R <- 100
overlap <- TRUE
distr <- "pois"

alpha <- c(1,0.5,-1)
beta <- c(-0.25, -0.7, 0.8)
theta <- 0.1

set.seed(1234)

x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rnorm(n)

#x1 <- rnorm(n, 0, 1)
#x2 <- runif(n, 0, 1)
#x3 <- rbinom(n, 1, 0.5)
#x4 <- runif(n, -1, 1)

XX <- model.matrix(~x1+x2+x3+x4)
ZZ <- XX

X <- model.matrix(~x1+x2)
TrueX <- "1100"
if (overlap) {
    Z <- model.matrix(~x3+x2)
    TrueZ <- "0110"
} else {
    Z <- model.matrix(~x3+x4)
    TrueZ <- "0011"
}
TrueXZ <- paste(TrueX, TrueZ, sep="-")
lambda <- exp(drop(X %*% alpha))
phi <- plogis(drop(Z %*% beta))

## define model set here and lookup matrices/distances

tmp <- lapply(1:4, function(k) combn(4, k))
mc1 <- matrix(0L, 16, 4)
r <- 1
for (i in 1:length(tmp)) {
    for (j in 1:ncol(tmp[[i]])) {
        r <- r + 1
        mc1[r, tmp[[i]][,j]] <- 1L
    }
}
colnames(mc1) <- paste0("x", 1:4)
rownames(mc1) <- apply(mc1, 1, function(z) paste0(z, collapse=""))

df <- expand.grid(z=rownames(mc1), x=rownames(mc1))
df <- df[,c("x","z")]
df$xz <- interaction(df$x, df$z, sep="-", drop=TRUE)
True <- which(df$x==TrueX & df$z==TrueZ)

mc2 <- cbind(mc1[as.character(df$x),], mc1[as.character(df$z),])
colnames(mc2) <- paste0(rep(c("x","z"), each=4), 1:4)
rownames(mc2) <- as.character(df$xz)

d1 <- matrix(NA, 16, 16)
dimnames(d1) <- list(rownames(mc1), rownames(mc1))
for (i in 1:16) {
    for (j in 1:16) {
        d1[i,j] <- sum(abs(mc1[i,] - mc1[j,]))
    }
}
d2 <- matrix(NA, 256, 256)
dimnames(d2) <- list(rownames(mc2), rownames(mc2))
for (i in 1:256) {
    for (j in 1:256) {
        d2[i,j] <- sum(abs(mc2[i,] - mc2[j,]))
    }
}


## loop for replicates /start
Run <- 1
allres <- list()
for (Run in 1:R) {
cat("R:", Run, "\t", date(), "\n");flush.console()

A <- rbinom(n, 1, 1-phi)
if (distr == "pois")
    Y <- rpois(n, lambda * A)
if (distr == "negbin")
    Y <- stats:::rnbinom(n, size=theta, mu=lambda*A)

## fit true model
#mTrue <- zi.fit(Y, X, Z, distr=distr, type=c("ML","CL","PL"))

## fit all models
out_ML <- list()
out_CL <- list()
out_CLPL <- list()
Zi <- ZZ[,1,drop=FALSE]
for (i in 1:16) {
    Xi <- XX[,c(TRUE,as.logical(mc1[i,])),drop=FALSE]
    out_CL[[rownames(mc1)[i]]] <- zi.fit(Y, Xi, Zi, distr=distr, type="CL")
    for (j in 1:16) {
#        cat("R:", Run, "\tCL:", i, "\tPL:", j, "\n");flush.console()
        Zi <- ZZ[,c(TRUE,as.logical(mc1[j,])),drop=FALSE]
        id <- paste(rownames(mc1)[i], rownames(mc1)[j], sep="-")
        out_ML[[id]] <- zi.fit(Y, Xi, Zi, distr=distr, type="ML")$ML
        out_CLPL[[id]] <- zi.fit(Y, Xi, Zi, distr=distr, type="PL", 
            fit=out_CL[[rownames(mc1)[i]]])$CLPL
    }
}

## identify path and best fit

aic <- function(x) -2*x$loglik + 2*length(x$coef)
#aic <- function(x) -2*x$loglik
#aic <- function(x, alpha=0.5) {
#    a <- -2*x$loglik + 2*length(x$coef)
#    b <- -2*x$loglik + log(n)*length(x$coef)
#    alpha * a + (1-alpha) * b
#}


res <- list()
res[["truth_ML"]] <- structure(aic(out_ML[[TrueXZ]]), names=TrueXZ)
res[["truth_CL"]] <- structure(aic(out_CL[[TrueX]]$CL), 
    names=paste(TrueX, "xxxx", sep="-"))
res[["truth_CLPL"]] <- structure(aic(out_CLPL[[TrueXZ]]), names=TrueXZ)

#### x-first forwards, ML

bestx <- "0000"
bestz <- "0000"
bestid <- paste(bestx, bestz, sep="-")
bestaic <- aic(out_ML[[bestid]])

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][(i+1):16]
    steps <- paste(names(steps)[steps == 1], bestz, sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestx <- substr(bestid, 1, 4)
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}

xzmods1 <- xmods
names(xzmods1) <- paste(names(xmods), bestz, sep="-")
zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][(i+1):16]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}
xzmods2 <- zmods
names(xzmods2) <- paste(bestx, names(zmods), sep="-")

res[["xfirst_forwards_ML"]] <- rev(sort(c(xzmods1, xzmods2)))

#### z-first forwards, ML

bestx <- "0000"
bestz <- "0000"
bestid <- paste(bestx, bestz, sep="-")
bestaic <- aic(out_ML[[bestid]])

zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][(i+1):16]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

xzmods1 <- zmods
names(xzmods1) <- paste(bestx, names(zmods), sep="-")

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][(i+1):16]
    steps <- paste(names(steps)[steps == 1], bestz, sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestx <- substr(bestid, 1, 4)
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}
xzmods2 <- xmods
names(xzmods2) <- paste(names(xmods), bestz, sep="-")

res[["zfirst_forwards_ML"]] <- rev(sort(c(xzmods1, xzmods2)))

#### x-z forwards, ML

bestxz <- "0000-0000"
bestaic <- aic(out_ML[[bestxz]])

xzmods <- structure(rep(NA, 256), names=rownames(d2))
xzmods[bestxz] <- bestaic

for (counter in 1:8) {
    i <- which(rownames(d2) == bestxz)
    steps <- d2[, bestxz][(i+1):256]
    steps <- names(steps)[steps == 1]
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestxz <- steps[which.min(daic)]
    bestaic <- unname(aicvals[which.min(daic)])
    xzmods[bestxz] <- bestaic
}

res[["xz_forwards_ML"]] <- rev(sort(xzmods))

#### x-first backwards, ML

bestx <- "1111"
bestz <- "1111"
bestid <- paste(bestx, bestz, sep="-")
bestaic <- aic(out_ML[[bestid]])

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][1:(i-1)]
    steps <- paste(names(steps)[steps == 1], bestz, sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestx <- substr(bestid, 1, 4)
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}

xzmods1 <- xmods
names(xzmods1) <- paste(names(xmods), bestz, sep="-")
zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][1:(i-1)]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}
xzmods2 <- zmods
names(xzmods2) <- paste(bestx, names(zmods), sep="-")

res[["xfirst_backwards_ML"]] <- rev(sort(c(xzmods1, xzmods2)))

#### z-first backwards, ML

bestx <- "1111"
bestz <- "1111"
bestid <- paste(bestx, bestz, sep="-")
bestaic <- aic(out_ML[[bestid]])

zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][1:(i-1)]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

xzmods1 <- zmods
names(xzmods1) <- paste(bestx, names(zmods), sep="-")

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][1:(i-1)]
    steps <- paste(names(steps)[steps == 1], bestz, sep="-")
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestx <- substr(bestid, 1, 4)
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}
xzmods2 <- xmods
names(xzmods2) <- paste(names(xmods), bestz, sep="-")

res[["zfirst_backwards_ML"]] <- rev(sort(c(xzmods1, xzmods2)))

#### x-z backwards, ML

bestxz <- "1111-1111"
bestaic <- aic(out_ML[[bestxz]])

xzmods <- structure(rep(NA, 256), names=rownames(d2))
xzmods[bestxz] <- bestaic

for (counter in 1:8) {
    i <- which(rownames(d2) == bestxz)
    steps <- d2[, bestxz][1:(i-1)]
    steps <- names(steps)[steps == 1]
    aicvals <- sapply(out_ML[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestxz <- steps[which.min(daic)]
    bestaic <- unname(aicvals[which.min(daic)])
    xzmods[bestxz] <- bestaic
}

res[["xz_backwards_ML"]] <- rev(sort(xzmods))

#### x-first forwards, CL

bestx <- "0000"
bestaic <- aic(out_CL[[bestx]]$CL)

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][(i+1):16]
    steps <- names(steps)[steps == 1]
    aicvals <- sapply(out_CL[steps], function(z) aic(z$CL))
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestx <- steps[which.min(daic)]
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}

xzmods1 <- xmods
names(xzmods1) <- paste(names(xmods), "xxxx", sep="-")

bestz <- "0000"
bestaic <- aic(out_CLPL[[paste(bestx, bestz, sep="-")]])
zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][(i+1):16]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_CLPL[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}
xzmods2 <- zmods
names(xzmods2) <- paste(bestx, names(zmods), sep="-")

res[["xfirst_forwards_CL"]] <- c(rev(sort(xzmods1)), rev(sort(xzmods2)))

#### x-first forwards, CL

bestx <- "1111"
bestaic <- aic(out_CL[[bestx]]$CL)

xmods <- structure(rep(NA, 16), names=rownames(d1))
xmods[bestx] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestx)
    steps <- d1[, bestx][1:(i-1)]
    steps <- names(steps)[steps == 1]
    aicvals <- sapply(out_CL[steps], function(z) aic(z$CL))
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestx <- steps[which.min(daic)]
    bestaic <- unname(aicvals[which.min(daic)])
    xmods[bestx] <- bestaic
}

xzmods1 <- xmods
names(xzmods1) <- paste(names(xmods), "xxxx", sep="-")

bestz <- "1111"
bestaic <- aic(out_CLPL[[paste(bestx, bestz, sep="-")]])
zmods <- structure(rep(NA, 16), names=rownames(d1))
zmods[bestz] <- bestaic

for (counter in 1:4) {
    i <- which(rownames(d1) == bestz)
    steps <- d1[, bestz][1:(i-1)]
    steps <- paste(bestx, names(steps)[steps == 1], sep="-")
    aicvals <- sapply(out_CLPL[steps], aic)
    daic <- aicvals - bestaic
    delta <- min(daic)
    if (delta >= 0)
        break
    bestid <- steps[which.min(daic)]
    bestz <- substr(bestid, 6, 9)
    bestaic <- unname(aicvals[which.min(daic)])
    zmods[bestz] <- bestaic
}
xzmods2 <- zmods
names(xzmods2) <- paste(bestx, names(zmods), sep="-")

res[["xfirst_backwards_CL"]] <- c(rev(sort(xzmods1)), rev(sort(xzmods2)))

allres[[Run]] <- res
}

save.image("c:/Dropbox/pkg/zi-clpl/zi2.Rdata")

## ----------------------------------------------------------------

load("c:/Dropbox/pkg/zi-clpl/zi2.Rdata")

xx <- t(sapply(allres, function(zz) sapply(zz, function(z) names(z)[length(z)])))
xx <- xx[,-c(2,3)]
#xxx <- apply(xx, 2, function(id) data.frame(100*table(id)/R))
id <- sort(unique(as.character(xx)))
fr <- matrix(0, length(id), ncol(xx))
colnames(fr) <- colnames(xx)
rownames(fr) <- id
for (i in 1:ncol(xx)) {
    tmp <- table(xx[,i])
    fr[names(tmp),i] <- tmp
}
rr <- data.frame(freq=fr[TrueXZ,]/R, nmodels=colSums(fr>0), 
    #H=apply(fr/100, 2, function(z) sum(z^2)),
    best=apply(fr, 2, function(z) ifelse(rownames(fr)[which.max(z)]=="1100-0110",1,0)))
rr
plot(rr, pch=19, col=c(1, rep(2,6), 4, 4))

## timings based on 1 run

tml <- colMeans(do.call(rbind, lapply(out_ML, "[[", "time")))
tcl <- colMeans(do.call(rbind, lapply(out_CLPL, "[[", "time")))

tcl / tml

res
summary(sapply(out_ML, aic) - aic(out_ML[[TrueXZ]]))
summary(sapply(out_CLPL, aic) - aic(out_CLPL[[TrueXZ]]))

data.frame(best=sapply(res, function(z) names(z)[length(z)]))


## package results
## save
