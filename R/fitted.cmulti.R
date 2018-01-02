fitted.cmulti <-
function(object, ...){
    type <- object$type
    if (type %in% c("mix", "fmix"))
        stop("fitted values are not available for finite mixture model")
    X <- model.matrix(object)
    poisson("log")$linkinv(drop(X %*% coef(object)))
}
