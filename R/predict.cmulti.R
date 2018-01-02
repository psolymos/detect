predict.cmulti <-
function(object, newdata = NULL, type = c("link", "response"), ...)
{
    if (object$type %in% c("mix", "fmix"))
        stop("predicted values are not available for finite mixture model")
    if (is.null(newdata)) {
        rval <- fitted(object)
        ## fitted gives 'response', need to use linkfun for 'link'
        if (type == "link")
            rval <- poisson("log")$linkfun(rval)
    } else {
        X <- model.matrix(formula(object$formula, lhs=c(FALSE, FALSE)), newdata)
        rval <- drop(X %*% coef(object))
        if (type == "response") {
            rval <- poisson("log")$linkinv(rval)
        }
    }
    rval
}
