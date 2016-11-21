model.matrix.cmulti <-
function (object, ...)
{
    if (!is.null(object$x))
        rval <- object$x
    else if (!is.null(object$model))
        rval <- model.matrix(object$formula, object$model)
    else stop("not enough information in fitted model to return model.matrix")
    return(rval)
}
