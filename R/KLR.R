#' @name KLR
#' 
#' @title Kernel Logistic Regression
#' 
#' @description This function fit a kernel logistic regression model to the data (\code{y}, \code{x}) using some pre-specified kernel. The return list contains the estimated kernel weights as well as the original data to perform predictions.
#'
#' @param y A \code{n x 1} column vector containing the responses (0-1).
#' @param x A \code{n x p} matrix containing the covariates.
#' @param kernel The kernel to use. Either \code{gaussian} (default) or \code{polynomial}.
#' @param lambda The regularization parameter.
#' @param sigma2 The scale in the \code{gaussian} and \code{polynomial} kernel. See details.
#' @param d The degree in the \code{polynomial} kernel.
#' @param threshold The convergence threshold.
#' @param max_iter The maximum number of iterations.
#' 
#' @details The \code{gaussian} kernel has the following form:
#' \deqn{exp(-||x-y||^2/sigma2).}
#' The \code{polynomial} kernel has the following form:
#' \deqn{(1+x'y/sigma2)^d.}
#'
#' @return A list containing:
#' \describe{
#' \item{\code{x}}{The original \code{x}.}
#' \item{\code{alpha}}{The vector of fitted weights.}
#' \item{\code{kernel}}{The kernel.}
#' \item{\code{sigma2}}{The scale parameter.}
#' \item{\code{d}}{The polynomial degree.}
#' }
#' @export
#' @seealso \link{predict.KLR}, \link{cv.KLR}, \link{contours.KLR}
KLR = function(
    y,
    x,
    kernel=c("gaussian", "polynomial")[1],
    lambda=0.01,
    sigma2=1.0,
    d=3,
    threshold=1.0e-6,
    max_iter=1e5
){
    # inputs check
    if(is.vector(y)) y = matrix(y, length(y), 1)
    if(lambda<1.0e-16) stop("lambda should be positive")
    if(sigma2<1.0e-16) stop("sigma2 should be positive")
    if(threshold<1.0e-16) stop("threshold should be positive")
    if(d < 1) stop("d should be larger than or equal to 1")
    
    # preparation
    n = nrow(x)
    p = ncol(x)
    if(is.null(colnames(x))) colnames(x) = paste("X",seq(p),sep="")
    vnames = colnames(x)
    
    # compute kernel
    KERNELS = c("gaussian", "polynomial")
    if(kernel == "gaussian"){
        D = as.matrix(stats::dist(x))
        K = exp( - D ^ 2 / sigma2 )
    }else if(kernel == "polynomial"){
        xs = scale(x, scale=T)
        D = xs %*% t(xs)
        K = ( 1 + D / sigma2 ) ^ d
    }else{
        stop(paste("only kernels", KERNELS, "are implemented"))
    }
    K = scale(t(scale(K, scale=F)), scale=F)
    
    # find stepsize
    mat = t(K) %*% K /4 + lambda * K
    step_size = 1./max(eigen(mat)$values)
    
    # obj function
    obj = function(alpha){
        lin_pred = K %*% alpha
        penalty = lambda * t(alpha) %*% K %*% alpha
        loss = mean(- y * lin_pred + log( 1. + exp(lin_pred) ))
        #print(c(loss, penalty, loss+penalty))
        return((loss+ penalty)[1,1])
    }
    
    # fit using gradient descent
    alpha = matrix(0, n, 1)
    obj_val_prev = obj(alpha)
    for(i in seq(max_iter)){
        # store previous
        alpha_prev = alpha
        # compute gradient
        lin_pred = y * K %*% alpha
        prob = 1. / (1. + exp(lin_pred))
        grad = - K %*% (y * prob) / n + lambda * K %*%  alpha 
        alpha = alpha - step_size * grad
        # check convergence
        obj_val = obj(alpha)
        if(abs(obj_val - obj_val_prev)/obj_val < threshold) break
        #if(max(abs(alpha - alpha_prev)) < threshold) break
        obj_val_prev = obj_val
    }
    
    # return
    out = list(x=x, alpha=alpha, kernel=kernel, sigma2=sigma2, d=d)
    class(out) = "KLR"
    return(out)
}



#' @name predict.KLR
#' 
#' @title Predict using a KLR fit
#' 
#' @description Prediction
#'
#' @param object An object of class \code{KLR}.
#' @param newx The \code{m x p} matrix of observations at which to perform prediction.
#' @param ... Extra arguments (not used).
#'
#' @return The \code{m x 1} vector of predicted probabilities.
#' @method predict KLR
#' @export 
#' @seealso \link{KLR}
predict.KLR = function(object, newx=object$x, ...){
    # construct kernel
    m = nrow(newx)
    n = nrow(object$x)
    p = ncol(object$x)
    if(object$kernel == "gaussian"){
        if(m==n){
            if(all(object$x == newx)){
                D = as.matrix(stats::dist(object$x))
            }else{
                D = matrix(pdist::pdist(object$x, newx)@dist, n, m, T)
            }
        }else{
            D = matrix(pdist::pdist(object$x, newx)@dist, n, m, T)
        }
        K = exp( - D ^ 2 / object$sigma2 )
    }else if(object$kernel == "polynomial"){
        object$x = scale(object$x, scale=T)
        newx = (newx - matrix(attr(object$x, 'scaled:center'), m, p, T)) /
            matrix(attr(object$x, 'scaled:scale'), m, p, T)
        D = object$x %*% t(newx)
        K = ( 1 + D / object$sigma2) ^ object$d
    }
    
    # compute predictors
    f = t(K) %*% object$alpha
    p = 1. / (1. + exp(-f))
    
    return(p)
}



#' @name contours.KLR
#' 
#' @title Produce level curve for a KLR object.
#'
#' @param object An object of class \code{KLR}.
#' @param dims Dimensions for which t0 produce contours. Other dimensions are set the mean.
#' @param res Resolution of the grid.
#' @param levels Levels at which to produce level curves.
#'
#' @return A list containing the desired curves. Each list has a \code{level} attribute stating the respective level as well as \code{x} and \code{y} attributes defining the curve.
#' @export
contours.KLR = function(
    object,
    dims = 1:2,
    res = 100,
    levels = c(0.5)
){
    if (res<11) stop("you should use res >10")
    if(!(length(dims) == 2)) stop("only 2D contours are possible")
    
    # get ranges
    x = object$x
    xm = colMeans(x)
    xrange = range(x[,dims[1]])
    yrange = range(x[,dims[2]])
    
    # create mesh grid
    xx = seq(xrange[1], xrange[2], length.out = res)
    yy = seq(yrange[1], yrange[2], length.out = res)
    newx_df = expand.grid(x=xx,y=yy)
    newx = matrix(xm, res*res, ncol(x), T)
    newx[,dims[1]] = newx_df$x
    newx[,dims[2]] = newx_df$y
    
    # get predictions
    preds = matrix(predict.KLR(object, newx), res, res)
    
    # produce contour
    curves = grDevices::contourLines(
            x=xx, y=yy, z=preds, levels=levels
    )
    
    return(curves)
}








