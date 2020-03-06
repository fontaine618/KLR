#' @name KLR
#' 
#' @title Kernel Logistic Regression
#' 
#' @description This function does
#'
#' @param y A \code{n x 1} column vector containing the responses (0-1).
#' @param x A \code{n x p} matrix containing the covariates.
#' @param kernel The kernel to use. Either \code{gaussian} (default) or \code{polynomial}.
#' @param lambda The regularization parameter.
#' @param sigma2 The scale in the \code{gaussian} kernel.
#' @param d The degree in the \code{polynomial} kernel.
#' @param threshold The convergence threshold.
#' @param max_iter The maximum number of iterations.
#' @param standardize 
#'
#' @return
#' @export
#'
#' @examples
KLR = function(
    y,
    x,
    kernel=c("gaussian", "polynomial")[1],
    lambda=0.01,
    sigma2=1.0,
    d=3,
    threshold=1.0e-6,
    max_iter=1e5,
    standardize=FALSE
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
    if(standardize) x = scale(x)
    KERNELS = c("gaussian", "polynomial")
    if(kernel == "gaussian"){
        D = as.matrix(dist(x))
        K = exp( - D ^ 2 / sigma2 )
    }else if(kernel == "polynomial"){
        D = x %*% t(x)
        K = ( 1 + D ) ^ d
    }else{
        stop(paste("only kernels", KERNELS, "are implemented"))
    }
    
    # find stepsize
    mat = t(K) %*% K /4 + lambda * K
    step_size = 1./max(eigen(mat)$values)
    
    # fit with gradient descent
    alpha = matrix(0, n, 1)
    for(i in seq(max_iter)){
        # store previous
        alpha_prev = alpha
        # compute gradient
        lin_pred = y * K %*% alpha
        prob = 1. / (1. + exp(lin_pred))
        grad = - K %*% (y * prob) /n + lambda * K %*% alpha 
        alpha = alpha - step_size * grad
        # check convergence
        if(max(abs(alpha - alpha_prev)) < threshold) break
    }
    
    # return
    if(standardize) x = x * 
        matrix(attr(x, 'scaled:scale'), n, p, T) + 
        matrix(attr(x, 'scaled:center'), n, p, T)
    out = list(x=x, alpha=alpha, kernel=kernel, sigma2=sigma2, d=d, standardize=standardize)
    class(out) = "KLR"
    return(out)
}



#' @name predict.KLR
#' 
#' @title Predict using a KLR fit
#' 
#' @description 
#'
#' @param KLRobj 
#' @param newx 
#'
#' @return
#' @export
#'
#' @examples
predict.KLR = function(KLRobj, newx){
    # construct kernel
    m = nrow(newx)
    n = nrow(KLRobj$x)
    p = ncol(KLRobj$x)
    if(KLRobj$standardize) {
        KLRobj$x = scale(KLRobj$x)
        newx = (newx - matrix(attr(KLRobj$x, 'scaled:center'), n, p, T)) / 
            matrix(attr(KLRobj$x, 'scaled:scale'), n, p, T)
    }
    if(KLRobj$kernel == "gaussian"){
        D = matrix(0, n, m)
        for(i in seq(n)){
            for(j in seq(m)){
                D[i, j] = norm(KLRobj$x[i, ] - newx[j, ], "2")
            }
        }
        K = exp( - D ^ 2 / KLRobj$sigma2 )
    }else if(KLRobj$kernel == "polynomial"){
        D = KLRobj$x %*% t(newx)
        K = ( 1 + D ) ^ KLRobj$d
    }
    
    # compute predictors
    f = t(K) %*% KLRobj$alpha
    p = 1. / (1. + exp(-f))
    
    return(p)
}












