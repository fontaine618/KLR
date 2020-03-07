#' @name cv.KLR
#' 
#' @title Cross-validation on Kernel Logistic Regression
#' 
#' @description This function performs cross-validation to select values of tuning parameters.
#' 
#' @param y A \code{n x 1} column vector containing the responses (0-1).
#' @param x A \code{n x p} matrix containing the covariates.
#' @param n_folds Number of folds in the CV.
#' @param kernel The kernel to use. Either \code{gaussian} (default) or \code{polynomial}.
#' @param lambda The regularization parameter(s).
#' @param sigma2 The scale(s) in the \code{gaussian} and \code{polynomial} kernel. See details.
#' @param d The degree(s) in the \code{polynomial} kernel.
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
#' \item{\code{mpe}}{The mean prediction error across all folds for all combinations of parameters.}
#' \item{\code{lambda_min}}{The selected value of \code{lambda}.}
#' \item{\code{sigma2_min}}{The selected value of \code{sigma2}.}
#' \item{\code{d_min}}{The selected value of \code{d}.}
#' }
#' @export
#' @seealso \link{KLR}
cv.KLR = function(
    y,
    x,
    n_folds=5,
    kernel=c("gaussian", "polynomial")[1],
    lambda=c(1.0,0.1,0.01,0.001),
    sigma2=c(5.0,2.0,1.0,0.5),
    d=3,
    threshold=1.0e-6,
    max_iter=1e5
){
    # parameter grid
    parms = expand.grid(lambda = lambda, sigma2=sigma2, d=d)
    
    # prepare folds
    N = nrow(x)
    ids = sample(N)
    folds = cut(ids,breaks=n_folds,labels=FALSE)
    
    
    # perform CV
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(cl, list(
        "n_folds", "y", "x", "kernel", "folds",
        "threshold", "max_iter", "parms"
    ), envir=environment())
    parallel::clusterEvalQ(cl, "library(KLR)")
    mpe = t(parallel::parSapply(cl, 
        seq(nrow(parms)), 
        function(i){
        lam = parms$lambda[i]
        sig = parms$sigma2[i]
        dd = parms$d[i]
        mp = mean(sapply(seq(n_folds), function(k){
            id_in = folds == k
            fit = KLR::KLR(
                y[!id_in, , drop=F], x[!id_in, , drop=F], kernel,
                lam, sig, dd, threshold, max_iter
            )
            pred = predict.KLR(fit, x[id_in, , drop=F]) > 0.5
            1-mean(pred == y[id_in,,drop=F])
        }, simplify="array"))
        return(c(lambda=lam, sigma2=sig, d=dd, mpe=mp))
    }))
    parallel::stopCluster(cl)
    
    # mpe = t(sapply(seq(nrow(parms)), function(i){
    #     lam = parms$lambda[i]
    #     sig = parms$sigma2[i]
    #     dd = parms$d[i]
    #     mp = mean(sapply(seq(n_folds), function(k){
    #         id_in = folds == k
    #         fit = KLR(
    #             y[!id_in, , drop=F], x[!id_in, , drop=F], kernel, 
    #             lam, sig, dd, threshold, max_iter
    #         )
    #         pred = predict(fit, x[id_in, , drop=F]) > 0.5
    #         1-mean(pred == y[id_in,,drop=F])
    #     }, simplify="array"))
    #     return(c(lambda=lam, sigma2=sig, d=dd, mpe=mp))
    # }, simplify="array"))
    
    # get min
    i = which.min(mpe[,4])
    
    #return
    list(mpe = mpe, lambda_min = parms$lambda[i], sigma2_min = parms$sigma2[i], d_min = parms$d[i])
}