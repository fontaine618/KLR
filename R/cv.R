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
    parms = expand.grid(lambda = lambda, sigma2=sigma2)
    
    # prepare folds
    N = nrow(x)
    ids = sample(N)
    folds = cut(ids,breaks=n_folds,labels=FALSE)
    
    # perform CV
    sapply(parms, function(parm){
        lam = parm[1]
        sig = parm[2]
        sapply(seq(n_folds), function(k){
            id_in = folds == k
            fit = KLR(
                y[!id_in, , drop=F], x[!id_in, , drop=F], kernel, 
                lam, sig, d, threshold, max_iter
            )
        }, simplify="array")
    }, simplify="array")
}