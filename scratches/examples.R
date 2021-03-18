library(KLR)

set.seed(1)
n = 200
p = 3
y = runif(n) > 0.5
y = y + 1
mu = matrix(rnorm(2*p), 2, p)
x = mu[y, ]
x = x + matrix(rnorm(n*p), n, p)

object = KLR(y, x, "gaussian", lambda=1.0, sigma2=0.5)
contours = contours.KLR(object, dims=1:2, res=100, level=0.5)
plot(x[, 1], x[, 2], col=y)
sapply(contours, function(c) lines(c$x, c$y, type="l"))

cv.KLR(y, x, n_folds=10, lambda=c( 0.1, 1.0, 10.), sigma2=c(1.0,2.0, 3.0, 4., 5.))

object$n_iter


kernel = "polynomial"
lambda = 0.1
sigma2 = 1.0
d = 3
threshold = 1.0e-6
max_iter = 100000

alpha = rnorm(n)/1000
beta0 = 1.
obj(alpha, beta0)
gradients(alpha, beta0)

lr = 0.00001
