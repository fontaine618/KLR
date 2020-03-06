data_train = read.csv("~/Dropbox/H2020/STATS601/601-HW4/data/classification_dat.txt", header=F, sep=" ")
colnames(data_train) = c("X1", "X2", "class")
x = as.matrix(data_train[, 1:2])
y = as.matrix(data_train[, 3])
newx = x


kernel = "polynomial"
sigma2 = 2000.0
d = 3
lambda = 0.005
threshold=1.0e-6
max_iter=1e5
standardize=FALSE


KLRobj = KLR(
    y=y,
    x=x,
    kernel=kernel,
    lambda=lambda,
    sigma2=sigma2,
    d=d,
    threshold=threshold,
    max_iter=max_iter,
    standardize=standardize
)

preds = predict.KLR(KLRobj, newx)
mean((preds > .5) == y)

cbind(
    preds > 0.5,
    y
)

KLRobj$alpha
