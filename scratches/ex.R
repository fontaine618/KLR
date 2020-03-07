library(ggplot2)

data_train = read.csv("~/Dropbox/H2020/STATS601/601-HW4/data/classification_dat.txt", header=F, sep=" ")
colnames(data_train) = c("X1", "X2", "class")
x = as.matrix(data_train[, 1:2])
y = as.matrix(data_train[, 3])
newx = x


kernel = "gaussian"
sigma2 = 1.4
d = 3
lambda = 0.012
threshold=1.0e-6
max_iter=1e5

# kernel = "polynomial"
# sigma2 = 20.
# d = 3
# lambda = 0.01
# threshold=1.0e-6
# max_iter=1e5


KLRobj = KLR(
    y=y,
    x=x,
    kernel=kernel,
    lambda=lambda,
    sigma2=sigma2,
    d=d,
    threshold=threshold,
    max_iter=max_iter
)

preds = predict(KLRobj)
1 - mean((preds > .5) == y)



res = 100
dims = 1:2
levels = c(0.5, 0.51)
lvl = 0.5


contours(KLRobj, 1:2, 100, 0.5)


n = 100
X1 = seq(-15,5,length.out = n)
X2 = seq(-20,5,length.out = n)
newx = as.matrix(expand.grid(X1=X1,X2=X2))
preds = predict(KLRobj, newx)
hist(preds)

KLR_df = data.frame(
    X1=newx[,1],
    X2=newx[,2],
    pred=preds,
    model="KLR"
)

data_train$class = factor(data_train$class)
ggplot() +
    geom_point(
        data=data_train, 
        aes(x=X1, y=X2, shape=class, color=class)
    ) +
    geom_contour(
        data=KLR_df, 
        aes(x=X1, y=X2, z=pred, linetype=model),
        breaks=c(0.25,0.5,0.75),
        color='black'
    )







data_train = read.csv("~/Dropbox/H2020/STATS601/601-HW4/data/classification_dat.txt", header=F, sep=" ")
colnames(data_train) = c("X1", "X2", "class")
x = as.matrix(data_train[, 1:2])
y = as.matrix(data_train[, 3])
newx = x



kernel = "gaussian"
d = 3
lambda=c(1.0,0.1,0.01,0.001)
sigma2=c(5.0,2.0,1.0,0.5)
threshold=1.0e-6
max_iter=1e5
n_folds=3

cv.KLR(
    y=y,
    x=x,
    n_folds=n_folds,
    kernel=kernel,
    lambda=lambda,
    sigma2=sigma2,
    d=d,
    threshold=threshold,
    max_iter=max_iter
)
