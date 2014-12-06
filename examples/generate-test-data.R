#!/usr/bin/env Rscript
args <- commandArgs()
outfile <- args[6]
targetfile <- file(args[7])

n.samples <- 1000
n.labeled <- 200
n.features <- 100


n.good.features <- 1 #number of features that predict y
n.feature.predictors <- 5 #number of features that predict the features that predict y

n.random.features <- n.features - n.good.features
n.unlabeled.samples <- 800

X.random <- matrix(rnorm(n.random.features*n.samples, 0, 1), n.random.features, n.samples)

X.good.features <- matrix(0, n.good.features, n.samples)


#features that predict the good features
feature.predictors <- sample(n.random.features, n.feature.predictors)
#weights for each of the feature predictors
alpha <- matrix(runif(n.feature.predictors, -1, 1), 1, n.feature.predictors)

for(i in 1:n.good.features) {
	X.good.features[i,] <- alpha %*% X.random[feature.predictors,]
}

#weights for each of the good features when predicting y
beta <- matrix(runif(n.good.features, -1, 1), 1, n.good.features)

y.total <- beta %*% X.good.features


X.total <- rbind(X.random, X.good.features)

rnamesX <- c()
for(i in 1:n.features) {
	rnamesX[i] <- paste("f", i, sep = "")
}
rownames(X.total) <- rnamesX
cnamesX <- c()
for(j in 1:n.samples) {
	cnamesX[j] <- paste("s", j, sep = "")
}
colnames(X.total) <- cnamesX

stopifnot(dim(X.total) == c(n.features, n.samples))

X <- X.total[,1:n.labeled]
Xu <- X.total[,(n.labeled+1):dim(X.total)[[2]]]

Y <- matrix(y.total[1:n.labeled])

cnamesY <- c()
for(k in 1:dim(Y)[[2]]) {
	cnamesY <- paste("prob", k, sep = "")
}
colnames(Y) <- cnamesY
rownames(Y) <- colnames(X)

target <- rownames(X)[dim(X)[[1]]]
problem <- colnames(Y)[[1]]


save(X, Xu, Y, file = outfile)
cat(problem, "\t", target, "\n", file = targetfile)
