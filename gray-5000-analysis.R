#use Ando algorithm to train predictor on y_aux, use to predict y
#ando_prediction <- function(X, y, y_aux,  h) {
#	source("multitask.R")
#	u <- joint_min(X, y_aux, h)

#ridge regression is a standard algorithm that will be used for comparison.
ridge.regression <- function(X, y, lambda) {
	library(MASS)
	X <- t(X)
	a <- dim(X)[[1]]
	b <- dim(X)[[2]]
	cat("Dimension of X: ", dim(X), "\n")
	cat("Length of y: ", length(y), "\n")

	x1 <- ginv(t(X) %*% X + lambda * diag(b))
	x2 <- t(X) %*% y
	#print(x2)
	beta <- x1 %*% x2
	cat("Length of Beta: ", length(beta), "\n")
	return(beta)
}
ridge.objective <- function(X, y, beta) {
	z <- t(beta) %*% X
	loss <- mean((y-z)^2)
	#reg <- sum(beta^2)
	loss
}
ando.objective <- function(X, y, w, v, theta) {
	z <- t(w) %*% X + t(v) %*% theta %*% X
	return(mean(y-z)^2)
}
	
remove.missing.data <- function(x, y) {
	m <- 1
	while(m < length(y)) {
		if(is.na(y[[m]])) {
			x <- x[,-m]
			y <- y[-m]
		}
		else {
			m <- m + 1
		}
	}
	list(X = x, y = y)
}

ando.run <- function(X.lab, y, X.unlab, feature.min, feature.max, n.unlabeled, use.unlabeled, inhibitor, gene) {
	source("multitask.R")
	source("ando.R")
	
	y.inhib <- y[,inhibitor]
	na.removed <- remove.missing.data(X.lab, y.inhib)
	X.labeled <- na.removed$X
	y.labeled <- na.removed$y
	#na.removed.unlabeled <- remove.missing.data(X.unlab, y.inhib)
	X.unlabeled <- X.unlab

	labeled.size <- dim(X.labeled)[[2]]

	#cat("dimension of x.unlabeled: ", dim(X.unlabeled), "\n")
	#cat("dimension of X.labeled: ", dim(X.labeled), "\n")
	
	y.gene.labeled <- X.labeled[gene,]
	y.gene.unlabeled <- X.unlabeled[gene,]
	X.labeled[-gene,]
	X.unlabeled[-gene,]

	#cat("dimension of X.labeled: ", dim(X.labeled), "\n")
	X.train <- X.labeled[feature.min:feature.max, 1:floor(labeled.size/2)]
	X.test <- X.labeled[feature.min:feature.max, ceiling(labeled.size/2):labeled.size]
	y.test <- y.labeled[ceiling(labeled.size/2):labeled.size]
	y.train <- y.labeled[1:floor(labeled.size/2)] 
	y.gene.train <- y.gene.labeled[1:floor(labeled.size/2)]
	X.unlabeled <- X.unlabeled[feature.min:feature.max,]


	X.all <- cbind(X.train, X.unlabeled)
	X.all <- X.all[,1:n.unlabeled]
	y.gene <- c(y.gene.labeled, y.gene.unlabeled)
	y.gene <- y.gene[1:n.unlabeled]

	ando.X <- list()
	ando.Y <- list()
	if(use.unlabeled) {
		ando.X <- list(X.train, X.all)
		ando.Y <- list(y.train, y.gene)
	}
	else {
		ando.X <- list(X.train)
		ando.Y <- list(y.train)
	}
	#test <- list(X.test = X.test, y.test = y.test)
	#list(ando.X = ando.X, ando.Y =  ando.Y, test = test)
	joint.min.out <- joint.min(ando.X, ando.Y, 10, 3)
	W.hat <- joint.min.out$W.hat
	V.hat <- joint.min.out$V.hat
	Theta.hat <- joint.min.out$Theta.hat
	cat("Objective for ", inhibitor, " with ", n.unlabeled, " samples for Auxiliary problem ", rownames(X.lab)[[gene]], ":", ando.objective(X.test, y.test, W.hat[,1], V.hat[,1], Theta.hat), "\n")

}
gray.analyze <- function(filename) {
	load(filename)
	sink("gray-5000.out")
	ando.run(X, Y, Xu, 0, 500, 0, 0, "Erlotinib", 195)
	ando.run(X, Y, Xu, 0, 500, 100, 1, "Erlotinib", 195)
	ando.run(X, Y, Xu, 0, 500, 200, 1, "Erlotinib", 195)
	ando.run(X, Y, Xu, 0, 500, 200, 1, "Erlotinib", 100)

	ando.run(X, Y, Xu, 2500, 3000, 0,  0, "Lapatinib", 2722)
	ando.run(X, Y, Xu, 2500, 3000, 100, 1, "Lapatinib", 2722)
	ando.run(X, Y, Xu, 2500, 3000, 200, 1, "Lapatinib", 2722)
	ando.run(X, Y, Xu, 2500, 3000, 100, 1, "Lapatinib", 2600)
	ando.run(X, Y, Xu, 2500, 3000, 200, 1, "Lapatinib", 2600)

	ando.run(X, Y, Xu, 0, 500, 200, 1, "Lapatinib", 195)
	ando.run(X, Y, Xu, 0, 500, 100, 1, "Lapatinib", 195)
}
