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

	beta <- ginv(t(X) %*% X + lambda * diag(b)) %*% t(X) %*% y
	cat("Length of Beta: ", length(beta), "\n")
	return(beta)
}
ridge.objective <- function(X, y, beta) {
	z <- t(beta) %*% X
	loss <- mean((y-z)^2)
	reg <- sum(beta^2)
	loss + lambda * reg
}
run.ando <- function(X.all, y.egfr) {
	m <- length(y.egfr)
  nproblems <- 2  #number of prediction problems to use. The i-th prediction problem is the i-th sample in X.all
  nfeatures <- 100 #in case I don't want to use all of the features due to computation time. 
  print(m)
  X.list <- list()
  y.list <- list()
  #build list of X.all repeated m times
  for(i in 1:nproblems) {
  	x.all.i <- matrix(unlist(X.all[,i]))
    x.all.i.manual <- matrix(0, nfeatures, 1)
    for(p in 1:nfeatures) {
    	x.all.i.manual[p,1] <- X.all[p,i]
    }
    X.list[[i]] <- x.all.i.manual
    y.list[[i]] <- y.egfr[[i]]
    #print(dim(x.all.i))
  }
   
	source("multitask.R")
  source("ando.R")
  min.out <- joint_min(X.list, y.list, 5, 2)
  W.hat <- min.out[["W.hat"]]
  V.hat <- min.out[["V.hat"]]
  cat("Dimension of W.hat: ", dim(W.hat), "\n")
  cat("Dimention of V.hat: ", dim(V.hat), "\n")
        
  Theta.hat <- min.out[["Theta.hat"]]
  cat("Dimension of Theta.hat: ", dim(Theta.hat), "\n")
  cat("Length of X.list: ", length(X.list), "\n")
  cat("Length of y.list: ", length(y.list), "\n");
  cat("Dimention of X.list[[1]]: ", dim(X.list[[1]]), "\n")
  cat("Dimension of X.list[[2]]: ", dim(X.list[[2]]), "\n")
  t.X.list <- lapply(X.list, t)
  #t.y.list <- lapply(y.list, t)
  cat("Dimension of t.X.list[[1]]: ", dim(t.X.list[[1]]), "\n")
  mydata <- list(X.list = t.X.list, y.list = y.list) #ando test code uses NxF data matrices
	ando.test(mydata, W.hat, V.hat, Theta.hat)
        
	#cat("Ando Objective: ", f.obj(t.X, t.y, W.hat, V.hat, Theta.hat))
  #ando_predictions <- t(W.hat) %*% X.all + t(V.hat) %*% Theta.hat %*% X.all
}


gray.analysis <- function(filename) {
	load(filename)
	X.labeled <- X
	X.all <- cbind(X, Xu)
	y.egfr <- X.all["EGFR",]
	X.all <- X.all[-195,]
	y.sens <- Y[,"Erlotinib"]

       
  run.ando(X.all, y.egfr) #ando joint predictor
  beta <- ridge.regression(X.labeled, y.sens, 1)
	cat("Ridge objective: ", ridge.objective(X.labeled, y.sens, beta), 1)

        
  #print(dim(ando_predictions))
              
	
}



