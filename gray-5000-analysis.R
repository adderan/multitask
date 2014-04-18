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
ridge.objective <- function(X, y, beta, lambda) {
	z <- t(beta) %*% X
	loss <- mean((y-z)^2)
	reg <- sum(beta^2)
	loss + lambda * reg
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
run.ando <- function(X.all, y.egfr, y.sens) {
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
        
	#make new weight vector for labeled data
	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)

}


gray.analysis <- function(filename) {
	load(filename)

	max.features <- 1000
	
	y.sens <- Y[,"Erlotinib"]
	na.removed <- remove.missing.data(X, y.sens)  #remove data that is missing from y list
	X <- na.removed$X
	y.sens <- na.removed$y

	labeled.samples <- length(y.sens)

	#test data
	X.test <- X[1:max.features, 21:42]
	X.unlabeled <- Xu[1:max.features,]
	y.test <- y.sens[21:42]

	#training data
	X.train <- X[1:max.features, 1:21]
	y.train <- y.sens[1:21]
  

	#data for ando algorithm
	X.all <- cbind(X.train, X.unlabeled)
	y.egfr <- X.all["EGFR",]
	X.all <- X.all[-195,]
		
   
  ando.out <- run.ando(X.all, y.egfr, y.train) #ando joint predictor
	Theta.hat <- ando.out$Theta.hat
	W.hat <- ando.out$W.hat
	V.hat <- ando.out$V.hat
	u <- W.hat + t(Theta.hat) %*% V.hat
	W.labeled <- w_min(X.train, y.train, u, Theta.hat)
	ando.predictor <- w.labeled + t(Theta.hat) %*% V.hat
	
	cat("Ando objective: ", f.obj1(X.test, y.test, w.labeled, V.hat, Theta.hat), "\n")

  #beta <- ridge.regression(X.train, y.train, 1)
	#cat("Ridge objective: ", ridge.objective(X.test, y.test, beta, 1), "\n")

        
  #print(dim(ando_predictions))
              
	
}



