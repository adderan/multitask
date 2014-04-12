#use Ando algorithm to train predictor on y_aux, use to predict y
#ando_prediction <- function(X, y, y_aux,  h) {
#	source("multitask.R")
#	u <- joint_min(X, y_aux, h)
	
#ridge.regression <- function(
gray.analysis <- function(filename) {
	load(filename)
	X.all <- cbind(X, Xu)
	y.egfr <- X.all["EGFR",]
	X.all <- X.all[-195,]
	y.sens <- Y[,"Erlotinib"]

        m <- length(y.egfr)
        nproblems <- 2  #number of prediction problems to use. The i-th prediction problem is the i-th sample in X.all
        nfeatures <- 100 #in case I don't want to use all of the features due to computation time. 
        print(m)
        X <- list()
        y <- list()
        #build list of X.all repeated m times
        for(i in 1:nproblems) {
          x.all.i <- matrix(unlist(X.all[,i]))
          x.all.i.manual <- matrix(0, nfeatures, 1)
          for(p in 1:nfeatures) {
            x.all.i.manual[p,1] <- X.all[p,i]
          }
          X[[i]] <- x.all.i.manual
          y[[i]] <- y.sens[[i]]
          #print(dim(x.all.i))
        }
        
        #ando joint predictor
        source("multitask.R")
        source("ando.R")
        min.out <- joint_min(X, y, 5, 2)
        W.hat <- min.out[["W.hat"]]
        V.hat <- min.out[["V.hat"]]
        cat("Dimension of W.hat: ", dim(W.hat), "\n")
        cat("Dimention of V.hat: ", dim(V.hat), "\n")
        
        Theta.hat <- min.out[["Theta.hat"]]
        cat("Dimension of Theta.hat: ", dim(Theta.hat), "\n")
        cat("Length of X: ", length(X), "\n")
        cat("Length of y: ", length(y), "\n");
        cat("Dimention of X[[1]]: ", dim(X[[1]]), "\n")
        cat("Dimension of X[[2]]: ", dim(X[[2]]), "\n")
        t.X <- lapply(X, t)
        t.y <- lapply(y, t)
        cat("Dimension of t.X[[1]]: ", dim(t.X[[1]]), "\n")
        mydata <- list(X.list = t.X, y.list = t.y) #ando test code uses NxF data matrices
	ando.test(mydata, W.hat, V.hat, Theta.hat)
        

        ando_predictions <- t(W.hat) %*% X.all + t(V.hat) %*% Theta.hat %*% X.all
        
        print(dim(ando_predictions))
        
        
	
}



