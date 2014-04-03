#use Ando algorithm to train predictor on y_aux, use to predict y
#ando_prediction <- function(X, y, y_aux,  h) {
#	source("multitask.R")
#	u <- joint_min(X, y_aux, h)
	
#ridge.regression <- function(
main <- function() {
	load("gray-5000.RData")
	X.all <- cbind(X, Xu)
	y.egfr <- X.all["EGFR",]
	X.all <- X.all[-195,]
	y.sens <- Y[,"Erlotinib"]

        #ando joint predictor
        source("multitask.R")
        source("ando.R")
        min.out <- joint_min(X.all, y.egfr, 2)
        W.hat <- min.out[["W.hat"]]
        V.hat <- min.out[["V.hat"]]
        Theta.hat <- min.out[["Theta.hat"]]
	pert.sol(X.all, y.egfr, W.hat, V.hat, Theta.hat)
        

        ando_predictions <- t(W.hat) %*% X.all + t(V.hat) %*% Theta.hat %*% X.all
        
        print(dim(ando_predictions))
        
        
	
}



