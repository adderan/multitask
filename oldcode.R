#no longer used
run.ando <- function(X.all, y.egfr, y.sens) {
  m <- length(y.egfr)
  nproblems <- dim(X.all)[[2]]  #number of prediction problems to use. The i-th prediction problem is the i-th sample in X.all
  p <- dim(X.all)[[1]]
  print(m)
  X.list <- list()
  y.list <- list()
  #build list of X.all repeated m times
  for(i in 1:nproblems) {
  	x.all.i <- matrix(unlist(X.all[,i]))
    x.all.i.manual <- matrix(0, p, 1)
    for(k in 1:p) {
    	x.all.i.manual[k,1] <- X.all[k,i]
    }
    X.list[[i]] <- x.all.i.manual
    y.list[[i]] <- y.egfr[[i]]
    #print(dim(x.all.i))
  }
   
	source("multitask.R")
  source("ando.R")
  min.out <- joint.min(X.list, y.list, 5, 2)
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
  #list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)
  Theta.hat
}

	
gray.analysis <- function(filename) { 
	source("multitask.R")
	source("ando.R")
	load(filename)

	max.features <- 500
	
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
	X.unlabeled <- X.unlabeled[, 1:100]
	X.all <- cbind(X.train, X.unlabeled)
	y.egfr <- X.all["EGFR",]
	X.all <- X.all[-195,]
	X.train <- X.train[-195,]  #remove erlotinib from training data for Ando, but not ridge regression. 
	X.test <- X.test[-195,]
		
   
  	#Theta.hat <- run.ando(X.all, y.egfr, y.train) #ando joint predictor
	#cat("Dimension of Theta.hat: ", dim(Theta.hat), "\n")
	#cat("Dimension of X: ", dim(X.train.ando), "\n")

	cat("dimension of X.all: ", dim(X.all), "\n")
	cat("length of y.egfr", length(y.egfr), "\n")
	cat("dimension of x.train: ", dim(X.train), "\n");
	cat("dimension of y.train: ", length(y.train), "\n");
	ando.X <- list(X.all, X.train);
	ando.Y <- list(y.egfr, y.train);
	joint.min.out <- joint.min(ando.X, ando.Y, 10, 3)

	Theta.hat <- joint.min.out[["Theta.hat"]]
	W.hat <- joint.min.out$W.hat
	V.hat <- joint.min.out$V.hat
	#data <- list(X.list = ando.X, y.list = ando.Y)
	#ando.test(t(ando.X), ando.Y, W.hat, V.hat, Theta.hat) 

	#labeled.predictor <- optimize.labeled(X.train, y.train, Theta.hat, 1)
	#W.hat <- labeled.predictor[["W.hat"]]
	#V.hat <- labeled.predictor[["V.hat"]]
	
	cat("Ando objective with unlabeled data: ", ando.objective(X.test, y.test, W.hat[,2], V.hat[,2], Theta.hat), "\n")

	#ando without unlabeled data
	ando.X <- list(X.train, X.train)
	ando.Y <- list(y.egfr, y.train)
	joint.min.out <- joint.min(ando.X, ando.Y, 10, 3)
	Theta.hat <- joint.min.out$Theta.hat
	W.hat <- joint.min.out$W.hat
	V.hat <- joint.min.out$V.hat
	cat("Ando objective without unlabeled data: ", ando.objective(X.test, y.test, W.hat[,2], V.hat[,2], Theta.hat), "\n")

  	beta <- ridge.regression(X.train, y.train, 1)
	cat("Training Ridge Objective: ", ridge.objective(X.train, y.train, beta), "\n")
	cat("Ridge objective: ", ridge.objective(X.test, y.test, beta), "\n")

        
  	#print(dim(ando_predictions))
             	
}


gray.analysis <- function(filename) {
	source("multitask.R")
	source("ando.R")
	load(filename)
	data <- create.ando.data(X, Y, Xu, 0, 500, 10, "Erlotinib", 195)
	joint.min.out <- joint.min(list(data$ando.X[[2]]), list(data$ando.Y[[2]]), 10, 3)
	W.hat <- joint.min.out$W.hat
	V.hat <- joint.min.out$V.hat
	Theta.hat <- joint.min.out$Theta.hat
	#print(W.hat)
	#cat("Dimension of W.hat: ", dim(W.hat), "\n")
	#cat("dimension of test data: ", dim(data$test$X.test), "\n")
	cat("Erlotinib, No unlabeled data: ", ando.objective(data$test$X.test, data$test$y.test, W.hat[,1], V.hat[,1], Theta.hat), "\n")
	data <- create.ando.data(X, Y, Xu, 0, 500, 10, "Lapatinib", 195)
	joint.min.out <- joint.min(list(data$ando.X[[2]]), list(data$ando.Y[[2]]), 10, 3)
	cat("Lapatinib, No unlabeled data: " , ando.objective(data$test$X.test, data$test$y.test, joint.min.out$W.hat[,1], joint.min.out$V.hat[,1], joint.min.out$Theta.hat), "\n")


	points <- 2
	stepsize <- 50
	for(i in 1:points) {
		erlotinibegfr <- create.ando.data(X, Y, Xu, 0, 500, stepsize*i, "Erlotinib", 195)
		lapatinibegfr <- create.ando.data(X, Y, Xu, 0, 500, stepsize*i, "Lapatinib", 195)
		lapatiniberbb2 <- create.ando.data(X, Y, Xu, 2500, 3000, stepsize*i, "Lapatinib", 2722)
		erlot.control <- create.ando.data(X, Y, Xu, 0, 500, stepsize*i, "Erlotinib", 100)
		lapat.control <- create.ando.data(X, Y, Xu, 0. 500, stepsize*i, "Lapatinib", 100)
		joint.min.erlotinibegfr <- joint.min(erlotinibegfr$ando.X, erlotinibegfr$ando.Y, 10, 3)
		joint.min.lapatinibegfr <- joint.min(lapatinibegfr$ando.X, lapatinibegfr$ando.Y, 10, 3)
		joint.min.lapatiniberbb2 <- joint.min(lapatiniberbb2$ando.X, lapatiniberbb2$ando.Y, 10, 3)
		joint.min.erlot.control <- joint.min(erlot.control$ando.X, erlot.control$ando.Y, 10, 3)
		joint.min.lapat.control <- joint.min(lapat.control$ando.X, lapat.control$ando.Y, 10, 3)

		cat("Erlotinib/EGFR with ", i*stepsize, " unlabeled samples: ", ando.objective(erlotinibegfr$test$X.test, erlotinibegfr$test$y.test, joint.min.erlotinibegfr$W.hat[,2], joint.min.erlotinibegfr$V.hat[,2], joint.min.erlotinibegfr$Theta.hat), "\n")
		cat("Lapatinib/EGFR with ", i*stepsize, " unlabeled samples: ", ando.objective(lapatinibegfr$test$X.test, lapatinibegfr$test$y.test, joint.min.lapatinibegfr$W.hat[,2], joint.min.lapatinibegfr$V.hat[,2], joint.min.lapatinibegfr$Theta.hat), "\n")	
		cat("Lapatinib/ERBB2 with ", i*stepsize, " unlabeled samples: ", ando.objective(lapatiniberbb2$test$X.test, lapatiniberbb2$test$y.test, joint.min.lapatiniberbb2$W.hat[,2], joint.min.lapatiniberbb2$V.hat[,2], joint.min.lapatiniberbb2$Theta.hat), "\n")
		cat("Erlotinib/ABCC12 (control) with ", i*stepsize, " unlabeled samples: ", ando.objective(erlot.control$test$X.test, erlot.control$test$y.test, joint.min.erlot.control$W.hat[,2], joint.min.erlot.control$V.hat[,2], joint.min.erlot.control$Theta.hat), "\n")
		cat("Lapatinib/ABCC12 (control) with ", i*stepsize, " unlabeled samples: ", ando.objective(lapat.control$test$X.test, lapat.control$test$y.test, joint.min.lapat.control$W.hat[,2], joint.min.lapat.control$V.hat[,2], joint.min.lapat.control$Theta.hat), "\n")

		cat("\n")
	}

}


