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
#generate some test data. 
generate.data <- function(f, g, n, epsilon) { #f is total # of features in data, g is number of prediction problems, g is number of true predictors,  n is number of feature vectors 	
	x <- matrix(runif(f*n), f, n)
	w <- norm(c(rep(1, times=g), rep(0, times=(f-g))))
	y <- matrix(runif(g*n), g, n)
	x <- normalize_matrix(x)
	y <- normalize_matrix(y)

	for(k in 1:(g/2)) {
		for(i in 1:n) {
			y[k, i] <- w%*%(x[, i]) + rnorm(1, mean = 0, sd = epsilon)
		}


	}
	return(list(x, y))
}

#normalizes each feature vector in a matrix of feature vectors. Only useful for normalizing the test data. 
normalize.matrix <- function(matr) {
	n <- dim(matr)[[2]]
	for(i in 1:n) {
		matr[ ,i] <- norm(matr[ ,i])
	}
	return(matr)
}
		
#normalizes a vector to one.		
norm <- function(x) {
	norm = 0;
	for(i in 1:length(x)) {
		norm <- norm + x[i]*x[i]
	}
	normed.vector <- x * (1/sqrt(norm))
	return(normed.vector)
}

w.gradient.descent <- function(x, y, v, theta, restarts, step, iters) {
	f <- dim(x)[[1]]
	w.hat <- c(rep(0, f))
	minobj <- 10000000
	for(r in 1:restarts) {
		new.w <- w.find.local.min.numeric(x, y, v, theta, step, iters)
		new.w.obj <- f.obj1(t(x), y, new.w, v, theta)
		if(new.w.obj < minobj) {
			w.hat <- new.w
			minobj <- new.w.obj
		}
	}
	return(w.hat)
}
		
w.find.local.min.numeric <- function(x, y, v, theta, step, iters) {
	f <- dim(x)[[1]]
	w <- c(runif(f, -1, 1))
	gmin <- f.obj1(t(x), y, w, v, theta)
	iter <- 0
	while(iter < iters) {
		wstep <- w - dgdw.numeric(x, y, w, v, theta, step)
		gnew <- f.obj1(t(x), y, wstep, v, theta)
		if(gnew > gmin) {
			step <- step /2;
		}
		else {
			w <- wstep;
			gmin <- gnew;
			iter <- iter + 1;
		}
	}
	return(w)
	
}
dgdw.numeric <- function(x, y, w, v, theta, step) {
	f <- dim(x)[[1]]
	dg <- c(rep(0, f))
	for(k in 1:f) {
		wstep <- w
		wstep[k] <- wstep[k] + step
		g0 <- f.obj1(t(x), y, w, v, theta)
		g1 <- f.obj1(t(x), y, wstep, v, theta)
		dg[k] <- (g1 - g0)/step
	}
	return(dg)
}
	
ando.test.output <- function(data, h) {
	#X <- data$X.list
	#y <- data$y.list
	m <- length(data$X.list)
	p <- dim(data$X.list[[1]])[[1]]
	print(m)
	print(p)
	u <- matrix(0, p, m)
	mat <- array(runif(h*p), dim=c(h,p))
	#theta <- qr.Q(qr(mat))
	#print(t(theta) %*% theta)
	#print(dim(theta))
	lambda <- c(rep(1, times=m))
	theta <- matrix(runif(h*p), h, p)
	V.hat <- theta %*% u
	W.hat <- w.min.matrix(data$X.list, data$y.list, u, theta)
	u <- W.hat + t(theta) %*% V.hat
	#theta <- theta_min(u, p, m, h, lambda)
	#V.hat <- theta %*% u
	#W.hat <- w_min_matrix(data$X.list, data$y.list, u, theta)
	

	return(list(W.hat = W.hat, V.hat = V.hat, Theta.hat = theta))
}

#calculate derivative of vT*theta*x for v minimization
dhdVq <- function(x, q, theta) {
	d <- 0
	for(k in 1:length(x)) {
		d <- d + x[k]*theta[q,k]
	}
	#cat("Dimension of dhdVq: ", length(d), "\n")
	return(d)
}


#compute minimum v given w and theta. Used to minimize v for target problem.
v.min <- function(x, y, w, theta) {
	p <- dim(x)[[1]]
	n <- dim(x)[[2]]
	h <- dim(theta)[[1]]
	
	cat("p = ", p, "\n")
	print(x)

	coefficients <- matrix(0, h, h) 

	for(q in 1:h) {
		#coefficients.q <- c(length=h)
		for(j in 1:h) {
				coefficient.qj <- 0
				for(i in 1:n) {
					coefficient.qji <- 0
					for(k in 1:p) {
						coefficient.qji <- coefficient.qj + theta[j,k]*x[k,i]
					}
					coefficient.qji <- coefficient.qj * dhdVq(x[,i], q, theta)
					coefficient.qj <- coefficient.qj + coefficient.qji
				}
				if(q==j) {
					coefficient.qj <- coefficient.qj + 2
				}
				coefficients[q,j] <- coefficient.qj
		}
		#sum <- c(length=p)
		#for(i in 1:n) {
		#	sum <- sum + x[,i]*dhdVq(x[,i], q, theta)
		#}
		#coefficients[q,] <- theta %*% sum
		#coefficients[q,q] <- coefficients[q,q] + 2
	}
	values <- c()
	for(q in 1:h) {
		value.q <- 0
		for(i in 1:n) {
			value.q <- (y[i] - t(w)%*%x[,i])*dhdVq(x[,i], q, theta)
		}
		values[q] <- value.q
	}
	v <- solve(coefficients, values)

}
v.gradient.descent <- function (x, y, w, theta, iters, restarts) {
	n <- dim(x)[[2]]
	p <- dim(x)[[1]]
	h <- dim(theta)[[1]]
	vmin <- c()
	objmin <- 1000000
	for(r in 1:restarts) {
		v <- c(runif(h, -1, 1))
		vobj <- f.obj1(t(x), y, w, v, theta)
		stepsize <- 0.1
		iter <- 0
		while(iter < iters) {
			vtest <- v - v.prime.numeric(x, y, w, v, theta, stepsize)
			vtestobj <- f.obj1(t(x), y, w, vtest, theta)
			if(vtestobj > vobj) stepsize <- stepsize/2
			else {
				v <- vtest
				vobj <- f.obj1(t(x), y, w, v, theta)
				iter <- iter + 1
			}
			#print(v.prime.numeric(x, y, w, v, theta, 0.000001))
			print(vtestobj)
		}
		finalobj <- f.obj1(t(x), y, w, v, theta)
		print(finalobj)
		if(finalobj < objmin) {
			vmin <- v
			objmin <- finalobj
		}
	}
	cat("Final ando training objective: ", objmin, "\n")
	vmin
}
v.prime.numeric <- function(x, y, w, v, theta, step) {
	source("ando.R")
	h <- dim(theta)[[1]]
	vprime <- c(rep(0, h))
	for(j in 1:h) {
		vstep <- v
		vstep[j] <- vstep[j] + step
		vstepobj <- f.obj1(t(x), y, w, vstep, theta)
		vobj <- f.obj1(t(x), y, w, v, theta)
		vprime[j] <- (vstepobj - vobj)/h
	}
	return(vprime)
}

v.prime <- function(x, y, w, v, theta) {
	vprime <- c(length = length(v))
	n <- dim(x)[[2]]
	for(q in 1:length(v)) {
		vprime.q <- 0
		for(i in 1:n) {
			vprime.q <- vprime.q + (t(w) %*% x[,i] + t(v) %*% theta %*% x[,i] - y[[i]])*dhdVq(x[,i],q,theta)
		}
		vprime[q] <- vprime.q
	}
	vprime
}
	
optimize.labeled <- function(x, y, theta, iters) {
	f <- dim(x)[[1]]
	n <- dim(x)[[2]]
	h <- dim(theta)[[1]]
	print(h)
	v <- c(runif(h, -1, 1))
	w <- c(runif(f, -1, 1))
	cat("Dimension of theta: ", dim(theta), "\n")
	cat("Dimension of x: ", dim(x), "\n")
	cat("Length of v:", h, "\n")
	for(i in 1:iters) {
		w <- w.min(x, y, v, theta)
		v <- v.gradient.descent(x, y, w, theta, 100, 20)
	}
	list(W.hat = w, V.hat = v)
}


