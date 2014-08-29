source("../multitask/ando.R")
source("../multitask/w-min-glm.R")
library(glmnet)

cache <- NULL

joint.min <- function(X, y, h, iters) { #x is matrix of feature vectors, y is matrix of output vectors, f is number of features, m is number of learning problems. 
	m <- length(y)
	f <- dim(X[[1]])[[1]]
	#print(f)
	X.t <- lapply(X, t)

	#lambda <- c(rep(1, times=m)) #this is the factor multiplying the regularizer in the regularized loss function. 
	Theta.hat <- matrix(runif(h*f), nrow=h, ncol=f); #start with arbitrary theta
	u <- matrix(rep(0, times=f*m), f, m) #start with zero u. 
	V.hat <- Theta.hat %*% u #initialize v = theta * u
    #alternately minimize the w-vector and theta.
	for(i in 1:iters) {
		V.hat <- Theta.hat %*% u
		W.hat <- w.min.matrix(X, y, u, Theta.hat, use.cache = FALSE)
		u <- W.hat + t(Theta.hat) %*% V.hat
		Theta.hat <- theta.min(u, f, m, h)
		cat("total objective: ", f.obj(X.t, y, W.hat, V.hat, Theta.hat), "\n")
  	} 
	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)
	#Theta.hat
}

aso.train <- function(x, y, h, iters, lambda = 1, recompute.cache = FALSE, use.cache = TRUE, cache.not.preloaded=TRUE) {

	#if a cache exists for this value of lambda, load it
	cachefile <- cache.filename(lambda)
	if(file.exists(cachefile) && use.cache && cache.not.preloaded) {
		load(cachefile)
		cache <- loaded.cache
	}
	else if(!file.exists(cachefile) && use.cache) {
		cache <- list()
	}
	n.problems <- length(x)
	n.features <- dim(x[[1]])[[1]]
	stopifnot(n.problems == length(y))

	u <- matrix(rep(0, times = n.features*n.problems), n.features, n.problems)
	colnames(u) <- names(x)
	rownames(u) <- rownames(x[[1]])

	Theta.hat <- matrix(runif(h*n.features), h, n.features)
	W.hat <- c()
	V.hat <- c()
	for(iter in 1:iters) {
		W.hat <- w.min.all.problems(x, y, u, Theta.hat, use.cache)
		V.hat <- c()
		for(l in 1:n.problems) {
			problem.name <- names(x)[[l]]
			V.hat[,problem.name] <- theta * u[,problem.name]
		}
		for(l in 1:n.problems) {
			problem.name <- names(x)[[l]]
			u[,problem.name] <- W.hat[,problem.name] + t(theta) %*% V.hat[,problem.name]
		}
		Theta.hat <- theta.min(u, f, m, h)
	}

	#write the updated cache
	loaded.cache <- cache
	save(loaded.cache, file=cachefile)
	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)
}
aso.predict <- function(aso.trained.model, new.x) {
	n.samples <- dim(new.x)[[2]]
	n.features <- dim(new.x)[[1]]

	w <- aso.trained.model$W.hat[,"primary"]
	v <- aso.trained.model$V.hat[,"primary"]
	y.pred <- t(w) %*% new.x + t(v) %*% theta %*% new.x
	stopifnot(length(y.pred) == n.samples)
	return(y.pred)
}


cache.filename <- function(lambda) {
	filename <- "cache/"
	filename <- paste(filename, "lambda", sep="")
	filename <- paste(filename, lambda, sep="")
	return(filename)
}
w.min.cache <- function(x, y, y_description, v, theta) {
	if(!(y_description %in% names(cache))) {
		add.to.cache(x, y, y_description)
	}
	else {
		cat("Retrieving ", y_description, " from the cache.\n")
	}
	w.precomputed <- cache[[y_description]]
	w.new <- v %*% theta
	w <- w.precomputed - w.new
	return(w)
}
add.to.cache <- function(x, y, y_description, lambda) {
	cat("Computing ", y_description, " from scratch and adding it to the cache.\n")
	fit <- glmnet(t(x), y, alpha = 0)
	w.precomputed <- predict(fit, type= "coef", s = lambda)
	cache[[y_description]] <- w.precomputed
}

#find the minimum w-vectors for each prediction problem, with a given theta. Returns the matrix. Assumes data is FxN 
w.min.all.problems <- function(X, y, u, theta, use.cache=FALSE) {
	h <- dim(theta)[[1]]  #number of dimensions for the lower dimensional map.
	m <- length(X) #number of prediction problems
 	f <- dim(X[[1]])[[1]] #number of features
	W.hat <- matrix(0, f, m)  
	for(l in 1:m) {
        #select the data for the l-th prediction problem
		problem.name <- names(X)[[l]]
		X_l <- X[[problem.name]]
		y_l <- y[[problem.name]]
		u_l <- u[, problem.name]
		v_l <- theta %*% u_l
		#print(l)

		#w.min.out <- w.min(X_l, y_l, v_l, theta)  #in test code, X[[l]] is n*p matrix, this code uses p*n

		if(!use.cache) {
			w.min.out.glm <- w.min.glm(X_l, y_l, v_l, theta)
			W.hat[, problem.name] <- w.min.out.glm
		}
		if(use.cache) {
			w.min.out.cache <- w.min.cache(X_l, y_l, problem.name, v_l, theta)
			W.hat[, problem.name] <- w.min.out.cache
		}

	}
	return(W.hat)
}
#compute the singular value decomposition of a matrix containing the minimum u vectors for each prediction problem, as calculated by w_min. Return a new theta from the first h rows of the SVD. 
theta.min <- function(u, f, m, h) {
	print("theta.min")
	U <- matrix(1, f, m)
	#for(l in 1:m) {
	#	U[ , l] <- u[ , l]
	#}
	svd_of_U <- svd(u, nu = h, nv = 0)
 	v1 <- svd_of_U$u
	#save(U, file="testmatrix.RData")
	new.theta <- matrix(0, nrow=h, ncol=f)

	for(k in 1:h) {
		new.theta[k,] <- t(v1)[k,]
	}
	return(new.theta)
}

#finds minimum w-vector for one prediction problem. x and y should be the l-th data matrix and the lth output
w.min <- function(x, y, v, theta) {
	print("w-min")
	n <- dim(x)[[2]]  #the number of samples for this prediction problem
	f <- dim(x)[[1]]  #the number of features. Features are the rows of the data matrix. 
	#v <- theta %*% u    #this is v subscript l

        #Compute the derivative of the loss function with respect to each feature. Set each derivative equal to zero and solve system of equations to find minimum w-vector. 
	coefficients <- matrix(0, f, f) #coefficient matrix of the linear system. Each row contains the coefficients for one equation. The i-th row contains the coefficients of the derivative of the loss function with rexpect to x_i  
	for(q in 1:f) {
		for(k in 1:f) {
			sum <- 0
			for(i in 1:n) {
				sum <- sum + x[k, i]*x[q, i]
			}
			if(k == q) sum <- sum + 2
			coefficients[q, k] = sum
		}
	}
	values <- c() #vector of the right-hand sides of each linear equation (the part not dependent on the x_i's.
	for(q in 1:f) {
		val <- 0
		for(i in 1:n) {
			val = val + y[[i]]*x[q, i]
                        #cat("Length of x[,i]: ", length(x[,i]), "\n")
                        #cat("Dimension of theta: ", dim(theta), "\n")
			thetatimesxi <- theta %*% x[ , i]
			vtimesthetatimesxi <- t(v) %*% thetatimesxi

			val = val - vtimesthetatimesxi*x[q, i]
		}
		values[[q]] = val
	}
	w <- solve(coefficients, values)   #this is the minimum w for problem l

	return(w)
}

