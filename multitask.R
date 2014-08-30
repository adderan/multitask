source("../multitask/ando.R")
source("../multitask/w-min-glm.R")
library(glmnet)

cache <- NULL
cache.filename <- ""

aso.train <- function(x, y, h, iters, lambda = 1, use.cache = TRUE) {
	if(use.cache && is.null(cache)) {
		stop("Must load a cache first if use.cache is TRUE.")
	}

	n.problems <- length(x)  #number of total problems to minimize over, including the primary problem and all auxiliary problems
	n.features <- dim(x[[1]])[[1]] #must use same set of features for all problems
	stopifnot(n.problems == length(y))
	
	u <- matrix(rep(0, times = n.features*n.problems), n.features, n.problems) #each column will contain the sum w + v*theta
	colnames(u) <- names(x)
	rownames(u) <- rownames(x[[1]])

	Theta.hat <- matrix(runif(h*n.features), h, n.features)  #the structural matrix, shared between all problems, that will be minimized.
	W.hat <- matrix(0, n.features, n.problems)  #the high-dimensional weight vector
	V.hat <- matrix(0, h, n.problems)  #the low-dimensional weight vector
	
	colnames(V.hat) <- names(x)
	for(iter in 1:iters) {
		W.hat <- w.min.all.problems(x, y, u, Theta.hat, use.cache, lambda)
		for(l in 1:n.problems) {
			problem.name <- names(x)[[l]]
			V.hat[,problem.name] <- Theta.hat %*% u[,problem.name]
		}
		for(l in 1:n.problems) {
			problem.name <- names(x)[[l]]
			u[,problem.name] <- W.hat[,problem.name] + t(Theta.hat) %*% V.hat[,problem.name]
		}
		Theta.hat <- theta.min(u, n.features, n.problems, h)
	}

	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)
}
load.cache <- function(filename) {
	if(file.exists(filename)) {
		print("Loading the specified cache.")
		load(filename)
		cache <<- loaded.cache
	}
	else {
		print("Cache does not exist yet. Creating it.")
		cache <<- list()
	}
	cache.filename <<- filename

}
save.cache <- function() {
	loaded.cache <- cache
	save(loaded.cache, file=cache.filename)
}
	
aso.predict <- function(aso.trained.model, new.x, primary.problem) {
	cat("Predicting for problem: ", primary.problem, "\n")
	n.samples <- dim(new.x)[[2]]
	n.features <- dim(new.x)[[1]]

	w <- aso.trained.model$W.hat[,primary.problem]
	v <- aso.trained.model$V.hat[,primary.problem]
	theta <- aso.trained.model$Theta.hat
	y.pred <- t(w) %*% new.x + t(v) %*% theta %*% new.x
	stopifnot(length(y.pred) == n.samples)
	return(y.pred)
}


cache.filename <- function(lambda, cache.name) {
	filename <- cache.name
	filename <- paste(filename, "-lambda", sep="")
	filename <- paste(filename, lambda, sep="")
	filename <- paste(filename, ".RData", sep="")
	return(filename)
}

#add v * theta to the precomputed weight vector to use as minimum w-vector for this value of theta
w.min.cache <- function(x, y, problem.name, v, theta, lambda) {
	if(!(problem.name %in% names(cache))) {
		add.to.cache(x, y, problem.name, lambda)
	}
	else {
		#cat("Retrieving ", problem.name, " from the cache.\n")
	}
	w.precomputed <- cache[[problem.name]]
	w.new <- t(v) %*% theta
	w <- w.precomputed - w.new
	return(w)
}
add.to.cache <- function(x, y, problem.name, lambda) {
	cat("Computing ", problem.name, " from scratch and adding it to the cache.\n")
	fit <- glmnet(t(x), y, alpha = 0)
	w.precomputed <- predict(fit, type= "coef", s = lambda)
	w.precomputed <- w.precomputed[-1]
	w.precomputed <- as.vector(w.precomputed)

	#must use different operator to change the global cache
	cache[[problem.name]] <<- w.precomputed
}

#find the minimum w-vectors for each prediction problem, with a given theta. Returns the matrix. Assumes data is FxN 
w.min.all.problems <- function(X, y, u, theta, use.cache=FALSE, lambda) {
	h <- dim(theta)[[1]]  #number of dimensions for the lower dimensional map.
	m <- length(X) #number of prediction problems
 	f <- dim(X[[1]])[[1]] #number of features

	W.hat <- matrix(0, f, m)
	colnames(W.hat) <- names(X)
	rownames(W.hat) <- rownames(X[[1]])

	for(l in 1:m) {
        #select the data for the l-th prediction problem
		problem.name <- names(X)[[l]]
		X_l <- X[[problem.name]]
		y_l <- y[[problem.name]]
		u_l <- u[, problem.name]
		v_l <- theta %*% u_l

		#w.min.out <- w.min(X_l, y_l, v_l, theta)  #in test code, X[[l]] is n*p matrix, this code uses p*n

		if(!use.cache) {
			w.min.out.glm <- w.min.glm(X_l, y_l, v_l, theta)
			W.hat[, problem.name] <- w.min.out.glm
		}
		if(use.cache) {
			w.min.out.cache <- w.min.cache(X_l, y_l, problem.name, v_l, theta, lambda)
			W.hat[, problem.name] <- w.min.out.cache
		}

	}
	return(W.hat)
}
#compute the singular value decomposition of a matrix containing the minimum u vectors for each prediction problem, as calculated by w_min. Return a new theta from the first h rows of the SVD. 
theta.min <- function(u, f, m, h) {
	U <- matrix(1, f, m)
	svd_of_U <- svd(u, nu = h, nv = 0)
 	v1 <- svd_of_U$u
	new.theta <- matrix(0, nrow=h, ncol=f)

	for(k in 1:h) {
		new.theta[k,] <- t(v1)[k,]
	}
	return(new.theta)
}

#finds minimum w-vector analytically for one prediction problem. x and y should be the l-th data matrix and the lth output. This is slow and should only be used for testing. 
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
			thetatimesxi <- theta %*% x[ , i]
			vtimesthetatimesxi <- t(v) %*% thetatimesxi

			val = val - vtimesthetatimesxi*x[q, i]
		}
		values[[q]] = val
	}
	w <- solve(coefficients, values)   #this is the minimum w for problem l

	return(w)
}

