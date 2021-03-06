##implementation of Ando and Zhang's Alternating Structure Optimization.
library(glmnet)

#Trains ASO over a list of training data and response vectors.
#The parameters x and y should be lists, such that x[[problem.name]] is
#a matrix of training data and y[[problem.name]] is the corresponding response
#vector. The features should be the same for all problems. The structural
#matrix will be optimized over all problems, and the W and V vectors will
#be optimized for the primary problem after the optimal Theta has been found. An
#object containing the optimal vectors is returned.
aso.train <- function(x, y, h = 10, iters = 3, lambda = 1, ANALYTIC = TRUE, primary) {

	#number of total problems to minimize over, 
	#including the primary problem and all auxiliary problems
	n.problems <- length(x)  
	n.features <- dim(x[[1]])[[1]] #must use same set of features for all problems
	stopifnot(n.problems == length(y))
	
	#The combined weight vector (W + V*theta), represented as a matrix
	#(n.features x n.problems), where each column is the weight vector 
	#for one of the problems.
	u <- matrix(rep(0, times = n.features*n.problems), n.features, n.problems)
	 
	colnames(u) <- names(x) #set the colnames to the problem names
	rownames(u) <- rownames(x[[1]])

	#the structural matrix, shared between all problems, initially random.
	Theta.hat <- matrix(runif(h*n.features), h, n.features)  
	W.hat <- matrix(0, n.features, n.problems)  #the high-dimensional weight vector
	V.hat <- matrix(0, h, n.problems)  #the low-dimensional weight vector
	
	colnames(V.hat) <- names(x)
	
	#optimize theta and W alternately
	for(iter in 1:iters) {
		W.hat <- w.min.all.problems(x, y, u, Theta.hat, ANALYTIC, lambda)
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
	w.primary.out <- w.min.glm(x[[primary]], y[[primary]], lambda)
	v.primary.out <- v.min.glm(x[[primary]], y[[primary]], Theta.hat, lambda)
	w.primary <- w.primary.out$w
	v.primary <- v.primary.out$v
	intercept <- -1*w.primary.out$intercept - v.primary.out$intercept
	#v.primary <- V.hat[,primary]

	list(Theta.hat = Theta.hat, W.hat = w.primary, V.hat = v.primary, intercept = intercept)
}
aso.obj <- function(x, y, w, v, theta) {
	y.pred <- t(w) %*% x + t(v) %*% theta %*% x
	mse <- mean((y.pred - y)^2)
	return(mse)
}

#Takes the optimal weight vectors and Theta matrix
#produced by aso.train, and a new data matrix. Predicts
#the response for each sample in the new matrix.
aso.predict <- function(aso.trained.model, new.x) {
	#cat("Predicting for problem: ", primary.problem, "\n")
	n.samples <- dim(new.x)[[2]]
	n.features <- dim(new.x)[[1]]

	w <- aso.trained.model$W.hat
	v <- aso.trained.model$V.hat
	theta <- aso.trained.model$Theta.hat

	y.pred <- t(w) %*% new.x + t(v) %*% theta %*% new.x + aso.trained.model$intercept
	stopifnot(length(y.pred) == n.samples)
	return(y.pred)
}

#find the minimum w-vectors for each prediction problem, with a given theta. Returns the matrix. Assumes data is FxN 
w.min.all.problems <- function(X, y, u, theta, analytic, lambda) {
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

		if(!analytic) {
			library(glmnet)
			w.min.out.glm <- w.min.glm(X_l, y_l, lambda)
			W.hat[, problem.name] <- w.min.out.glm
		}
		#if(use.cache) {
		#	w.min.out.cache <- w.min.cache(X_l, y_l, problem.name, v_l, theta, lambda)
		#	W.hat[, problem.name] <- w.min.out.cache
		#}
		if(analytic) {
			w.min.out <- w.min(X_l, y_l)
			W.hat[, problem.name] <- w.min.out
		}

	}
	return(W.hat)
}
#compute the singular value decomposition of a matrix containing the minimum u vectors for each prediction problem, as calculated by w_min. Return a new theta from the first h rows of the SVD. 
theta.min <- function(u, f, m, h) {
	svd_of_U <- svd(u, nu = h, nv = 0)
 	v1 <- svd_of_U$u
	new.theta <- matrix(0, nrow=h, ncol=f)
	for(k in 1:h) {
		new.theta[k,] <- t(v1)[k,]
	}
	return(new.theta)
}
w.min <- function(x, y) {
	n <- dim(x)[[1]]
	I <- diag(n)
	w <- solve(x %*% t(x) + I, x %*% y)
	return(w)
}
v.min <- function(x, y, theta) {
	n <- dim(x)[[1]]
	I <- diag(n)
	xtheta <- theta %*% x
	h <- dim(xtheta)[[1]]
	#print(xtheta, file = sys.stderr)
	v <- try(solve(xtheta %*% t(xtheta), xtheta %*% y))
	if(class(v) == "try-error") {return(rep(NA,times=h)) }
	return(v)
}
	

#finds minimum w-vector analytically for one prediction problem. x and y should be the l-th data matrix and the lth output. 
w.min.old <- function(x, y, v, theta) {
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

w.min.glm <- function(X, y, lambda) {
	#print(length(v))
	#print(dim(theta))
	#print(dim(X))
	fit <- glmnet(x = t(X), y = y, alpha = 0)

	w <- predict(fit, type="coef", s = lambda)
	intercept <- w[1]
	w <- w[-1]
	#print(rownames(w)[1:10])
	#print(rownames(X)[1:10])
	
	return(list(w = w, intercept = intercept))
}
v.min.glm <- function(x, y, theta, lambda) {
	xtheta <- theta %*% x
	fit <- glmnet(x = t(xtheta), y = y, alpha = 0)
	v <- predict(fit, type = "coef", s = lambda)
	intercept <- v[1]
	v <- v[-1]
	return(list(v = v, intercept = intercept))
}
	
