#generate some test data. 
generateData <- function(f, g, n, epsilon) { #f is total # of features in data, g is number of prediction problems, g is number of true predictors,  n is number of feature vectors 	
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
normalize_matrix <- function(matr) {
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
	normed_vector <- x * (1/sqrt(norm))
	return(normed_vector)
}


joint_min <- function(X, y, h, iters) { #x is matrix of feature vectors, y is matrix of output vectors, f is number of features, m is number of learning problems. 
	m <- length(y)
	f <- dim(X[[1]])[[1]]
	print(f)


	lambda <- c(rep(1, times=m)) #this is the factor multiplying the regularizer in the regularized loss function. 
	Theta.hat <- matrix(runif(h*f), nrow=h, ncol=f); #start with arbitrary theta
	u <- matrix(rep(0, times=f*m), f, m) #start with zero u. 
	V.hat <- Theta.hat %*% u #initialize v = theta * u
        #alternately minimize the w-vector and theta.
	for(i in 1:iters) {
		V.hat <- Theta.hat %*% u
		W.hat <- w_min_matrix(X, y, u, Theta.hat)
		u <- W.hat + t(Theta.hat) %*% V.hat
		Theta.hat <- theta_min(u, f, m, h, lambda)
        } 
	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)

}

#finds minimum w-vector for one prediction problem. x and y should be the l-th data matrix and the lth output
w_min <- function(x, y, u, theta) {
	n <- dim(x)[[2]]  #the number of samples for this prediction problem
	f <- dim(x)[[1]]  #the number of features. Features are the rows of the data matrix. 
	v <- theta %*% u    #this is v subscript l

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
#find the minimum w-vectors for each prediction problem, with a given theta. Returns the matrix. Assumes data is FxN 
w_min_matrix <- function(X, y, u, theta) {
	h <- dim(theta)[[1]]  #number of dimensions for the lower dimensional map.
	m <- length(X) #number of prediction problems
 	f <- dim(X[[1]])[[1]] #number of features
	W.hat <- matrix(0, f, m)  
	for(l in 1:m) {
                #select the data for the l-th prediction problem
		X_l <- X[[l]]
		y_l <- y[[l]]
		u_l <- u[, l]
		print(l)
		w_min_out <- w_min(X_l, y_l, u_l, theta)  #in test code, X[[l]] is n*p matrix, this code uses p*n
		print(dim(W.hat))
		W.hat[ , l] <- w_min_out
	}
	return(W.hat)
}
ando_test_output <- function(data, h) {
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
	W.hat <- w_min_matrix(data$X.list, data$y.list, u, theta)
	u <- W.hat + t(theta) %*% V.hat
	#theta <- theta_min(u, p, m, h, lambda)
	#V.hat <- theta %*% u
	#W.hat <- w_min_matrix(data$X.list, data$y.list, u, theta)
	

	return(list(W.hat = W.hat, V.hat = V.hat, Theta.hat = theta))
}

#compute the singular value decomposition of a matrix containing the minimum u vectors for each prediction problem, as calculated by w_min. Return a new theta from the first h rows of the SVD. 
theta_min <- function(u, f, m, h, lambda) {
	U <- matrix(1, f, m)
	for(l in 1:m) {
		U[ , l] <- sqrt(lambda[[l]])*u[ , l]
	}
	svd_of_U <- svd(U, nu = h)
 	v1 <- svd_of_U$u
	newTheta <- matrix(0, nrow=h, ncol=f)

	for(k in 1:h) {
		newTheta[k,] <- t(v1)[k,]
	}
	return(newTheta)
}


