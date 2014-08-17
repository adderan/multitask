source("../multitask/ando.R")
source("../multitask/gd.R")

joint.min <- function(X, y, h, iters) { #x is matrix of feature vectors, y is matrix of output vectors, f is number of features, m is number of learning problems. 
	m <- length(y)
	f <- dim(X[[1]])[[1]]
	#print(f)


	lambda <- c(rep(1, times=m)) #this is the factor multiplying the regularizer in the regularized loss function. 
	Theta.hat <- matrix(runif(h*f), nrow=h, ncol=f); #start with arbitrary theta
	u <- matrix(rep(0, times=f*m), f, m) #start with zero u. 
	V.hat <- Theta.hat %*% u #initialize v = theta * u
    #alternately minimize the w-vector and theta.
	for(i in 1:iters) {
		V.hat <- Theta.hat %*% u
		W.hat <- w.min.matrix(X, y, u, Theta.hat)
		u <- W.hat + t(Theta.hat) %*% V.hat
		Theta.hat <- theta.min(u, f, m, h, lambda)
  	} 
	list(W.hat = W.hat, V.hat = V.hat, Theta.hat = Theta.hat)
	#Theta.hat
}


w.prime <- function(x, y, w, v, theta) {
	n <- dim(x)[[2]]
	f <- dim(x)[[1]]
	gprime <- c()
	for(q in 1:f) {
		gprimeq <- 0
		for(i in 1:n) {
			isum <- 0
			isum <- isum + y[[i]]
			isum <- isum - (t(w) %*% x[ , i])
			isum <- isum - (t(v) %*% theta) %*% x[ , i]
			isum <- isum * (-x[q, i])
			gprimeq <- gprimeq + isum
		}
		gprimeq <- gprimeq + 2*w[[q]]
		gprime[[q]] <- gprimeq
	}
	return(gprime)
}

#find the minimum w-vectors for each prediction problem, with a given theta. Returns the matrix. Assumes data is FxN 
w.min.matrix <- function(X, y, u, theta) {
	h <- dim(theta)[[1]]  #number of dimensions for the lower dimensional map.
	m <- length(X) #number of prediction problems
 	f <- dim(X[[1]])[[1]] #number of features
	W.hat <- matrix(0, f, m)  
	for(l in 1:m) {
                #select the data for the l-th prediction problem
		X_l <- X[[l]]
		y_l <- y[[l]]
		u_l <- u[, l]
		v_l <- theta %*% u_l
		#print(l)


		w.min.out <- w.min(X_l, y_l, v_l, theta)  #in test code, X[[l]] is n*p matrix, this code uses p*n
		w.min.out.gd <- w.min.analytic.gd(X_l, y_l, v_l, theta, w.prime, f.obj1, 2)
		cat("exact solution objective: ", f.obj1(t(X_l), y_l, w.min.out, v_l, theta), "\n")
		cat("gd objective: ", f.obj1(t(X_l), y_l, w.min.out.gd, v_l, theta), "\n")

		#w.min.out <- w.gradient.descent(X_l, y_l, v_l, theta, 1, 1, 10); 
		#print(dim(W.hat))
		W.hat[ , l] <- w.min.out
	}
	return(W.hat)
}
#compute the singular value decomposition of a matrix containing the minimum u vectors for each prediction problem, as calculated by w_min. Return a new theta from the first h rows of the SVD. 
theta.min <- function(u, f, m, h, lambda) {
	print("theta.min")
	U <- matrix(1, f, m)
	for(l in 1:m) {
		U[ , l] <- sqrt(lambda[[l]])*u[ , l]
	}
	svd_of_U <- svd(U, nu = h, nv = 0)
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

