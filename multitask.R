
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

#normalizes each feature vector in matrix of feature vectors
normalize_matrix <- function(matr) {
	n <- dim(matr)[[2]]
	for(i in 1:n) {
		matr[ ,i] <- norm(matr[ ,i])
	}
	return(matr)
}
		
		
norm <- function(x) {
	norm = 0;
	for(i in 1:length(x)) {
		norm <- norm + x[i]*x[i]
	}
	normed_vector <- x * (1/sqrt(norm))
	return(normed_vector)
}


joint_min <- function(X, y, h) { #x is matrix of feature vectors, y is matrix of output vectors, f is number of features, m is number of learning problems. 
	m <- length(y)
	f <- dim(X[[1]])[[1]]
	print(f)


	lambda <- c(rep(1, times=m))
	Theta.hat <- matrix(runif(h*f), nrow=h, ncol=f);
	u <- matrix(rep(0, times=f*m), f, m)
	V.hat <- Theta.hat %*% u
	for(i in 1:4) {
		V.hat <- Theta.hat %*% u
		W.hat <- w_min_matrix(X, y, u, Theta.hat)
		u <- W.hat + t(Theta.hat) %*% V.hat
		
		Theta.hat <- theta_min(x, y, u, Theta.hat, lambda)
	}
	source("ando.R")
	f.obj(X, y, W.hat, V.hat, Theta.hat) 

}
#finds minimum w vector for prediction problem l. x and y should be the lth data matrix and the lth output
w_min <- function(x, y, u, theta) {
	#n <- dim(x)[[2]]
	n <- dim(x)[[2]]
	f <- dim(x)[[1]]
	#newU <- matrix(0, nrow = f, ncol = m)
	newU <- c(rep(0, times=f))
	#wmatrix <- matrix(0, nrow = f, ncol = m)
	
	#for(l in 1:m) {
		#x <- X[[l]] #choose data for l-th prediction problem
		#n <- dim(x)[[2]]
	v <- theta %*% u    #this is v subscript l
		#print(length(v))
	coefficients <- matrix(0, f, f)
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
	values <- c()
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
	#print(w)

	#newU[, l] <- w + t(theta) %*% v

	#return(list(newU, wmatrix))
	return(w)
}
w_min_matrix <- function(X, y, u, theta) {
	h <- dim(theta)[[1]]
	m <- length(X)
	f <- dim(X[[1]])[[2]]
	W.hat <- matrix(0, f, m)
	#V.hat <- matrix(0, h, m)
	for(l in 1:m) {
		X_l <- X[[l]]
		y_l <- y[[l]]
		u_l <- u[, l]
		print(l)
		w_min_out <- w_min(t(X_l), y_l, u_l, theta)  #in test code, X[[l]] is n*p matrix, this code uses p*n
		print(dim(W.hat))
		W.hat[ , l] <- w_min_out
	}
	return(W.hat)
}
ando_test_output <- function(data, h) {
	#X <- data$X.list
	#y <- data$y.list
	m <- length(data$X.list)
	p <- dim(data$X.list[[1]])[[2]]
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


