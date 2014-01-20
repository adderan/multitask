
generateData <- function(f, g, n, epsilon) { #f is total # of features in data, g is number of prediction problems, g is number of true predictors,  n is number of feature vectors 	
	x <- matrix(runif(f*n), f, n)
	w <- c(rep(1, times=g), rep(0, times=(f-g)))
	y <- matrix(runif(g*n), g, n)

	for(k in 1:(g/2)) {
		for(i in 1:n) {
			y[k, i] <- w%*%(x[, i]) + rnorm(1, mean = 0, sd = epsilon)
		}
		#y[k, ] <- norm(y[k, ])


	}
	return(list(x, y))
	#print(w)
}
normalize_data <- function(matr) {
	n <- dim(matr)[[2]]
	for(i in 1:n) {
		
norm <- function(x) {
	norm = 0;
	for(i in 1:length(x)) {
		norm <- norm + x[i]*x[i]
	}
	normed_vector <- x * (1/sqrt(norm))
	return(normed_vector)
}


joint_min <- function(x, y, h) { #x is matrix of feature vectors, y is matrix of output vectors, f is number of features, m is number of learning problems. 
	f <- dim(x)[[1]]
	m <- dim(y)[[1]]
	n <- dim(x)[[2]]


	lambda <- c(rep(1, times=m))
	theta <- matrix(runif(h*f), nrow=h, ncol=f);
	u <- matrix(rep(0, times=f*m), f, m)
	for(i in 1:1000) {
		u <- w_min(x, y, u, theta)
		theta <- theta_min(x, y, u, theta, lambda)
	}
}

w_min <- function(x, y, u, theta) {
	n <- dim(x)[[2]]
	f <- dim(x)[[1]]
	m <- dim(y)[[1]] 
	newU <- matrix(f, m)
	#w <- matrix(rep(0, times=n*f), f, m)
	#v <- matrix(f, m)
	for(l in 1:m) {
		v <- theta%*%u[,l]
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
				val = val + y[l, i]*x[q, i]
				#should v be v transpose?
				vtimestheta <- v%*%theta
				vtimesthetatimesxi <- vtimestheta%*%x[ ,i]
				val = val - vtimesthetatimesxi*x[q, i]
			}
			values[[q]] = val
		}
		w <- solve(coefficients, values)
		#should theta be theta transpose?
		newU[, l] <- w %+% v%*%theta
	}
	return(newU)
}

theta_min <- function(x, y, u, theta, lambda) {
	f <- dim(x)[[1]]
	m <- dim(u)[[2]]
	h <- dim(theta)[[1]]

	U <- matrix(f, m)
	for(l in 1:m) {
		U[ , l] <- lambda[l]*u[ , l]
	}
	svd_of_U <- svd(U)
  v <- svd_of_u$u
	newTheta <- matrix(0, nrow=h, ncol=f)

	for(k in 1:h) {
		newTheta[k,] <- v[k,]
	}
	return(newTheta)
}
tests <- function() {
	data <- generateData(100, 10, 100, 0.1)
	print(data[[1]])
	print(data[[2]])

}

	
