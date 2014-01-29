
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
# l indexes into m, i indexes into n
w_min <- function(x, y, u, theta) {
	n <- dim(x)[[2]]
	f <- dim(x)[[1]]
	m <- dim(y)[[1]] 
	newU <- matrix(0, nrow = f, ncol = m)
	wmatrix <- matrix(0, nrow = f, ncol = m)
	
	for(l in 1:m) {
		v <- theta %*% cbind(u[,l])     #this is v subscript l
		print(length(v))
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
				#vtimestheta <- rbind(v)%*%theta
				#vtimesthetatimesxi <- vtimestheta%*%x[ ,i]
				thetatimesxi <- theta %*% x[ , i]
				vtimesthetatimesxi <- t(v) %*% thetatimesxi

				val = val - vtimesthetatimesxi*x[q, i]
			}
			values[[q]] = val
		}
		w <- solve(coefficients, values)   #this is the minimum w for problem l
		#print(w)

		wmatrix[,l] <- w
		#should theta be theta transpose?
		newU[, l] <- w + t(theta) %*% v
	}
	

	return(list(newU, wmatrix))
}

theta_min <- function(x, y, u, theta, lambda) {
	f <- dim(x)[[1]]
	m <- dim(u)[[2]]
	h <- dim(theta)[[1]]

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

