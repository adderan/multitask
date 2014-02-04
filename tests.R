
#calculates the derivative of the error plus regularizer. Used to confirm results of minimization with respect to w. Should be zero after minimization. 
g_prime <- function(x, y, wmatrix, l, v, theta) {
	n <- dim(x)[[2]]
	f <- dim(x)[[1]]
	w <- wmatrix[ , l]
	gprime <- c()
	for(q in 1:f) {
		gprimeq <- 0
		for(i in 1:n) {
			isum <- 0
			isum <- isum + y[l, i]
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
# finds the value of the regularized error for the w minimization. 
find_w_error <- function(x, y, wmatrix, l, u, theta) {
	n <- dim(x)[[2]]
	g <- 0
	v <- theta %*% u[, l]
	w <- wmatrix[,l]
	for(i in 1:n) {
		isum <- 0
		isum <- isum + y[l, i]
		isum <- isum - t(w) %*% x[ , i]
		isum <- isum - t(v) %*% theta %*% x[ , i]
		g <- g + isum*isum
	}
	g <- g + t(w) %*% w
	return(g)
}
find_theta_error <- function(x, y, u, theta, lambda) {
	n <- dim(x)[[2]]
	m <- dim(y)[[1]]
	error <- 0
	base_error <- 0
	for(l in 1:m) {
		v <- theta %*% u[, l]
		#isum <- 0
		#for(i in 1:n) {
		#	square_difference <- 0
		#	square_difference <- square_difference + y[l, i]
		#	square_difference <- square_difference - t(u[ , l]) %*% x[ , l]
		#	square_difference <- square_difference * square_difference
		#	isum <- isum + square_difference
		#}
		regularizer <- lambda[[l]] * vector_magnitude(u[, l] - t(theta) %*% v)
		baseline <- vector_magnitude(t(theta) %*% v)
		base_error <- base_error + baseline
	
		error <- error + regularizer
	}
	#print(base_error)
	#print(error)
	return(error)
}
vector_magnitude <- function(x) {
	magnitude <- 0
	n <- length(x)
	for(i in 1:n) {
		magnitude <- magnitude + x[[i]]*x[[i]]
	}
	return(magnitude)
}
graph_theta_error <- function(x, y, u, theta, lambda) {
	min_error <- find_theta_error(x, y, u, theta, lambda)
	error_dist <- c()
	points <- 10
	h <- dim(theta)[[1]]
	p <- dim(theta)[[2]]

	for(i in 1:points) {
		new_theta <- matrix(rnorm(h*p, mean = 0, sd = 1), h, p)
		#new_theta <- normalize_matrix(new_theta)
		norm_theta <- matrix(1, h, p)
		for(i in 1:h) {
			norm_theta[i, ] <- new_theta[i, ] * sqrt(vector_magnitude(theta[i,]))/sqrt(vector_magnitude(new_theta[i,]))
		}
		 
		permute_theta <- theta[sample(nrow(theta)),]	 
		er <- find_theta_error(x, y, u, norm_theta, lambda)
		error_dist[[i]] <- er
		#print(er)
		#print(norm_theta)
	}
	error_dist <- c(error_dist, min_error)
	label <- c(rep("", times = points), "min")
	dotchart(error_dist, labels = label)
}
		
#graphs the regularized error for the suspected miniumum wmin as well as surrounding w vectors to verify that wmin is the minimum. 
graph_error_dist <- function(x, y, wmatrix, l, v, theta) {
	w <- wmatrix [ , l]
	min_error <- find_w_error(x, y, w, l, v, theta)
	error_dist <- c()
	permuted_error_dist <- c()
	points <- 1000
	for(j in 1:points) {
		delta_w <- c(rnorm(dim(x)[[2]], mean=0, sd = 0.001))
		new_w <- w + delta_w
		append(error_dist, find_w_error(x, y, new_w, l, v, theta))
		delta_w <- c(rnorm(dim(x)[[2]], mean=0, sd=0.005))
		new_w <- w + delta_w
		append(error_dist, find_w_error(x, y, new_w, l, v, theta))
		permuted_error_dist[[j]] <- find_w_error(x, y, sample(w), l, v, theta)
	}
	colors = rep("blue", times = length(error_dist) + 1)
	colors[length(error_dist) + 1] = "red" 
	error_dist <- c(error_dist, min_error)
	permuted_error_dist <- c(permuted_error_dist, min_error)
	dotchart(error_dist, color = colors)
	#dotchart(permuted_error_dist, color = colors)
	print(min_error)


}
	
tests <- function() {
	source("multitask.R")

	data <- generateData(100, 10, 100, 0.1)
	x <- data[[1]]
	y <- data[[2]]
	f <- dim(x)[[1]]
	n <- dim(x)[[2]]
	m <- dim(y)[[1]]
	h <- 20

	l <- 1
	u <- matrix(rep(0, times=f*m), nrow = f, ncol = m)
	theta <- matrix(runif(h*f), nrow = h, ncol = f)
	v <- theta %*% u[ , l]
	w_min_out <- w_min(x, y, u, theta)
	w <- w_min_out[[2]]
	u <- w_min_out[[1]]
	lambda <- c(rep(1, times=m))

	print(find_w_error(x, y, w, l, u, theta))
	gprime <- g_prime(x, y, w, l, v, theta)
	print(gprime)
	#graph_error_dist(x, y, w, l, v, theta)

	
	#min_theta_out <- theta_min(x, y, u, theta, lambda)
	#print(u)
	#print(min_theta_out)
	#print("Random thetas")
	#graph_theta_error(x, y, u, min_theta_out, lambda)
	#print("Real min theta")
	#print(find_theta_error(x, y, u, min_theta_out, lambda))
	
}
joint_min_test <- function() {
	source("multitask.R")
	data <- generateData(100, 10, 100, 0.1)
	x <- data[[1]]
	y <- data[[2]]
	h <- 20
	joint_min(x, y, h)
}







