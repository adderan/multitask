
function <- generateData(f, g, n, epsilon) { #f is total # of features in data, g is number of true predictors, n is number of feature vectors 	
	data <- matrix(runif(f*n), f, n)
	w <- c(rep(1, times=g), rep(0, times=(f-g)))
	y <- matrix(runif(g*n), g, n)

	for(i in 1:(g/2)) {
		for(j in 1:n) {
			y[i,j] = w%*%data[, j]/norm(w) + rnorm(1, mean = 0, sd = epsilon)
		}
	}
	#print(w)
}
norm(x) {
	norm = 0;
	for(i in 1:length(x)) {
		norm += x[i]*x[i]
	}
	return sqrt(norm)
}


joint_min <- function(x, y) { #x is matrix of feature vectors, y is matrix of output vectors, m is number of features. 
	lambda <- c(rep(1, times=n))
	#theta <- matrix(runif(h*p), nrow=h, ncol=p);
	u <- c(rep(0, times=n))
	n <- dim(x)[[2]]
	m <- dim(y)[[1]] 
	for(l in 1:m) {
		v <- theta%*%u


}




	
	

#y <- w%*%x
#print(y)


