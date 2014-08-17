w.min.sgd <- function(X, y, v, theta, obj, time.allowed) {
	p <- dim(X)[[1]]
	start <- proc.time()

	w <- runif(p, -1, 1)
	while(time < time.allowed) {
		time <- proc.time() - start
		w <- adjust.w
	}

	return(index)
}

adjust.w <- function(w) {
	step <- 0.01

	w.new <- w
	index <- sample(1, length(w))
	w.new[[index]] <- w.new[[index]] + step
	return(w.new)
}

	

w.min.analytic.gd <- function(X, y, v, theta, grad, obj, time.allowed) {
	start <- proc.time()
	time <- 0

	alpha <- 0.001
	p <- dim(X)[[1]]
	w <- runif(p, -1, 1)
	print(obj(t(X), y, w, v, theta))
	while(time < time.allowed) {
		time <- proc.time() - start
		w <- w + grad(X, y, w, v, theta) * alpha
		print(obj(t(X), y, w, v, theta))
	}
	return(w)
}

