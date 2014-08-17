w.min.sgd <- function(X, y, v, theta, time.allowed) {
	alpha <- 0.1

	p <- dim(X)[[1]]
	start <- proc.time()
	time <- 0

	w <- runif(p, -1, 1)
	obj0 <- s.obj(t(X), y, w, v, theta)
	print(obj0)
	prev.obj <- obj0
	obj.new <- obj0
	while(time < time.allowed) {
		prev.obj <- obj.new
		time <- proc.time() - start
		w <- w - alpha * numerical.gradient(X, y, w, v, theta, s.obj)
		obj.new <- s.obj(t(X), y, w, v, theta)
		print(obj.new)
		if(obj.new >= prev.obj) {
			alpha <- alpha/10
		}
		else {
			alpha <- alpha * 2
		}
	}
	return(w)

}
s.obj <- function(x, y, w, v, theta) {
	sample.size <- 5
	n <- dim(x)[[1]]
	start <- sample(n,1)
	stop <- min(start + sample.size, n)

	short.x <- x[start:stop,]
	short.y <- y[start:stop]

	z <- short.x %*% w + short.x %*% t(theta) %*% v
	loss <- mean((short.y - z)^2)
	reg <- sum(w^2)
	return(loss + reg)
}


numerical.gradient <- function(X, y, w, v, theta, obj) {
	step <- 0.01

	obj0 <- obj(t(X), y, w, v, theta)
	gradient <- c()
	for(i in 1:length(w)) {
		new.w <- w
		new.w[i] <- new.w[i] + step
		objnew <- obj(t(X), y, new.w, v, theta)
		gradient[i] <- (objnew - obj0)/step
	}
	return(gradient)


	
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
