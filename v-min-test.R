

v.min.test <- function() {
	source("ando.R")
	source("multitask.R")

	mydata <- test1.data()
	x <- (mydata$X.list)[[1]]   ##NxP data matrix
	y <- (mydata$y.list)[[1]]

	print(dim(x))
	print(dim(y))
	
	n <- dim(x)[[1]]
	p <- dim(x)[[2]]
	h <- 10
	
	Theta.hat <- matrix(runif(h*p, 0, 1), h, p)
	V.hat <- c(runif(h, -1, 1))
	W.hat <- w.min(t(x), y, V.hat, Theta.hat)   ##length is P
	V.hat <- v.min(t(x), y, W.hat, Theta.hat)   ##length is h
	cat("length of W.hat = ", length(W.hat), "\n")
	tol <- 1e-5
	
	v1 <- V.hat
	v2 <- V.hat
	for(k in 1:h) v1[k] <- V.hat[k] + tol
	for(k in 1:h) v2[k] <- V.hat[k] - tol
	cat("Length of V.hat: ", length(V.hat), "\n")
	print(V.hat)
	cat("Dimension of Theta.hat: ", dim(Theta.hat), "\n")
	obj0 <- f.obj1(x, y, W.hat, V.hat, Theta.hat)
	obj1 <- f.obj1(x, y, W.hat, v1, Theta.hat)
	obj2 <- f.obj1(x, y, W.hat, v2, Theta.hat)
	
	cat("Optimum objective = ", obj0, "\n")
	cat("Positive perturbation objective = ", obj1, "\n")
	cat("Negative perturbation objective = ", obj2, "\n")
	
	vprime <- v.prime(t(x), y, W.hat, V.hat, Theta.hat)
	print(vprime)

}
