svd.approximate.lower.rank <- function(A, k, h) {
	m <- dim(A)[[1]]
	n <- dim(A)[[2]]

	omega <- matrix(runif(n*k, -1, 1), n, k)
	
	Y <- A %*% omega

	qrout <- qr(Y)

	Q <- qr.Q(qrout)

	B <- t(Q) %*% A

	svdout <- svd(B, nu = h)
	uhat <- svdout$u

	U <- Q %*% uhat

	return(U)
}


svd.test <- function() {
	m <- 10
	n <- 20

	A <- matrix(runif(n*m, -1 ,1), m, n)
	return(svd.approximate.lower.rank(A, 5))
}


