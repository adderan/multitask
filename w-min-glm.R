library(glmnet)

w.min.glm <- function(X, y, v, theta) {
	#print(length(v))
	#print(dim(theta))
	#print(dim(X))
	y.adjusted <- y - t(v) %*% theta %*% X
	fit <- glmnet(x = t(X), y = y.adjusted, alpha = 0)

	w <- predict(fit, type="coef", s = 1)
	w <- w[-1]
	#print(rownames(w)[1:10])
	#print(rownames(X)[1:10])
	
	return(w)
}


	
