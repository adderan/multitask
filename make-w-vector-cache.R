#X is list of data matrices, y is list of corresponding response vectors

make.w.vector.cache <- function(X, y, lambda, batch.name) {
	library(glmnet)
	n.problems <- length(X)
	stopifnot(n.problems == length(y))

	w.vector.cache <- list()
	for(i in 1:n.problems) {
		problem.name <- names(X)[[i]]
		stopifnot(problem.name == names(y)[[i]])
		X_i <- t(X[[i]])
		y_i <- y[[i]]
		fit <- glmnet(X_i, y_i, alpha = 0)
		w.vector <- predict(fit, type= "coef", s = lambda)
		w.vector <- w.vector[-1]
		w.vector.cache[[problem.name]] <- w.vector
	}
	filename <- cache.filename(batch.name, lambda)
	save(w.vector.cache, file=filename)
	return(w.vector.cache)
}

cache.filename <- function(batch.name, lambda) {
	filename <- "cache/"
	filename <- paste(filename, batch.name, sep="")
	filename <- paste(filename, "_lambda", sep="")
	filename <- paste(filename, lambda, sep="")
	return(filename)
}
		
