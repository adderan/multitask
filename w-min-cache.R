
w.min.cache(x, y, y_description, v, theta, w.vector.cache) {
	w.hat.unshifted <- w.vector.cache[[y_description]]
	stopifnot(!is.null(w.hat))
	w.hat <- w.hat.unshifted - 
}

