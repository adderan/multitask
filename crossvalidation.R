train <- function(X, Y) {
	load("/soe/sokolov/ssl/gray/gray-5000.RData");
	source("multitask.R")

	size <- dim(Xu)
	for(i in 1:size[[1]]) {
		for(j in size[[1]]) {
			if(i == j) next
			ijremoved = Xu[-i,]
			ijremoved = ijremoved[-j,]
			n <- 100
			X <- list(rep(ijremoved, times=2*n))
			y <- list(n)
			for(nij in 1:n) {
				

			
			
			

}

