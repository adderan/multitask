train <- function(h) {
	load("gray-5000.RData")
	source("multitask.R")

	size <- dim(Xu)
	i <- 5
	j <- 6
	#for(i in 1:size[[1]]) {
		#for(j in size[[1]]) {
			if(i == j) next
			ijremoved = Xu[-i,]
			ijremoved = ijremoved[-j,]
			print(dim(ijremoved))
			n <- 1
			X <- replicate(2*n, ijremoved, simplify=FALSE)
			print(length(X))
			print(dim(X[[1]]))
			y <- list(2*n)
			for(nij in 1:n) { 
				y[[nij]] <- Xu[i, nij]
				y[[nij + n]] <- Xu[j, nij]
			}
			min <- joint_min(X, y, h)
			print(min)

						
}
