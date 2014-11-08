

#deletes the row corresponding to a gene from the matrices 
#X.unlabeled and X.labeled and returns the matrices and the row that was removed
remove.auxiliary.gene <- function(X.train, X.test, X.unlabeled, gene) {
	aux <- X.unlabeled[gene,]
	index <- which(rownames(X.unlabeled) == gene)
	X.unlabeled.reduced <- X.unlabeled[-index,]
	X.train.reduced <- X.train[-index,]
	X.test.reduced <- X.test[-index,]
	return(list(X.unlabeled = X.unlabeled.reduced, X.train = X.train.reduced, X.test = X.test.reduced, aux = aux))
}

#accepts a dataset consisting of a matrix of (n features, m samples) and a 
#response matrix of (m samples, l problems) and splits it in half by samples.
partition.data <- function(X, Y) {
	m <- dim(X)[[2]]
	split <- floor(m/2)
	X.train <- X[,1:split]
	Y.labels <- Y[1:split,]

	X.test <- X[,(split+1):m]
	Y.answers <- Y[(split+1):m,]
	return(list(X.train = X.train, Y.labels = Y.labels, X.test = X.test, Y.answers = Y.answers))
}

#accepts a matrix of (n features, m samples) and a matrix of (m features, l problems), 
#selects a particular problem from the response matrix, and removes all the samples from the response matrix and the 
#data matrix where the response is NA
remove.bad.data <- function(X, Y, problem) {
	drug.response <- Y[,problem]
	good.samples <- !is.na(drug.response)
	X.good <- X[,good.samples]
	drug.response.good <- drug.response[good.samples]

	return(list(X = X.good, drug.response = drug.response.good))
}

#reads a list of drug targets from a file
read.drug.targets <- function(filename, X) {
	target.table <- read.table(filename, fill = TRUE)

	n.drugs <- dim(target.table)[[1]]

	drug.targets <- list()
	for(i in 1:n.drugs) {
		drug.name <- as.character(target.table[i,1])
		drug.targets[[drug.name]] <- list()
		target.list <- as.character(target.table[i,2])
		target.list <- strsplit(target.list, ",")[[1]]
		n.targets <- length(target.list)
		if(n.targets == 0) next
		for(t in 1:n.targets) {
			target <- target.list[[t]]
			if(!is.null(target) && length(target) != 0 && target %in% rownames(X)) {
				insert.location <- length(drug.targets[[drug.name]]) + 1
				drug.targets[[drug.name]][[insert.location]] <- target
			}
		}
	}
	return(drug.targets)
}
