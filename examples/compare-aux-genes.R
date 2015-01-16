#!/usr/bin/env Rscript
source("aso.R")
source("examples/ssl-gray-setup.R")
args <- commandArgs()
infile <- args[6]
drug <- args[7]
target <- args[8]
n.aux <- args[9]
load(infile)
partitioned.data <- partition.data(X, Y)
Y.labels <- partitioned.data$Y.labels
Y.answers <- partitioned.data$Y.answers

bad.data.removed <- remove.bad.data(partitioned.data$X.train, Y.labels, drug)
drug.response <- bad.data.removed$drug.response
X.train <- bad.data.removed$X

	
bad.test.data.removed <- remove.bad.data(partitioned.data$X.test, Y.answers, drug)
X.test <- bad.test.data.removed$X
drug.answers <- bad.test.data.removed$drug.response

n.genes <- dim(X.train)[[1]]

scores <- list()
#start with the auxiliary problem as the drug target gene
aux.gene <- target
for(i in 1:n.aux) {
	Y.aux <- Xu[aux.gene,]
	Xu.reduced <- Xu
	Xu.reduced[aux.gene,] <- 0
	
	aso.X <- list()
	aso.Y <- list()
	
	aso.X[[drug]] <- X.train
	aso.Y[[drug]] <- drug.response
	aso.X[[aux.gene]] <- Xu.reduced
	aso.Y[[aux.gene]] <- Y.aux
	
	aso.model <- aso.train(aso.X, aso.Y, primary = drug)
	aso.predictions <- aso.predict(aso.model, X.test)
	aso.score <- cor(as.vector(aso.predictions), drug.answers, method = "spearman")
	scores[[aux.gene]] <- aso.score
	cat(aux.gene, "\t", aso.score, "\n")
	#set auxiliary problem to a random gene
	aux.gene <- rownames(X.train)[[sample(n.genes,1)]]
}
#print(scores)
