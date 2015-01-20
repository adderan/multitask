#!/usr/bin/env Rscript

source("../aso.R")
source("ssl-gray-setup.R")

args <- commandArgs()
datafile <- args[6]
#drug.aux <- args[7]
#auxiliary <- args[8]
load(datafile)


partitioned.data <- partition.data(X, Y)
Y.labels <- partitioned.data$Y.labels
Y.answers <- partitioned.data$Y.answers
cat("h\tDrug\tScore\n") 
for(d in 7:length(args)) {
	drug.target <- args[d]
	drug <- unlist(strsplit(drug.target, split = "[.]"))[[1]]
	auxiliary <- unlist(strsplit(drug.target, split = "[.]"))[[2]]
	bad.data.removed <- remove.bad.data(partitioned.data$X.train, Y.labels, drug)
	drug.response <- bad.data.removed$drug.response
	X.train <- bad.data.removed$X

	
	bad.test.data.removed <- remove.bad.data(partitioned.data$X.test, Y.answers, drug)
	X.test <- bad.test.data.removed$X
	drug.answers <- bad.test.data.removed$drug.response

	Y.aux <- Xu[auxiliary,]
	Xu.reduced <- Xu
	Xu.reduced[auxiliary,] <- 0

	aso.X <- list()
	aso.Y <- list()
	aso.X[[drug]] <- X.train
	aso.Y[[drug]] <- drug.response
	aso.X[[auxiliary]] <- Xu.reduced
	aso.Y[[auxiliary]] <- Y.aux


	h0 <- 5
	hstep <- 5
	for(i in 1:2) {
		aso.model <- aso.train(aso.X, aso.Y, primary = drug, h = h0, ANALYTIC = TRUE)
		aso.predictions <- aso.predict(aso.model, X.test)
		score <- cor(as.vector(aso.predictions), drug.answers, method = "spearman")
		cat(h0, "\t", drug.target, "\t", score, "\n")
		h0 <- h0 + hstep
	}
}
