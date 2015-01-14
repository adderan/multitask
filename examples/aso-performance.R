#!/usr/bin/env Rscript

source("aso.R")
source("examples/ssl-gray-setup.R")
suppressMessages(library(glmnet))

args <- commandArgs()
infile <- args[6]
targetfile <- args[7]

load(infile)

drug.targets <- read.drug.targets(targetfile, X)

partitioned.data <- partition.data(X, Y)
Y.labels <- partitioned.data$Y.labels
Y.answers <- partitioned.data$Y.answers


n.drugs <- dim(Y)[[2]]
n.samples <- dim(Y)[[1]]
cat("Drug\tAuxiliary_Problem\tASO_base_score\tASO_Score\tglmnet_score\n")

h0 <- 10
analytic <- TRUE
for(d in 1:n.drugs) {
	drug <- colnames(Y)[[d]]
	if(drug %in% names(drug.targets) && length(drug.targets[[drug]]) > 0) {
		target <- drug.targets[[drug]][[1]]
	}
	else {
		next
	}
	#print(partitioned.data$Y.labels)

	bad.data.removed <- remove.bad.data(partitioned.data$X.train, Y.labels, drug)
	drug.response <- bad.data.removed$drug.response
	X.train <- bad.data.removed$X

	
	bad.test.data.removed <- remove.bad.data(partitioned.data$X.test, Y.answers, drug)
	X.test <- bad.test.data.removed$X
	drug.answers <- bad.test.data.removed$drug.response

	#print("Predicting with glmnet.")
	glmnet.model <- glmnet(t(X.train), drug.response)
	glmnet.predictions <- predict(glmnet.model, t(X.test), s = 0.1)

	#print("Predicting with ASO, no auxiliary.")
	aso.base.X <- list()
	aso.base.Y <- list()
	aso.base.X[[drug]] <- X.train
	aso.base.Y[[drug]] <- drug.response
	aso.base.model <- aso.train(aso.base.X, aso.base.Y, primary = drug, h = h0, ANALYTIC = analytic)
	aso.base.predictions <- aso.predict(aso.base.model, X.test)
	aso.base.score <- cor(drug.answers, as.vector(aso.base.predictions), method = "spearman")

	aso.X <- list()
	aso.Y <- list()
	aso.X[[drug]] <- X.train
	aso.Y[[drug]] <- drug.response

	Y.aux <- Xu[target,]
	Xu.reduced <- Xu
	Xu.reduced[target,] <- 0
	aso.X[[target]] <- Xu.reduced
	aso.Y[[target]] <- Y.aux
	#for (i in 1:5) {
	#	Y.aux <- Xu[i,]
	#	Xu.reduced <- Xu
	#	Xu.reduced[target,] <- 0
	#	aso.X[[as.character(i)]] <- Xu.reduced
	#	aso.Y[[as.character(i)]] <- Y.aux
	#}
	aso.model <- aso.train(aso.X, aso.Y, primary = drug, h = h0, ANALYTIC = analytic)

	aso.predictions <- aso.predict(aso.model, X.test)
	aso.score <- cor(drug.answers, as.vector(aso.predictions), method="spearman")
	glmnet.score <- cor(drug.answers, as.vector(glmnet.predictions), method="spearman")
	
	cat(drug, "\t", target, "\t", aso.base.score, aso.score, glmnet.score, "\n")
}

