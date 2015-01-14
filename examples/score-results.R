#!/usr/bin/env Rscript

input <- file("stdin")
scores <- read.table(file = input, header = TRUE)
print(scores)
n.drugs <- dim(scores)[[1]]
cat("Data for ", n.drugs, " drugs.\n")
cat("ASO base average: ", sum(scores[,3])/n.drugs, "\n")
cat("ASO auxiliary average: ", sum(scores[,4])/n.drugs, "\n")
cat("glmnet average: ", sum(scores[,5])/n.drugs, "\n")
