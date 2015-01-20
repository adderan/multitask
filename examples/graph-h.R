#!/usr/bin/env Rscript
library(reshape2)
library(ggplot2)
args <- commandArgs()
data <- read.table(args[6], header = TRUE)

#data <- melt(read.table(args[7], header = TRUE), id = "h")
#for(i in 8:length(args)) {
#	h.graph <- melt(read.table(args[i], header = TRUE), id = "h")
#	data <- rbind(data, h.graph)
#}
#colnames(data) <- c("h", "Drug", "Score")
print(data)
ggplot(data, aes(x = h, y = Score, colour = Drug)) + geom_line()
ggsave(file = args[7])
