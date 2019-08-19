require(ggplot2); require(scales); require(reshape2)
d=read.csv('E2.csv')
head(d)
melt(d)