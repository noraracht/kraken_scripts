#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require("phangorn")
require("ape")

tree = read.tree(args[1])
ref = read.tree(args[2])
cat(wRF.dist(tree,ref,normalize=T))
cat("\t")
cat(RF.dist(tree,ref))
cat("\n")
