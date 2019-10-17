#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ape)

print(args[1])
print(args[2])

tree = read.tree( args[1] )

svg_filename = args[2]

svg(filename=svg_filename)
plot(tree,no.margin=TRUE,edge.width=2)
dev.off()
