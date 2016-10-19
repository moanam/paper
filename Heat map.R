rm(list=ls())

# Required packages
library(caTools)
library(gplots)
require(lattice)
require(gridExtra)

# Paths to input/output files
file.path = "C:/Users/Moana/Documents/Uni/2016/Publication/"
tax.path = paste(file.path, "Fasta file/outseqs_by_sample_Moana.fasta", sep="")
# abs.path = paste(file.path, "Data/Absolute_abundances/particleTrajs.txt", sep="")
out.path = paste(file.path, "Output/", sep="")
source(paste(file.path, "fig2_v1_functions.R", sep=""))

mydata = read.csv("C:/Users/Moana/Documents/Uni/2016/Publication/paper/Fulltaxonomy.csv", header = T)
tax = mydata[,c(2,35)]
colnames(tax) = c("OTU", "Taxonomy")

map = mydata[,c(3,15,16,13)]