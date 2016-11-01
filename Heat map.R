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
colnames(map) = c("Sample", "RevBarcode", "FwdBarcode", "Timepoint")

########################################################################################################################
# Separate data by replicate
dat.M1 = subset(mydata, Replicate == "A")
dat.M2 = subset(mydata, Replicate == "B")


dat.M1 = dat.M1[,c(1,3,12)]
dat.M2 = dat.M2[,c(1,3,12)]

####This sums the OTU abundances from different water masses from the same time
library(dplyr)
df <- group_by(dat.M1, OTU, Time)
dat.M1 <- summarise(df, Abundance = sum(value))

df <- group_by(dat.M2, OTU, Time)
dat.M2 <- summarise(df, Abundance = sum(value))

###################################################
#rearrange to make times columns
library(reshape2)
dat.M1 = dcast(dat.M1, OTU ~ Time) #Rearrages so that each time has a column
dat.M2 = dcast(dat.M2, OTU ~ Time)

##########################
#Make OTU row names
dat.M1 <- data.frame(dat.M1[,-1], row.names=dat.M1[,1])
dat.M2 <- data.frame(dat.M2[,-1], row.names=dat.M2[,1])

#Change column names
colnames(dat.M1) <- c("M1.0","M1.126", "M1.182", "M1.325", "M1.406", "M1.448")
colnames(dat.M2) <- c("M2.0","M2.126", "M2.182", "M2.325", "M2.406", "M2.448")
