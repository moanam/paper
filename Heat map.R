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

####This sums the OTU abundances from different water masses from the same time

library(dplyr)

df <- group_by(dat.M1, OTU, Time)

df.summary <- summarise(df, Abundance = sum(value))

library(dplyr)
res <- dat.M1 %>%
  group_by(OTU, Time) %>% 
  mutate(value=sum(value))
dat.M3 = as.data.frame(res)

dat.M4 = aggregate(OTU ~ Time, data=dat.M1, sum)

dat.M4 = aggregate(dat.M1, OTU)



head(dat.M1)

dat.M3 <- data.frame(dat.M1[,-1], row.names=dat.M1[,1]) #This is currently producing duplicate row names as I have each OTU for each time point