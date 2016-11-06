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
source(paste(file.path, "paper/fig2_v1_functions.R", sep=""))

mydata = read.csv("C:/Users/Moana/Documents/Uni/2016/Publication/paper/Fulltaxonomy.csv", header = T)
tax = mydata[,c(1,34)]
colnames(tax) = c("OTU", "Taxonomy")

library(splitstackshape)

tax = cSplit(tax, 'Taxonomy', sep="_", type.convert=FALSE)
tax$Taxonomy_1 <- NULL
tax$Taxonomy_2 <- NULL
tax$Taxonomy_3 <- NULL
names(tax)[names(tax)=="Taxonomy_4"] <- "Taxonomy"

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

############################################

reps = c("M1", "M2")

# Calculate relative abundance data for each replicate
for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.rel = count2rel(temp)
  assign(paste("dat.", replicate, ".rel", sep=""), temp.rel)
}

############################################
#Lines 62-77

########################################################################################################################
#times = unique(abs$Timepoint)
times = c(0,126, 182, 325, 406, 448)

# Calculate smoothed total abundance trajectory
K = 1.270e+05
P0 = 7.969e-04
r = 2.102e-01

total.ab.traj.med.smooth = (K * P0 * exp(r * times))/(1 + P0 * (exp(r * times)-1))
########################################################################################################################
# Convert relative abundances to absolute abundances for each replicate
for (i in 1:length(reps)){
  replicate = reps[i]
  temp.rel = get(paste("dat.", replicate, ".rel", sep=""))
  temp.abs = rel2abs(temp.rel, total.ab.traj.med.smooth)
  
  assign(paste("dat.", replicate, ".abs", sep=""), temp.abs)
}

##############################################

# Find top OTUs from each replicate individually (by relative abundance)
# Include OTUs present at >1% relative abundance at any timepoint

include_M1 = c(1:6)
include_M2 = c(1:6)

for (i in 1:length(reps)){
  replicate = reps[i]
  include = get(paste("include_", replicate, sep=""))
  temp.rel = get(paste("dat.", replicate, ".rel", sep=""))[ , include]
  temp.subset.names = names(which(apply(temp.rel, 1, max) > 0.01))
  assign(paste(replicate, ".rel.subset.names", sep=""), temp.subset.names)
}

########################################################################################################################
# Calculate smoothed trajectories for all OTUs
for (i in 1:length(reps)){
  replicate = reps[i]
  temp.abs.subset.names = get(paste(replicate, ".rel.subset.names", sep=""))
  temp.abs.subset = get(paste("dat.", replicate, ".abs", sep=""))
  
  smoothResults = medSmoothOTUs(temp.abs.subset.names, dat.M1.abs, dat.M2.abs)
  #   smoothResults = smoothOTUs(temp.abs.subset.names, temp.abs.subset)
  temp.abs.subset.smooth = smoothResults[[1]]
  temp.abs.subset.smooth.norm = smoothResults[[2]]
  
  removeNA <- function(dat){
    return(dat[which(is.na(apply(dat, 1, sum)) == FALSE), ])
  }
  
  temp.abs.subset.smooth.norm.rm = removeNA(temp.abs.subset.smooth.norm)
  temp.abs.subset.smooth.rm = temp.abs.subset.smooth[row.names(temp.abs.subset.smooth.norm.rm), ]
  
  temp.abs.subset.smooth.rm.ordered = temp.abs.subset.smooth.rm[rev(order(apply(temp.abs.subset.smooth.rm, 1, calc_COM))), ]
  temp.abs.subset.smooth.norm.rm.ordered = temp.abs.subset.smooth.norm.rm[rev(order(apply(temp.abs.subset.smooth.norm.rm, 1, calc_COM))), ]
  
  assign(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""), temp.abs.subset.smooth.norm.rm.ordered)
}
########################################################################################################################
library(devtools)
source_url('https://gist.github.com/menugget/7689145/raw/dac746aa322ca4160a5fe66c70fec16ebe26faf9/image.scale.2.r')
#image.scale function

################
# Plot heatmap
ncolors=599
scaleblackorangered <- colorRampPalette(c("black", "orange", "orangered","red"))(n=ncolors)

for (i in 1:length(reps)){
  replicate = reps[i]
  mat2plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""))
  mat2max = get(paste("dat.", replicate, ".abs", sep=""))[row.names(mat2plot), ]
  col_breaks = c(seq(min(mat2plot, na.rm=TRUE), max(mat2plot, na.rm=TRUE), length=600))
  
  # Prepare for plotting taxonomic distinctions
#  tax$Taxonomy <- as.character(tax$Taxonomy)
  #famNames = unlist(lapply(row.names(mat2plot), function(OTU) strsplit(tax[which(tax$OTU == OTU), "Taxonomy"], "|", fixed=TRUE)[[1]][4]))
  famNames <- tax[['Taxonomy']] #This takes all the taxonomy. I think mat2plot has the OTUs in my figure
  
  library(data.table)
  OTUs_in_fig = setDT(mat2plot, keep.rownames = TRUE)[] #Makes what I think are the oTUs a column
  OTUnames = subset(OTUs_in_fig, select=c("rn")) #Isolate the OTUs. But OTUS also gives this!
  
  ##############!!!!!!!!!!!!!!!!!!!!!!!!!!!#########################!!!!!!!!!!!!!!!!!!!!
  #Now: Work out how to add the matching Taxonomy from tax to OTUnames, in the order of OTUnames. Make THIS new thing famNames
  
  famNames = unlist(lapply(famNames, function(fam) ifelse((is.na(fam) | fam == "unclassified"), yes="UATL", no=fam)))
  famNames.collapsed = as.factor(famNames)
  levels(famNames.collapsed) = c(levels(famNames.collapsed), "Other")
  famNamesToKeep = c('Alteromonadales','Flavobacteriales','Oceanospirillales','Sphingobacteriales','Vibrionales','Rhodobacterales')
  
  famNames.collapsed[!(famNames.collapsed %in% famNamesToKeep)] = "Other"
  famNames.collapsed = factor(famNames.collapsed) # clean up unused levels  
  
  # Prepare for plotting bar plot
  OTUs = row.names(mat2max)
  
  max.dat <- as.data.frame(matrix(nrow=length(OTUs), ncol=3))
  for (i in 1:length(OTUs)){
    OTU = OTUs[i]
    
    trajs <- NULL
    for (j in 1:length(reps)){
      rep = reps[j]
      traj = unlist(get(paste("dat.", rep, ".abs", sep=""))[OTU, ])
      trajs = rbind(trajs, traj)
    }
    
    traj.max = apply(trajs, 1, max, na.rm=TRUE)
    traj.max.mean.log = log10(mean(traj.max))
    traj.max.log.sd = sd(log10(traj.max))
    
    max.dat[i, ] = c(OTU, traj.max.mean.log, traj.max.log.sd)
  }
  
  colnames(max.dat) = c("OTU", "max.mean", "max.sd")
  
  traj.max.mean.log.sub = as.numeric(max.dat$max.mean)-1
  
  # Plot heat map
  pdf(file=paste(out.path, "fig2a_", replicate, "_heatmap.pdf", sep=""), height=20, width=14.2)
  colorScheme = c("slategrey", "#a6cee3", "#1f78b4", "mediumpurple1", "#b2df8a", "#33a02c", "black")
  legendOrder = c("Vibrionales", "Alteromonadales", "Oceanospirillales", "Rhodobacterales", "Flavobacteriales", "Sphingobacteriales", "Other")
  
  temp = data.frame(cbind(legendOrder, colorScheme))
  row.names(temp) = temp$legendOrder
  temp.ordered = temp[levels(famNames.collapsed), ]
  colorScheme.ordered = temp.ordered$colorScheme
  
  par(mar=c(6,5,0,1), mgp=c(7,2,0), oma=c(2,15.5,2,0))
  layout(matrix(c(1,2),1,2), widths=c(1.2,7))
  image(as.matrix(t(as.numeric(famNames.collapsed))), xaxt='n', yaxt='n', bty='n',
        breaks=0:length(levels(famNames.collapsed)), 
        col=as.vector(colorScheme.ordered))
  mtext(text="Taxon", side=2, line=1.5, outer=FALSE, cex=4, font=2)
  par(mar=c(6,0,0,2))
  image(x=times, y=seq(1:length(row.names(mat2plot))), z=as.matrix(t(mat2plot)),
        col=scaleblackorangered, breaks=col_breaks,
        xaxt='n', yaxt="n", xlab="", ylab="")
  axis(1, at=seq(0,448,by=64), cex.axis=3)
  mtext(text="Time (days since 28 Jan 2014)", side=1, line=6, outer=FALSE, cex=4, font=2)
  
  legend(x=-330, y=49.5, legendOrder, fill=colorScheme, 
         bty='o', box.lwd=2, xpd=NA, cex=2, title=expression(bold("Orders")))
  dev.off()
  
  # Plot bar plot of maximum values
  at.labels = c(1,2,3,4,5,6)
  pdf(file=paste(out.path, "fig2a_", replicate, "_barPlot.pdf", sep=""), height=30, width=6)
  par(oma=c(0,0,9,0), mar=c(0,2,2,2), xaxs="i")
  mp = barplot(traj.max.mean.log.sub, horiz=TRUE, col="dodgerblue4", las=2, names.arg="", xaxt="n", ylim=c(1,83), space=0, width=1, xlim=c(0, 4))
  plotCI(traj.max.mean.log.sub, mp, gap=0, sfrac=0.002, add=TRUE, cex=0, 
         ui=traj.max.mean.log.sub+as.numeric(max.dat$max.sd), li=traj.max.mean.log.sub-as.numeric(max.dat$max.sd), 
         err="x", col="dodgerblue4", lwd=2)
  axis(side=3, at=c(0,1,2,3,4,5), labels=sapply(at.labels, function(x) as.expression(substitute(list(10^x), list(x=x)))), cex.axis=2.5, line=1)
  abline(v=c(0,1,2,3,4,5), col="dimgrey", lwd=2, lty=2)
  dev.off()
  
  pdf(file=paste(out.path, "fig2a_color_scale_", replicate, ".pdf", sep=""), height=1, width=5)
  par(mar=c(1,1.5,1,1.5), oma=c(0,0,2,0), mgp=c(3,0.5,0), tck=-0.3)
  image.scale(as.matrix(t(mat2plot)), col=scaleblackorangered, breaks=col_breaks)
  #   axis(side=1, labels=FALSE, tck=0)
  axis(side=3, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=2)
  dev.off()
}
