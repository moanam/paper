# GENERATE HEATMAP & BAR PLOT OF OTU ABSOLUTE ABUNDANCE TRAJECTORIES
# MANOSHI SEN DATTA
# 27 OCTOBER 2016
########################################################################################################################
# Required packages
library(caTools)
library(gplots)
require(lattice)
require(gridExtra)
########################################################################################################################
# Paths to input/output files
file.path = "C:/Users/Moana/Documents/Uni/2016/Publication/"
dat.path = paste(file.path, "bulkPartSea.otus.97.otu_table.headers.txt", sep="")
map.path = paste(file.path, "map_file_all.txt", sep="")
# tax.path = paste(file.path, "Data/16S/output/otu/bulkPartSea.otus.97.nonchimeras.labeled.rdp.cleaned.fasta", sep="")
# abs.path = paste(file.path, "Data/Absolute_abundances/particleTrajs.txt", sep="")
# out.path = paste(file.path, "draft/mainTextFigures/fig2/", sep="")
source(paste(file.path, "fig2_v1_functions.R", sep=""))
# source(paste(file.path, "draft/mainTextFigures/fig2/image.scale.2.r", sep=""))
########################################################################################################################
# Read in variables from files
dat = read.table(dat.path, header=TRUE, stringsAsFactors=FALSE) #Reads in data
row.names(dat) = dat$OTUId #Makes the row names the OTU ids
dat = dat[, -1] #Reomoves column 1 (OTU ID)
map = read.table(map.path, header=FALSE, stringsAsFactors=FALSE) #Reads in data
colnames(map) = c("Sample", "RevBarcode", "FwdBarcode", "SampleType", "Timepoint", "SizeFrac", "Desc") #Changes column names from V1 etc
# tax = read.table(tax.path, header=FALSE, stringsAsFactors=FALSE, sep="\t") #Had to comment out this bit as I was not given the data
# colnames(tax) = c("OTU", "Taxonomy")
# abs = read.table(abs.path, header=TRUE, stringsAsFactors=FALSE)

reps = c("M1", "M2", "M3")
# times = unique(abs$Timepoint) #Not given
########################################################################################################################
# Remove OTUs that correspond to non-bacterial taxa (chloroplast, bacilliariophyta, etc.)
# I can't run this bit as it relies on tax
# OTUs2remove <- NULL
# 
# for (i in 1:nrow(dat)){
#   OTU = row.names(dat)[i]
#   taxonomy = strsplit(tax[which(tax$OTU == OTU), "Taxonomy"], "|", fixed=TRUE)[[1]]
#   
#   if ("Chloroplast" %in% taxonomy){
#     OTUs2remove = c(OTUs2remove, OTU)
#   }
# }
# 
# dat = dat[-which(row.names(dat) %in% OTUs2remove), ]
########################################################################################################################
# Separate data by replicate
for (i in 1:length(reps)){
  replicate = reps[i]
  temp = dat[, map[which(map$Desc == replicate & map$SampleType == "Particles"), "Sample"]]
  assign(paste("dat.", replicate, sep=""), temp)
}
#Makes a file called "temp" (temporary?) with only one replicate and then re-names it. "temp" is left with the last iteration
########################################################################################################################
# Calculate relative abundance data for each replicate
for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.rel = count2rel(temp)
  assign(paste("dat.", replicate, ".rel", sep=""), temp.rel)
}
########################################################################################################################
# Calculate smoothed total abundance trajectory
K = 1.270e+05
P0 = 7.969e-04
r = 2.102e-01

# total.ab.traj.med.smooth = (K * P0 * exp(r * times))/(1 + P0 * (exp(r * times)-1)) # Didn't work because I don't have "times"
########################################################################################################################
# Convert relative abundances to absolute abundances for each replicate
# Missing above
# for (i in 1:length(reps)){
#   replicate = reps[i]
#   temp.rel = get(paste("dat.", replicate, ".rel", sep=""))
#   temp.abs = rel2abs(temp.rel, total.ab.traj.med.smooth)
#   
#   assign(paste("dat.", replicate, ".abs", sep=""), temp.abs)
# }
########################################################################################################################
# Find top OTUs from each replicate individually (by relative abundance)
# Include OTUs present at >1% relative abundance at any timepoint, excluding timepoints at which the total 
# abundance was below the limit of detection for qPCR (not reliable)

include_M1 = c(3:16) #This excludes the early time points, as explained above
include_M2 = c(3:16)
include_M3 = c(3:16)

for (i in 1:length(reps)){
  replicate = reps[i]
  include = get(paste("include_", replicate, sep=""))
  temp.rel = get(paste("dat.", replicate, ".rel", sep=""))[ , include]
  temp.subset.names = names(which(apply(temp.rel, 1, max) > 0.01))
  assign(paste(replicate, ".rel.subset.names", sep=""), temp.subset.names)
} # Makes a lift of the most abundant OTUs in each replicate
########################################################################################################################
# # Calculate smoothed trajectories for all OTUs
# Didn't work due to not having some data
# for (i in 1:length(reps)){
#   replicate = reps[i]
#   temp.abs.subset.names = get(paste(replicate, ".rel.subset.names", sep=""))
#   temp.abs.subset = get(paste("dat.", replicate, ".abs", sep=""))
#   
#   smoothResults = medSmoothOTUs(temp.abs.subset.names, dat.M1.abs, dat.M2.abs, dat.M3.abs)
# #   smoothResults = smoothOTUs(temp.abs.subset.names, temp.abs.subset)
#   temp.abs.subset.smooth = smoothResults[[1]]
#   temp.abs.subset.smooth.norm = smoothResults[[2]]
#   
#   removeNA <- function(dat){
#     return(dat[which(is.na(apply(dat, 1, sum)) == FALSE), ])
#   }
# 
#   temp.abs.subset.smooth.norm.rm = removeNA(temp.abs.subset.smooth.norm)
#   temp.abs.subset.smooth.rm = temp.abs.subset.smooth[row.names(temp.abs.subset.smooth.norm.rm), ]
# 
#   temp.abs.subset.smooth.rm.ordered = temp.abs.subset.smooth.rm[rev(order(apply(temp.abs.subset.smooth.rm, 1, calc_COM))), ]
#   temp.abs.subset.smooth.norm.rm.ordered = temp.abs.subset.smooth.norm.rm[rev(order(apply(temp.abs.subset.smooth.norm.rm, 1, calc_COM))), ]
#   
#   assign(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""), temp.abs.subset.smooth.norm.rm.ordered)
# }
########################################################################################################################
# Plot heatmap
ncolors=599
scaleblackorangered <- colorRampPalette(c("black", "orange", "orangered","red"))(n=ncolors)

for (i in 1:length(reps)){
  replicate = reps[i]
  mat2plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""))
  mat2max = get(paste("dat.", replicate, ".abs", sep=""))[row.names(mat2plot), ]
  col_breaks = c(seq(min(mat2plot, na.rm=TRUE), max(mat2plot, na.rm=TRUE), length=600))
  
  # Prepare for plotting taxonomic distinctions
  famNames = unlist(lapply(row.names(mat2plot), function(OTU) strsplit(tax[which(tax$OTU == OTU), "Taxonomy"], "|", fixed=TRUE)[[1]][4]))
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
  axis(1, at=seq(0,144,by=24), cex.axis=3)
  mtext(text="Colonization time (hours)", side=1, line=6, outer=FALSE, cex=4, font=2)
  
  legend(x=-84, y=49.5, legendOrder, fill=colorScheme, 
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
