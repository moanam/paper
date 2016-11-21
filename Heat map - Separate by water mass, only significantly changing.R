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
source(paste(file.path, "paper/fig2_v1_functions-4 replicates.R", sep=""))

mydata = read.csv("C:/Users/Moana/Documents/Uni/2016/Publication/paper/Fulltaxonomy.csv", header = T)
tax = mydata[,c(1,34)]
colnames(tax) = c("OTU", "Taxonomy")

library(splitstackshape)

tax = cSplit(tax, 'Taxonomy', sep="_", type.convert=FALSE)
tax$Taxonomy_1 <- NULL
tax$Taxonomy_2 <- NULL
tax$Taxonomy_3 <- NULL
names(tax)[names(tax)=="Taxonomy_4"] <- "Taxonomy"

map = mydata[,c(3,15,16,12)]
colnames(map) = c("Sample", "RevBarcode", "FwdBarcode", "Timepoint")

#Significantly changing OTUs
Significance = read.delim("C:/Users/Moana/Documents/Uni/2016/Publication/Andrew Gray/results 29-04-2016 HB.tab")
Changing = subset(Significance, tsignificant=="1")

#Subsets to only significant OTUs
selectedRows <- (mydata$OTU %in% Changing$otu)
Changingdf <- mydata[selectedRows,]

########################################################################################################################
# Separate data by water mass
dat.M1 = subset(Changingdf, Water_mass_PCA == "NW")
dat.M2 = subset(Changingdf, Water_mass_PCA == "STW")
dat.M3 = subset(Changingdf, Water_mass_PCA == "FRONT")
dat.M4 = subset(Changingdf, Water_mass_PCA == "SAW")

reps = c("M1", "M2", "M3", "M4")

for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.dat = temp[,c(1,3,12)]
  assign(paste("dat.", replicate, sep=""), temp.dat)
}

####This sums the OTU abundances from the replicates and different samples in the same water mass from the same time
library(dplyr)
for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  df <- group_by(temp, OTU, Time)
  temp.dat <- summarise(df, Abundance = sum(value))
  assign(paste("dat.", replicate, sep=""), temp.dat)
}

###################################################
#rearrange to make times columns
library(reshape2)
#Rearrages so that each time has a column

for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.dat = dcast(temp, OTU ~ Time)
  assign(paste("dat.", replicate, sep=""), temp.dat)
}

#FRONT is missing the second time point (the front was really narrow then)

##########################
#Make OTU row names

for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.dat = data.frame(temp[,-1], row.names=temp[,1])
  assign(paste("dat.", replicate, sep=""), temp.dat)
}

#Change column names
colnames(dat.M1) <- c("M1.0","M1.126", "M1.182", "M1.325", "M1.406", "M1.448")
colnames(dat.M2) <- c("M2.0","M2.126", "M2.182", "M2.325", "M2.406", "M2.448")
colnames(dat.M3) <- c("M3.0","M3.182", "M3.325", "M3.406", "M3.448")
colnames(dat.M4) <- c("M4.0","M4.126", "M4.182", "M4.325", "M4.406", "M4.448")

############################################

# Calculate relative abundance data for each replicate
for (i in 1:length(reps)){
  replicate = reps[i]
  temp = get(paste("dat.", replicate, sep=""))
  temp.rel = count2rel(temp)
  assign(paste("dat.", replicate, ".rel", sep=""), temp.rel)
}

########################################################################################################################

times = c(0,126, 182, 325, 406, 448)

# ########################################################################################################################
# Convert relative abundances to absolute abundances for each replicate
#Workaround to just name what we had "absolute abundances", because I'm not using the trajectory
dat.M1.abs=dat.M1
dat.M2.abs=dat.M2
dat.M3.abs=dat.M3
dat.M4.abs=dat.M4

##############################################

# Find top OTUs from each replicate individually (by relative abundance)
# Include OTUs present at >1% relative abundance at any timepoint

include_M1 = c(1:6)
include_M2 = c(1:6)
include_M3 = c(1:5)
include_M4 = c(1:6)

for (i in 1:length(reps)){
  replicate = reps[i]
  include = get(paste("include_", replicate, sep=""))
  temp.rel = get(paste("dat.", replicate, ".rel", sep=""))[ , include]
  temp.subset.names = names(which(apply(temp.rel, 1, max) > 0.0))
  assign(paste(replicate, ".rel.subset.names", sep=""), temp.subset.names)
}

########################################################################################################################
# Calculate smoothed trajectories for all OTUs and do the heatmap taxonomy. This used to be two different bits but something weird happnened.

library(devtools)
source_url('https://gist.github.com/menugget/7689145/raw/dac746aa322ca4160a5fe66c70fec16ebe26faf9/image.scale.2.r')
#image.scale function

################
#Heatmap taxonomy
library(data.table)

###############

for (i in 1:length(reps)){
  replicate = reps[i]
  temp.abs.subset.names = get(paste(replicate, ".rel.subset.names", sep="")) #Names of all the oTUs to be used in the heat map for that figure. Includes the changing OTUs across everything, but only those actually found in the water mass
  temp.abs.subset = get(paste("dat.", replicate, ".abs", sep="")) #OTUs and abundances (table)
  
  smoothResults = medSmoothOTUs(temp.abs.subset.names, dat.M1.abs, dat.M2.abs, dat.M3.abs, dat.M4.abs) #medSmoothOTUs is Manoshi's function to calculate median, smoothed trajectories for all OTUs. Maybe this is why the plots are so similar and I should just call the one replicate I'm plotting each time?
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
  
  OTUs_in_fig = setDT(temp.abs.subset.smooth.norm.rm.ordered, keep.rownames = TRUE)[] #Makes what I think are the oTUs a column
  OTUnamesused = subset(OTUs_in_fig, select=c("rn")) #Isolate the OTUs. But OTUS also gives this!
  
  OTUnamesused$Taxonomy <- tax$Taxonomy[match(OTUnamesused$rn,tax$OTU)]
  assign(paste(replicate, ".OTUnames", sep=""), OTUnamesused)
  
}
########################################################################################################################
# library(devtools)
# source_url('https://gist.github.com/menugget/7689145/raw/dac746aa322ca4160a5fe66c70fec16ebe26faf9/image.scale.2.r')
# #image.scale function
# 
# ################
# #Heatmap taxonomy
# library(data.table)
# for (i in 1:length(reps)){
#   mat3plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep="")) #This only has the OTUs as row names for M1-M3
#   OTUs_in_fig = setDT(mat3plot, keep.rownames = TRUE)[] #Makes what I think are the oTUs a column
#   OTUnamesused = subset(OTUs_in_fig, select=c("rn")) #Isolate the OTUs. But OTUS also gives this!
# 
#   OTUnamesused$Taxonomy <- tax$Taxonomy[match(OTUnamesused$rn,tax$OTU)]
#   assign(paste(replicate, ".OTUnames", sep=""), OTUnamesused)
#   
# }

################
# Plot heatmap
ncolors=599
scaleblackorangered <- colorRampPalette(c("black", "orange", "orangered","red"))(n=ncolors)

for (i in 1:length(reps)){
  replicate = reps[i]
  mat2plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""))
  mat2max = get(paste("dat.", replicate, ".abs", sep=""))#[row.names(mat2plot), ]
  col_breaks = c(seq(min(mat2plot, na.rm=TRUE), max(mat2plot, na.rm=TRUE), length=600))
  
  # Prepare for plotting taxonomic distinctions
  OTUnames <- get(paste(replicate, ".OTUnames", sep=""))
  famNames <- OTUnames[['Taxonomy']]
  famNames = unlist(lapply(famNames, function(fam) ifelse((is.na(fam) | fam == "unclassified"), yes="UATL", no=fam)))
  famNames.collapsed = as.factor(famNames)
  levels(famNames.collapsed) = c(levels(famNames.collapsed), "Other")
  #famNamesToKeep = c('Alteromonadales','Flavobacteriales','Oceanospirillales','Sphingobacteriales','Vibrionales','Rhodobacterales')
  famNamesToKeep = c("Flavobacteriales", "Rhodobacterales", "Oceanospirillales", "Alteromonadales", "Rhodospirillales", "SAR11 clade", "Other")
  
  famNames.collapsed[!(famNames.collapsed %in% famNamesToKeep)] = "Other"
  famNames.collapsed = factor(famNames.collapsed) # clean up unused levels  
  
  # Prepare for plotting bar plot
  OTUs = row.names(mat2max)
  
  max.dat <- as.data.frame(matrix(nrow=length(OTUs), ncol=3)) #OTUs, mean, sd
  
  #trajs <- NULL
  
  #Get only the OTUs that I want for the bar graph. I should make this in to a for loop.
#   dat.M1.abs2 = setDT(dat.M1.abs, keep.rownames = TRUE)[]#Makes the row names a column so I can subest only the OTUs that I want
#   dat.M1.abs3 = dat.M1.abs2[dat.M1.abs2$rn %in% mat3plot$rn, ]#Subsets only the most abundant OTUs
#   rownames(dat.M1.abs3) <- dat.M1.abs3$rn#Makes the oTU names the row names
#   dat.M1.abs3$rn <- NULL#Remove the OTU column
#   
#   dat.M2.abs2 = setDT(dat.M2.abs, keep.rownames = TRUE)[]#Makes the row names a column so I can subest only the OTUs that I want
#   dat.M2.abs3 = dat.M2.abs2[dat.M2.abs2$rn %in% mat3plot$rn, ]#Subsets only the most abundant OTUs
#   rownames(dat.M1.abs3) <- dat.M2.abs3$rn#Makes the oTU names the row names
#   dat.M2.abs3$rn <- NULL#Remove the OTU column
#   
  
  for (i in 1:3){#length(OTUs)
#     temp.OTU=get(paste("dat.", replicate, ".abs", sep=""))
#     temp.OTU=setDT(temp.OTU, keep.rownames = TRUE)[]
    OTU = Changing$otu[i]
    # OTU = OTUs[i] #I should just change THIS to the oTUs that I want!!
    
    trajs <- NULL
    for (j in 1:3){#length(reps)
      rep = reps[j]
      traj = unlist(get(paste("dat.", rep, ".abs", sep=""))[OTU, ]) #Get gives the name of an object. Type dat.M1.abs[OTU, ] to see where this gets its thing from. This gives the trajectory of each OTU (it overwrites the previous one)
      trajs = rbind(trajs, traj)
    }
    
    traj.max = apply(trajs, 1, max, na.rm=TRUE) #trajs = the object used. 1 = indicates rows. max = the function applied. Gives the maximum value in each row - maximum abundance in each water mass, for each OTU
    traj.max.mean.log = log10(mean(traj.max)) #Calculates the log10 of the mean of the maximum values - log(1)=0, log(0.5)=-0.3, log(0)=-inf (probably from the picking and rareification)
    traj.max.log.sd = sd(log10(traj.max)) #Calclutes the sd of the above. Gives NA, I think if it was only found in one replicate
    
    max.dat[i, ] = c(OTU, traj.max.mean.log, traj.max.log.sd) #This gives the log of the means of the 4 maximum abundances (1 in each water mass), so of course it is the same for each water mass!
  }
  
  colnames(max.dat) = c("OTU", "max.mean", "max.sd")
  
  traj.max.mean.log.sub = as.numeric(max.dat$max.mean)-1
  
  # Plot heat map
  pdf(file=paste(out.path, "fig2a_", replicate, "_heatmap.pdf", sep=""), height=20, width=14.2)
  colorScheme = c("slategrey", "#a6cee3", "#1f78b4", "mediumpurple1", "#b2df8a", "#33a02c", "black")
  legendOrder = c("Flavobacteriales", "Rhodobacterales", "Oceanospirillales", "Alteromonadales", "Rhodospirillales", "SAR11 clade", "Other")
  #legendOrder = c("Vibrionales", "Alteromonadales", "Oceanospirillales", "Rhodobacterales", "Flavobacteriales", "Sphingobacteriales", "Other")
  
    #Legend
  temp = data.frame(cbind(legendOrder, colorScheme))
  row.names(temp) = temp$legendOrder
  temp.ordered = temp[levels(famNames.collapsed), ]
  colorScheme.ordered = temp.ordered$colorScheme
  
    #Graph layout
  par(mar=c(6,5,0,1), mgp=c(7,2,0), oma=c(2,15.5,2,0)) #Par=Parameters. Mar=margin, (bottom, left, top, right) gives the number of lines of margin. mgp=margin line for the axis title, axis labels and axis line. Oma=(bottom, left, top, right) outer margins in lines of text
  layout(matrix(c(1,2),1,2), widths=c(1.2,7)) #matrix=Divides the graph in to 1 row and 2 columns. Figure 1 gets column 1 and figure 2 gets column 2. widths= relative witdths of the columns
    #This makes the taxonomy bar!
  image(as.matrix(t(as.numeric(famNames.collapsed))), xaxt='n', yaxt='n', bty='n',
        breaks=0:length(levels(famNames.collapsed)), 
        col=as.vector(colorScheme.ordered)) #Colours to use
  mtext(text="Taxon", side=2, line=1.5, outer=FALSE, cex=4, font=2)
  par(mar=c(6,0,0,2))
    #This makes the heat map
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
  mp = barplot(traj.max.mean.log.sub, horiz=TRUE, col="dodgerblue4", las=2, names.arg="", xaxt="n", ylim=c(1,54), space=0, width=1, xlim=c(0, 4))
  plotCI(traj.max.mean.log.sub, mp, gap=0, sfrac=0.002, add=TRUE, cex=0, 
         ui=traj.max.mean.log.sub+as.numeric(max.dat$max.sd), li=traj.max.mean.log.sub-as.numeric(max.dat$max.sd), 
         err="x", col="dodgerblue4", lwd=2)
  axis(side=3, at=c(0,1,2,3,4,5), labels=sapply(at.labels, function(x) as.expression(substitute(list(10^x), list(x=x)))), cex.axis=2.5, line=1)
  abline(v=c(0,1,2,3,4,5), col="dimgrey", lwd=2, lty=2)
  dev.off()
  
    #Colour scale
  pdf(file=paste(out.path, "fig2a_color_scale_", replicate, ".pdf", sep=""), height=1, width=5)
  par(mar=c(1,1.5,1,1.5), oma=c(0,0,2,0), mgp=c(3,0.5,0), tck=-0.3)
  image.scale(as.matrix(t(mat2plot)), col=scaleblackorangered, breaks=col_breaks)
  #   axis(side=1, labels=FALSE, tck=0)
  axis(side=3, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=2)
  dev.off()
}
