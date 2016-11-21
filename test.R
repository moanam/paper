
# famNames.collapsed is the file that has the order of the families and determines the order of the colours in the bar
# mat2plot is the matrix of heat map data, with the OTUs in order (41, atm)
mat2plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep=""))
# .abs.subset.smooth.norm.ordered.rm has the OTU names as row names. I can probably rearrange this.
#OTUnames links the OTU with the family! But only for the last replicate. I have a mistake in the code! This means that all the lagends have been for the last replicate.