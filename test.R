
# famNames.collapsed is the file that has the order of the families and determines the order of the colours in the bar
# mat2plot is the matrix of heat map data, with the OTUs in order (41, atm)
mat2plot = get(paste(replicate, ".abs.subset.smooth.norm.ordered.rm", sep="")) # .abs.subset.smooth.norm.ordered.rm are almost identical
#.abs.subset.smooth.norm.ordered.rm is made from "dat.", replicate, ".abs", which are nice and different, like we'd want. So the reason the heat maps all look similar is somewhere in lines 150-172

# .abs.subset.smooth.norm.ordered.rm has the OTU names as row names. I can probably rearrange this.
#M.OTUnames links the OTU with the family