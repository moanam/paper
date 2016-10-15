rm(list=ls())

```{r Load_libraries_and_data}
library(phyloseq)
library(ggplot2)

otutable_biom_file <- paste("C:/Users/Moana/Documents/Uni/2015/ECOL491/All data/data/", "merged_otu_table.biom", sep = "")
map_file <- paste("C:/Users/Moana/Documents/Uni/2015/ECOL491/All data/data/", "Horizontal_Transect.txt", sep = "")

```


# Now import the .biom-formatted otu_table-tax_table file.
```{r}
biom_otu_tax <- import_biom(otutable_biom_file)

# Add sample data to the dataset using merge
bmsd <- import_qiime_sample_data(map_file)
class(bmsd)
dim(bmsd)#dimensions of bmsd (what is the structure of the map file?) Should match the number of samples and variables
sample_data(bmsd)
biom_otu_tax

#Order levels
bmsd$Month <- factor(bmsd$Month, levels = c("January", "June", "July", "December", "March", "April"))
bmsd$Water_mass_PCA <- factor(bmsd$Water_mass_PCA, levels = c("NW", "STW", "FRONT", "SAW"))

#Merge into phyloseq format
KN_phyloseq <- merge_phyloseq(biom_otu_tax, bmsd)
KN_phyloseq
```

#####Create average result for multiple rarefaction by transforming data using (divide by 10) and check counts per sample
```{r Create average result for multiple rarefaction by transforming data using (divide by 10), results='markup'}
KN_phyloseq = transform_sample_counts(KN_phyloseq, function(x) x/10)
sample_sums(KN_phyloseq)
```

##### Round and confirm count number
```{r Round and confirm count number, results='markup'}
KN_phyloseq = transform_sample_counts(KN_phyloseq, round)
sample_sums(KN_phyloseq)
#check that all numbers are integers
otu_table_df<- as.data.frame(otu_table(KN_phyloseq))
```

Check that all OTUs have representative counts  
For here taxa = OTU  
__Commands interpretation:__  
_Total number of taxa in dataset:_ sum(taxa_sums(KN_phyloseq) > 0)   

_Any taxa with no hits:_ any(taxa_sums(KN_phyloseq)== 0)

```{r identify taxa with only zeros, results='markup', echo=TRUE}
sum(taxa_sums(KN_phyloseq) > 0)
any(taxa_sums(KN_phyloseq)== 0)
sum(taxa_sums(KN_phyloseq) == 0)
any(taxa_sums(KN_phyloseq) > 1)
sum(taxa_sums(KN_phyloseq) > 1)
any(taxa_sums(KN_phyloseq) < 1)
sum(taxa_sums(KN_phyloseq) < 1)
```

#####Prune taxa with less than 1 count and check taxa numbers again
```{r  Save original file and create new file with only present (no zeroes) taxa, results='markup', echo=TRUE}

#Create new file with only present (no zeroes) taxa

KN_phyloseq = prune_taxa(taxa_sums(KN_phyloseq) > 1, KN_phyloseq)
any(sample_sums(KN_phyloseq) == 0)
any(sample_sums(KN_phyloseq) > 0)
sum(taxa_sums(KN_phyloseq) > 0)
any(sample_sums(KN_phyloseq) < 1)
sum(taxa_sums(KN_phyloseq) < 1)
```

#plot across variable
```{r}
sample_variables(KN_phyloseq)
```

#rename Ranks
```{r}
colnames(tax_table(KN_phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

```

#Order levels
```{r}
sample_data(KN_phyloseq)$Water_mass_PCA = factor(sample_data(KN_phyloseq)$Water_mass_PCA, levels = c("NW", "STW", "FRONT", "SAW"))
sample_data(KN_phyloseq)$Month = factor(sample_data(KN_phyloseq)$Month, levels = c("January", "June", "July", "December", "March", "April"))
```

####################################################################################################################

```{r Load_libraries_and_data}

Jan_otutable_biom_file <- paste("C:/Users/Moana/Documents/Uni/2015/ECOL491/All data/data/", "Jan_8000_merged_otu_table.biom", sep = "")
Jan_map_file <- paste("C:/Users/Moana/Documents/Uni/2015/ECOL491/All data/data/", "JanHorizontal_Transect.txt", sep = "")

```


# Now import the .biom-formatted otu_table-tax_table file.
```{r}
Jan_biom_otu_tax <- import_biom(Jan_otutable_biom_file)

# Add sample data to the dataset using merge
Jan_bmsd <- import_qiime_sample_data(Jan_map_file)
class(Jan_bmsd)
dim(Jan_bmsd)#dimensions of Jan_bmsd (what is the structure of the map file?) Should match the number of samples and variables
sample_data(Jan_bmsd)
Jan_biom_otu_tax

#Order levels
Jan_bmsd$Month <- factor(Jan_bmsd$Month, levels = c("January", "June", "July", "December", "March", "April"))
Jan_bmsd$Water_mass_PCA <- factor(Jan_bmsd$Water_mass_PCA, levels = c("NW", "STW", "FRONT", "SAW"))

#Merge into phyloseq format
Jan_phyloseq <- merge_phyloseq(Jan_biom_otu_tax, Jan_bmsd)
Jan_phyloseq
```

#####Create average result for multiple rarefaction by transforming data using (divide by 10) and check counts per sample
```{r Create average result for multiple rarefaction by transforming data using (divide by 10), results='markup'}
Jan_phyloseq = transform_sample_counts(Jan_phyloseq, function(x) x/10)
sample_sums(Jan_phyloseq)
```

##### Round and confirm count number
```{r Round and confirm count number, results='markup'}
Jan_phyloseq = transform_sample_counts(Jan_phyloseq, round)
sample_sums(Jan_phyloseq)
#check that all numbers are integers
otu_table_df<- as.data.frame(otu_table(Jan_phyloseq))
```

Check that all OTUs have representative counts  
For here taxa = OTU  
__Commands interpretation:__  
_Total number of taxa in dataset:_ sum(taxa_sums(Jan_phyloseq) > 0)   

_Any taxa with no hits:_ any(taxa_sums(Jan_phyloseq)== 0)

```{r identify taxa with only zeros, results='markup', echo=TRUE}
sum(taxa_sums(Jan_phyloseq) > 0)
any(taxa_sums(Jan_phyloseq)== 0)
sum(taxa_sums(Jan_phyloseq) == 0)
any(taxa_sums(Jan_phyloseq) > 1)
sum(taxa_sums(Jan_phyloseq) > 1)
any(taxa_sums(Jan_phyloseq) < 1)
sum(taxa_sums(Jan_phyloseq) < 1)
```

#####Prune taxa with less than 1 count and check taxa numbers again
```{r  Save original file and create new file with only present (no zeroes) taxa, results='markup', echo=TRUE}

#Create new file with only present (no zeroes) taxa

Jan_phyloseq = prune_taxa(taxa_sums(Jan_phyloseq) > 1, Jan_phyloseq)
any(sample_sums(Jan_phyloseq) == 0)
any(sample_sums(Jan_phyloseq) > 0)
sum(taxa_sums(Jan_phyloseq) > 0)
any(sample_sums(Jan_phyloseq) < 1)
sum(taxa_sums(Jan_phyloseq) < 1)
```

#plot across variable
```{r}
sample_variables(Jan_phyloseq)
```

#rename Ranks
```{r}
colnames(tax_table(Jan_phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

```

#Order levels
```{r}
sample_data(Jan_phyloseq)$Water_mass_PCA = factor(sample_data(Jan_phyloseq)$Water_mass_PCA, levels = c("NW", "STW", "FRONT", "SAW"))
sample_data(Jan_phyloseq)$Month = factor(sample_data(Jan_phyloseq)$Month, levels = c("January", "June", "July", "December", "March", "April"))
```

####################################################################################################################

#Merge files
```{r}
complete_phyloseq = merge_phyloseq(KN_phyloseq, Jan_phyloseq)
```


#subset samples by experiment and calculate alpha
```{r}
Exp1_phyloseq = complete_phyloseq
#Exp1_phyloseq = subset_samples(complete_phyloseq, Exp1 == "yes")
alpha_summary_Exp1 <- estimate_richness(Exp1_phyloseq, measures = c("Observed", "Shannon")) #gives a tabel of important numbers for when you are writing up

```

##########################################################################################################################

#Make OTU data frame
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85,1)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data
melt.phylum = melt(pOTU.phylum,id.vars="OTU") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"
agg.phylum=aggregate(.~OTU+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ OTU) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_OTUMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_OTUMetadataLayout2.csv")

```

#Make Genus data frame
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85,87, 91)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

pOTU.phylum = subset(pOTU.phylum, Phylum == "D_1__Actinobacteria") #Subset to only one phylum
pOTU.phylum$Phylum <- NULL

melt.phylum = melt(pOTU.phylum,id.vars="Genus") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Genus', sep="_", type.convert=FALSE)
melt.phylum$Genus_1 <- NULL
melt.phylum$Genus_2 <- NULL
melt.phylum$Genus_3 <- NULL
melt.phylum$Genus_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Genus_4"] <- "Genus"

agg.phylum=aggregate(.~Genus+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Genus) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_ActinobacteriaMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_ActinobacteriaMetadataLayout2.csv")

```

#Make Phylum data frame
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85,87)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

melt.phylum = melt(pOTU.phylum,id.vars="Phylum") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Phylum', sep="_", type.convert=FALSE)
melt.phylum$Phylum_1 <- NULL
melt.phylum$Phylum_2 <- NULL
melt.phylum$Phylum_3 <- NULL
melt.phylum$Phylum_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Phylum_4"] <- "Phylum"

agg.phylum=aggregate(.~Phylum+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Phylum) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_PhylaMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_PhylaMetadataLayout2.csv")

```

#Make Genus data frame (all phyla)
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85, 91)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

melt.phylum = melt(pOTU.phylum,id.vars="Genus") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Genus', sep="_", type.convert=FALSE)
melt.phylum$Genus_1 <- NULL
melt.phylum$Genus_2 <- NULL
melt.phylum$Genus_3 <- NULL
melt.phylum$Genus_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Genus_4"] <- "Genus"

agg.phylum=aggregate(.~Genus+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Genus) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_GeneraMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_GeneraMetadataLayout2.csv")

```


#Make Order data frame (all phyla)
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85, 89)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

melt.phylum = melt(pOTU.phylum,id.vars="Order") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Order', sep="_", type.convert=FALSE)
melt.phylum$Order_1 <- NULL
melt.phylum$Order_2 <- NULL
melt.phylum$Order_3 <- NULL
melt.phylum$Order_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Order_4"] <- "Order"

agg.phylum=aggregate(.~Order+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Order) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_OrderMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_OrderMetadataLayout2.csv")

```



#Make Family data frame (all phyla)
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85, 90)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

melt.phylum = melt(pOTU.phylum,id.vars="Family") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Family', sep="_", type.convert=FALSE)
melt.phylum$Family_1 <- NULL
melt.phylum$Family_2 <- NULL
melt.phylum$Family_3 <- NULL
melt.phylum$Family_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Family_4"] <- "Family"

agg.phylum=aggregate(.~Family+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Family) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_FamilyMetadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_FamilyMetadataLayout2.csv")

```



#Make Order data frame
```{r}
library(ggplot2)
library(plyr)
library(scales)
library(reshape)
library(RColorBrewer)
library(grid)
library(phyloseq)

speciesList <- tapply(sample_names(Exp1_phyloseq), get_variable(Exp1_phyloseq, "SampleID"), c) #Makes columns of water mass types and says which samles are found there
speciesPhyseq <- lapply(speciesList, prune_samples, Exp1_phyloseq) #only samples that have these species will be kept
speciesOTUtable <- lapply(speciesPhyseq,otu_table) #go through data and extract the information that matches that data
speciesAvg <- lapply(speciesOTUtable,rowMeans) #what does each species contribute to the 100%? percent transformation
pooledOTUtable = t(do.call(rbind,speciesAvg)) #Take the new average and bind it to original document
#This is as far back as I could trace the x
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable) #converts data to a dataframe (note that this is subsetted!!)
TT = tax_table(Exp1_phyloseq) #attaches taxonomy to those species
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Exp1_phyloseq)) #attaches OTUs to data frame
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:85])) #numbers show in which columns treatment data is located
pOTU.phylum =pOTU[,c(2:85,89, 90)] #need to change this and the above numbers manually. The 6 refers to the column of taxonomy data

pOTU.phylum = subset(pOTU.phylum, Order == "D_3__SAR11 clade") #Subset to only one phylum
pOTU.phylum$Order <- NULL

melt.phylum = melt(pOTU.phylum,id.vars="Family") #melt all the OTUs that are one particular thing into a rank. Still has many observations for each OTU in each water mass
colnames(melt.phylum)[2]="SampleID"

library(splitstackshape)

melt.phylum = cSplit(melt.phylum, 'Family', sep="_", type.convert=FALSE)
melt.phylum$Family_1 <- NULL
melt.phylum$Family_2 <- NULL
melt.phylum$Family_3 <- NULL
melt.phylum$Family_5 <- NULL
names(melt.phylum)[names(melt.phylum)=="Family_4"] <- "Family"

agg.phylum=aggregate(.~Family+SampleID,melt.phylum,sum) #Turns it into percentages
agg.phylum$SampleID = gsub( "X", "", paste(agg.phylum$SampleID)) ##########THIS IS MY FIX!!!!!!#######################

library(reshape2)
data = dcast(agg.phylum, SampleID ~ Family) #Rearrages so that each species has a column


PhylumCount = data
Metadata = rbind(bmsd, Jan_bmsd)
NewDF <- merge(PhylumCount, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF$Chlorophylla <- NULL
NewDF$BA <- NULL
NewDF$BP <- NULL
NewDF <- na.omit(NewDF)

#write.csv(NewDF, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_SAR11Metadata.csv")


PhylumCount2 = agg.phylum
NewDF2 <- merge(PhylumCount2, Metadata, by="SampleID", all.x=TRUE, all.y=TRUE)
NewDF2$Chlorophylla <- NULL
NewDF2$BA <- NULL
NewDF2$BP <- NULL
NewDF2 <- na.omit(NewDF2)

#write.csv(NewDF2, file = "C://Users//Moana//Documents//Uni//2015//ECOL491//All data//data//S_SAR11MetadataLayout2.csv")

```

Oder 89 SAR11 clade
Family 90
