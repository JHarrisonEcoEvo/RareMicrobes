#input otu table, sample-well key, and taxonomy file
tax=read.delim("data/16S_sequenceProcessing/combinedTaxonomy97otus.txt", header=F)
otus=read.delim("data/16S_sequenceProcessing/otuTable97otus.txt", header=T)
sampNames=read.csv("data/Sample_Well_Key_includesQuants.csv")
sampNames$pl_well=paste(sampNames$plate, sampNames$well, sep=".")

#now we need to pick which OTUs to keep
#I am first removing any OTUs that the UNITE database said were plants. 
#then I will subset the data to just keep only those OTUs that one or the other database
#had >80% confidence in placement to phylum. In previous work I have found that <80% placement to
#phylum often match nontarget taxa on NCBI...so I remove those too.

fungiOnly= tax[tax[,4]!="d:Plantae",]

#this keeps any row that for either database there is more than 7 characters in the cell
#describing the taxonomic hypothesis.
#this gets rid of anything that was not identified to phylum by one or the other db

fungiOnly = fungiOnly[which(nchar(as.character(fungiOnly[,3])) > 10 | nchar(as.character(fungiOnly[,4])) > 10),]

######################################################################################################
#NOTE: I strongly recommend doing some spot checks of the OTUs that made it through. Find some that were
#not id'd as fungi by both databases and go to the NCBI web interface and paste them in there. 
#if they show up as some other eukaryote then it may be worth scrubbing the data more
######################################################################################################

#overwrite our OTU table so that it includes only the taxon of interest
otus=otus[otus$X.OTU.ID %in% fungiOnly$V1,]

#remove empty columns (no taxa of interest in a sample)
otus=otus[,which(colSums(otus[,2:length(otus)])!=0)]

otus=t(otus)
colnames(otus)=as.character(otus[1,])
otus=otus[-1,]
row.names(otus)=gsub("X(\\d+\\.[a-h]\\d+).*","\\1",row.names(otus))

#check to see if anything in the row.names is not in the key
which(!(row.names(otus) %in% sampNames$pl_well))

otus=data.frame(otus)

#replace well based row names with the sample names
otus$samps=row.names(otus)
newotus=merge(otus, sampNames, by.x="samps", by.y="pl_well")

#make sure to edit the range in the call to the newotus object
write.csv(data.frame(newotus$sampleName, newotus[,grep("Otu",names(newotus))]),file="data/bacteriaOTUtable97otus.csv", row.names = F)
