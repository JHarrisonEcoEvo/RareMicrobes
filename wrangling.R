rm(list=ls())
plantorder = read.csv("./data/planting_order.csv")
table(plantorder$treament)

#Add sla and leaf area
sla = read.csv("./data/sla.csv")

alldat <- merge(sla, plantorder, by.x = "Plant", by.y="rns", all=T)
dim(alldat)
dim(plantorder)

alldat$leafarea_avg = rowMeans(data.frame(alldat[,5:7], na.rm=T))

alldat$sla_1 = alldat$a1/alldat$m1
alldat$sla_2 = alldat$a2/alldat$m2
alldat$sla_3 = alldat$a3/alldat$m3
alldat$sla_4 = alldat$a4/alldat$m4

alldat$sla_avg = rowMeans(alldat[,27:30],na.rm=T)

#Add CN data
cn <- read.csv("./data/CNmeasurements.csv")
key <- read.csv("./data/CNkey.csv")
names(key)[which(names(key)=="Sample")] <- "Sample2"

#first replace duplicate samples with their mean (tech.replicates)
key$well <- paste("T", key$Tray, "-", key$Row., sep="")
cn$Sample <- as.character(cn$Sample)
dat <- merge(cn, key, by.x="Sample", by.y = "well", all = T)
dat$percN <- dat$wt...N/dat$Weight.mg.
dat$percC <- dat$wt...C/dat$Weight.mg.

dupes <- dat[grep("*[ab]", dat$Sample2),]

dupes$samp <- gsub("(\\d+)[ab]","\\1", dupes$Sample2)

d15N<- NA
whichDupe <- NA
l<- 1
c13<- NA
percC<- NA
percN<- NA
for(i in unique(dupes$samp)){
  d15N[l] <- mean(dupes[which(dupes$samp==i),2])
  whichDupe[l]<- i
  c13[l] <- mean(dupes[which(dupes$samp==i),4])
  percN[l] <- mean(dupes[which(dupes$samp==i),13])
  percC[l] <- mean(dupes[which(dupes$samp==i),14])
  l<-l+1
}

derepDupes <- data.frame(whichDupe, NA, d15N, NA, c13, NA, percN, percC ) #NAs are so name vector matches, see blw

#Remove dupes from orig data frame and add the derep data back
dat = dat[-grep("*[ab]", as.character(dat$Sample2)),] #remove dupes
dat = subset(dat, dat$SampleType != "Extra")
dat2 = dat[,c(9,10, 2:5, 13, 14)]
names(derepDupes) <- names(dat2) #the NAs are so the name vector matches the correct length
cndat <- rbind(dat2, derepDupes) #add dupes back in

cndat<- cndat[which(cndat$Sample2 !="other"),]

#merge data
str(alldat$Plant)
cndat$Sample2 <- as.numeric(as.character(cndat$Sample2))

dim(alldat)
alldat2 <- merge(cndat, alldat, by.x = "Sample2", by.y = "Plant", all=T)
dim(alldat2)

#Calculate NDFA
for(i in 1:length(alldat2[,1])){
  alldat2$ndfa[i] <- 100*((5.14 - alldat2$d15N....vs..air.[i])/(5.14 - min(na.omit(as.numeric(alldat2$d15N....vs..air.)))))
  #see the Func Ecol supplemental for this formula from Hogsberg. The 5.14 was the 'other' the nonfixer forb im the cndata.
}

dim(cndat)
dim(alldat2)
dim(alldat)

#End of season measurements
eos <- read.csv("./data/AsleBoxMeasurementsEndofSeason.csv")
alldat2 <- merge(eos, alldat2, by.x = "Plant.Number", by.y = "Sample2", all=T)

#Add swainsonine results
swain <- read.csv("./data/SwainsonineMg_JHarrison.csv")
alldat2 <- merge(alldat2, swain, by.x = "Plant.Number", by.y = "sample_num", all=T)
dim(alldat2)

#Add culturing information

dat <- read.csv("data/Culturing_xp_results.csv")

alldat2 <- merge(dat, alldat2, 
                by.x = "plant", 
                by.y = "Plant.Number",
                all = T)
dim(alldat2)

#################
#Remove those plants where treatment appeared to be unsuccessful. 
#Use this to remove swainsonine positive plants with reduced a. fulva from consideration
#these are likely those plants for which treatment was unsuccessful

#Recode the #Value! data as zero as these had very small decimal values associated with them. 
alldat2$X..Swainsonine <- as.character(alldat2$X..Swainsonine)
alldat2$X..Swainsonine[alldat2$X..Swainsonine == "#VALUE!"] <- 0
alldat2$X..Swainsonine <- as.numeric(as.character(alldat2$X..Swainsonine))

#thees had swainsonine when they were not supposed to
alldat2$treatment_failed <- "no"
alldat2$treatment_failed[which(alldat2$X..Swainsonine > 0.001 & alldat2$treament %in% c("untreated_neg","treated_neg"))] <- "yes"

#Check out any plants that A. fulva was cultured from but that had their seed coats removed and didn't have swainsonine
#All but one of these were identified as failures via the preceeding code
alldat2$treatment_failed[which(alldat2$Undifilum > 1 & alldat2$treament %in% c("untreated_neg","treated_neg"))] <- "yes"

#remove the NA rows caused by technical replicates, inoculum, and cultures
alldat2 <- alldat2[grep("^\\d+$",alldat2$plant),]

write.csv(alldat2, file="./data/trait_and_treatment_data.csv", row.names = F)

rm(list=ls())
dat <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum", stringsAsFactors = T, header = T)
traits <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = T)
names(dat) <- gsub("ITS_\\w+_\\w+_Harrison_asle_(\\d+).*", "\\1", names(dat))

negs <- traits$plant[grep("neg", traits$treament)]

negs_afulva_counts <- dat[3331,names(dat) %in% negs]
names(negs_afulva_counts)[which(negs_afulva_counts > 0)]
#178, 52, 108, 61, 163, 43, 58, 156, 86

traits$plant[traits$treatment_failed == "yes"]

#lets look at the positive plants and see how many had afulva. 
pos <- traits$plant[grep("plus", traits$treament)]
pos_afulva_counts <- dat[3331,names(dat) %in% pos]

#combine pcr duplicates
newdat <- vector()
k <- 1
nms <- gsub("(\\d+)\\.\\d+", "\\1",names(pos_afulva_counts))

for(i in unique(nms)){
  if(length(grep(i, names(pos_afulva_counts))) > 1){
    newdat[k] <- sum(pos_afulva_counts[grep(i, names(pos_afulva_counts))])
    names(newdat)[k] <- i
    k <- k + 1
  }else{next}
}
#lets combine the count data, culture data, and swainsonine, to see which
#positives are duds

toconsider <- traits[grep("plus",traits$treament),]

culture_plus <- as.vector(na.omit(toconsider$plant[toconsider$Undifilum > 0]))
swain_plus <- as.vector(na.omit(toconsider$plant[toconsider$X..Swainsonine > 0]))
seq_plus <- names(newdat)[newdat > 0]

output <- toconsider$plant %in% as.numeric(c(culture_plus, swain_plus, seq_plus))
table(output)

#if we cut all the suspect ones, what would we have left?
bad <- toconsider$plant[!output]
testdf <- traits[traits$treatment_failed == "no",]
testdf <- testdf[!(testdf$plant %in% bad), ]
table(testdf$treament)
#still quite a lot. 

#lets write that to our traits csv
traits$conservative_treatmentFailed <- "yes"

traits$conservative_treatmentFailed[traits$plant %in% testdf$plant] <- "no"
table(traits$conservative_treatmentFailed)

write.csv(traits, file = "./data/trait_and_treatment_data.csv", row.names = F)

#look for contaminants

rm(list=ls())
dat <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum", stringsAsFactors = F, header = T)

blanks <- dat[,grep("lank", names(dat))]
contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)

forcounts <- dat[-contams, 2:length(dat)]

sum(rowSums(forcounts[c(1:(length(forcounts[,1])-3),length(forcounts[,1])), ])) / (length(forcounts[1,])/2)

length(dat$OTUID)

sum(rowSums(forcounts))
sum(rowSums(blanks)[contams])

write.csv(as.character(dat$OTUID[contams]), file = "./data/fungal_contams.csv", row.names = F)

rm(list=ls())
dat <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum", stringsAsFactors = T, header = T)

blanks <- dat[,grep("lank", names(dat))]
contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)

forcounts <- dat[-contams, 2:length(dat)]
sum(rowSums(forcounts))

sum(rowSums(forcounts[1:(length(forcounts[,1])-4), ])) / (length(forcounts[1,])/2)
length(dat$OTUID)
sum(rowSums(blanks)[contams])

dat$OTUID[contams]
sum(rowSums(blanks)[contams])

#last element was ISD, so cut that.
write.csv(as.character(dat$OTUID[contams])[1:2], file = "./data/bacterial_contams.csv", row.names = F)

