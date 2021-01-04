#Nov 3, 2020. 
#J. Harrison. 
#Univ. of Wyoming. 
rm(list=ls())
library(CNVRG)
dat <- read.table("./data/otutableJOINED_ISD_ASLE_duds_sum", 
                  fill = T, header = T, stringsAsFactors = F)

metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)

metadat <- metadat[order(metadat$treament),]
treatment <- metadat$treament

#the order of my OTU file needs to follow this:
#metadat$plant

#sanity check
#cbind(gsub("X16S_\\w+_\\w+_Harrison_asle_(\\d+)asle", "\\1", names(dat)), names(dat))
newnames <- gsub("ITS_\\w+_\\w+_Harrison_asle_(\\d+).*", "\\1", names(dat))

#combine PCR duplicates
newdat <- data.frame(matrix(nrow = dim(dat)[1]))
k <- 1
for_names <- vector()
for(i in unique(newnames)){
  if(length(grep(paste("^",i,"$", sep = ""), newnames)) > 1){ #this weeds out the non-plant stuff
    newdat[,k] <- rowSums(dat[,grep(paste("^",i,"$", sep = ""), newnames)])
    for_names <- c(for_names, i)
    k <- k + 1
  }else{next}
}

names(newdat) <- for_names

#now need to make the data follow the order of the metadata
#sanity check, compare these
#names(newdat)[match(metadat$plant, as.numeric(names(newdat)))]
#metadat$plant
#dim(newdat)

#These sequenced plants are not in the metadata
names(newdat)[!as.numeric(names(newdat)) %in% metadat$plant]
#we good

#next we need to remove those plants that are somewhat suspect. This step is optional.
metadat <- metadat[metadat$conservative_treatmentFailed == "no",]
newdat <- newdat[,names(newdat) %in% metadat$plant]

#need to get newdat into the same order as the metadata
metadat <- metadat[!duplicated(metadat$plant),]

newdat <- newdat[,na.omit(match(metadat$plant, as.numeric(names(newdat))))]

#trim metadata to remove the stuff that didnt get sequenced
metadat <- metadat[metadat$plant %in% names(newdat),]

#doublecheck that the metadat plant follows the names of the data
table(metadat$plant == names(newdat))
cbind(metadat$plant[metadat$plant != names(newdat)], names(newdat)[metadat$plant != names(newdat)])

#Next need to get the newdat df into the format needed for CNVRG and determine start and end indices

# transpose and put the names in place
t_newdat <- as.data.frame(t(newdat))
colnames(t_newdat) <- dat$OTUID
dim(t_newdat)
t_newdat <- data.frame(metadat$plant, t_newdat)
#metadat$treament

#remove taxa with zero counts
t_newdat <- t_newdat[,colSums(t_newdat[,2:length(t_newdat)]) > 0]

indexer <- function(x){
  starts <- vector()
  ends <- vector()
  k <- 1
  for(i in unique(x)){
    #Extract the indices for the treatment
    indices <- which(x == i)
    #Make a sequence from the min to the max of the indices.
    #We will use this to test that the indices are in order and that the data are formatted properly.
    test_indices <- seq(min(indices), max(indices), by = 1)
    if(any((indices == test_indices) == F)){
      print("ERROR: it does not appear that all the replicates for a treatment group are adjacent.")
    }else{
      starts[k] <- min(indices)
      ends[k] <- max(indices)
      k <- k + 1
    }
  }
  if(length(starts) != length(unique(x))){
    print("ERROR: there was a problem trying to calculate starting and ending indices for each treatment group. Check data formatting.")
    return("FAILED")
  }
  if(length(ends) != length(unique(x))){
    print("ERROR: there was a problem trying to calculate starting and ending indices for each treatment group. Check data formatting.")
    return("FAILED")
  }
  return(list(starts = starts,
              ends = ends))
}

#Note, modeling using the unedited count data failed when using VB.
modelOut <- varHMC(
  countData = t_newdat+1,
  starts = indexer(metadat$treament)$starts,
  ends = indexer(metadat$treament)$ends,
  algorithm = "NUTS",
  chains = 2,
  burn = 500,
  samples = 1500,
  thinning_rate = 2,
  cores = 16,
  params_to_save = c("pi", "p")
)

diffs <- diff_abund(model_output = modelOut, countData = counts)

save.image(file = "./CNVRG_JOINED_ITS_conservative.Rdata")


