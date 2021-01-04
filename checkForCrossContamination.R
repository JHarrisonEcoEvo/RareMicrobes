rm(list=ls())
dat <- read.table("./data/coligoTableASLE", header = T, stringsAsFactors = F)

#bring in the key that shows which coligo was in which well
demux <- read.csv("./data/NovaSeq2_DemuxJH.csv", header = T, stringsAsFactors = F)

#clip off the ISD
dat <- dat[-97,]

# head(demux)
# head(dat)[1:5,1:5]

samps <- demux[grep("asle", demux$samplename), c(1,2,4)]
samps$combo <- apply( samps[ , ] , 1 , paste , collapse = "-" )

coligohits <- data.frame(matrix(ncol = 2))
k <- 1
consideredSamps <- vector()
for(i in samps$combo ) {
  #Should have two hits because of PCR duplicates.
  #These should be in different plates, but in the same well
  consideredSamps[k] <- i
  i <- gsub("\\w+-\\w+-(\\w+_\\w+_\\d+)", "\\1",i)
  x <- dat[, c(1, grep(paste(i, "$", sep = ""), names(dat)))]
  
  if(is.null(dim(x))){
    next
  }
  hits <- x[x[, 2] > 1 | x[, 3] > 1, ]
  
  if (length(hits) > 3) {
    print(i)
    print("the preceeding sample had too many hits")
    next
  }
  #coligo position
  position <- unique(demux$wellposition[demux$samplename == i])
    
  coligohits[k, 1] <-
    sum(hits[grep(position, hits$OTUID), 2:length(hits)])
  coligohits[k, 2] <-
    sum(colSums(hits[-grep(position, hits$OTUID), 2:length(hits)]))
 
  k <- k + 1
}

coligohits$perc_contam <- coligohits$X2/rowSums(coligohits)
mean(coligohits$perc_contam)
sum(coligohits$X2) / sum(rowSums(coligohits))

#write concerning samples to a text file
samps <- samps[-grep("146$", samps)]
write.csv(data.frame(consideredSamps[coligohits$perc_contam > 0.01], coligohits[coligohits$perc_contam > 0.01, ]),
          file = "suspectSamples.csv", row.names = F)

