library(rjags)

#input otu table, sample-well key, and taxonomy file
otus=read.delim("./data/ITSdata/otuTableZotus_rawITS.txt", header=T)
#otus <- read.csv("./data/fungiOTUtableZotus.csv", header=T)
sampNames=read.csv("./data/Sample_Well_Key_includesQuants (1).csv")
sampNames$pl_well=paste(sampNames$plate, sampNames$well, sep=".")
treatments <- read.csv("./data/planting_order (1).csv")

names(otus) <- gsub("X(\\d+\\.[a-h]\\d+)_.*","\\1",names(otus))

#check to see if anything is not in the key
#which(!(names(otus) %in% sampNames$pl_well))

#Wrangle a bit and get treatment info in there and everything ordered proper
# otus <- data.frame(t(otus), stringsAsFactors = FALSE)
# names(otus) <- as.character(as.vector(otus[1,]))
# otus <- otus[-1,]

#replace well based row names with the sample names
otus$samps <- row.names(otus)
newotus <- merge(otus, treatments, by.x="newotus.sampleName", by.y="rns")
#newotus <- merge(newotus, treatments, by.x="sampleName", by.y="rns")
dim(newotus)


newotus <- newotus[order(newotus$treament),]
treat <- newotus$treament

for(i in 3:100){
  newotus[,i] <- as.numeric(newotus[,i])
}
newotus <- newotus[,3:100]

#Define operating range for the model 

starts <-  c(1,9, 15, 21, 25, 64, 101, 138)
ends <- c(8, 14, 20, 24, 63, 100, 137, 173)

#####################
#Model specification#
#####################

community.model.level <- "model{
  for(i in 1:N){
    for(j in start[i]:end[i]){
      datamatrix[j,] ~ dmulti(p[j,], nreads[j])
      p[j,1:notus] ~ ddirch(pi[i,]*theta[i])   
    }

    pi[i,1:notus] ~ ddirch(alpha*hypertheta)
    theta[i] ~ dunif(1.0E-3, max(nreads)) 
    #changing the max value here makes little difference.
  }

  alpha ~ ddirch(hyperalpha)
  hypertheta ~ dunif(1.0E-3, max(nreads))

  for(k in 1:notus){
    hyperalpha[k] <-1/notus
  }
}"


#compile model
sim.mod.jags <- jags.model(textConnection(community.model.level),
                           data = list(
                             datamatrix = 1+newotus,
                             notus = dim(newotus)[2],               
                             nreads = rowSums(newotus),  
                             N = 8,
                             start = starts,
                             end = ends
                           ),
                            n.chains=2, 
                           n.adapt=0)


#adapt model
iter_needed <- 0
y=FALSE
while(y==FALSE){
  y <-  adapt(sim.mod.jags, 
              n.iter=1000,
              end.adaptation=FALSE)
  iter_needed <- 1000 + iter_needed
  if(iter_needed > 8000){break}
}

#burn in model.
burn<-5000
update(sim.mod.jags,
       n.iter=burn)

#here we extract the MCMC samples
sim.mod.sam <- jags.samples(model=sim.mod.jags,
                            variable.names=c(
                              "pi",
                              "alpha"
                            ),
                            n.iter=4000,
                            thin=4)

save.image("./BayesianRelativeAbundanceModelingFUNGI.RData")