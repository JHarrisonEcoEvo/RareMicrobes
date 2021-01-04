library(rjags)

#input otu table, sample-well key, and taxonomy file
otus=read.csv("./data/16s_data/bacteriaOTUtableZotus.csv", header=T)
sampNames=read.csv("./data/Sample_Well_Key_includesQuants (1).csv")
sampNames$pl_well=paste(sampNames$plate, sampNames$well, sep=".")
treatments <- read.csv("./data/planting_order.csv")

newotus <- merge(otus, sampNames, by.x="samps", by.y="pl_well")
newotus <- merge(newotus, treatments, by.x="sampleName.y", by.y="rns")
dim(newotus)
dim(otus)

newotus <- newotus[order(newotus$treament),]

for(i in 3:153){
  newotus[,i] <- as.numeric(as.character(newotus[,i]))
}
newotus <- newotus[,3:153]
dim(newotus)
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
                           n.adapt=0,
                           n.chains=2)


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
                              "p",
                              "alpha"
                            ),
                            n.iter=4000,
                            thin=4)

save.image("./BayesianRelativeAbundanceModelingBACTERIA_noplant_with_p.RData")