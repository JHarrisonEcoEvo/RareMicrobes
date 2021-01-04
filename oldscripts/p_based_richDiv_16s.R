#J. Harrison
#Nov 13, 2018
set.seed(666)
library(rjags)

#This script plots max likelihood richness and diversity estimates for each
#replicate for each treatment group.
#The underlying data are the p terms from multinomial in modeling scripts
#not the pi terms. 

ml <- read.csv("./data/ML_bacteria_withShannonSimpson.csv")

dat <- read.csv("./data/bacteriaOTUtableZotus.csv")
dat$samps <- as.character(dat$samps)

samps <- read.csv("./data/planting_order.csv")
dim(dat)

#Merge datasets
dat <- merge(dat, 
      samps,
      by.x = "samps",
      by.y = "rns",
      all.x = T)
sum(colSums(dat[,2:112])) #still 33k

dat <- dat[order(dat$treament),]
dim(dat)

#clean up a little
dat <- dat[-(173:174),]

dat <- dat[,-1]
treatment <- dat$treament

#ML should be in the same order as dat bc that is how we modeled it. 
out <- data.frame(ml, dat)
#####################################

########
# Plot #
########

#Function from Mage that adds transparency to colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf(width = 7,
    height = 8,
    file = "./visuals/Supplemental_16s_RichDivVsTreatmentPlot.pdf")
par(oma = c(10,2,4,2),
    mfrow = c(2,1),
    mar = c(2,3,0,0))

#plot diversity

 stripchart(out$shannon ~ out$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("gray", alpha=0.8),
           cex = 2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Shannon-Weiner",
           xlim=c(0.5,8.5),
           ylim = c(3.2,4.2),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(out$shannon ~ out$treament,
        las = 2,
        add = T,
       # yaxt = "n",
       xaxt="n",
            col = c(NA,
                NA,
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA),
        outline = F,
        frame = F
        )
# axis(side = 2, 
#      at = c(0,1,2,3),
#      labels = c(0,1,2,3)
# )
text("a)",
     x = -0.5,
     y = 4.5,
     xpd = NA, 
     cex = 2.5)

stripchart(out$simpson ~ out$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("gray", alpha=0.8),
           cex=2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Simpson's",
           xlim = c(0.5,8.5),
           ylim = c(0.88,1),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(out$simpson ~ out$treament,
        las = 2,
        add = T,
       # yaxt = "n",
        col = c(NA,
                NA,
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA), 
        outline = F,
        frame = F,
        names = c("Control -",
                  "Control +",
                  "Treated,Control -",
                  "Treated,Control +",
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +")
        )
# axis(side = 2, 
#      at = c(0,25,50),
#      labels = c(0,25,50)
# )
text("b)",
     x = -0.5,
     y = 1.01,
     xpd = NA,
     cex = 2.5)
dev.off()

out2 <- out[which(out$treament %in% c("treated_neg", "treated_plus", "untreated_neg", "untreated_plus")),]
out2$treament <- droplevels(out2$treament)

pdf(width = 4,
    height = 8,
    file = "./visuals/16S_RichDivVsTreatmentPlot.pdf")
par(oma = c(10,2,4,2),
    mfrow = c(3,1),
    mar = c(2,3,0,0))

#Plot richness
stripchart(out2[,rich] ~ out2$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("gray", alpha=0.8),
           cex=2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Richness",
           xlim=c(0.5,4.5),
           ylim = c(0,100),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(out2[,rich] ~ out2$treament,
        las = 2,
        add = T,
       # yaxt = "n",
         cex.axis = 1.5,
        xaxt = "n",
            col = c(add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA),
        outline = F,
        frame = F
        )

text("d)",
     x = -0.4,
     y = 103,
     xpd = NA, 
     cex = 2.5)

mtext("Bacteria",
    side = 3,
    cex = 2,
    line = 2, 
    xpd = NA)
#plot diversity

stripchart(exp(out2[,rich + 1]) ~ out2$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("gray", alpha=0.8),
           cex=2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Shannon-Weiner",
           xlim=c(0.5,4.5),
           ylim = c(0,12),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(exp(out2[,rich + 1]) ~ out2$treament,
        las = 2,
        add = T,
       # yaxt = "n",
         cex.axis = 1.5,
       xaxt="n",
            col = c(add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA),
        outline = F,
        frame = F
        )
# axis(side = 2, 
#      at = c(0,1,2,3),
#      labels = c(0,1,2,3)
# )
text("e)",
     x = -0.4,
     y = 12.5,
     xpd = NA, 
     cex = 2.5)

stripchart(1/(out2[,rich + 2]) ~ out2$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("gray", alpha=0.8),
           cex=2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Simpson's",
           xlim=c(0.5,4.5),
           ylim = c(0,10),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(1/(out2[,rich + 2]) ~ out2$treament,
        las = 2,
        add = T,
       # yaxt = "n",
            col = c(add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA),
        outline = F,
        frame = F,
         cex.axis = 1.5,
        names = c(
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +")
        )
# axis(side = 2, 
#      at = c(0,25,50),
#      labels = c(0,25,50)
# )
text("f)",
     x = -0.4,
     y = 10.5,
     xpd = NA,
     cex = 2.5)
dev.off()

#Rarefaction analysis

##Compare ML estimates to proportions in data to see how different they are. 
# Rarefying leads to pretty different richness estimates, though overall patterns seem same within 
# a host. It does not show stark differences between hosts like the multinomial sampling does. 
# 
# rowSums(larvae[,7:length(larvae)])
# 
# r_larv <- vegan::rrarefy(larvae[which(rowSums(larvae[,7:length(larvae)]) > 999)
#   ,7:length(larvae)], 1000)
# rich_r <- NA
# for(i in 1:dim(r_larv)[1]){
#   rich_r[i] <- length(which(r_larv[i,] > 0))
# }
# 
# par(oma = c(8,2,2,2))
# boxplot(rich_r~treats[which(rowSums(larvae[,7:length(larvae)]) > 999)],
#         las = 2)


############################################################
# Analysis of differential richness among treatment groups #
############################################################

################
# Define model #
################

meanModeler <- "model{
  #model likelihood by sampling group
  for(j in 1:groups){
     for(i in start[j]:end[j]){
        response[i] ~ dnorm(mu[j], tau[j])
     }
  
  #priors, note that we allow unequal variance among groups
    mu[j] ~ dnorm(0, 0.0001)
    tau[j] <- pow(sigma[j], -2)
    sigma[j] ~ dunif(0, 100)
  }
}"

####################
# Define functions # 
####################

modelRun <- function(x, 
                     start, 
                     end){
    sim.mod.jags <- jags.model(textConnection(meanModeler),
                           data = list(
                             response = x,
                            groups = 8,
                             start = start,
                             end = end
                           ), 
                           n.chains = 2, 
                           n.adapt = 0)


#adapt
  iter_needed <- 0
  y=FALSE
  while(y==FALSE){
    y <-  adapt(sim.mod.jags, 
              n.iter=1000,
              end.adaptation=FALSE)
    iter_needed <- 1000 + iter_needed
    if(iter_needed > 5000){break}
  }

#burn
  update(sim.mod.jags, 
       n.iter = 10000)

#sample
#if one were interested in modeling variance here, then one could do that. 
#I dont see any reason to compare between and within variance as one does in an Anova
#to determine if means differ because we are directly modeling the means themselves. 

  sim.mod.sam <- jags.samples(model = sim.mod.jags, 
                            variable.names = c(
                              "mu"
                              ), 
                            n.iter = 4000,
                            thin = 4)
  return(sim.mod.sam)
}

diagnostics <- function(x) {
  #Diagnostic statistics.
  #I do not calculate these for delta, because delta depends on mu. If mus are
  #good, then deltas are good.
  
  #Calculate the Gelman-Rubin statistic
  gr <- gelman.diag(as.mcmc.list(x),
                    multivariate = F)
  print(gr)
  
  #Calculate temporal autocorrelation in the chains
  #this function determines if the mean between the first and last parts of the chain differ
  #if they don't differ then burn in was long enough.
  #The values output are z scores (diff in means divided by SE)
  #The "frac" options  determine the fractions of the chains to compare
  
  gk <- geweke.diag(as.mcmc.list(x),
                    frac1 = 0.1,
                    frac2 = 0.5)
  print(gk)
  
  print("Parameters that didn't have a long enough burn in.")
  print(names(which(2 * pnorm(-abs(
    gk[[1]]$z
  )) < 0.05)))
}

resultify <- function(x) {
  resultsTot <- list()
  for (i in 1:8) {
    resultsTot[[i]] <-
      rowMeans(apply(x[i, , ], 2, quantile, probs = c(0.125, 0.5, 0.975)))
  }
  
  outT <- matrix(nrow = 8,
                 ncol = 8)
  
  for (i in 1:8) {
    for (j in 1:8) {
      diffs <- x[i, , ] - x[j, , ]
      outT[i, j] <-  length(which(diffs < 0)) / length(diffs)
    }
  }
  
  return(list(resultsTot,
              outT)
         )
}

#Analyze

#treats
out$treament
start <- c(1,9,15,21,25,64,101,138)
end <- c(8,14,20,24,63,100,137,173)

modelout <- modelRun(out[,rich],
         start = start,
         end = end)
diagnostics(modelout$mu)
resultify(modelout$mu)
write.csv(resultify(modelout$mu)[[2]],
          file = "./richVsTreatment_16s_results.csv")

modelout <- modelRun(exp(out[,rich + 1]),
         start = start,
         end = end)
diagnostics(modelout$mu)
resultify(modelout$mu)
write.csv(resultify(modelout$mu)[[2]],
          file = "./shannonVsTreatment_16S_results.csv")

modelout <- modelRun(1/out[,rich + 2],
         start = start,
         end = end)
diagnostics(modelout$mu)
resultify(modelout$mu)
write.csv(resultify(modelout$mu)[[2]],
          file = "./simpsonVsTreatment_16S_results.csv")