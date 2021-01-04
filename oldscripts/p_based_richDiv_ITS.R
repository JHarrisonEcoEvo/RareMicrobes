#J. Harrison
#Nov 13, 2018
set.seed(666)
library(rjags)

#This script plots max likelihood richness and diversity estimates for each
#replicate for each treatment group.
#The underlying data are the p terms from multinomial in modeling scripts
#not the pi terms. 

ml <- read.csv("./data/withOUTAlfuML_otuProportionTableFungi_last3fieldsRichDiv.csv")
#ml <- read.csv("./data/withAlfuML_otuProportionTableFungi_last3fieldsRichDiv.csv")

dat <- read.csv("./data/fungiOTUtableZotus.csv")
dat$newotus.sampleName <- as.character(dat$newotus.sampleName)

samps <- read.csv("./data/planting_order.csv")
dim(dat)

dat <- merge(dat, 
      samps,
      by.x = "newotus.sampleName",
      by.y = "rns",
      all.x = T)

dat <- dat[order(dat$treament),]
dim(dat)
dat <- dat[-(174:192),]

#ML should be in the same order as dat bc that is how we modeled it. 
out <- data.frame(ml, dat)
#####################################

########
# Plot #
########

#Figure out indices to plot. Just so the code works for A. fulva + or - datasets. 
rich <- which(names(out) == "newotus.sampleName") - 3

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
    file = "./visuals/Supplemental_RichDivVsTreatmentPlot.pdf")
par(oma = c(10,2,4,2),
    mfrow = c(3,1),
    mar = c(2,3,0,0))

#Plot richness
stripchart(out[,rich] ~ out$treament,
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
           xlim=c(0.5,8.5),
           ylim = c(0,24),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(out[,rich] ~ out$treament,
        las = 2,
        add = T,
        yaxt = "n",
        xaxt = "n",
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

axis(side = 2, 
     at = c(seq(0,24, by = 4)),
     labels = c(seq(0,24, by = 4)),
     las = 2)
mtext("Fungi",
    side = 3,
    cex = 2,
    line = 2,
    xpd = NA)
text("a)",
     x = -0.5,
     y = 25.5,
     xpd = NA, 
     cex = 2.5)
text("7",
     x = 3,
     y = 25.5,
     xpd = NA)
text("7",
     x = 6,
     y = 25.5,
     xpd = NA)
text("3,6",
     x = 7,
     y = 25.5,
     xpd = NA)
#plot diversity

stripchart(exp(out[,rich + 1]) ~ out$treament,
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
           xlim=c(0.5,8.5),
           ylim = c(0,6),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(exp(out[,rich + 1]) ~ out$treament,
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
text("b)",
     x = -0.5,
     y = 6.5,
     xpd = NA, 
     cex = 2.5)
text("2,6,8",
     x = 1,
     y = 5,
     xpd = NA)
text("1,3,5,7",
     x = 2,
     y = 5,
     xpd = NA)
text("2,6,8",
     x = 3,
     y = 5,
     xpd = NA)
text("2,6,8",
     x = 5,
     y = 5,
     xpd = NA)
text("1,3,5,7",
     x = 6,
     y = 5,
     xpd = NA)
text("2,6,8",
     x = 7,
     y = 5,
     xpd = NA)
text("1,3,5,7",
     x = 8,
     y = 5,
     xpd = NA)


stripchart(1/(out[,rich + 2]) ~ out$treament,
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
           xlim=c(0.5,8.5),
           ylim = c(0,35),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(1/(out[,rich + 2]) ~ out$treament,
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
text("c)",
     x = -0.5,
     y = 37,
     xpd = NA,
     cex = 2.5)
text("7",
     x = 6,
     y = 30,
     xpd = NA)
text("6",
     x = 7,
     y = 30,
     xpd = NA)

dev.off()

out2 <- out[which(out$treament %in% c("treated_neg", "treated_plus", "untreated_neg", "untreated_plus")),]
out2$treament <- droplevels(out2$treament)

pdf(width = 4,
    height = 8,
    file = "./visuals/ITS_RichDivVsTreatmentPlot.pdf")
par(oma = c(10,3,4,2),
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
           ylim = c(0,24),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(out2[,rich] ~ out2$treament,
        las = 2,
        add = T,
        yaxt = "n",
        xaxt = "n",
        cex.axis = 1.5,
        col = c(add.alpha("light gray", alpha=0.6),
                add.alpha("light gray", alpha=0.6),
                NA,
                NA),
        outline = F,
        frame = F
        )
axis(side = 2, 
     at = c(seq(0,24, by = 4)),
     labels = c(seq(0,24, by = 4)),
    cex = 1.5,
    cex.axis = 1.5,
    las = 2)

mtext("Fungi",
    side = 3,
    cex = 2,
    line = 2,
    xpd = NA)

text("a)",
     x = -0.4,
     y = 24.5,
     xpd = NA, 
     cex = 2.5)
text("ab",
     x = 1,
     y = 26,
     xpd = NA,
     cex = 1.5)
text("a",
     x = 2,
     y = 26,
     xpd = NA,
     cex = 1.5)
text("b",
     x = 3,
     y = 26,
     xpd = NA,
     cex = 1.5)
text("ab",
     x = 4,
     y = 26,
     xpd = NA,
     cex = 1.5)
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
           ylim = c(0,6),
           #,ylim=c(0,0.6)
           frame.plot=F,
           xpd = NA
)
boxplot(exp(out2[,rich + 1]) ~ out2$treament,
        las = 2,
        add = T,
       # yaxt = "n",
       xaxt="n",
       cex.axis = 1.5,
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
text("b)",
     x = -0.4,
     y = 6.5,
     xpd = NA, 
     cex = 2.5)
text("a",
     x = 1,
     y = 6.5,
     xpd = NA,
     cex = 1.5)
text("b",
     x = 2,
     y = 6.5,
     xpd = NA,
          cex = 1.5)
text("a",
     x = 3,
     y = 6.5,
     xpd = NA,
          cex = 1.5)
text("b",
     x = 4,
     y = 6.5,
     xpd = NA,
          cex = 1.5)
stripchart(1/(out2[,rich + 2]) ~ out2$treament,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("light gray", alpha=0.6),
           cex=2.5,
           cex.lab = 1.5,
           xaxt="n",
           yaxt="n",
           ylab="Simpson's",
           xlim=c(0.5,4.5),
           ylim = c(0,35),
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
        names = c( "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +")
        )
# axis(side = 2, 
#      at = c(0,25,50),
#      labels = c(0,25,50)
# )
text("c)",
     x = -0.4,
     y = 37,
     xpd = NA,
     cex = 2.5)
text("ab",
     x = 1,
     y = 37,
     xpd = NA,
          cex = 1.5)
text("a",
     x = 2,
     y = 37,
     xpd = NA,
          cex = 1.5)
text("b",
     x = 3,
     y = 37,
     xpd = NA,
          cex = 1.5)
text("ab",
     x = 4,
     y = 37,
     xpd = NA,
          cex = 1.5)
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
          file = "./richVsTreatment_ITS_results.csv")

modelout <- modelRun(exp(out[,rich + 1]),
         start = start,
         end = end)
diagnostics(modelout$mu)
resultify(modelout$mu)
write.csv(resultify(modelout$mu)[[2]],
          file = "./shannonVsTreatment_ITS_results.csv")

modelout <- modelRun(1/out[,rich + 2],
         start = start,
         end = end)
diagnostics(modelout$mu)
resultify(modelout$mu)
write.csv(resultify(modelout$mu)[[2]],
          file = "./simpsonVsTreatment_ITS_results.csv")