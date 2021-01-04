rm(list=ls())
library(rjags)
dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv")
metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x  %in% afulva$V1

dat <- merge(dat, metadat, by.x = "sample", by.y="plant")
afulvaHits <- which(names(dat) %in% afulva$V1)
afulva <- rowSums(dat[, afulvaHits]) / dat$ISD
controls <- grep("ontrol",dat$treament)

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F) 
#make sure none of these are in contams. Nope
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x %in% levtaur$V1

levtaurHits <- which(names(dat) %in% levtaur$V1)
levtaur <- rowSums(dat[, levtaurHits]) / dat$ISD

#less abundant taxa analysis
newdat <- dat[,-c(levtaurHits, afulvaHits)]
#make sure there arent any taxa that are way more abundant than the others
plot(rev(sort(colSums(newdat[,3:which(names(newdat)=="Otu997")]) ))) #not too bad

bytaxon <- newdat[,3:which(names(newdat)=="Otu997")] 

output <- data.frame(matrix(ncol = dim(bytaxon)[2], nrow = 8))

for(i in 1:dim(bytaxon)[2]){
  output[,i] <- aggregate(bytaxon[,i] ~ newdat$treament, FUN = mean)[,2]
}
ord <- output[,rev(order(colSums(output)))]
ord[,1:10]
mean(colMeans(ord))
mean(unlist(ord))

#make isd transformed data
isd_tr <- newdat[,3:which(names(newdat)=="Otu997")] /  newdat$ISD

indexer <- function(x){
  start <- c(min(which(x =="treated_neg")),
             min(which(x =="treated_plus")),
             min(which(x=="untreated_neg")),
             min(which(x =="untreated_plus"))
  )
  end <- c(max(which(x =="treated_neg")),
           max(which(x =="treated_plus")),
           max(which(x =="untreated_neg")),
           max(which(x =="untreated_plus"))
  ) 
  return(list(start, end))
}

meanModeler <- "model{
#model likelihood by sampling group
for(j in 1:groups){
for(i in start[j]:end[j]){
response[i] ~ dnorm(mu[j], tau[j])
}
#priors
mu[j] ~ dnorm(0, 0.0001)
tau[j] <- pow(sigma[j], -2)
sigma[j] ~ dunif(0, 100)
}
#calculate differences between sampling groups
delta12 <- mu[1]-mu[2]
delta13 <- mu[1]-mu[3]
delta14 <- mu[1]-mu[4]
delta23 <- mu[2]-mu[3]
delta24 <- mu[2]-mu[4]
delta34 <- mu[3]-mu[4]
}"

#run model function
modelRun <- function(x, start, end){
  sim.mod.jags <- jags.model(textConnection(meanModeler),
                             data = list(
                               response = x,
                               groups = 4,
                               start = start,
                               end = end
                             ), 
                             n.chains=2, 
                             n.adapt=0)
  
  
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
  burn<-5000
  update(sim.mod.jags, 
         n.iter=burn)
  
  #sample
  #if one were interested in modeling variance here, then one could do that. 
  #I dont see any reason to compare between and within variance as one does in an Anova
  #to determine if means differ because we are directly modeling the means themselves. 
  
  sim.mod.sam <- jags.samples(model=sim.mod.jags, 
                              variable.names=c(
                                "mu",
                                "delta12",
                                "delta13",
                                "delta14",
                                "delta23",
                                "delta24",
                                "delta34",
                                "sigma",
                                "tau"
                              ), 
                              n.iter=4000,
                              thin=4)
  
  #Diagnostic statistics. 
  #I do not calculate these for delta, because delta depends on mu. If mus are
  #good, then deltas are good. 
  
  #Calculate the Gelman-Rubin statistic
  
  gr <- gelman.diag(as.mcmc.list(sim.mod.sam$mu), 
                    multivariate=F)
  print(gr)
  
  
  #Calculate temporal autocorrelation in the chains
  #this function determines if the mean between the first and last parts of the chain differ
  #if they don't differ then burn in was long enough. 
  #The values output are z scores (diff in means divided by SE)
  #The "frac" options  determine the fractions of the chains to compare
  
  gk <- geweke.diag(as.mcmc.list(sim.mod.sam$mu), 
                    frac1=0.1, frac2=0.5) 
  print(gk)
  
  print("Parameters that didn't have a long enough burn in.")
  print(names(which(2*pnorm(-abs(gk[[1]]$z))<0.05)))
  
  ###################
  #Get them Results!#
  ###################
  
  #Get estimates of mus and deltas
  
  calcCertainty <- function(x){
    plot(density(x))
    abline(v = 0)
    
    #calculate estimate of certainty for the threshold lovers
    p <- 1-length(which(x > 0))/length(x)
    return(p)
  }
  print("Comparisons")
  print(calcCertainty(sim.mod.sam$delta12))
  print(calcCertainty(sim.mod.sam$delta13))
  print(calcCertainty(sim.mod.sam$delta14))
  print(calcCertainty(sim.mod.sam$delta23))
  print(calcCertainty(sim.mod.sam$delta24))
  print(calcCertainty(sim.mod.sam$delta34))
  
  return(list(quantiles = list(rowMeans(apply(sim.mod.sam$mu[1,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               #Treatment group 2
                               rowMeans(apply(sim.mod.sam$mu[2,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               
                               #Treatment group 3
                               rowMeans(apply(sim.mod.sam$mu[3,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975))),
                               
                               #Treatment group 4
                               rowMeans(apply(sim.mod.sam$mu[4,,1:2], 2, quantile, 
                                              probs = c(0.025,0.5, 0.975)))),
              means = list(mean(apply(sim.mod.sam$mu[1,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[2,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[3,,1:2], 2, FUN=mean)),
                           mean(apply(sim.mod.sam$mu[4,,1:2], 2, FUN=mean))),
              tau = sim.mod.sam$tau
  ))
}

relabund <- rowSums(newdat[-grep("ontrol",newdat$treament),3:which(names(newdat)=="Otu997")])
treats <- newdat$treament[-grep("ontrol",newdat$treament)]
relabund <- relabund[order(treats)]
start <- indexer(treats[order(treats)])[[1]]
end <- indexer(treats[order(treats)])[[2]]

out <- modelRun(relabund, start, end)

#isd
isd <- rowSums(isd_tr)
treats <- newdat$treament[-grep("ontrol",newdat$treament)]
isd <- isd[order(treats)]
start <- indexer(treats[order(treats)])[[1]]
end <- indexer(treats[order(treats)])[[2]]

out <- modelRun(isd, start, end)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#plot
pdf(width = 8, height = 10, file = "./visuals/rarefungi_vs_treatment.pdf")
#par(mfrow = c(2,1), oma = c(4,4,1,1))

# stripchart(relabund ~ treats[order(treats)],
#            vertical = TRUE,
#            data = dat,
#            method = "jitter",
#            pch = 20,
#            col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
#            cex=1.5,
#            xaxt="n",
#            cex.lab = 1.5,
#            yaxt="n",
#            bty = "n",
#            xpd = NA,
#            ylab= "",
#            xlim=c(0.5,4.5)
#          #  ,ylim=c(0,1)
#            ,frame.plot=F
# )
# 
# boxplot(relabund ~ treats[order(treats)],
#         outline = F, las = 2,
#         ylab = "Summed fungal rel. abund.",
#         xlab = "",
#         xpd = NA,
#         add = T,
#         names = c(
#                   "Treated -",
#                   "Treated +",
#                   "Untreated -",
#                   "Untreated +"),
#         col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
#         cex.lab = 1.5)
# text(x = 1, y = 0.8, "a", xpd = NA)
# text(x = 2, y = 0.8, "b", xpd = NA)
# text(x = 3, y = 0.8, "a", xpd = NA)
# text(x = 4, y = 0.8, "ab", xpd = NA)
# text(x = -0.5, y = 0.9, "a)", cex = 2, xpd = NA)


stripchart(log10(isd) ~ treats[order(treats)],
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
           cex=2,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           bty = "n",
           xpd = NA,
           ylab= "",
           xlim=c(0.5,4.5)
           #,ylim=c(0,1)
           ,frame.plot=F
)

boxplot( log10(isd) ~ treats[order(treats)],
        outline = F, las = 2,
        ylab = "Summed fungal rel. abund. / ISD",
        xlab = "",
        add = T,
        cex=2,
        xpd = NA, 
        frame.plot = F,
        names = c(
          "Treated -",
          "Treated +",
          "Untreated -",
          "Untreated +"),
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
        cex.lab = 1.5)
text(x = 1, y = 17, "a", xpd = NA)
text(x = 2, y = 17, "a", xpd = NA)
text(x = 3, y = 17, "a", xpd = NA)
text(x = 4, y = 17, "a", xpd = NA)
text(x = -0.5, y = 19, "b)", cex = 2, xpd = NA)
dev.off()

#bacteria
dat <- read.csv("./data/ml_p_table_16s_conservative.csv")
metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)

dat <- dat[,-c(300,588)]

dat <- merge(dat, metadat, by.x = "sample", by.y="plant")

controls <- grep("ontrol",dat$treament)

dat <- dat[-controls,]
dat <- dat[order(dat$treament),]

aggregate(rowSums(dat[,3:642])~ dat$treament, FUN=mean)

dat[,3:642] <- dat[,3:642] / dat$ISD

modelRun(unlist(rowSums(dat[,3:642])), indexer(dat$treament)[[1]], indexer(dat$treament)[[2]])

pdf(width = 8, height = 10, file = "./visuals/bacteria_vs_treatment.pdf")

stripchart(rowSums(dat[,3:642]) ~ dat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
           cex=2,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           bty = "n",
           xpd = NA,
           ylab= "",
           xlim=c(0.5,4.5)
           #  ,ylim=c(0,1)
           ,frame.plot=F
)

boxplot(rowSums(dat[,3:642]) ~ dat$treament,
        outline = F, las = 2,
        ylab = "Summed fungal rel. abund.",
        xlab = "",
        xpd = NA,
        add = T,
        cex=2,
        bty = "n",
        frame.plot = F,
        names = c(
          "Treated -",
          "Treated +",
          "Untreated -",
          "Untreated +"),
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
        cex.lab = 1.5)

dev.off()
