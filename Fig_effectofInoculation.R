#############################
#Make non supplemental plot #
#############################
rm(list=ls())

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

newdat <- read.csv("./data/trait_and_treatment_data.csv")
newdat <- newdat[newdat$treatment_failed=="no",]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

newdat <- newdat[-grep("ontrol", newdat$treament),]
newdat$treament <- droplevels(newdat$treament)

pdf(width = 8, height = 10, file = "./visuals/Fig_inoculationEffect.pdf")
par(oma = c(9,4,2,0),
    mar = c(1,5,1,3),
    mfrow = c(2,2))


#plot morphospecies
stripchart(newdat$morphospecies_not_counting_Undifilum ~ newdat$treament,
           vertical = TRUE,
           data = newdat,
           method = "jitter",
           pch = 20,
           bty = "n",
           cex.lab = 1.5,
           col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.6)),
           cex=2,
           xaxt="n",
           yaxt="n",
           ylab="culture richness",
           xlim=c(0.5,4.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)

boxplot(newdat$morphospecies_not_counting_Undifilum ~ newdat$treament, 
        las = 2, 
        xaxt="n",
        frame.plot=F,
        col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.5)),
        add = T,
        outline = F)
mtext("a)", side = 2,
      line = 3.5, 
      padj = -6,
      cex = 2,
      las = 2,
      xpd = NA)

#Plot percentage colonization
stripchart(newdat$notAlfuInfected/newdat$num_leaf_segments ~ newdat$treament,
           vertical = TRUE,
           data = newdat,
           method = "jitter",
           pch = 20,
           col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.6)),
           cex=2,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           ylab="% leaves infected",
           xlim=c(0.5,4.5)
           #,ylim=c(0,0.6)
           ,frame.plot=F
)

boxplot(newdat$notAlfuInfected/newdat$num_leaf_segments ~ newdat$treament, 
        las = 2, 
        xaxt="n",
        col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.5)),
        frame.plot=F,
        add = T,
        outline = F)


mtext("b)", side = 2,
      line = 3.5, 
      padj = -6,
      
      cex = 2,
      las = 2)

#############
# Read data #
#pull in ml estimates, inocula taxa, and treatment info
dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv", stringsAsFactors = F, header = T)
#convert to absolute abundances, or not. Doing both for plot
dat_relabund <- data.frame(dat$sample, dat[,3:length(dat)])
dat_isd <- data.frame(dat$sample, dat[,3:length(dat)] / dat$ISD)

dim(dat)

#sanitycheck, using different object for the transform
#abs_dat[1,1:3]
#dat[1,1:5] / dat$ISD[1]

inocula <- read.table("./data/linking_ITS_inoculumOTUs_to_NOVAseqOTUS", stringsAsFactors = F, header = F)
traits <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F, header = T)
traits <- traits[!duplicated(traits$plant),]
traits <- traits[traits$treatment_failed=="no",]

#merge ml with treatment
merge_dat_relabund <- merge(dat_relabund, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)
merge_dat_isd <- merge(dat_isd, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)

#identify inoculum taxa
inoculum <- unique(c(inocula$V2, "Otu91"))

dat3 <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum",stringsAsFactors = F, header = T)

names(dat3) <- gsub("ITS_.*_(\\d+)", "\\1", names(dat3))
names(dat3) <- gsub("(\\d+)_.*", "\\1", names(dat3))
#combine PCR duplicates
combodat <- data.frame(matrix(nrow = dim(dat3)[1]))
k <- 1
for_names <- vector()
for(i in unique(names(dat3))){
  if(length(grep(paste("^",i,"$", sep = ""), names(dat3))) > 1){
    combodat[,k] <- rowSums(dat3[,grep(paste("^",i,"$", sep = ""), names(dat3))])
    for_names <- c(for_names, i)
    k <- k + 1
  }else{next}
}

names(combodat) <- for_names
inoculated <- combodat[dat3$OTUID %in% inoculum,]

otus <- dat3$OTUID[dat3$OTUID %in% inoculum]
inoculated <- data.frame(t(inoculated))
names(inoculated) <- otus
inoculated$plants <- row.names(inoculated)
merged_inoculated <- merge(inoculated, traits, by.x = "plants", by.y = "plant")

prevalence <- data.frame()

k <- 1
for(j in unique(merged_inoculated$treament)){
  df <- merged_inoculated[merged_inoculated$treament == j, 2:5]
  
  for(i in 1:length(df)){
    prevalence[i,k] <- length(which(df[,i] > 0)) / length(df[,i])
  }
  
  names(prevalence)[k] <- j
  k <- k + 1
}

#zeros are cases where the taxon was not present at all. 
#what about comparing when a taxon is more prevalent in one gropu than the other
table(prevalence[,2] > prevalence[,3]) #should have more trues than falses 
table(prevalence[,1] > prevalence[,5]) #


#make sure none of these are in contams. Nope
contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)

#note that 91 was also in the inoculum, but didnt show up in the miseq data. 
#This merge causes us to lose some of the stuff in inoculum object, because that object contains duds, ISD.

inoculum_df_isd <- merge_dat_isd[,c(1,which(names(merge_dat_isd) == "treament") , which(names(merge_dat_isd) %in% inoculum))]
inoculum_df_relab <- merge_dat_relabund[,c(1,which(names(merge_dat_relabund) == "treament") , which(names(merge_dat_relabund) %in% inoculum))]

inoculum_df_isd$sums <- rowSums(inoculum_df_isd[,3:length(inoculum_df_isd)])
inoculum_df_relab$sums <- rowSums(inoculum_df_relab[,3:length(inoculum_df_relab)])

inoculum_df_isd <- inoculum_df_isd[-grep("ontrol", inoculum_df_isd$treament),]
inoculum_df_relab <- inoculum_df_relab[-grep("ontrol", inoculum_df_relab$treament),]

summary(aov(inoculum_df_isd$sums ~ inoculum_df_isd$treament))
TukeyHSD(aov(inoculum_df_isd$sums ~ inoculum_df_isd$treament))
TukeyHSD(aov(log10(inoculum_df_isd$sums) ~ inoculum_df_isd$treament))

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


inoculum_df_isd <- inoculum_df_isd[order(inoculum_df_isd$treament),]
start <- indexer(inoculum_df_isd$treament)[[1]]
end <- indexer(inoculum_df_isd$treament)[[2]]

out <- modelRun(log10(inoculum_df_isd$sums), start, end)

stripchart(log10(inoculum_df_isd$sums) ~ inoculum_df_isd$treament,
           vertical = TRUE,
           method = "jitter",
           pch = 20,
           bty = "n",
           cex.lab = 1.5,
           col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.6)),
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="Ratio of summed fungal\nproportions to ISD",
           xlim=c(0.5,4.5),
           xpd = NA,
           #,ylim=c(0,0.6)
           frame.plot=F
)
boxplot(log10(inoculum_df_isd$sums) ~ inoculum_df_isd$treament,  outline = F, las = 2, ylab = "", xlab = "",
        names = c("","", "", ""),
        add = T,
        frame.plot = F,
        col = c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.5)))

mtext("c)", side = 2,
      line = 3.5, 
      padj = -6,
      
      cex = 2,
      las = 2)


##########################################
# determine efficacy of inoculum FOR 16S #
##########################################
#pull in ml estimates, inocula taxa, and treatment info
dat <- read.csv("./data/ml_p_table_16s_conservative.csv", stringsAsFactors = F, header = T)
#convert to absolute abundances
dat_ra <- data.frame(dat$sample, dat[,3:length(dat)])
dat_isd <- data.frame(dat$sample, dat[,3:length(dat)] / dat$ISD)
#sanitycheck, using different object for the transform
#abs_dat[1,1:3]
#dat[1,1:5] / dat$ISD[1]

inocula <- read.table("./data/linking_16S_inoculumOTUs_to_NOVAseqOTUS", stringsAsFactors = F, header = F)
length(unique(inocula$V1)) #206 taxa with 13823 matches?
traits <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F, header = T)
traits <- traits[!duplicated(traits$plant),]
traits <- traits[traits$treatment_failed == "no",]
#merge ml with treatment
merge_dat_ra <- merge(dat_ra, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)
merge_dat_isd <- merge(dat_isd, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)

dat3 <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum",stringsAsFactors = F, header = T)
dat3 <- dat3[1:657,c(1,grep("inoculum", names(dat3)))]

dat3 <- dat3[rowSums(dat3[,2:length(dat3)]) > 0,]

bacts <- dat3$OTUID
dat3 <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum",stringsAsFactors = F, header = T)

names(dat3) <- gsub("X16S_.*_(\\d+)", "\\1", names(dat3))
names(dat3) <- gsub("(\\d+)asle", "\\1", names(dat3))
#combine PCR duplicates
combodat <- data.frame(matrix(nrow = dim(dat3)[1]))
k <- 1
for_names <- vector()
for(i in unique(names(dat3))){
  if(length(grep(paste("^",i,"$", sep = ""), names(dat3))) > 1){
    combodat[,k] <- rowSums(dat3[,grep(paste("^",i,"$", sep = ""), names(dat3))])
    for_names <- c(for_names, i)
    k <- k + 1
  }else{next}
}

names(combodat) <- for_names

#do prevalence analysis
#dat3 is in same order as combodat
#identify inoculum taxa
inoculum <- unique(c(inocula$V2, bacts))
inoculum <- gsub("(centroid=Otu\\d+);.*", "\\1", inoculum)
inoculum <- gsub("=","\\.", inoculum)
inoculated <- combodat[dat3$OTUID %in% gsub("\\.","=", inoculum),]

otus <- dat3$OTUID[dat3$OTUID %in% gsub("\\.","=", inoculum)]
inoculated <- data.frame(t(inoculated))
names(inoculated) <- otus
inoculated$plants <- row.names(inoculated)
merged_inoculated <- merge(inoculated, traits, by.x = "plants", by.y = "plant")

prevalence <- data.frame()

k <- 1
for(j in unique(merged_inoculated$treament)){
  df <- merged_inoculated[merged_inoculated$treament == j, 2:56]
  
  for(i in 1:length(df)){
    prevalence[i,k] <- length(which(df[,i] > 0)) / length(df[,i])
  }
  
  names(prevalence)[k] <- j
  k <- k + 1
}

#identify inoculum taxa
inoculum <- unique(c(inocula$V2, bacts))
inoculum <- gsub("(centroid=Otu\\d+);.*", "\\1", inoculum)
inoculum <- gsub("=","\\.", inoculum)
names(merge_dat_isd) <- gsub("(centroid.Otu\\d+);.*", "\\1", names(merge_dat_isd))
names(merge_dat_ra) <- gsub("(centroid.Otu\\d+);.*", "\\1", names(merge_dat_ra))

inoculum_df_ra <- merge_dat_ra[,c(1,which(names(merge_dat_ra) == "treament") , which(names(merge_dat_ra) %in% inoculum))]
inoculum_df_isd <- merge_dat_isd[,c(1,which(names(merge_dat_isd) == "treament") , which(names(merge_dat_isd) %in% inoculum))]

inoculum_df_ra$sums <- rowSums(inoculum_df_ra[,3:length(inoculum_df_ra)])
inoculum_df_isd$sums <- rowSums(inoculum_df_isd[,3:length(inoculum_df_isd)])

#make a whole panel of boxplots
#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

inoculum_df_isd <- inoculum_df_isd[-grep("ontrol", inoculum_df_isd$treament),]

TukeyHSD(aov(log10(inoculum_df_isd$sums)~ inoculum_df_isd$treament))c

inoculum_df_isd <- inoculum_df_isd[order(inoculum_df_isd$treament),]
start <- indexer(inoculum_df_isd$treament)[[1]]
end <- indexer(inoculum_df_isd$treament)[[2]]

out <- modelRun(log10(inoculum_df_isd$sums), start, end)

stripchart(log10(inoculum_df_isd$sums)~ inoculum_df_isd$treament,
           vertical = TRUE,
           method = "jitter",
           pch = 20,
           bty = "n",
           cex.lab = 1.5,
           col= c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.6)),
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
          # ylim = c(-3.5,-1.5),
           xlim=c(0.5,4.5),
           xpd = NA,
           #,ylim=c(0,0.6)
           frame.plot=F
)
boxplot(log10(inoculum_df_isd$sums) ~ inoculum_df_isd$treament,  outline = F, las = 2, xlab = "",
        add = T,
        names = c("","","",""),
        ylab = "",
        frame.plot = F,
        col = c(add.alpha(c("cyan2", "cyan4", "darkgoldenrod2", "darkgoldenrod4"),alpha=0.5)),
        main = "")

text(c(expression(paste("Inoc. treated, ", italic("A. fulva")," -",
                           sep = "")),
          expression(paste("Inoc. treated, ", italic("A. fulva")," +",
                           sep = "")),
          expression(paste("No inoculum, ", italic("A. fulva")," -",
                           sep = "")),
          expression(paste("No inoculum, ", italic("A. fulva")," +",
                           sep = ""))),
     srt = 70, 
     xpd = NA,
     pos = 2,
     y = c(-3.7,-3.7,-3.7,-3.7),
     x = c(1.5,2.5,3.5,4.5))

text(c(expression(paste("Inoc. treated, ", italic("A. fulva")," -",
                        sep = "")),
       expression(paste("Inoc. treated, ", italic("A. fulva")," +",
                        sep = "")),
       expression(paste("No inoculum, ", italic("A. fulva")," -",
                        sep = "")),
       expression(paste("No inoculum, ", italic("A. fulva")," +",
                        sep = ""))),
     srt = 70, 
     xpd = NA,
     pos = 2,
     y = c(-3.7,-3.7,-3.7,-3.7),
     x = c(-7,-6,-5,-4))

text("Ratio of summed bacterial\nproportions to ISD (log10)",
                y = -2.5, x = -2.3, xpd = NA, srt = 90, cex = 1.5)

mtext("d)", side = 2,
      line = 3.5, 
      padj = -6,
      
      cex = 2,
      las = 2)

dev.off()

