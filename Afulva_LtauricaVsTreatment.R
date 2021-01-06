rm(list=ls())
dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv")
metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
metadat <- metadat[metadat$treatment_failed == "no",]
dim(metadat)

afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x  %in% afulva$V1

dat <- merge(dat, metadat, by.x = "sample", by.y="plant")
dim(dat)

afulvaHits <- which(names(dat) %in% afulva$V1)
afulva <- rowSums(dat[, afulvaHits]) / dat$ISD
aggregate(log10(afulva)~dat$treament, FUN = mean)
TukeyHSD(aov(log10(afulva)~dat$treament)) #this is too conservative. Need to use Bayesian estimation of means. 
#Or compare the pis between groups. 

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

#make isd transformed data
isd_tr <- newdat[,3:which(names(newdat)=="Otu997")] /  newdat$ISD

pdf(width = 8, height = 10, file = "./visuals/rarefungi_vs_treatment.pdf")
par(mfrow = c(2,1), oma = c(4,3,1,1))
boxplot(rowSums(newdat[,3:which(names(newdat)=="Otu997")]) ~ newdat$treament,
        outline = F, las = 2,
        ylab = "Summed fungal relative abundances",
        xlab = "",
        xpd = NA, 
        names = c("",
                  "",
                  "",
                  "",
                  "",
                  "",
                  "",
                  ""),
        col = c("white", "white", "lightgray", "lightgray","lightgray", "lightgray", "white", "white"),
        cex.lab = 1.5)

boxplot(rowSums(isd_tr) ~ newdat$treament,
        outline = F, las = 2,
        ylab = "Summed fungal rel. abund. / ISD",
        xlab = "",
        xpd = NA, 
        names = c("Control -",
                  "Control +",
                  "Treated, Control -",
                  "Treated, Control +",
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +"),
        col = c("white", "white", "lightgray", "lightgray","lightgray", "lightgray", "white", "white"),
        cex.lab = 1.5)
dev.off()
######
# plot
######
aggregate(afulva~dat$treament, FUN = mean)
TukeyHSD(aov(afulva~dat$treament))

pdf(width = 8, height = 6, file = "./visuals/A.fulva_boxplot_w_controls.pdf")
par(oma = c(4,4.5,0,0), mar = c(4,4.5,2,2))
boxplot(afulva~dat$treament, outline = F, las = 2,
        ylab = expression(paste("Ratio of ", italic('A. fulva')," to ISD", sep = "")),
        xlab = "",
        xpd = NA, 
        names = c("Control -",
                  "Control +",
                  "Treated, Control -",
                  "Treated, Control +",
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +"),
        col = c("white", "white", "lightgray", "lightgray","lightgray", "lightgray", "white", "white"),
        cex.lab = 1.5)
dev.off()

pdf(width = 8, height = 6, file = "./visuals/A.fulva_boxplot.pdf")
par(oma = c(4,2,0,0))

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

stripchart(afulva[-controls] ~ dat$treament[-controls],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(c("darkgray", "darkgray", "lightgray", "lightgray"), alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           ylab=expression(paste("Ratio of ", italic('A. fulva')," to ISD", sep = "")),
           xlim=c(0.5,4.5)
           ,ylim=c(0,0.35)
           ,frame.plot=F
)

boxplot(afulva[-controls]~dat$treament[-controls], 
        outline = F, las = 2,
        add = T,
        ylab = "",
        xlab = "",
        names = c(
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +"),
        col = add.alpha(c("lightgray", "lightgray", "white", "white"), alpha=0.5),
        cex.lab = 1.5, xpd = NA)
dev.off()

##############
# L. taurica #
##############

pdf(width = 8, height = 6, file = "./visuals/L.taurica_boxplot_w_controls.pdf")
par(oma = c(4,4.5,0,0), mar = c(4,4.5,2,2))
boxplot(levtaur~dat$treament, outline = F, las = 2,
        ylab = expression(paste("Ratio of ", italic('L.taurica')," to ISD", sep = "")),
        xlab = "",
        names = c("Control -",
                  "Control +",
                  "Treated, Control -",
                  "Treated, Control +",
                  "Treated -",
                  "Treated +",
                  "Untreated -",
                  "Untreated +"),
        col = c("white", "white", "lightgray", "lightgray","lightgray", "lightgray", "white", "white"),
        cex.lab = 1.5, xpd = NA)
dev.off()

controls <- grep("ontrol",dat$treament)

pdf(width = 8, height = 6, file = "./visuals/L.taurica_boxplot.pdf")
par(oma = c(4,2.5,0,0))

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

stripchart(levtaur[-controls] ~ dat$treament[-controls],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(c("darkgray", "darkgray", "lightgray", "lightgray"), alpha=0.5),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           xpd = NA,
           ylab=expression(paste("Ratio of ", italic('L.taurica')," to ISD", sep = "")),
           xlim=c(0.5,4.5)
           ,ylim=c(0,100)
           ,frame.plot=F
)

boxplot(
  levtaur[-controls]~dat$treament[-controls], 
  outline = F, 
  las = 2,
  add = T,
  ylab = "",
  xlab = "",
  names = c(
    "Treated -",
    "Treated +",
    "Untreated -",
    "Untreated +"),
  col = add.alpha(c("lightgray", "lightgray", "white", "white"), alpha=0.5),
  cex.lab = 1.5, xpd = NA)
dev.off()

##############
# Plot of both
##############

pdf(width = 8, height = 10, file = "./visuals/afulva_levtaur_boxplot.pdf")
par(oma = c(10,4,1,1), mar = c(1,4,1,1), mfrow = c(2,2))

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

stripchart(log10(afulva[-controls]) ~ dat$treament[-controls],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.6),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           xpd = NA,
           ylab=expression(paste("Ratio of ", italic('A. fulva')," to ISD", sep = "")),
           xlim=c(0.5,4.5),
           cex.lab = 1.5
           ,ylim=c(-3,1) # c(0,0.1)  for rel. abund
           ,frame.plot=F
)

boxplot(log10(afulva[-controls])~dat$treament[-controls], 
        outline = F, las = 2,
        add = T,
        ylab = "",
        xlab = "",
        frame.plot=F,
        names = c(
          "",
          "",
          "",
          ""),
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
        cex.lab = 1.5, xpd = NA)

##Plot culture data
newdat <- read.csv("./data/trait_and_treatment_data.csv")
newdat <- newdat[newdat$treatment_failed == "no",]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

newdat <- newdat[-grep("ontrol", newdat$treament),]
newdat$treament <- droplevels( newdat$treament)

#Plot A. fulva colonization

barplot(height = aggregate(newdat$Undifilum/newdat$num_leaf_segments ~ 
                             newdat$treament, FUN = mean)[,2],
        ylim = c(0,0.25),
        las = 2,
        ylab = "",
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.6),
        border = add.alpha("white", alpha = 0.001),
        names = "")

text(y = 0.12, x = -1.1, srt = 90, xpd = NA, expression(paste( 
                                         italic('A. fulva')," infection rate",sep="")),cex = 1.5)

#Then I calculate the standard deviation for each specie and condition :
stdev <- aggregate(newdat$Undifilum/newdat$num_leaf_segments ~ 
                     newdat$treament, FUN= sd)
lens <- aggregate(newdat$Undifilum/newdat$num_leaf_segments ~ 
            newdat$treament, FUN= length)
means <- aggregate(newdat$Undifilum/newdat$num_leaf_segments ~ 
                    newdat$treament, FUN= mean)
ci <- stdev[,2]* 1.96 / lens[,2]

colz <- c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4")
loc <- c(.7,1.9,3.1,4.3)
for(i in 1:4){
  segments(x0 = loc[i], x1 = loc[i],
      y0 = means[i,2] - ci[i],
       y1 = means[i,2] + ci[i],
       col = add.alpha(colz[i], alpha = 0.6),lwd = 2)
  segments(x0 = loc[i]-.1, x1 = loc[i]+.1,
           y0 = means[i,2] - ci[i],
           y1 = means[i,2] - ci[i],
           col = add.alpha(colz[i], alpha = 0.6),lwd = 2)
  segments(x0 = loc[i]+.1, x1 = loc[i]-.1,
           y0 = means[i,2] + ci[i],
           y1 = means[i,2] + ci[i],
           col = add.alpha(colz[i], alpha = 0.6),lwd = 2)
}

axis(side = 1, labels = c("","","",""), at = loc, mgp = c(3, -1, 0.2))
# 
# stripchart(newdat$Undifilum/newdat$num_leaf_segments ~ newdat$treament,
#            vertical = TRUE,
#            data = dat,
#            method = "jitter",
#            pch = 20,
#            col = add.alpha(c("darkgray", "darkgray", "lightgray", "lightgray"), alpha=0.5),
#            cex=2,
#            xaxt="n",
#            cex.lab = 1.5,
#            yaxt="n",
#            bty = "n",
#            xpd = NA,
#            ylab= substitute(paste("% ", italic("A. fulva"), " infected", sep = "")),
#            xlim=c(0.5,4.5)
#            ,ylim=c(0,1)
#            ,frame.plot=F
# )
# 
# boxplot(newdat$Undifilum/newdat$num_leaf_segments ~ newdat$treament, 
#         las = 2, 
#         col = add.alpha(c("darkgray", "darkgray", "lightgray", "lightgray"), alpha=0.5),
#         add = T,
#         frame.plot=F,
#         outline = F, 
#         ylim=c(0,1),
#         names = c(
#                   "",
#                   "",
#                   "",
#                   ""
#         ),
#         cex = 1.5)

##Levtaur
stripchart(levtaur[-controls] ~ dat$treament[-controls],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.6),
           cex = 1.5,
           xaxt="n",
           yaxt="n",
           xpd = NA,
           cex.lab = 1.5,
           ylab=expression(paste("Ratio of ", italic('L.taurica')," to ISD", sep = "")),
           xlim=c(0.5,4.5)
           ,ylim=c(0,160)
           ,frame.plot=F
)

boxplot(
  levtaur[-controls]~dat$treament[-controls], 
  outline = F, 
  las = 2,
  add = T,
  ylab = "",
  xlab = "",
  frame.plot=F,
  names = c("", "", "", ""),
  col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
  cex.lab = 1.5, xpd = NA)

#############
#Mildew plot#
#############

dat <- read.csv("data/trait_and_treatment_data.csv")
dat <- dat[dat$treatment_failed == "no",]

dat$lvs_with_fungus[which(is.na(dat$lvs_with_fungus))] <- 0
#Remove dead plants
dat <- dat[-which(is.na(dat$X.LeavesCollection)),]

dat <- dat[-which(is.na(dat$lvs_with_fungus / dat$X.LeavesCollection )),]

dat <- dat[-grep("ontrol", dat$treament),]
dat$treament <- droplevels( dat$treament)

stripchart(dat$lvs_with_fungus / dat$X.LeavesCollection ~ dat$treament,
           vertical = TRUE,
           data = dat,
           method = "jitter",
           pch = 20,
           col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.6),
           cex=1.5,
           xaxt="n",
           cex.lab = 1.5,
           yaxt="n",
           bty = "n",
           xpd = NA,
           ylab= "Mildew infection rate",
           xlim=c(0.5,4.5)
           ,ylim=c(0,1)
           ,frame.plot=F
)

boxplot(dat$lvs_with_fungus / dat$X.LeavesCollection ~ dat$treament,
        las = 2,
        add = T,
        #xaxt = "n",
        frame.plot=F,
        col=add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.5),
        outline = F,
        names = c(
                  "",
                  "",
                  "",
                  "")
)


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
     y = c(-0.1,-0.1, -0.1, -0.1),
     x = c(1.2,2.2, 3.2, 4.2))


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
     y = c(-0.1,-0.1, -0.1, -0.1),
     x = c(-4.3,-3.3,-2.3,-1.3))

dev.off()

#compare Pis

rm(list=ls())
load("CNVRG_JOINED_ITS_conservative.Rdata")
library(CNVRG)
afulva <- read.table("./editedSeqdata/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F)

afulvaHits <- which(names(t_newdat) %in% afulva$V1)
isdOut <- isd_transform(model_output = modelOut,
                        isd_index=2294,
                        countData=t_newdat)


afu <- isdOut$pi[,,afulvaHits-1]
afu_ests <- apply(afu, c(1,2), sum )
colMeans(afu_ests)  #point estimates
apply(afu_ests, 2, FUN=function(x){mean(x) + 1.96*(sd(x)/sqrt(length(x)))})

table((afu_ests[,6] - afu_ests[,5]) > 0)[1] / sum(table((afu_ests[,6] - afu_ests[,5]) > 0))
#This seems wrong, but what could be the problem?

#Lets check with ml

isdOut <- isd_transform(model_output = modelOut,
                        isd_index=which(names(t_newdat) == "ISD"),
                        countData=t_newdat,
                        format = "ml")
afulvaHits <- which(names(isdOut) %in% afulva$V1)
rowSums(isdOut[,afulvaHits])

ests <- extract_point_estimate(modelOut = modelOut, countData = t_newdat, treatments = 8)
afulvaHits <- which(names(ests$pointEstimates_p) %in% afulva$V1)

metadat <- metadat[metadat$treatment_failed == "no",]
#Without dividing by the ISD, the relative abundance of A. fulva goes up, as expected.
aggregate(rowSums(ests$pointEstimates_p[,afulvaHits]) ~ metadat$treament, FUN = mean)

aggregate(rowSums(ests$pointEstimates_p[,afulvaHits]) ~ metadat$treament, FUN=
  function(x){mean(x) + 1.96*(sd(x)/sqrt(length(x)))})

aggregate(rowSums(ests$pointEstimates_p[,afulvaHits]) ~ metadat$treament, FUN=
            function(x){mean(x) - 1.96*(sd(x)/sqrt(length(x)))})

#with ISD the pattern changes a tad, but still 8 is way more than 7, when using ISD
aggregate((rowSums(ests$pointEstimates_p[,afulvaHits]) /
            ests$pointEstimates_p$ISD) ~ metadat$treament, FUN = mean)

TukeyHSD(aov(rowSums(log(ests$pointEstimates_p[,afulvaHits]) /
  ests$pointEstimates_p$ISD) ~ metadat$treament))



aggregate(rowSums(ests$pointEstimates_p[,afulvaHits]) /
            ests$pointEstimates_p$ISD ~ metadat$treament,
          FUN = function(x){mean(x) - 1.96*(sd(x)/sqrt(length(x)))})


pis <- rstan::extract(modelOut, "pi")
afulvaHits <- which(names(t_newdat) %in% afulva$V1)
avs <- apply(pis$pi[,,afulvaHits - 1], c(1,2), sum)
colMeans(avs)
apply(avs, 2, quantile, probs=c(0.05, 0.95))

#Is the problem that my pi estimates are quite uncertain, due to the many zeros in the data?
#Moreover, how many of the plus plants actually had A. fulva?

#First, check that the estimates are indeed wide, for the pis.
isdOut <- isd_transform(model_output = modelOut,
                        isd_index=which(names(t_newdat) == "ISD"),
                        countData=t_newdat)
afulvaHits <- which(names(t_newdat) %in% afulva$V1)

#note the dimensions are NOT the same, need to subtract one from all indices
dim(t_newdat)
str(isdOut$pi)

afu <- isdOut$pi[,,c(afulvaHits-1)]
afu_ests <- apply(afu, c(1,2), sum )
apply(afu_ests, 2, mean)
apply(afu_ests, 2, quantile, probs=c(0.05, 0.95))

#trying without apply to make sure I am not jacking somehting up
outdf <- data.frame(matrix(nrow=1000, ncol = 8))
for(i in 1:8){
  for(j in 1:1000){
    outdf[j,i] <- sum(afu[j,i,1:57])
  }
}
colMeans(outdf)


#They aren't particularly wide, so scratch that idea.

#Determine how many plants had A. fulva

aggregate(rowSums(ests$pointEstimates_p[,afulvaHits]) /
            ests$pointEstimates_p$ISD ~ metadat$treament, FUN = function(x){length(which(x>0.05))})
#There are fewer plants in the plus categories, probably because of the filtering I did.

#This still doesn't explain anything though.Maybe the ISD doesn't work perfectly for A. fulva?
#I dont think that is a critical issue, and it shouldn't explain the discrepency between the pi
#and p estimates. Also the controls look as expected...

#If we assume everything is correct, then what could explain why the pis and ps dont match?

#is there anyway I could have the treatment order messed up?
#I reloaded the data and checked my metadat order and it seemed good. not shown

#count how many diffs there were for a. fulva and L. taur.

dfs <- diffs[,which(names(diffs) %in% afulva$V1)]
dfs[63,] #recall that the output is proportion greater than zero for the first treatment minus the second.
#This suggests that every taxon went up from 7 to 8.
dfs[45,] #same story

#This is further evidence that A.fulva goes up when expected...but what is up with the pi estimates
#Is summing them together causing problems, some how?

test <- ests$pointEstimates_p[,afulvaHits]
outag <- list()
for(i in 1:57){
 outag[[i]] <- aggregate(test[,i] ~ metadat$treament, FUN = mean)
}

result <- sapply(outag, FUN=function(x){x[8,2] > x[7,2]})
table(result)

result <- sapply(outag, FUN=function(x){x[6,2] > x[5,2]})
table(result)

#with pis
test <- ests$pointEstimates_pi[,afulvaHits]

testout <- NA
for(i in 1:57){
#testout[i] <- test[8,i] > test[7,i]
testout[i] <- test[6,i] > test[5,i]
}

#So individually, all pis go up...but when I summed them earlier they don't?
rowSums(test)  #this looks ok.
rowSums(test) / ests$pointEstimates_pi$ISD #but this does not!

test <- test / ests$pointEstimates_pi$ISD
testout <- NA
for(i in 1:57){
  #testout[i] <- test[8,i] > test[7,i]
  testout[i] <- test[6,i] > test[5,i]
}

#So division by ISD causes there to be no effect of treatment for non-control samples, but
#if we don't divide by the ISD all is well. Moreover, division by the ISD does not cause problems
#when we use the max likelihood estimates for p.

#I could imagine pi values being more constrained then p values, perhaps that is the issue.
#The p values seem potentially more conservative.
#I could also imagine the estimates for the ISD to be very uncertain.
#Checked this and doesnt seem to be the case: summary(pis$pi[,7,2293])

#paste p estimates versus treatments. Trying to determine what treatments the real
#high a. fulva plants are in.
cbind(c(rowSums(ests$pointEstimates_p[,afulvaHits]) /
            ests$pointEstimates_p$ISD), metadat$treament)

summary(aov((rowSums(ests$pointEstimates_p[,afulvaHits]) /
      ests$pointEstimates_p$ISD) ~ metadat$treament))

TukeyHSD(aov(rowSums(ests$pointEstimates_p[,afulvaHits]) /
              ests$pointEstimates_p$ISD ~ metadat$treament))

#Make a table that has all the info in it, doing by hand bc easier.


##############
# L. taurica #
##############

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F)
levtaurHits <- which(names(ests$pointEstimates_p) %in% levtaur$V1)

dfs <- diffs[,which(names(diffs) %in% levtaur$V1)]


TukeyHSD(aov(rowSums(ests$pointEstimates_p[,levtaurHits]) /
               ests$pointEstimates_p$ISD ~ metadat$treament))

rowSums(ests$pointEstimates_pi[,levtaurHits]) /
               ests$pointEstimates_pi$ISD

aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) ~ metadat$treament, FUN=
            function(x){mean(x) + 1.96*(sd(x)/sqrt(length(x)))})

aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) ~ metadat$treament, FUN=
            function(x){mean(x) - 1.96*(sd(x)/sqrt(length(x)))})

aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) ~ metadat$treament, FUN=mean)

#First, check that the estimates are indeed wide, for the pis.
isdOut <- isd_transform(model_output = modelOut,
                        isd_index=which(names(t_newdat) == "ISD"),
                        countData=t_newdat)
levtaurHits <- which(names(t_newdat) %in% levtaur$V1)

afu <- isdOut$pi[,,c(levtaurHits-1)]
afu_ests <- apply(afu, c(1,2), sum )
apply(afu_ests, 2, mean)
apply(afu_ests, 2, quantile, probs=c(0.05, 0.95))


aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) / ests$pointEstimates_p$ISD ~
            metadat$treament, FUN=
            function(x){mean(x) + 1.96*(sd(x)/sqrt(length(x)))})

aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) / ests$pointEstimates_p$ISD ~
            metadat$treament, FUN=
            function(x){mean(x) - 1.96*(sd(x)/sqrt(length(x)))})

aggregate(rowSums(ests$pointEstimates_p[,levtaurHits]) / ests$pointEstimates_p$ISD ~
            metadat$treament, FUN=mean)

################
#See if there are any samples without either A. fulva or L. taur
#################

rm(list=ls())
dat <- read.table("./data/otutableJOINED_ISD_ASLE_duds_sum", header = T, stringsAsFactors = F)
afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x  %in% afulva$V1

afulvaHits <- which(dat$OTUID %in% afulva$V1)

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F) 
#make sure none of these are in contams. Nope
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x %in% levtaur$V1

levtaurHits <- which(dat$OTUID %in% levtaur$V1)

notpresent <- 1 + which(rowSums(dat[c(afulvaHits, levtaurHits),2:length(dat)]) < 100)
k <- 1
out <- NA
for(i in notpresent){
out[k]<-vegan::diversity(dat[,i])
k <- k + 1
} 
dat2 <- dat[,-notpresent]
k <- 1
out2 <- NA
for(i in 2:length(dat2)){
  out2[k]<-vegan::diversity(dat2[,i])
  k <- k + 1
}  
mean(out2)
mean(out)

########################################################
#BAyesin version
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
              tau = sim.mod.sam$tau,
              mu = sim.mod.sam$mu
  ))
}

dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv")
metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
metadat <- metadat[metadat$treatment_failed == "no",]
dim(metadat)

afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
#contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
#contam$x  %in% afulva$V1

dat <- merge(dat, metadat, by.x = "sample", by.y="plant")
dim(dat)

afulvaHits <- which(names(dat) %in% afulva$V1)
dat <- dat[order(dat$treament),]
afulva <- rowSums(dat[, afulvaHits]) / dat$ISD
aggregate(log10(afulva)~dat$treament, FUN = mean)
TukeyHSD(aov(log10(afulva)~dat$treament))

library(rjags)

start <- indexer(dat$treament)[[1]]
end <- indexer(dat$treament)[[2]]
out <- modelRun(log10(afulva), start, end)


levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F)
levhits <- which(names(dat) %in% levtaur$V1)
lev <- rowSums(dat[, levhits]) / dat$ISD
out <- modelRun(lev, start, end)

