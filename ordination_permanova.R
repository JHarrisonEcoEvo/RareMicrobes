rm(list=ls())
library(vegan)
library(ecodist)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

dat2 <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv", stringsAsFactors = F)

#Bring in trait data
#Bring in trait data
dat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
dat <- dat[dat$treatment_failed == "no",]

dat <- merge(dat, dat2, by.x = "plant", by.y = "sample")

dat <- dat[match(dat2$sample,dat$plant),]
# 
afulva <- read.table("./data//afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F)
#do for lev taur. Not changing object names for ease
levtaur <- read.table("./data//levtaurMatchesITS", stringsAsFactors=F, header = F)

afulvaHits <- which(names(dat) %in% afulva$V1)
levtaurHits <- which(names(dat) %in% levtaur$V1)

#Divide by ISD
dat[,65:length(names(dat))]  <- dat[,65:length(names(dat))] / dat$ISD

#identify controls
controls <- grep("ontrol", dat$treament)

#change colors
dat$colors <- NA
dat$colors[dat$treament == "untreated_neg"] <- "deeppink"
dat$colors[dat$treament == "untreated_plus"] <- "deeppink4"
dat$colors[dat$treament == "treated_neg"] <- "deepskyblue"
dat$colors[dat$treament == "treated_plus"] <- "deepskyblue4"

pdf(width = 10, height = 6, file = "./visuals/ordination.pdf")
par(mfrow = c(1,3), oma=c(2,2,4,5))
maiUse <- par()$mai
ord <- cmdscale(d = distance(dat[-controls,65:(length(names(dat))-3)], method = "bray-curtis"),
                k = 2)
plot(ord, type = "n", xlab = "PCoA 1", ylab = "PCoA 2", main = substitute(paste("All taxa")))
points(ord, pch=16, col=add.alpha(dat$colors[-controls],0.5), cex=2)
ordispider(ord, dat$treament[-controls], col=unique(add.alpha(dat$colors[-controls],0.5)))

#note that manhattan, kulcsynksi and bray-curtis are all the same for proportion data. 

#omit the big guns

dat2 <- dat[,-c(afulvaHits, levtaurHits)]
ord <- cmdscale(d = distance(dat2[-controls,65:(length(names(dat2))-3)], method = "bray-curtis"),
                k = 2)
plot(ord, type = "n", xlab = "PCoA 1", ylab = "PCoA 2", 
     main = substitute(paste("All taxa except ", italic('A. fulva')," and ",italic('L. taurica'))))
points(ord, pch=16, col=add.alpha(dat$colors[-controls],0.5), cex=2)
ordispider(ord, dat$treament[-controls], col=unique(add.alpha(dat$colors[-controls],0.5)))

par(mai=c(0,0.5,0,0))
plot(NULL)
legend("center",legend =  c(expression(paste("No inoculum,", italic('A. fulva -'))),
                            expression(paste("No inoculum,", italic('A. fulva +'))),
                            expression(paste("Inoc. treated,", italic('A. fulva -'))),
                            expression(paste("Inoc. treated,", italic('A. fulva +')))), 
       col=unique(add.alpha(dat$colors[-controls],0.5)),
       pch = 16,
       xpd = NA,
       bty = "n",
       cex = 3)

dat$colors <- NA
dat$colors[dat$treament == "untreated_neg"] <- "chartreuse3"
dat$colors[dat$treament == "untreated_plus"] <- "brown3"
dat$colors[dat$treament == "treated_neg"] <- "chartreuse"
dat$colors[dat$treament == "treated_plus"] <- "brown1"
# The following ordinations seem of little use
# par(mai=maiUse)
# 
# ord <- cmdscale(d = distance(rowSums(dat[-controls,afulvaHits]), method = "bray-curtis"),
#                 k = 2)
# plot(ord, type = "n", xlab = "PCoA 1", ylab = "PCoA 2", main = substitute(paste(italic('A. fulva'),"*")))
# points(ord, pch=16, col=add.alpha(dat$colors[-controls],0.5), cex=2)
# ordispider(ord, dat$treament[-controls], col=unique(add.alpha(dat$colors[-controls],0.5)))
# 
# ord <- cmdscale(d = distance(rowSums(dat[-controls,levtaurHits]), method = "bray-curtis"),
#                 k = 2)
# plot(ord, type = "n", xlab = "PCoA 1", ylab = "PCoA 2", main = substitute(italic('L. taurica')))
# points(ord, pch=16, col=add.alpha(dat$colors[-controls],0.5), cex=2)
# ordispider(ord, dat$treament[-controls], col=unique(add.alpha(dat$colors[-controls],0.5)))
# 
# par(mai=c(0,0.5,0,0))
# 
# plot(NULL)
# legend("center",legend =  c(expression(paste("No inoculum,", italic('A. fulva -'))),
#                             expression(paste("No inoculum,", italic('A. fulva +'))),
#                             expression(paste("Inoc. treated,", italic('A. fulva -'))),
#                             expression(paste("Inoc. treated,", italic('A. fulva +')))), 
#        col=unique(add.alpha(dat$colors[-controls],0.5)),
#        pch = 16,
#        xpd = NA,
#        bty = "n",
#        cex = 3)
dev.off()

set.seed(666)
# I dont think we want to permute within the treatment group because we are trying to see if the treatment group itself matters. 
# We would permute across some other aspect of experimental design that we were worryed about, like some other blocking factor
#or date sampled, or somethign.
adonis(dat[-controls,65:(length(names(dat2)))] ~ dat$treament[-controls])#, strata = dat$treament[-controls])

adonis(dat[-controls,levtaurHits] ~ dat$treament[-controls])#, strata = dat$treament[-controls])
adonis(dat[-controls,afulvaHits] ~ dat$treament[-controls])#, strata = dat$treament[-controls])
adonis(dat2[-controls,65:(length(names(dat2))-3)] ~ dat2$treament[-controls])#, strata = dat2$treament[-controls])

#BACTERIA

rm(list=ls())
library(vegan)
library(ecodist)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

dat2 <- read.csv("./data/ml_p_table_16s_conservative.csv", stringsAsFactors = F)

#Bring in trait data
#Bring in trait data
dat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
dat <- dat[dat$treatment_failed == "no",]

dat <- merge(dat, dat2, by.x = "plant", by.y = "sample")

dat <- dat[match(dat2$sample,dat$plant),]
controls <- grep("ontrol", dat$treament)

#Divide by ISD
dat[,65:length(names(dat))]  <- dat[,65:length(names(dat))] / dat$ISD

#change colors
dat$colors <- NA
dat$colors[dat$treament == "untreated_neg"] <- "deeppink"
dat$colors[dat$treament == "untreated_plus"] <- "deeppink4"
dat$colors[dat$treament == "treated_neg"] <- "deepskyblue"
dat$colors[dat$treament == "treated_plus"] <- "deepskyblue4"

bacts <- dat[-controls,65:705]
forplot <- na.omit(data.frame(dat$colors[-controls], dat[-controls,65:705]))
# big <- which(colSums(na.omit(bacts)) > 11) #If desired, can omit abundant taxa
# bacts <- bacts[,-which(names(bacts) %in% names(big))]
ord2 <- cmdscale(d = distance(na.omit(bacts), method = "bray-curtis"),
                 k = 2)
#sort(colSums(bacts[,-big]))

pdf(width = 10, height = 6, file = "./visuals/ordinationBacteria.pdf")
par(mfrow =c(1,2))

plot(ord2, type = "n", xlab = "PCoA 1", ylab = "PCoA 2", main = substitute(paste("All taxa")))
points(ord2[,1], ord2[,2], pch=16, col=add.alpha(forplot$dat.colors..controls.,0.5), cex=2)

ordispider(ord2, na.omit(dat$treament[-controls]), col=unique(add.alpha(na.omit(dat$colors[-controls]),0.5)))

par(mai=c(0,0,0,0))

plot(NULL)
legend("center",legend =  c(expression(paste("No inoculum,", italic('A. fulva -'))),
                            expression(paste("No inoculum,", italic('A. fulva +'))),
                            expression(paste("Inoc. treated,", italic('A. fulva -'))),
                            expression(paste("Inoc. treated,", italic('A. fulva +')))), 
       col=unique(add.alpha(dat$colors[-controls],0.5)),
       pch = 16,
       xpd = NA,
       bty = "n",
       cex = 2)
dev.off()

adonis(dat[-controls,65:(length(names(dat2))-2)] ~ dat$treament[-controls])#, strata = dat$treament[-controls])

       