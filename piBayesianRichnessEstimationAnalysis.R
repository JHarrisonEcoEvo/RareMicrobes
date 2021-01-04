# 
# load("./CNVRG_JOINED_ITS_conservative.Rdata") #This has to be changed to bacteria for the bacteria plot, obviously.
# #This script should be run interactively as certain parts of it require the input to be changed from fungi to bacteria. 
# 
# isdOut <- isd_transform(model_output = modelOut, isd_index=which(names(t_newdat) == "ISD"), countData=t_newdat)
# 
# #remove Afulva. for ITS. This is not relevant for bacteria.
# afulva <- read.table("./editedSeqdata/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
# 
# levtaur <- read.table("./editedSeqdata/levtaurMatchesITS", stringsAsFactors=F, header = F) 
# levtaurHits <- which(names(t_newdat) %in% c(levtaur$V1, as.character(dat[contams, 1])))
# 
# dat <- read.table("~/synced//otutableJOINED_ISD_ASLE_AFULVA_sum", stringsAsFactors = T, header = T)
# blanks <- dat[,grep("lank", names(dat))]
# contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)
# 
# afulvaHits <- which(names(t_newdat) %in% c(afulva$V1, as.character(dat[contams, 1])))
# afulvaHits <- unique(afulvaHits)
# 
# outSh <- diversity_calc(model_output = isdOut$pi[,,-afulvaHits], countData=t_newdat, params = "pi", 
#                       entropy_measure = "shannon", equivalents = T)
# simps <- diversity_calc(model_output = isdOut$pi[,,-afulvaHits], countData=t_newdat, params = "pi", 
#                         entropy_measure = "simpson", equivalents = T)
# 
# simps <- diversity_calc(model_output = isdOut$pi[,,-c(afulvaHits,levtaurHits)], countData=t_newdat, params = "pi", 
#                         entropy_measure = "simpson", equivalents = T)
# 
# c(quantile(simps$entropy_pi[[5]], c(0.025, 0.975)),mean(simps$entropy_pi[[5]]))
# c(quantile(outSh$entropy_pi[[5]], c(0.025, 0.975)),mean(outSh$entropy_pi[[5]]))
# 
# #write output, so that one only needs to load these to make the plots, instead of the huge fitted model object.
# save(simps , file = "outITSdeletable_ISDtransformedNoafulva.simps")
# save(outSh , file = "outITSdeletable_ISDtransformedNoafulva.sshan")


#for bacteria:
rm(list=ls())
load("./CNVRG_16S_conservative.RData") #This has to be changed to bacteria for the bacteria plot, obviously.
isdOut <- isd_transform(model_output = modelOut, isd_index=which(names(t_newdat) == "ISD"), countData=t_newdat)

#Quick check to see if anything is interesting with archae. nope
#afulva <- read.table("./data/tax_gg16S.sintax", stringsAsFactors=F, header = F, fill = T)
#archs <- gsub("centroid=(Otu\\d+);seqs.*", "centroid.\\1",afulva[grep("k:Arch", afulva$V2),1])
#afulvaHits <- which(names(t_newdat) %in% archs)
dat <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum", stringsAsFactors = T, header = T)
blanks <- dat[,grep("lank", names(dat))]
contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)
contams <- dat[contams,1][1:2]
contams <- gsub("=", "\\.", contams)
contamHits <- which(names(t_newdat) %in% contams)

outSh <- diversity_calc(model_output = isdOut$pi[,,-contamHits], countData=t_newdat, params = "pi", 
                        entropy_measure = "shannon", equivalents = T)

simps <- diversity_calc(model_output = isdOut$pi[,,-contamHits], countData=t_newdat, params = "pi", 
                        entropy_measure = "simpson", equivalents = T)

save(simps , file = "out16Sdeletable_ISDtransformedNoafulva.simps")
save(outSh , file = "out16Sdeletable_ISDtransformedNoafulva.sshan")

#Example for how to calculate certainty of differences in diversity estimates among treatment groups 
# quantile(out$entropy_pi[[5]] - out$entropy_pi[[6]], c(0.025, 0.975))
# quantile(out$entropy_pi[[5]] - out$entropy_pi[[7]], c(0.025, 0.975))
# quantile(out$entropy_pi[[6]] - out$entropy_pi[[8]], c(0.025, 0.975))
# quantile(out$entropy_pi[[7]] - out$entropy_pi[[8]], c(0.025, 0.975))

###############
#Fig for paper#
###############

#Function to make credible interval lines
limner <- function(coords, index, colz = "black") {
  #add vertical lines
  lines(
    x = c(index, index),
    y = c(coords[[1]], coords[[2]]),
    lty = 1,
    col = colz,
    lwd = 2
  )
  #add horizontal end lines
  lines(
    x = c(index - 0.03,
          index + 0.03),
    y = c(coords[[1]],
          coords[[1]]),
    col = colz,
    lty = 1,
    lwd = 2
  )
  lines(
    x = c(index - 0.03,
          index + 0.03),
    y = c(coords[[2]],
          coords[[2]]),
    col = colz,
    lty = 1,
    lwd = 2
  )
  points(index, coords[[3]], 
         pch = 23,
         col = colz,
         bg = colz,
         cex = 2)
}

load("~/Desktop/teton/outITSdeletable_ISDtransformedNoafulva.simps")
load("~/Desktop/teton/outITSdeletable_ISDtransformedNoafulva.sshan")

pdf(width=5,height=6, file = "./Supplemental_RichDiv_vs_treatment_withControl_ITS.pdf")
par(mfrow = c(2,1),
    oma=c(6,0,3.1,0),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(-0.1,10),
     ylab = "Shannon's",
     xlab = "",
     las = 2,
     xaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))

coords <- c(quantile(outSh$entropy_pi[[5]], c(0.025, 0.975)), mean(outSh$entropy_pi[[5]]))
limner(log2(coords), 1)

coords <- c(quantile(outSh$entropy_pi[[6]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[6]]))
limner(log2(coords), 1.5)
coords <- c(quantile(outSh$entropy_pi[[7]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[7]]))
limner(log2(coords), 2)
coords <- c(quantile(outSh$entropy_pi[[8]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[8]]))
limner(log2(coords), 2.5)

#controls
coords <- c(quantile(outSh$entropy_pi[[1]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[1]]))
limner(log2(coords), 3)
coords <- c(quantile(outSh$entropy_pi[[2]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[2]]))
limner(log2(coords), 3.5)
coords <- c(quantile(outSh$entropy_pi[[3]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[3]]))
limner(log2(coords), 4)
coords <- c(quantile(outSh$entropy_pi[[4]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[4]]))
limner(log2(coords), 4.5)

plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(0.9,2),
     ylab = "Simpson's",
     xlab = "",
     las = 2,
     xaxt = "n",
     yaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))
axis(side=2,
     at = seq(1,2, by = .2),
     labels = c(seq(1,2, by = .2)),
     las = 1
)
axis(side=1,
     at = seq(1,4.5, by = 0.5),
     labels = F
)
text(seq(1,4.5, by = 0.5),
     par("usr")[3]-0.3,
     #pos = 2,
     xpd = NA,
     srt = 45,
     labels = c( expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),

                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = ""))
     )
)

text(3.7,
     labels = "Controls",
     xpd= NA,
     par("usr")[3]-0.6
     )
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 0.8,
     xright = 1.8,
     ybottom = 0.9,
     ytop = 3.45)
text("Inoculum\ntreated",
     x = 1.3,
     y = 3.6)

rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 2.8,
     xright = 3.8,
     ybottom = 0.9,
     ytop = 3.45)
text("Inoculum\ntreated",
     x = 3.3,
     y = 3.6)


coords <- c(quantile(simps$entropy_pi[[5]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[5]]))
limner(coords, 1)

coords <- c(quantile(simps$entropy_pi[[6]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[6]]))
limner(coords, 1.5)
coords <- c(quantile(simps$entropy_pi[[7]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[7]]))
limner(coords, 2)
coords <- c(quantile(simps$entropy_pi[[8]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[8]]))
limner(coords, 2.5)

#controls
coords <- c(quantile(simps$entropy_pi[[1]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[1]]))
limner(coords, 3)
coords <- c(quantile(simps$entropy_pi[[2]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[2]]))
limner(coords, 3.5)
coords <- c(quantile(simps$entropy_pi[[3]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[3]]))
limner(coords, 4)
coords <- c(quantile(simps$entropy_pi[[4]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[4]]))
limner(coords, 4.5)

dev.off()

####################
#No control plotting
####################

pdf(width=5,height=6, file = "./RichDiv_vs_treatment_ITS_Conservative.pdf")
par(mfrow = c(2,1),
    oma=c(6,0,2.5,2),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(6,10),
     ylab = "Shannon's (log 2 scale)",
     xlab = "",
     las = 2,
     xaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))

cz <- c(add.alpha(c("cyan4", "cyan4", "darkgoldenrod2", "darkgoldenrod2"),alpha=0.6))
coords <- c(quantile(outSh$entropy_pi[[5]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[5]]))
limner(log2(coords), 1, colz = cz[1])

coords <- c(quantile(outSh$entropy_pi[[6]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[6]]))
limner(log2(coords), 1.5, colz = cz[2])
coords <- c(quantile(outSh$entropy_pi[[7]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[7]]))
limner(log2(coords), 2, colz = cz[3])
coords <- c(quantile(outSh$entropy_pi[[8]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[8]]))
limner(log2(coords), 2.5, colz = cz[4])

plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(1,1.3),
     ylab = "Simpson's",
     xlab = "",
     las = 2,
     xaxt = "n",
     yaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))
axis(side=2,
     at = seq(1,1.3, by = .1),
     labels = c(seq(1,1.3, by = .1)),
     las = 1
)
axis(side=1,
     at = seq(1,2.5, by = 0.5),
     labels = F
)
text(seq(1,2.5, by = 0.5),
     par("usr")[3]-0.05,
     #pos = 2,
     xpd = NA,
     srt = 45,
     labels = c( expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = ""))
     )
)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 0.8,
     xright = 1.8,
     ybottom = 1,
     ytop = 1.66)
text("Inoculum\ntreated",
     x = 1.3,
     y = 1.7)

coords <- c(quantile(simps$entropy_pi[[5]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[5]]))
limner(coords, 1, cz[1])

coords <- c(quantile(simps$entropy_pi[[6]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[6]]))
limner(coords, 1.5, cz[2])
coords <- c(quantile(simps$entropy_pi[[7]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[7]]))
limner(coords, 2, cz[3])
coords <- c(quantile(simps$entropy_pi[[8]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[8]]))
limner(coords, 2.5, cz[4])

dev.off()

###############
#Plot bacteria#
###############

#Function to make credible interval lines
limner <- function(coords, index, colz = "black") {
  #add vertical lines
  lines(
    x = c(index, index),
    y = c(coords[[1]], coords[[2]]),
    lty = 1,
    col = colz,
    lwd = 2
  )
  #add horizontal end lines
  lines(
    x = c(index - 0.03,
          index + 0.03),
    y = c(coords[[1]],
          coords[[1]]),
    lty = 1,
    col = colz,
    lwd = 2
  )
  lines(
    x = c(index - 0.03,
          index + 0.03),
    y = c(coords[[2]],
          coords[[2]]),
    lty = 1,
    col = colz,
    lwd = 2
  )
  points(index, coords[[3]], 
         pch = 23,
         col = colz,
         bg = colz,
         cex = 2)
}
load("~/Desktop/teton/out16Sdeletable_ISDtransformedNoafulva.simps")
load("~/Desktop/teton/out16Sdeletable_ISDtransformedNoafulva.sshan")
pdf(width=5,height=6, file = "./BacteriaRichDiv_vs_treatment_16S_Conservative.pdf")
par(mfrow = c(2,1),
    oma=c(6,0,2.5,2),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(1.5,3.5),
     ylab = "Shannon's",
     xlab = "",
     las = 2,
     xaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))

cz <- c(add.alpha(c("cyan4", "cyan4", "darkgoldenrod2", "darkgoldenrod2"),alpha=0.6))

coords <- c(quantile(outSh$entropy_pi[[5]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[5]]))
limner((coords), 1, cz[1])


coords <- c(quantile(outSh$entropy_pi[[6]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[6]]))
limner((coords), 1.5, cz[2])


coords <- c(quantile(outSh$entropy_pi[[7]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[7]]))
limner((coords), 2, cz[3])
coords <- c(quantile(outSh$entropy_pi[[8]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[8]]))
limner((coords), 2.5, cz[4])

plot(NULL,
     xlim = c(0.8,2.7),
     ylim = c(1.5,4.5),
     ylab = "Simpson's",
     xlab = "",
     las = 2,
     xaxt = "n",
    # yaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))
# axis(side=2,
#      at = seq(1.5,1.3, by = .1),
#      labels = c(seq(1,1.3, by = .1)),
#      las = 1
# )
axis(side=1,
     at = seq(1,2.5, by = 0.5),
     labels = F
)
text(seq(1,2.5, by = 0.5),
     par("usr")[3]-0.4,
     #pos = 2,
     xpd = NA,
     srt = 45,
     labels = c( expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = ""))
     )
)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 0.8,
     xright = 1.8,
     ybottom = 1.5,
     ytop = 8)
text("Inoculum\ntreated",
     x = 1.3,
     y = 8.7)

coords <- c(quantile(simps$entropy_pi[[5]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[5]]))
limner(coords, 1, cz[1])

coords <- c(quantile(simps$entropy_pi[[6]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[6]]))
limner(coords, 1.5, cz[2])
coords <- c(quantile(simps$entropy_pi[[7]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[7]]))
limner(coords, 2, cz[3])
coords <- c(quantile(simps$entropy_pi[[8]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[8]]))
limner(coords, 2.5, cz[4])

dev.off()


pdf(width=5,height=6, file = "./visuals/Supplemental_RichDiv_vs_treatment_withControl_16S.pdf")
par(mfrow = c(2,1),
    oma=c(6,0,2.5,0),
    mar = c(1,5,0,0),
    xpd = NA)
plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(0.5,2.0),
     ylab = "Shannon's",
     xlab = "",
     las = 2,
     xaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))

coords <- c(quantile(outSh$entropy_pi[[5]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[5]]))
limner(log2(coords), 1)

coords <- c(quantile(outSh$entropy_pi[[6]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[6]]))
limner(log2(coords), 1.5)
coords <- c(quantile(outSh$entropy_pi[[7]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[7]]))
limner(log2(coords), 2)
coords <- c(quantile(outSh$entropy_pi[[8]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[8]]))
limner(log2(coords), 2.5)

#controls
coords <- c(quantile(outSh$entropy_pi[[1]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[1]]))
limner(log2(coords), 3)
coords <- c(quantile(outSh$entropy_pi[[2]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[2]]))
limner(log2(coords), 3.5)
coords <- c(quantile(outSh$entropy_pi[[3]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[3]]))
limner(log2(coords), 4)
coords <- c(quantile(outSh$entropy_pi[[4]], c(0.025, 0.975)),
            mean(outSh$entropy_pi[[4]]))
limner(log2(coords), 4.5)

plot(NULL,
     xlim = c(0.8,4.7),
     ylim = c(1.5,5.5),
     ylab = "Simpson's",
     xlab = "",
     las = 2,
     xaxt = "n",
     #yaxt = "n",
     bty = "n",
     cex.lab = 1.3,
     mar=c(0,0,0,0))
# axis(side=2,
#      at = seq(1,2, by = .2),
#      labels = c(seq(1,2, by = .2)),
#      las = 1
# )
axis(side=1,
     at = seq(1,4.5, by = 0.5),
     labels = F
)
text(seq(1,4.5, by = 0.5),
     par("usr")[3]-0.8,
     #pos = 2,
     xpd = NA,
     srt = 45,
     labels = c( expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste("No ", italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; inoculum",
                                  sep = "")),
                 
                 expression(paste("No ", italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = "")),
                 expression(paste(italic("A. fulva"),
                                  #"; no inoculum",
                                  sep = ""))
     )
)

text(3.7,
     labels = "Controls",
     xpd= NA,
     par("usr")[3]-1.7
)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 0.8,
     xright = 1.8,
     ybottom = 1.5,
     ytop = 10.3)
text("Inoculum\ntreated",
     x = 1.3,
     y = 10.7)

rect(col=add.alpha("cyan", alpha = 0.2),border=NA,
     xleft = 2.8,
     xright = 3.8,
     ybottom = 1.5,
     ytop = 10.3)
text("Inoculum\ntreated",
     x = 3.3,
     y = 10.7)


coords <- c(quantile(simps$entropy_pi[[5]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[5]]))
limner(coords, 1)

coords <- c(quantile(simps$entropy_pi[[6]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[6]]))
limner(coords, 1.5)
coords <- c(quantile(simps$entropy_pi[[7]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[7]]))
limner(coords, 2)
coords <- c(quantile(simps$entropy_pi[[8]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[8]]))
limner(coords, 2.5)

#controls
coords <- c(quantile(simps$entropy_pi[[1]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[1]]))
limner(coords, 3)
coords <- c(quantile(simps$entropy_pi[[2]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[2]]))
limner(coords, 3.5)
coords <- c(quantile(simps$entropy_pi[[3]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[3]]))
limner(coords, 4)
coords <- c(quantile(simps$entropy_pi[[4]], c(0.025, 0.975)),
            mean(simps$entropy_pi[[4]]))
limner(coords, 4.5)

dev.off()
