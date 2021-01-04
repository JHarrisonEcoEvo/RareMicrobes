rm(list=ls())
dat2 <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv", stringsAsFactors = F)

dat <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum", stringsAsFactors = T, header = T)
blanks <- dat[,grep("lank", names(dat))]
contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)
contams <- dat[contams,1][1:2]
contams <- gsub("=", "\\.", contams)
contams <- gsub("centroid\\.(Otu\\d+)", "\\1", contams)
contamHits <- which(names(dat2) %in% contams)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#transform by the ISD
for(i in 1:length(dat2[,1])){
  dat2[i,3:length(dat2)]  <- dat2[i,3:length(dat2)] / dat2$ISD[i]
}

#remove contaminants
dat2 <- dat2[,-contamHits]
dat2[1:5,1:5]

#Remove A. fulva and make indices to remove L. taurica
afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
afulvaHits <- which(names(dat2) %in% afulva$V1)
dat2 <- dat2[,-afulvaHits]

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F) 
levtaurHits <- which(names(dat2) %in% levtaur$V1)
dat2 <- dat2[,-levtaurHits]


#calculate richness proxy
a <- vector()
for(i in 1:length(dat2[,1])){
  a[i] <- length(which(100*dat2[i,3:(length(dat2)-2)] > 0.03))
}

#merge with treatment data
dat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)

dat <- merge(dat, dat2, by.x = "plant", by.y = "sample")

dat <- dat[match(dat2$sample,dat$plant),]
table(dat$plant== dat2$sample)

controls <- grep("ontrol", dat$treament)

pdf(width=7, height = 13, file = "./visuals/richnessFungBact25_20.pdf")

par(oma = c(7,6,1,1), mfrow = c(2,1), mar = c(3,3,0,0))

aggregate(a[-controls] ~ dat$treament[-controls], FUN = mean)

stripchart(a[-controls] ~ dat$treament[-controls],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("lightgray", alpha=0.5),
           cex = 2.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
          ,ylim=c(0,100)
           ,frame.plot=F
)
boxplot(a[-controls] ~ dat$treament[-controls], las = 2, xlab = "", ylab = "fungi", outline = F,
        names = c("","", "", ""),
        add = T,
        frame.plot=F,
        col = c(add.alpha("lightgray", alpha=0.5),
                add.alpha("lightgray", alpha=0.5),
                add.alpha("white", alpha=0.5),
                add.alpha("white", alpha=0.5)))

text("fungi",
     xpd =NA,
     cex = 2,
     x = -0.5,
     srt = 90,
     y = 1000)

#bacteria
dat2 <- read.csv("./data/ml_p_table_16s_conservative.csv", stringsAsFactors = F)
dat <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum", stringsAsFactors = T, header = T)
blanks <- dat[,grep("lank", names(dat))]
contams <- which((rowSums(blanks) / rowSums(dat[,2:length(dat)])) > 0.05)
contams <- dat[contams,1][1:2]
contams <- gsub("=", "\\.", contams)
#contams <- gsub("centroid\\.(Otu\\d+)", "\\1", contams)
contamHits <- which(names(dat2) %in% contams)

for(i in 1:length(dat2[,1])){
  dat2[i,3:length(dat2)]  <- dat2[i,3:length(dat2)] / dat2$ISD[i]
}

dat2 <- dat2[,-contamHits]
dat2[1:5,1:5]

b <- vector()
for(i in 1:length(dat2[,1])){
  b[i] <- length(which(100*dat2[i,3:(length(dat2)-2)] > 20))
}

dat1 <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)

dat1 <- merge(dat1, dat2, by.x = "plant", by.y = "sample")

dat1 <- dat1[match(dat2$sample,dat1$plant),]
table(dat1$plant== dat2$sample)

controls1 <- grep("ontrol", dat1$treament)

stripchart(a[-controls1] ~ dat1$treament[-controls1],
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = add.alpha("lightgray", alpha=0.5),
           cex = 2.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           xlim=c(0.5,4.5)
           ,ylim=c(0,2000)
           ,frame.plot=F
)
boxplot(a[-controls1] ~ dat1$treament[-controls1], las = 2, xlab = "", ylab = "bacteria", outline = F,
        names = c(expression(paste("Inoc. treated, ", italic('A. fulva')," -")),
                  expression(paste("Inoc. treated, ", italic('A. fulva')," +")),
                  expression(paste("No inoculum, ", italic('A. fulva')," -")),
                  expression(paste("No inoculum, ", italic('A. fulva')," -"))),
        add = T,
        frame.plot=F,
        xpd = NA,
        col = c(add.alpha("lightgray", alpha=0.5),
                add.alpha("lightgray", alpha=0.5),
                add.alpha("white", alpha=0.5),
                add.alpha("white", alpha=0.5)))

text("Proxy for richness",
     xpd =NA,
     cex = 2,
     x = -1,
     srt = 90,
     y = 2100)

text("bacteria",
     xpd =NA,
     cex = 2,
     x = -0.5,
     srt = 90,
     y = 1000)


dev.off()

#Try using vegan, bc our method seems capricious

otus <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum", stringsAsFactors = F, header = T)

metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)

metadat <- metadat[order(metadat$treament),]
treatment <- metadat$treament

#the order of my OTU file needs to follow this:
#metadat$plant

#sanity check
#cbind(gsub("X16S_\\w+_\\w+_Harrison_asle_(\\d+)asle", "\\1", names(dat)), names(dat))
newnames <- gsub("ITS_\\w+_\\w+_Harrison_asle_(\\d+).*", "\\1", names(dat))

#combine PCR duplicates
newdat <- data.frame(matrix(nrow = dim(dat)[1]))
k <- 1
for_names <- vector()
for(i in unique(newnames)){
  if(length(grep(paste("^",i,"$", sep = ""), newnames)) > 1){ #this weeds out the non-plant stuff
    newdat[,k] <- rowSums(dat[,grep(paste("^",i,"$", sep = ""), newnames)])
    for_names <- c(for_names, i)
    k <- k + 1
  }else{next}
}

names(newdat) <- for_names

#now need to make the data follow the order of the metadata
#sanity check, compare these
#names(newdat)[match(metadat$plant, as.numeric(names(newdat)))]
#metadat$plant
#dim(newdat)

#These sequenced plants are not in the metadata
names(newdat)[!as.numeric(names(newdat)) %in% metadat$plant]
#we good

#next we need to remove those plants that are somewhat suspect. This step is optional.
metadat <- metadat[metadat$conservative_treatmentFailed == "no",]
newdat <- newdat[,names(newdat) %in% metadat$plant]

#need to get newdat into the same order as the metadata
metadat <- metadat[!duplicated(metadat$plant),]

newdat <- newdat[,na.omit(match(metadat$plant, as.numeric(names(newdat))))]

#trim metadata to remove the stuff that didnt get sequenced
metadat <- metadat[metadat$plant %in% names(newdat),]

#doublecheck that the metadat plant follows the names of the data
table(metadat$plant == names(newdat))
cbind(metadat$plant[metadat$plant != names(newdat)], names(newdat)[metadat$plant != names(newdat)])

#Next need to get the newdat df into the format needed for CNVRG and determine start and end indices

# transpose and put the names in place
t_newdat <- as.data.frame(t(newdat))
colnames(t_newdat) <- dat$OTUID
dim(t_newdat)
t_newdat <- data.frame(metadat$plant, t_newdat)
#metadat$treament

#remove taxa with zero counts
t_newdat <- t_newdat[,colSums(t_newdat[,2:length(t_newdat)]) > 0]

rich <- vector()
spec <- vector()
for(i in 1:length(t_newdat[,1])){
  rich[i] <- vegan::estimateR(t_newdat[i,2:length(t_newdat)])[4]
  spec[i] <- vegan::specpool(t_newdat[i,2:length(t_newdat)])[1]
}
boxplot(unlist(spec) ~ metadat$treament)
