tax <- read.table("./data/tax_uniteJOINED_fungiOnly.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
dat <- read.table("./data/otutableJOINED_ISD_ASLE_duds_sum",stringsAsFactors = F, header = T)
head(dat)

contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
dat <- dat[!(dat$OTUID %in% contam$x),] 

present <- tax[which(tax$V1 %in% dat$OTUID),]
dim(present)
tax <- strsplit(x = present$V4, split = ",")

table(sapply(tax,"[", 2))
table(sapply(tax,"[", 3))

sum(rowSums(dat[1:(length(dat[,1])-3),2:length(dat)]))

#################
#Do for bacteria#
#################

tax <- read.table("./data/tax_gg16S.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
dat <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum",stringsAsFactors = F, header = T)
head(dat)

contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)
dat <- dat[!(dat$OTUID %in% contam$x),] 

tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

present <- tax[which(tax$V1 %in% dat$OTUID),]
dim(present)
tax <- strsplit(x = present$V4, split = ",")

table(sapply(tax,"[", 2))

sum(rowSums(dat[1:(length(dat[,1])-4),2:length(dat)]))


dat$OTUID[which.max(rowSums(dat[1:(length(dat[,1])-4),2:length(dat)]))]
sort(rowSums(dat[1:(length(dat[,1])-4),2:length(dat)]))

dat$OTUID[119]
dat$OTUID[108]
dat$OTUID[188]
dat$OTUID[283]

tax <- read.table("./data/tax_gg16S.sintax",
                  fill = T, stringsAsFactors = F)
tax[grep("Otu14;", tax$V1),]
tax[grep("Otu15;", tax$V1),]
tax[grep("Otu17;", tax$V1),]
tax[grep("Otu20;", tax$V1),]


#############################
# Make biodiversity figures #
#############################
rm(list=ls())
metadata <- read.csv("./data/trait_and_treatment_data.csv",stringsAsFactors = F, header = T)
metadata <- metadata[metadata$treatment_failed == "no",]

dat <- read.csv("./data/ml_p_table_16s_conservative.csv",stringsAsFactors = F, header = T)
head(dat)
names(dat) <- gsub("\\.","=", names(dat))

#convert to ratios with ISD
dat[,3:length(dat)] <- dat[,3:length(dat)] /dat$ISD

#remove contaminants
contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)
dat <- dat[,!(names(dat) %in% contam$x)] 
dat <- dat[,-((length(dat)-3):length(dat))]
rownames(dat) <- dat$sample
dat <- dat[,-(1:2)]

#wrangle taxonomy
tax <- read.table("./data/tax_gg16S.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

present <- tax[which(tax$V1 %in% names(dat)),]
dim(present)
taxNames <- present$V1

#make sure taxonomy file is in the same order as the data
present <- present[match(names(dat), present$V1),]
present$V1 == names(dat)

tax <- strsplit(x = present$V4, split = ",")

names(tax) <- taxNames

#Calculate taxon proportions by treatment group, first merge data and metadata
dat$sample <- row.names(dat)
merge_dat <- merge(dat, metadata, by.x = "sample", by.y = "plant")

#can calculate mean per group like this
#aggregate(merge_dat$`centroid=Otu1000` ~ merge_dat$treament, FUN = mean)

#however we need to asssign taxonomy to taxa and lump them by phyla first
taxMeans <- data.frame(matrix(nrow = 8, ncol = 1))
k <- 1
for(i in names(table(sapply(tax,"[", 2)))){
 hits <-  tax[grep(i, sapply(tax,"[", 2))]
 if(length(which(names(merge_dat) %in% names(hits))) > 1){
  sumHits <- rowSums(merge_dat[,names(merge_dat) %in% names(hits)])
 }else{
   sumHits <- merge_dat[,names(merge_dat) %in% names(hits)]
 }
 meansHits <- aggregate(sumHits ~ merge_dat$treament, FUN = mean)
 taxMeans[,k] <- meansHits$sumHits
 names(taxMeans)[k] <- i
 k <- k + 1
}
#doublechecked by running aggregate that the treatment order is:
#treated -, treated +, untreated -, untreated +

#change colors and names. make suitable for main fig. 
#maybe color L. tuarica and A. fulva separately.
#the idea being to pair with rare taxa and culturing figs in a three part figure

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

names(taxMeans) <- gsub("p:", "", names(taxMeans))
names(taxMeans) <- gsub("\\[", "", names(taxMeans))
names(taxMeans) <- gsub("]", "", names(taxMeans))

pdf(width = 7, height = 10, file = "./visuals/raretaxaBarplot16s.pdf")
par(mfrow =c(2,1), mai = c(0.8,2,2.2,0))

barplot(height = unlist(taxMeans[5:8,]), las = 2,
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
        border = add.alpha("white", alpha = 0.001),
        names = "")
text("ratio with ISD", srt = 90, x = -21, y = 7, cex = 2, xpd = NA)

legend(x = "topright",
       legend = c(expression(paste("Inoc. treated, ", italic("A. fulva")," -",
                                    sep = "")),
                  expression(paste("Inoc. treated, ", italic("A. fulva")," +",
                                    sep = "")),
                  expression(paste("No inoculum, ", italic("A. fulva")," -",
                                    sep = "")),
                  expression(paste("No inoculum, ", italic("A. fulva")," +",
                                    sep = ""))),
       pch = rep(15,4),
       col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
       pt.cex = 2,
       cex = .8,
       bty = "n"
)

###################
#Do for rel. abund. 

dat <- read.csv("./data/ml_p_table_16s_conservative.csv",stringsAsFactors = F, header = T)
names(dat) <- gsub("\\.","=", names(dat))

#remove contaminants
contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)
dat <- dat[,!(names(dat) %in% contam$x)] 
dat <- dat[,-((length(dat)-3):length(dat))]
rownames(dat) <- dat$sample
dat <- dat[,-(1:2)]

#wrangle taxonomy
tax <- read.table("./data/tax_gg16S.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

present <- tax[which(tax$V1 %in% names(dat)),]
dim(present)
taxNames <- present$V1

#make sure taxonomy file is in the same order as the data
present <- present[match(names(dat), present$V1),]
present$V1 == names(dat)

tax <- strsplit(x = present$V4, split = ",")

names(tax) <- taxNames

#Calculate taxon proportions by treatment group, first merge data and metadata
dat$sample <- row.names(dat)
merge_dat <- merge(dat, metadata, by.x = "sample", by.y = "plant")

#can calculate mean per group like this
#aggregate(merge_dat$`centroid=Otu1000` ~ merge_dat$treament, FUN = mean)

#however we need to asssign taxonomy to taxa and lump them by phyla first
taxMeans <- data.frame(matrix(nrow = 8, ncol = 1))
k <- 1
for(i in names(table(sapply(tax,"[", 2)))){
  hits <-  tax[grep(i, sapply(tax,"[", 2))]
  if(length(which(names(merge_dat) %in% names(hits))) > 1){
    sumHits <- rowSums(merge_dat[,names(merge_dat) %in% names(hits)])
  }else{
    sumHits <- merge_dat[,names(merge_dat) %in% names(hits)]
  }
  meansHits <- aggregate(sumHits ~ merge_dat$treament, FUN = mean)
  taxMeans[,k] <- meansHits$sumHits
  names(taxMeans)[k] <- i
  k <- k + 1
}

names(taxMeans) <- gsub("p:", "", names(taxMeans))
names(taxMeans) <- gsub("\\[", "", names(taxMeans))
names(taxMeans) <- gsub("]", "", names(taxMeans))

par(mai = c(3,2,0,0))
barplot(height = unlist(taxMeans[5:8,]), las = 2,
        ylim = c(0,0.005),
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
        border = add.alpha("white", alpha = 0.001),
        names = "")
text("Proportional rel. abund.", srt = 90, x = -21, y = .002, cex = 2, xpd = NA)

locs <- seq(3.5,86.3, by = 86.3/18)
for(i in 1:length(locs)){
  text(names(taxMeans)[i], x = locs[i], y = -0.0004, xpd = NA, srt = 70, pos = 2)
}

text("a)", x = -30, cex =2, y = 0.012, xpd = NA)
text("b)", x = -30, cex =2, y = 0.006, xpd = NA)

dev.off()


#############################
# Make biodiversity figure for ITS #
#############################
rm(list=ls())
metadata <- read.csv("./data/trait_and_treatment_data.csv",stringsAsFactors = F, header = T)
metadata <- metadata[metadata$treatment_failed == "no",]
dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv",stringsAsFactors = F, header = T)
head(dat)
names(dat) <- gsub("\\.","=", names(dat))

#convert to ratios with ISD
dat[,3:length(dat)] <- dat[,3:length(dat)] /dat$ISD

#remove contaminants
contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
dat <- dat[,!(names(dat) %in% contam$x)] 
dat <- dat[,-((length(dat)-3):length(dat))]
rownames(dat) <- dat$sample
dat <- dat[,-(1:2)]

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F) 
levtaurHits <- which(gsub("centroid=","",names(dat)) %in% levtaur$V1)

afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
afulvaHits <- which(gsub("centroid=","",names(dat)) %in% afulva$V1)

dat <- dat[,-c(levtaurHits, afulvaHits)]

#wrangle taxonomy
tax <- read.table("./data/tax_uniteJOINED_fungiOnly.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

present <- tax[which(tax$V1 %in% names(dat)),]
dim(present)
taxNames <- present$V1

#make sure taxonomy file is in the same order as the data
present <- present[match(names(dat), present$V1),]
present$V1 == names(dat)

tax <- strsplit(x = present$V4, split = ",")

names(tax) <- taxNames

#Calculate taxon proportions by treatment group, first merge data and metadata
dat$sample <- row.names(dat)
merge_dat <- merge(dat, metadata, by.x = "sample", by.y = "plant")

#can calculate mean per group like this
#aggregate(merge_dat$`centroid=Otu1000` ~ merge_dat$treament, FUN = mean)

#however we need to asssign taxonomy to taxa and lump them by phyla first
taxMeans <- data.frame(matrix(nrow = 8, ncol = 1))
k <- 1
for(i in names(table(sapply(tax,"[", 2)))){
  hits <-  tax[grep(i, sapply(tax,"[", 2))]
  if(length(which(names(merge_dat) %in% names(hits))) > 1){
    sumHits <- rowSums(merge_dat[,names(merge_dat) %in% names(hits)])
  }else{
    sumHits <- merge_dat[,names(merge_dat) %in% names(hits)]
  }
  meansHits <- aggregate(sumHits ~ merge_dat$treament, FUN = mean)
  taxMeans[,k] <- meansHits$sumHits
  names(taxMeans)[k] <- i
  k <- k + 1
}
#doublechecked by running aggregate that the treatment order is:
#treated -, treated +, untreated -, untreated +

#change colors and names. make suitable for main fig. 
#maybe color L. tuarica and A. fulva separately.
#the idea being to pair with rare taxa and culturing figs in a three part figure

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

names(taxMeans) <- gsub("p:", "", names(taxMeans))
names(taxMeans) <- gsub("\\[", "", names(taxMeans))
names(taxMeans) <- gsub("]", "", names(taxMeans))

pdf(width = 7, height = 10, file = "./visuals/raretaxaBarplotITS.pdf")
par(mfrow =c(2,1), mai = c(0.8,2,2.2,0))

barplot(height = unlist(taxMeans[5:8,]), las = 2,
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
        border = add.alpha("white", alpha = 0.001),
        names = "")
text("ratio with ISD", srt = 90, x = -21, y = 7, cex = 2, xpd = NA)

legend(x = "topright",
       legend = c(expression(paste("Inoc. treated, ", italic("A. fulva")," -",
                                   sep = "")),
                  expression(paste("Inoc. treated, ", italic("A. fulva")," +",
                                   sep = "")),
                  expression(paste("No inoculum, ", italic("A. fulva")," -",
                                   sep = "")),
                  expression(paste("No inoculum, ", italic("A. fulva")," +",
                                   sep = ""))),
       pch = rep(15,4),
       col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
       pt.cex = 2,
       cex = .8,
       bty = "n"
)

###################
#Do for rel. abund. 

dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv",stringsAsFactors = F, header = T)
names(dat) <- gsub("\\.","=", names(dat))

#remove contaminants
contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
dat <- dat[,!(names(dat) %in% contam$x)] 
dat <- dat[,-((length(dat)-3):length(dat))]
rownames(dat) <- dat$sample
dat <- dat[,-(1:2)]

levtaur <- read.table("data/levtaurMatchesITS", stringsAsFactors=F, header = F) 
levtaurHits <- which(gsub("centroid=","",names(dat)) %in% levtaur$V1)

afulva <- read.table("data/afulvaMatchesITSCOMBO", stringsAsFactors=F, header = F) 
afulvaHits <- which(gsub("centroid=","",names(dat)) %in% afulva$V1)

dat <- dat[,-c(levtaurHits, afulvaHits)]

#wrangle taxonomy
tax <- read.table("./data/tax_uniteJOINED_fungiOnly.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

present <- tax[which(tax$V1 %in% names(dat)),]
dim(present)
taxNames <- present$V1

#make sure taxonomy file is in the same order as the data
present <- present[match(names(dat), present$V1),]
present$V1 == names(dat)

tax <- strsplit(x = present$V4, split = ",")

names(tax) <- taxNames

#Calculate taxon proportions by treatment group, first merge data and metadata
dat$sample <- row.names(dat)
merge_dat <- merge(dat, metadata, by.x = "sample", by.y = "plant")

#can calculate mean per group like this
#aggregate(merge_dat$`centroid=Otu1000` ~ merge_dat$treament, FUN = mean)

#however we need to asssign taxonomy to taxa and lump them by phyla first
taxMeans <- data.frame(matrix(nrow = 8, ncol = 1))
k <- 1
for(i in names(table(sapply(tax,"[", 2)))){
  hits <-  tax[grep(i, sapply(tax,"[", 2))]
  if(length(which(names(merge_dat) %in% names(hits))) > 1){
    sumHits <- rowSums(merge_dat[,names(merge_dat) %in% names(hits)])
  }else{
    sumHits <- merge_dat[,names(merge_dat) %in% names(hits)]
  }
  meansHits <- aggregate(sumHits ~ merge_dat$treament, FUN = mean)
  taxMeans[,k] <- meansHits$sumHits
  names(taxMeans)[k] <- i
  k <- k + 1
}

names(taxMeans) <- gsub("p:", "", names(taxMeans))
names(taxMeans) <- gsub("\\[", "", names(taxMeans))
names(taxMeans) <- gsub("]", "", names(taxMeans))

par(mai = c(3,2,0,0))
barplot(height = unlist(taxMeans[5:8,]), las = 2,
        #ylim = c(0,0.005),
        col = add.alpha(c("cyan2", "cyan4","darkgoldenrod1","darkgoldenrod4"), alpha = 0.8),
        border = add.alpha("white", alpha = 0.001),
        names = "")
text("Proportional rel. abund.", srt = 90, x = -21, y = .002, cex = 2, xpd = NA)

locs <- seq(3.5,86.3, by = 86.3/18)
for(i in 1:length(locs)){
  text(names(taxMeans)[i], x = locs[i], y = -0.0004, xpd = NA, srt = 70, pos = 2)
}

text("a)", x = -30, cex =2, y = 0.012, xpd = NA)
text("b)", x = -30, cex =2, y = 0.006, xpd = NA)

dev.off()
