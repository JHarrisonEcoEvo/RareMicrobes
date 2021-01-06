rm(list=ls())
dat <- read.csv("./data/ml_p_table_16s_conservative.csv")
metadat <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F)
metadat<- metadat[metadat$treatment_failed=="no",]
contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)
dat <- dat[,-which(names(dat) %in% gsub("=","\\.", contam$x))]

dat <- merge(dat, metadat, by.x = "sample", by.y="plant")
which.max(colSums(dat[,3:642]))

taxon <- dat[, 454] / dat$ISD

boxplot(taxon ~ dat$treament)
#do for other taxa
sort(colSums(dat[,3:642])) #centroid.Otu79   centroid.Otu66   centroid.Otu60 

taxon <- dat[, which(names(dat) == "centroid.Otu79")] / dat$ISD

boxplot(taxon ~ dat$treament, outline = F)


taxon <- dat[, which(names(dat) == "centroid.Otu66")] / dat$ISD

boxplot(taxon ~ dat$treament, outline = F)

taxon <- dat[, which(names(dat) == "centroid.Otu60")] / dat$ISD

boxplot(taxon ~ dat$treament, outline = F)
