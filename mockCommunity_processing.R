#J. Harrison
set.seed(666)
options(scipen=99)
rm(list=ls())
dat <- read.table("./data/otutableJOINED_ISD_ASLE_duds_sum",stringsAsFactors = F, header = T)
head(dat)

contam <- read.csv("./data/fungal_contams.csv", stringsAsFactors = F)
dat <- dat[!(dat$OTUID %in% contam$x),] 

hits <- dat[rowSums(dat[,grep("zymo", names(dat))]) > 0,1]

tax <- read.table("./data/tax_uniteJOINED_fungiOnly.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)

tax[tax$V1 %in% hits,]


#BACTERIA
set.seed(666)
options(scipen=99)
rm(list=ls())
dat <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum",stringsAsFactors = F, header = T)
head(dat)

contam <- read.csv("./data/bacterial_contams.csv", stringsAsFactors = F)
dat <- dat[!(dat$OTUID %in% contam$x),] 

hits <- dat[rowSums(dat[,grep("zymo", names(dat))]) > 0,1]

tax <- read.table("./data/tax_gg16S.sintax",
                  fill = T, stringsAsFactors = F)

tax$V1 <- gsub("(centroid=Otu\\d+);.*", "\\1",tax$V1)

tax[tax$V1 %in% hits,]
