rm(list=ls())
tax <- read.table("./data/tax_uniteJOINED_fungiOnly.sintax",
                  fill = T, stringsAsFactors = F)
head(tax)
dat <- read.table("./data/otutableJOINED_ISD_ASLE_AFULVA_sum",stringsAsFactors = F, header = T)
head(dat)

dat <- dat[,c(1,grep("inoculum", names(dat)))]

dat <- dat[rowSums(dat[,2:length(dat)]) > 0,]

dat$OTUID #all of these except OTU 91 show up in the miseq data.
# 
##########################################
# determine efficacy of inoculum FOR ITS #
##########################################
rm(list=ls())
#pull in ml estimates, inocula taxa, and treatment info
dat <- read.csv("./data/ml_p_table_COMBO_ITS_conservative.csv", stringsAsFactors = F, header = T)
#convert to absolute abundances, or not. Doing both for plot
dat2 <- data.frame(dat$sample, dat[,3:length(dat)])
dat <- data.frame(dat$sample, dat[,3:length(dat)] / dat$ISD)

dim(dat)

#sanitycheck, using different object for the transform
#abs_dat[1,1:3]
#dat[1,1:5] / dat$ISD[1]

inocula <- read.table("./data/linking_ITS_inoculumOTUs_to_NOVAseqOTUS", stringsAsFactors = F, header = F)
traits <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F, header = T)
traits <- traits[!duplicated(traits$plant),]
traits <- traits[traits$treatment_failed == "no",]

#merge ml with treatment
dim(dat)
merge_dat <- merge(dat, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)
merge_dat2 <- merge(dat2, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)

dim(merge_dat)

# #why dont the dimensions match? FIXED. Turns out dupes in the traits data.
# dat$sample  %in%  traits$plant
# merge_dat$sample %in% dat$sample
# duplicated(merge_dat$sample) #ah some dupes
# merge_dat[45:50,1:10]
# duplicated(traits$plant) #they appear to be in og data. 

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

inoculum_df <- merge_dat[,c(1,which(names(merge_dat) == "treament") , which(names(merge_dat) %in% inoculum))]
inoculum_df2 <- merge_dat2[,c(1,which(names(merge_dat2) == "treament") , which(names(merge_dat2) %in% inoculum))]

head(inoculum_df)
inoculum_df$sums <- rowSums(inoculum_df[,3:length(inoculum_df)])
inoculum_df2$sums <- rowSums(inoculum_df2[,3:length(inoculum_df2)])

dim(inoculum_df)



# #make boxplot
# boxplot(inoculum_df$sums ~ inoculum_df$treament, outline = F, ylim = c(0,500))
# 


#make a whole panel of boxplots
#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf(width = 6, height = 10, file = "./visuals/inoculumITS_by_ISD.pdf")

par(mfrow = c(4,2), oma = c(4,3,1,1))
for(i in 3:(length(inoculum_df)-1)){
boxplot(log10(inoculum_df[,i]) ~ inoculum_df$treament, outline = F, las = 2, ylab = "", xlab = "",
        names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                  "Treated -","Treated +", "Untreated -", "Untreated +"), xpd = NA,
        col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ))
  title(names(inoculum_df[i]),  line = -2, adj = 0.1)
#relative abundance plot
  boxplot(inoculum_df2[,i] ~ inoculum_df2$treament, outline = F, las = 2, ylab = "", xlab = "",
          names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                    "Treated -","Treated +", "Untreated -", "Untreated +"), xpd = NA,
          col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ))
  title(paste(names(inoculum_df2[i]), "rel. abund.", sep = " "),  line = -2, adj = 0.1)
  
}
boxplot(log10(inoculum_df$sums) ~ inoculum_df$treament,  outline = F, las = 2, ylab = "", xlab = "",
        names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                  "Treated -","Treated +", "Untreated -", "Untreated +"),
        col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ))
title("Summed",  line = -2, adj = 0.1)
boxplot(inoculum_df2$sums ~ inoculum_df2$treament,  outline = F, las = 2, ylab = "", xlab = "",
        names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                  "Treated -","Treated +", "Untreated -", "Untreated +"),
        col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ))
title("Summed (rel. abund.)",  line = -2, adj = 0.1)
text( "Proportion / ISD (log10)", xpd = NA, x = -7.8, y = 7.6, cex = 1.8 )
text( "Relative abundance", xpd = NA, x = 4, y = 7.6, cex = 1.8 )

dev.off()



diffs <- read.csv("./data/diffs_ITS.csv", stringsAsFactors = F)
diffs_df <- data.frame(diffs$comparison, diffs[,names(diffs) %in% inoculum])

unique(sort(traits$treament)) #so we care about 5 vs 7 and 6 vs 8 

diffs_df[diffs_df$diffs.comparison %in% c("treatment_5_vs_treatment_7","treatment_6_vs_treatment_8"),]


diff_5v7 <- vector()
diff_6v8 <- vector()
k <- 1
for(i in 3:length(inoculum_df)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)

  diff_5v7[k] <- avgs[5,2] - avgs[7,2]
  diff_6v8[k] <- avgs[6,2] - avgs[8,2]
  k <- k + 1
}

table(diff_5v7 < 0)
table(diff_5v7 < 0)[1] / sum(table(diff_5v7 < 0))
table(diff_6v8 < 0)


pdf(width = 10, height = 8, file = "./visuals/fungalInoculumGraph.pdf")
par(mfrow = c(2,1), mar = c(1,4,0,1), oma = c(2,6,0,0))
plot(NULL,
     xlim = c(0, 3.5),
     ylim = c(0,1),
     las = 2,
     bty = "n",
     ylab = "",
     xlab = "",
     xaxt = "n", 
     cex.lab = 1.2)
k <- 1
for(i in 3:(length(inoculum_df2)-1)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)
  print(i)
  points(k, avgs[8,2], col = add.alpha("darkgoldenrod2",0.6), pch = 16, cex = 2)
  points(k, avgs[6,2], col = add.alpha("cyan4",0.6), pch = 16, cex = 2)
  k <- k + 1
}
text(substitute(paste( italic('A. fulva')," present", sep = "")), srt = 0, xpd = NA, x = 3, y = 0.9, cex = 1.5)
legend(x = 0.2, y = 0.97, bty = "n", y.intersp = 1, legend = c("Inoculated", "Not inoculated"), pch = 16, col = c(add.alpha("cyan4",0.6),
                                                                                                                  add.alpha("darkgoldenrod2",0.6)), cex = 1.5)
axis(side = 1, at = c(0,1,2,3, 3.5), labels = c("","","", "", ""))
plot(NULL,
     xlim = c(0, 3.5),
     ylim = c(0,1),
     las = 2,
     ylab = "",
     xlab = "",
     xaxt = "n",
     bty = "n",
     cex.lab = 1.2)
k <- 1
for(i in 3:(length(inoculum_df2)-1)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)
  print(i)
  points(k, avgs[7,2], col = add.alpha("darkgoldenrod2",0.6), pch = 16, cex = 2)
  points(k, avgs[5,2], col = add.alpha("cyan4",0.6), pch = 16, cex = 2)
  k <- k + 1
}
text("Proportion / ISD", srt = 90, xpd = NA, x = -0.6, y = 1.1, cex = 2)
text(substitute(paste("No ", italic('A. fulva'), sep = "")), srt = 0, xpd = NA, x = 3, y = 0.9, cex = 1.5)
axis(side = 1, at = c(0,1,2,3, 3.5), labels = c("","","", "", ""))

dev.off()

##########################################
# determine efficacy of inoculum FOR 16S #
##########################################
rm(list=ls())
#pull in ml estimates, inocula taxa, and treatment info
dat <- read.csv("./data/ml_p_table_16s_conservative.csv", stringsAsFactors = F, header = T)
#convert to absolute abundances
dat2 <- data.frame(dat$sample, dat[,3:length(dat)])
dat <- data.frame(dat$sample, dat[,3:length(dat)] / dat$ISD)
#sanitycheck, using different object for the transform
#abs_dat[1,1:3]
#dat[1,1:5] / dat$ISD[1]

inocula <- read.table("./data/linking_16S_inoculumOTUs_to_NOVAseqOTUS", stringsAsFactors = F, header = F)
length(unique(inocula$V1)) #206 taxa with 13823 matches?
traits <- read.csv("./data/trait_and_treatment_data.csv", stringsAsFactors = F, header = T)
traits <- traits[!duplicated(traits$plant),]
traits <- traits[traits$treatment_failed == "no",]

#merge ml with treatment
dim(dat)
merge_dat <- merge(dat, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)
merge_dat2 <- merge(dat2, traits[,c(1,33)], by.x = "dat.sample", by.y = "plant", all.x = T)

dim(merge_dat)

dat3 <- read.table("./data/otutable16S_gg_cpDNAmtDNAISDduds_sum",stringsAsFactors = F, header = T)
dim(dat3)
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

#zeros are cases where the taxon was not present at all. 
#so we can compare number of zeros between comparable treated and untreated groups
table(prevalence[,2] == 0) #treated
table(prevalence[,3] == 0) #untreated, should have more zeros thn treated

table(prevalence[,1] == 0) #treated
table(prevalence[,5] == 0) #untreated

#what about comparing when a taxon is more prevalent in one gropu than the other
table(prevalence[,2] > prevalence[,3]) #should have more trues than falses 
table(prevalence[,1] > prevalence[,5]) #

#identify inoculum taxa
inoculum <- unique(c(inocula$V2, bacts))
inoculum <- gsub("(centroid=Otu\\d+);.*", "\\1", inoculum)
inoculum <- gsub("=","\\.", inoculum)
names(merge_dat) <- gsub("(centroid=Otu\\d+);.*", "\\1", names(merge_dat))

inoculum_df <- merge_dat[,c(1,which(names(merge_dat) == "treament") , which(names(merge_dat) %in% inoculum))]
inoculum_df2 <- merge_dat2[,c(1,which(names(merge_dat2) == "treament") , which(names(merge_dat2) %in% inoculum))]

head(inoculum_df)
inoculum_df$sums <- rowSums(inoculum_df[,3:length(inoculum_df)])
inoculum_df2$sums <- rowSums(inoculum_df2[,3:length(inoculum_df2)])

colSums(inoculum_df[,3:length(inoculum_df)])
#make a whole panel of boxplots
#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf(width = 8, height = 6, file = "./visuals/inoculum16S_by_ISD.pdf")
 par(mfrow = c(1,2),oma = c(4,1,1,1))
# for(i in 3:9){
#   boxplot(inoculum_df[,i] ~ inoculum_df$treament, outline = F, las = 2, ylab = "", xlab = "",
#           names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
#                     "Treated -","Treated +", "Untreated -", "Untreated +"), xpd = NA,
#           col = c(rep(add.alpha("gray", alpha = 0.8),2), rep("white", 2)))
#   title(names(inoculum_df[i]),  line = -2, adj = 0.9)
# }
boxplot(log10(inoculum_df$sums) ~ inoculum_df$treament,  outline = F, las = 2, xlab = "",
        names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                  "Treated -","Treated +", "Untreated -", "Untreated +"), ylab = "Proportion / ISD (log 10)",
        col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ),main = "Summed")
#title("Summed",  line = -1.5, adj = 0.9, cex = 0.5)

boxplot(inoculum_df2$sums ~ inoculum_df2$treament,  outline = F, las = 2, xlab = "",
        names = c("Control -", "Control +", "Treated, Control -", "Treated, Control + ",
                  "Treated -","Treated +", "Untreated -", "Untreated +"), ylab = "",
        col = c(rep("white", 2), rep(add.alpha("gray", alpha = 0.8), 4),rep("white", 2) ), main = "Summed (rel. abund.)")
#title("Summed (rel. abund.)",  line = 1, adj = 0.9, cex = 0.5)
text("Proportion", x = -2, y = 0.5, xpd = NA, srt = 90)
dev.off()

diffs <- read.csv("./data/diffs_16S.csv", stringsAsFactors = F)
inoculum <- gsub("centroid.(Otu\\d+)", "\\1", inoculum)
diffs_df <- data.frame(diffs$comparison, diffs[,names(diffs) %in% inoculum])

a <- diffs_df[as.character(diffs_df$diffs.comparison) %in% c("treatment_5_vs_treatment_7","treatment_6_vs_treatment_8"),]

diff_5v7 <- vector()
diff_6v8 <- vector()
k <- 1
for(i in 3:length(inoculum_df)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)

  diff_5v7[k] <- avgs[5,2] - avgs[7,2]
  diff_6v8[k] <- avgs[6,2] - avgs[8,2]
  k <- k + 1
}

table(diff_5v7 < 0)
table(diff_5v7 < 0)[1] / sum(table(diff_5v7 < 0))
table(diff_6v8 < 0)


########
# Plot##
#######
options(scipen = 99)
pdf(width = 10, height = 8, file = "./visuals/bacterialInoculumGraph.pdf")
par(mfrow = c(2,1), mar = c(1,4,0,1), oma = c(2,6,0,0))
plot(NULL,
     xlim = c(0, length(inoculum_df)+0.5),
     ylim = c(0,0.0003),
     las = 2,
     bty = "n",
     ylab = "",
     xlab = "",
     xaxt = "n", 
     cex.lab = 1.2)
k <- 1
for(i in 3:(length(inoculum_df2)-1)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)
  print(i)
  points(k, avgs[8,2], col = add.alpha("darkgoldenrod2",0.6), pch = 16, cex = 2)
  points(k, avgs[6,2], col = add.alpha("cyan4",0.6), pch = 16, cex = 2)
  k <- k + 1
}
text(substitute(paste( italic('A. fulva')," present", sep = "")), srt = 0, xpd = NA, x = 50, y = 0.00025, cex = 1.5)
legend(x = 1, y = 0.0003, bty = "n", y.intersp = 1, legend = c("Inoculated", "Not inoculated"), pch = 16, col = c(add.alpha("cyan4",0.6),
                                                            add.alpha("darkgoldenrod2",0.6)), cex = 1.5)
axis(side = 1, at = c(seq(0,50,by = 10),55), labels = rep("",7))

plot(NULL,
     xlim = c(0, length(inoculum_df)+0.5),
     ylim = c(0,0.0005),
     las = 2,
     ylab = "",
     xlab = "",
     xaxt = "n",
     bty = "n",
     cex.lab = 1.2)
k <- 1
for(i in 3:(length(inoculum_df2)-1)){
  avgs <- aggregate(inoculum_df2[,i] ~ inoculum_df2$treament, FUN = mean)
  print(i)
  points(k, avgs[7,2], col = add.alpha("darkgoldenrod2",0.6), pch = 16, cex = 2)
  points(k, avgs[5,2], col = add.alpha("cyan4",0.6), pch = 16, cex = 2)
  k <- k + 1
}
text("Proportion / ISD", srt = 90, xpd = NA, x = -11.5, y = 0.00055, cex = 2)
text(substitute(paste("No ", italic('A. fulva'), sep = "")), srt = 0, xpd = NA, x = 50, y = 0.00045, cex = 1.5)
axis(side = 1, at = c(seq(0,50,by = 10),55), labels = rep("",7))

dev.off()
