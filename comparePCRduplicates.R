dat <- read.table("~/Desktop/teton/otutable16S_gg_cpDNAmtDNAISDduds_sum", header = T, stringsAsFactors = F)
head(dat)


altnames <- gsub("X16S_\\w+_\\w+_Harrison_(asle_\\d+)asle", "\\1", names(dat))
altnames <- gsub("X16S_\\w+_\\w+_zymo_even_10xdilutionasle", "zymo", altnames)
altnames <- gsub("X16S_\\w+_\\w+_asle_blankasle", "blank", altnames)
altnames <- gsub("X16S_\\w+_\\w+_asle_inoculum1_6", "inoc1", altnames)
altnames <- gsub("X16S_\\w+_\\w+_asle_inoculum2asle", "inoc2", altnames)

altnames <- gsub("ITS_\\w+_\\w+_Harrison_(asle_\\d+)asle", "\\1", names(dat))
altnames <- gsub("ITS_\\w+_\\w+_zymo_even_10xdilutionasle", "zymo", altnames)
altnames <- gsub("ITS_\\w+_\\w+_asle_blankasle", "blank", altnames)
altnames <- gsub("ITS_\\w+_\\w+_asle_inoculum1_6", "inoc1", altnames)
altnames <- gsub("ITS_\\w+_\\w+_asle_inoculum2asle", "inoc2", altnames)

k <- 1
out <- list()
for(i in unique(altnames[2:length(altnames)])){
 indices <- grep(paste("^",i,"$", sep = ""), altnames)
 tocor <- dat[,indices]
 out[[k]] <-  cor.test(tocor[,1], tocor[,2])
 k <- k + 1
}

out[[181]]$estimate

hist(unlist(lapply( out , "[[" , "estimate" )))
mean(unlist(lapply( out , "[[" , "estimate" )))

#looking at the odd one that isn't matching well
which(unlist(lapply( out , "[[" , "estimate" )) ==  min(unlist(lapply( out , "[[" , "estimate" ))))
#this one is odd      out[[7]]
unique(altnames[2:length(altnames)])[7]
#asle_183

i = "asle_183"
indices <- grep(paste("^",i,"$", sep = ""), altnames)
tocor <- dat[,indices]
cor.test(tocor[,1], tocor[,2])
colSums(tocor)

#seems like there was a failure of this sample.

#what other problems?
which(unlist(lapply( out , "[[" , "estimate" )) < 0.8)

i = unique(altnames[2:length(altnames)])[1]
indices <- grep(paste("^",i,"$", sep = ""), altnames)
tocor <- dat[,indices]
cor.test(tocor[,1], tocor[,2])
colSums(tocor)
#except for 105 and 106 one of the mismatches has very low read counts
#see if they get better when converting to proportions
tocor[,1] <- tocor[,1] / sum(tocor[,1])
tocor[,2] <- tocor[,2] / sum(tocor[,2])
#doesn't help. Probably should have known.
      