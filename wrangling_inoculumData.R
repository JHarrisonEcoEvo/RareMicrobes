#For fungi
rm(list=ls())
tax <- read.table("./data/tax_inoculumNS_MS_dataITS.sintax", 
                  stringsAsFactors = F, fill = T)

otus <- read.table("./data/otutable_inoculum_ITS", stringsAsFactors = F)

otus$tax <- NA
for(i in 1:length(otus$V1)){
  otus$tax[i] <- tax[grep(paste(otus$V1[i],";", sep = ""), tax$V1),4]
}

#bring in actual sequences
seqs <- read.table("./data/all_ITS_OTUs_inoculum_NS_and_MiseqNAMES_OF_OTUS_NOTSAME_AS_NS.fasta", 
                   stringsAsFactors = F)

tmp <- scan("./data/all_ITS_OTUs_inoculum_NS_and_MiseqNAMES_OF_OTUS_NOTSAME_AS_NS.fasta",sep="\n",what="character")

indices <- grep(">", tmp)

#combine sequences
otus$seq <- NA

for(i in 2:(length(indices)-1)){
otus$seq[grep(paste(gsub(">(centroid=Otu.*);.*;.*", "\\1",
                         tmp[indices[i]]),"$",sep=""), 
          otus$V1)] <- paste(tmp[(indices[i]+1):(indices[i+1]-1)], sep = "", collapse = "")
}

#add the missing one
otus[11,]
otus$seq[11]<- paste(tmp[93:96], sep = "", collapse = "")

otus[grep("taurica", otus$tax),]
otus$V2

#add the Novaseq OTU
otus[1+length(otus$seq),] <- NA
otus$seq[length(otus$seq)] <- "AAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAGTGAATTAAACATGCTTGGTGTCTTCGCTTCGGCAAGGCCCTTGCTTAACTCACATCCCAACACCTGTGAACTGTAAGGCGCATGACTAGGTTCGCCCAAGTCATCGTCTGCCCTTTTTAACAAACAATTAATGTAACAAACGTAGTCTTATTATAACCTAATAAAACTTTCAACAACGGATCTCTTGGCTCTC"

#from taxonomy file from NovaSeq
otus$tax[length(otus$seq)] <- "d:Fungi,p:Ascomycota,c:Pezizomycetes,o:Pezizales,f:Pyronemataceae"
otus$V2[length(otus$seq)] <- 1

names(otus)[1:2] <- c("otu", "read_count")
otus <- otus[,-1]
write.csv(otus, row.names = F, file = "./data/inoculum_otu_table_ITS.csv")

##############
#For bacteria
#############

rm(list=ls())
tax <- read.table("./data/tax_inoculumNS_MS_data16S.sintax", 
                  stringsAsFactors = F, fill = T)

otus <- read.table("./data/otutable_inoculum16S", stringsAsFactors = F)

otus$tax <- NA
for(i in 1:length(otus$V1)){
  if(any(grepl(paste(otus$V1[i],";", sep = ""), tax$V1)) == T){
    otus$tax[i] <- tax[grep(paste(otus$V1[i],";", sep = ""), tax$V1),4]
  }else{
    otus$tax[i] <- "NA"
  }
}

#bring in actual sequences
seqs <- read.table("./data/all_16S_OTUs_inoculum_NS_and_Miseq_NAMES_OF_OTUS_NOTSAME_AS_NS.fasta", 
                   stringsAsFactors = F)

tmp <- scan("./data/all_16S_OTUs_inoculum_NS_and_Miseq_NAMES_OF_OTUS_NOTSAME_AS_NS.fasta",sep="\n",what="character")

indices <- grep(">", tmp)

#combine sequences
otus$seq <- NA

for(i in 2:(length(indices)-1)){
  otus$seq[grep(paste(gsub(">(centroid=Otu.*);.*;.*", "\\1",
                           tmp[indices[i]]),"$",sep=""), 
                otus$V1)] <- paste(tmp[(indices[i]+1):(indices[i+1]-1)], sep = "", collapse = "")
}

#add the missing one
otus[174,]
tmp[grep("centroid=Otu255", tmp):1303]

#add the Novaseq OTU
otus[174,4] <- "CTTGGTCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGATGCCTTACATGCAGACCAACTCGTGAATTTGTTTGAATACATAGGGATGGCACGGGTGTTTTTGCACCACGACCTCCCTTTGGGTAGGAGGGGGTGCGTCCCCCTCATGCCTGAACACAAACCCCGGCGTTCAATGCGCCAAGGAACATAAAAATCGATCAATGCGCCCCGTCGGCCCGGAGACGGTGCTCTGGCGGTGGTGCCTTGTCACATGATACAGAATGACTCTCGGCAACGGATATCTAGGCTCTTGCATCGATGAAGAACGCAGC"

tax[grep("centroid=Otu255", tax$V1),]
#from taxonomy file from NovaSeq
otus[174,3] <-  ""

names(otus)[1:2] <- c("otu", "read_count")
otus <- otus[,-1]
write.csv(otus, row.names = F, file = "./data/inoculum_otu_table_16S.csv")

