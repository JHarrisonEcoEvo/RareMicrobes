ord <- read.csv("~/Downloads/planting_order.csv")
h <- read.csv("~/Downloads/plantsHeather_used.csv")

dat = ord[-which(ord$rns %in% h$plant_id),]

dat = dat[dat$status=="good",]

table(dat$treament)

dat$harvest=NA

dat$harvest[sample(which(dat$treament == "treated_neg"), size=20)]="yes"
dat$harvest[sample(which(dat$treament == "treated_plus"), size=20)]="yes"

dat$harvest[sample(which(dat$treament == "untreated_neg"), size=20)]="yes"
dat$harvest[sample(which(dat$treament == "untreated_plus"), size=20)]="yes"


dat$harvest[grep("control", dat$treament)] = "yes"

write.csv(dat, file="toharvest.csv")