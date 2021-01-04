leaves = read.csv("./data/SizeEndofSeason.csv")
plantorder = read.csv("./data/planting_order.csv")


dat = merge(leaves, plantorder, by.x = "Plant.Number", by.y="rns")
dat$volumeEnd = dat[,2]*dat[,3]*dat[,4]
par(mar=c(10,2,2,2))
boxplot(dat$volumeEnd~dat$treament, las=2)

reg = aov(dat$volumeEnd~dat$treament)
TukeyHSD(reg)

table(dat$Flowering, dat$treament)

#need 14 more control neg to be treated, and 11 more control pos treated.

dat2 = dat[which(dat$stem_amount > 3),]
table(dat2$treats)


round(runif(11,1,length(dat[which(dat$treats=="controlNeg"),1])))

