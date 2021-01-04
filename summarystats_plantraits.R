dat = read.csv("~/Downloads/planting_order.csv", header=T)
dat$vol = dat[,6]*dat[,7]*dat[,8]
dat$perc_infected = dat$lvs_with_fungus/dat$X.LeavesCollection 
dat$perc_eaten = dat$lvs_with_herbivory/dat$X.LeavesCollection 

aggregate(dat$vol~dat$treament, FUN=mean)
plot(dat$treament, dat$vol)
reg = aov(dat$vol~dat$treament)
summary(reg)
TukeyHSD(reg)

table(dat$treament, dat$X.Fl_Stems_atCollectionTime)

aggregate(dat$X.LeavesCollection~dat$treament, FUN=mean)
plot(dat$treament, dat$X.LeavesCollection)
reg = aov(dat$X.LeavesCollection~dat$treament)
summary(reg)
TukeyHSD(reg)


aggregate(dat$perc_eaten ~dat$treament, FUN=mean)
plot(dat$treament, dat$perc_eaten)
reg = aov(dat$perc_eaten ~dat$treament)
summary(reg)
TukeyHSD(reg)

aggregate(dat$perc_infected ~dat$treament, FUN=mean)
plot(dat$treament, dat$lvs_with_fungus)
reg = aov(dat$perc_infected ~dat$treament)
summary(reg)
TukeyHSD(reg)

culture = read.csv("~/Downloads/Culturingxp.csv")
culture2 = merge(dat, culture, by.x="rns", by.y="plant")
#NOTE if black spore infections are included this pattern is reversed/obscured. 
#this is because many of the black spore infections are contaminants

culture2$percInfNonAlfu = culture2$notAlfuInfected/culture2$num_leaf_segments
plot(culture2$treament,culture2$percInfNonAlfu)
plot(culture2$treament,culture2$notAlfuInfected)
plot(culture2$treament,culture2$morphospecies_not_counting_Undifilum)

reg = aov(culture2$morphospecies_not_counting_Undifilum~culture2$treament)
summary(reg)
TukeyHSD(reg)

reg = aov(culture2$percInfNonAlfu~culture2$treament)
summary(reg)
TukeyHSD(reg)


