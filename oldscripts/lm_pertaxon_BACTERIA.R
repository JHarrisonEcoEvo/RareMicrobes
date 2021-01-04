library(rjags)
#library(vegan)


#####################
#Input and wrangling#
#####################

dat <- read.csv("~/rproj/asle_ddig/data/alldata.csv")

dat_mod <- data.frame(dat$Plant.Number,
                      dat$treament,
                      dat$vol/(100^3), 
                      dat$sla_avg, 
                      dat$percC, 
                      dat$percN, 
                      dat$ndfa, 
                      dat$X..Swainsonine,
                      dat[,99:242])

dat_mod_noNA <- na.omit(dat_mod)
dim(dat_mod_noNA)

#dat_mod_noNA <- merge(dat_mod_noNA, mlests, by.x = "dat.Plant.Number", by.y = "sampsAnalyzed")
dim(dat_mod_noNA)

#order
dat_mod_noNA <- dat_mod_noNA[order(dat_mod_noNA$dat.treament),]
dim(dat_mod_noNA)

dat_mod_noNA$dat.X..Swainsonine <- as.numeric(as.character(dat_mod_noNA$dat.X..Swainsonine))
dat_mod_noNA$dat.X..Swainsonine[which(is.na(dat_mod_noNA$dat.X..Swainsonine))] <- 0

for(i in 3:8){
  dat_mod_noNA[,i] <- as.vector(scale(dat_mod_noNA[,i], center = T, scale = T))
}

#need to redfine start and end points and count samples per treatment group
table(dat_mod_noNA$dat.treament) #not bad, only lost 16 data points by removing NAs

#define start and end points for model
starts <-  c(1,7, 13, 19, 23, 62, 96, 129)
ends <- c(6, 12, 18, 22, 61, 95, 128, 157)

#define model

community.model.level <- "model{
for(i in 1:N){                           
  for(j in start[i]:end[i]){
    datamatrix[j,] ~ dmulti(p[j,], nreads[j])
    p[j,1:notus] ~ ddirch(pi[j,]*theta[i]) 

    for(k in 1:notus){
      pi[j,k] <- 1/(1+exp(-(intercept[i,k] + beta1[i,k]*covariate1[j]
                  +	beta2[i,k]*covariate2[j]
                  +	beta3[i,k]*covariate3[j]
                  +	beta4[i,k]*covariate4[j]
                  +	beta5[i,k]*covariate5[j]
                  +	beta6[i,k]*covariate6[j]
                )))
    }
  }
  for(l in 1:notus){
    beta1[i,l] ~ dnorm(mu_b1[l], tau_b1[l])
    beta2[i,l] ~ dnorm(mu_b2[l], tau_b2[l])
    beta3[i,l] ~ dnorm(mu_b3[l], tau_b3[l])
    beta4[i,l] ~ dnorm(mu_b4[l], tau_b4[l])
    beta5[i,l] ~ dnorm(mu_b5[l], tau_b5[l])
    beta6[i,l] ~ dnorm(mu_b6[l], tau_b6[l])
    intercept[i,l] ~ dnorm(mu_i[l], tau_i[l])
  }

  #site specific conc. parameter of Dirichlet 
  theta[i] ~ dunif(hypertheta,2000)
}

#hyper prior for precision
hypertheta ~ dunif(1.0E-3, 2000)

#hyper priors for model terms
for(z in 1:notus){
  mu_b1[z] ~ dnorm(0,50)
  mu_b2[z] ~ dnorm(0,50)
  mu_b3[z] ~ dnorm(0,50)
  mu_b4[z] ~ dnorm(0,50)
  mu_b5[z] ~ dnorm(0,50)
  mu_b6[z] ~ dnorm(0,50)
  mu_i[z] ~ dnorm(0,50)

  tau_b1[z] ~ dunif(0, 200)
  tau_b2[z] ~ dunif(0, 200)
  tau_b3[z] ~ dunif(0, 200)
  tau_b4[z] ~ dunif(0, 200)
  tau_b5[z] ~ dunif(0, 200)
  tau_b6[z] ~ dunif(0, 200)
  tau_i[z] ~ dunif(0, 200)
}
}"

  #compile model
  sim.mod.jags <- jags.model(
    textConnection(community.model.level),
    data = list(
      datamatrix = dat_mod_noNA[,9:152],
      notus = 144, 
      nreads = as.vector(rowSums(dat_mod_noNA[,9:152])),
      covariate1 = dat_mod_noNA$dat.vol..100.3.,
      covariate2 = dat_mod_noNA$dat.sla_avg,
      covariate3 = dat_mod_noNA$dat.percN,
      covariate4 = dat_mod_noNA$dat.percC,
      covariate5 = dat_mod_noNA$dat.ndfa,
      covariate6 = dat_mod_noNA$dat.X..Swainsonine,
      start = starts,
      end = ends,
      N = 8
    ),
    n.chains = 2,
    n.adapt = 0
  )
  
  #adapt
  iter_needed <- 0
  y = FALSE
  while (y == FALSE) {
    y <-  adapt(sim.mod.jags,
                n.iter = 1000,
                end.adaptation = FALSE)
    iter_needed <- 1000 + iter_needed
    if (iter_needed > 5000) {
      break
    }
  }
  
  #burn
  burn <- 5000
  update(sim.mod.jags,
         n.iter = burn)
  
  #sample
  sim.mod.sam <- jags.samples(
    model = sim.mod.jags,
    variable.names = c(
      "mu",
      "beta1",
      "beta2",
      "beta3",
      "beta4",
      "beta5",
      "beta6",
      "mu_b1",
      "mu_b2",
      "mu_b3",
      "mu_b4",
      "mu_b5",
      "mu_b6"
    ),
    n.iter = 4000,
    thin = 4
  )
  
  #Diagnostic statistics.
  
  #Calculate the Gelman-Rubin statistic, spot check
  gr <- gelman.diag(as.mcmc.list(sim.mod.sam$beta4),
                    multivariate = F)
  gr
  
  gk <- geweke.diag(as.mcmc.list(sim.mod.sam$beta4),
                    frac1 = 0.1,
                    frac2 = 0.5)
  gk


save.image("./Bayesian_lm_traits_BACTERIA.RData")
