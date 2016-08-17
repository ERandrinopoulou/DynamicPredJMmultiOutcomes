####
# data = long formal data for the continuous longitudinal analysis
# data2 = long formal data for the ordinal longitudinal analysis
# data.id = data for the survival analysis
####
library(rjags)
library(R2WinBUGS)
library(JMbayes)
library(rms)
library(arm)

####
load("data.id1.Rdata")
load("data1.Rdata")
source("Functions.R")
source("PrepareData.R")

################################################
# Save model as txt
source("ModelJAGS.R")
filename <- file.path("JM_multiOutcomes_CR.txt")
write.model(model, filename)

################################################
# run model in jags
model.fit <- jags.model(file = "JM_multiOutcomes_CR.txt", data = Data, n.chains = con$n.chains, 
                         n.adapt = con$n.adapt, quiet = con$quiet)


update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
codaFit <- as.mcmc.list(res)

################################################
# results
bss <- do.call(rbind,codaFit)
colnames(bss)
n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, ii]
}

save(codaFit, file = "results.RData")
######### sims.list is a list with the results. Then you will need to obtain the mean/median/mode of each chain


