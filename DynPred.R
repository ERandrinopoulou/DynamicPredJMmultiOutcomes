###########################################################
# load results from MCMC - change the path
e1 <- new.env() 
load("results.RData", e1)
ls(e1) 
codaFit <- get('codaFit', e1) 
ls() 

parms <- c("betas", "betas2", "tau", "inv.D", "gammasR", "gammasD","b", "alphasR", "alphas2R", "alphasD", "alphas2D",
           "Bs.gammasR", "Bs.gammasD")

bss <- do.call(rbind,codaFit)
colnames(bss)
n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  iik <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, iik]
}

# name the parameters
betas1 <- sims.list[[1]][,1:ncX]
betas2 <- sims.list[[2]]
taus <- sims.list[[3]]
taus_b <- sims.list[[4]]
gR <- cbind(sims.list[[5]], sims.list[[8]], sims.list[[9]])
gD <- cbind(sims.list[[6]], sims.list[[10]], sims.list[[11]])
ksiR <- sims.list[[12]]
ksiD <- sims.list[[13]]



###########################################################
set.seed(2016)

# select a sample of the MCMC results - n=300
vec <- sample(length(taus), 300)     

### select a patient to predict
ii <- 46
# index of the patient for prediction
data.id$numb <- 1:dim(data.id)[1]
iii <- data.id$numb[data.id$id == ii]

# specify his time visits
ti <- data$Time[data$id == ii]

# select his survival design matrices
xWR <- x$WR[ii, ] 
xWD <- x$WD[ii, ] 

# select his offset for the ordinal dataset
offset22 <- as.vector(c(1, 1 + cumsum(table(data.id22$kk[data.id22$id == ii]))))


# create new parameters to save the results -  we are calculating the predictions at 12 time points
vvD <- matrix(0, length(vec), 12)
vvR <- matrix(0, length(vec), 12)

A <- list()
B <- list()

vector2 <- list()

# set number of observations
N <- length(data$y1)

# create overall random effects matrix
b1 <- data.matrix(ranef(lmeObject))
b2 <- (ranef(glmObject))
b2 <- as.matrix(b2$id)

b <- cbind(b1,b2)
b <- data.frame(b)

# set random effects matrix columns
raneff <- dim(b)[2]
nb <- ncZ + ncZ2

# offset for the design matrices of the ordinal outcome
qq <- table(X2[, 6]) # time variable
q <- cumsum(qq)

###########################################################
# set the cohort of the ordinal outcome for the predictive patient
coh <- rep(0, 4)
cohort.i <- data2.id$cohort[data2.id$id == ii][1]
cohort.i <- unclass(cohort.i) 
if (cohort.i > 1) {
  coh[cohort.i-1] <- 1
}
coh.mat <- c(1, coh)

###########################################################
# obtain predictions

for (k in 1:length(ti)) {
  
  vvD <- matrix(0, length(vec), 12)
  vvR <- matrix(0, length(vec), 12)
  DDD <- NULL
  RRR <- NULL
  
  eta.tR <- NULL
  eta.tD <- NULL
  
  vec2 <- NULL
  
  for (j in 1:length(vec)) {
    
    # set parameters
    gammasR <- gR[vec[j], 1:2]                                                                                                                                                                                                     
    gammasD <- gD[vec[j], 1:2]
    alphasR <- gR[vec[j], 3] 
    alphas2R <- gR[vec[j], 4] 
    alphasD <- gD[vec[j], 3]
    alphas2D <- gD[vec[j], 4] 
    
    Bs.gammasD <- xiD <- ksiD[vec[j],]
    Bs.gammasR <- xiR <- ksiR[vec[j],]
    betas11 <- betas1[vec[j],1:ncX]
    betas22 <- betas2[vec[j],1:ncX2]
    
    eta.tR <- (xWR%*%(gammasR)) 
    eta.tD <- (xWD%*%(gammasD)) 
    
    
    tauss <- taus[vec[j]]
    taus_bb <- taus_b[vec[j],]
    
    
    
    # random effects
    tau_b <- matrix(0,raneff,raneff)
    tau_b[1,] <- tau_b[,1] <- taus_bb[1:raneff]
    tau_b[2,] <- tau_b[,2] <- taus_bb[(raneff+1):(raneff+raneff)]
    tau_b[3,] <- tau_b[,3] <- taus_bb[(raneff*2+1):(raneff*2+raneff)]
    tau_b[4,] <- tau_b[,4] <- taus_bb[(raneff*3+1):(raneff*3+raneff)]
    tau_b[5,] <- tau_b[,5] <- taus_bb[(raneff*4+1):(raneff*4+raneff)]
    tau_b[6,] <- tau_b[,6] <- taus_bb[(raneff*5+1):(raneff*5+raneff)]
    
    tau_b[raneff,raneff] <- taus_bb[dim(taus_b)[2]]
    
    
    # draw the random effects from the posterior
    muCond1 <- rep(0, raneff)
    varCond1 <- solve(tau_b)
    
    bAll<-  b[iii,]
    raneff <- dim(b)[2]
    
    # proposal for new random effects
    b.propA <- c(rmvnorm(1, as.numeric(bAll), diag(0.1, raneff)) ) 
    
    bprop <-b
    b <- as.matrix(b)
    bprop[iii,1] <- b.propA[1]
    bprop[iii,2] <- b.propA[2]
    bprop[iii,3] <- b.propA[3]
    bprop[iii,4] <- b.propA[4]
    bprop[iii,5] <- b.propA[5]
    bprop[iii,6] <- b.propA[6]
    
    bprop <- as.matrix(bprop)
    b.prop <-  matrix(t(bprop), , 1, byrow = T)  
    
    
    p11 <-( -(tauss/2)*(t(y$y[1:k]-X[ , ][1:k,] %*% betas11-Z[ , ][1:k, ] %*% bprop[iii, (1:ncZ)])) %*% (y$y[1:k]-X[ , ][1:k, ] %*% betas11-Z[ , ][1:k, ] %*% bprop[iii, (1:ncZ)]) + #1st longitudinal outcome (continuous)
              sum(like(pii(X2[1:q[k], ], betas22, as.matrix(Z2[1:q[k], ]), bprop[iii, (ncZ+1):(ncZ+ncZ2)]), y2$y[1:q[k]])) + #2nd longitudinal outcome (ordinal)
              sum(likelihoodSurvCR(gammasR, gammasD, alphasR, alphasD, W2sR, W2sD, Xs, Xs2, Zs, Zs2, xiR, xiD, bprop))- #survival (competing risks)
              (1/2*(t(bprop[iii, ]-muCond1) %*% solve(varCond1) %*% (bprop[iii, ]-muCond1)))  ) #random effects
    p12 <- (-(tauss/2)*(t(y$y[1:k]-X[ , ][1:k, ] %*% betas11-Z[ , ][1:k, ] %*% b[iii,(1:ncZ)])) %*% (y$y[1:k]-X[ , ][1:k, ] %*% betas11-Z[ , ][1:k, ] %*% b[iii, (1:ncZ)]) +
              sum(like(pii(X2[1:q[k], ], betas22, as.matrix(Z2[1:q[k], ]), b[iii, (ncZ+1):(ncZ+ncZ2)]), y2$y[1:q[k]])) +
              sum(likelihoodSurvCR(gammasR, gammasD, alphasR, alphasD, W2sR, W2sD, Xs, Xs2, Zs, Zs2, xiR, xiD, b))-
              (1/2*(t(b[iii, ]-muCond1) %*% solve(varCond1) %*% (b[iii, ]-muCond1)))   )
    a <- (p11-p12)
    u <- log(runif(1))
    
    # accept/reject algorithm
    if (u <= a) {
      b[iii, ] <-  bprop[iii, ]
    } else {
      b[iii, ] <- b[iii, ]
    }
    
    #  calculate the CIF
    t1 <- ti[k]
    t2 <- max(Time)
    
    tt <- seq(t1,t2, length = 12)
    
    surva <- S(tt[1])
    
    
    if (surva == 0)     {
      surva <- 9.881313e-324
      vec2[j] <- vec[j]
    }
    
    # calculate the CIF/S
    try(vvD[j, ] <- c(CIFd(tt)/surva))
    
    try(vvR[j, ] <- c(CIFr(tt)/surva))
    
    # save it 
    DDD <- vvD[j, ]                                   
    RRR <- vvR[j, ]
    
    # correct for small differences betweeen predictions of next visit
    try(indD <- diff(DDD) < -1e-06)
    lD <- c(2:12) # 12 are the time points to calculate predictions
    pD <- lD[indD == T]
    if (length(pD) != 0) {
      for (uu in 1:length(pD))   { 
        DDD[pD[uu]] <- DDD[pD[uu]-1]      }
    }
    
    try(indR <- diff(RRR) < -1e-06)
    lR <- c(2:12)
    pR <- lR[indR == T]
    if (length(pR) != 0) {
      for (uu in 1:length(pR))   {  
        RRR[pR[uu]] <- RRR[pR[uu]-1] }
    }
    
    vvD[j, ] <- DDD
    vvR[j, ] <- RRR
    
  }
  
  # save which iterations where used from the MCMC
  vector2[[k]] <- vec2
  
  # save all iterations
  A[[k]] <- vvD
  B[[k]] <- vvR
print(k)
}

