crLongNEW <- function(data, statusVar, censLevel, nameStrata = "strata", nameStatus = "status2") {
  n <- nrow(data)
  status <- data[[statusVar]]
  unqLevs <- unique(status)
  unqLevs <- unqLevs[unqLevs != "alive"]
  ncr <- length(unqLevs)
  dataOut <- data[rep(seq_len(n), each = ncr), ]
  dataOut[[nameStrata]] <- rep(unqLevs, n)
  dataOut[[nameStatus]] <- as.numeric(dataOut[[statusVar]] ==
                                        dataOut[[nameStrata]])
  dataOut[[nameStrata]] <- factor(dataOut[[nameStrata]])
  dataOut
}

pii <- function(X, betas, Z, b) {
  exp(X%*%betas + Z%*%b)/(1 + exp(X%*%betas + Z%*%b)) #1/(1 + exp(-X%*%betas - Z%*%b))
}

like <- function(p, y){
  dbinom(y, 1, p, log = TRUE)
}


crback.setup <- function(y)
{
  yname <- as.character(substitute(y))
    y <- factor(y, exclude = NA)
  y <- unclass(y)
  ylevels <- levels(y)
  kint <- length(ylevels) - 1
  y <- as.integer(y - 1)
  reps <- ifelse(is.na(y), 1, ifelse(y > 1, kint - y + 1, kint))
  subs <- rep(1:length(y), reps)
  cuts <- vector("list", kint + 2)
  cuts[[1]] <- NA
  for(u in 0:kint)
    cuts[[u + 2]] <- if(u < 1) 1:kint else u:kint
  cuts <- unlist(cuts[ifelse(is.na(y), 1, y + 2)])
  y <- rep(y, reps)
  Y <- 1 * (y == cuts)
  labels <- c(paste(yname, "<=", ylevels[2:kint], sep = ""), "all")
  cohort <- factor(cuts, levels = 1:kint, labels = labels)
  list(y = Y, cohort = cohort, subs = subs, reps = reps)
}



likelihoodSurvCR <- function(gammasR, gammasD, alphasR, alphasD, W2sR, W2sD, Xs, Xs2, Zs, Zs2,xiR, xiD, b) {
  
  log.h0.sR <- matrix(0,1,K)
  log.h0.sD <- matrix(0,1,K)
  
  f.s <- matrix(0,1,K)
  f.s2 <- matrix(0,1,K)
  
  SurvLongR <- matrix(0,1,K)
  SurvLongD <- matrix(0,1,K)
  
  for (o in 1 : K) {
    log.h0.sR[o] <- xiR %*% W2sR[o,]
    log.h0.sD[o] <- xiD %*% W2sD[o,]
    
    f.s[o] <- betas11%*%Xs[o,] + Zs[o,] %*% (b[iii,1:ncZ])
    f.s2 <- Xs2[offset22[o]:(offset22[o+1]-1),]%*%betas22 + Zs2[offset22[o]:(offset22[o+1]-1),]%*%(b[iii,(ncZ + 1):(ncZ + ncZ2)])
    f.s2[o] <- sum(f.s2)
    
    SurvLongR[o] <- wk[o] * exp(log.h0.sR[o] + alphasR * f.s[o] + alphas2R * f.s2[o])
    SurvLongD[o] <- wk[o] * exp(log.h0.sD[o] + alphasD * f.s[o] + alphas2D * f.s2[o])
    
  }

  etaBaselineR <- xWR%*%gammasR
  etaBaselineD <- xWD%*%gammasD
  
  log.survivalR <- -exp(etaBaselineR) * P[iii] * sum(SurvLongR)
  log.survivalD <- -exp(etaBaselineD) * P[iii] * sum(SurvLongD)
  
  (log.survivalR + log.survivalD) 
  
}


log.basHazR <- function(r) {
  W2R <- splineDesign(rr, r, ord = con$ordSpline, outer.ok = TRUE)
  W2R%*%Bs.gammasR
}

log.basHazD <- function(r) {
  W2D <- splineDesign(rr, r, ord = con$ordSpline, outer.ok = TRUE)
  W2D%*%Bs.gammasD
}


ff <- function(d) {
  ns(d, knots = c(2.121834, 5.516769), Boundary.knots = c(0.01642711, 19.54004107)) 
}



m1 <- function(d) {
  XX <-  cbind(1, ff(d))
  ZZ <-   cbind(1,ff(d))
  c(XX%*%betas11 + ZZ%*%(b[iii,1:ncZ]))
}


m2 <- function(d) {
  if (length(d)>1)  {coh.mat <- matrix(rep(coh.mat, length(d)), ,5, byrow = T)} 
  bb <- cbind(coh.mat, d)#, Age2.i, Sex2.i)
  XX2 <- bb 
  ZZ2 <- cbind(1,d)
  XX2%*%betas22 + ZZ2%*%(b[iii,(ncZ+1):(ncZ+ncZ2)])
}



log.hR <- function(y) {
  c(matrix(rep(eta.tR, length(y)),,1) + matrix(alphasR%*%t(m1(y)),,1) + matrix(alphas2R%*% t(m2(y)),,1) + log.basHazR(y))
}

log.hD <- function(x) {
  c(matrix(rep(eta.tD, length(x)),,1) + matrix(alphasD%*%t(m1(x)),,1) + matrix(alphas2D%*% t(m2(x)),,1) + log.basHazR(x))
}


hazR <- function(s)  exp(log.hR(s))
hazD <- function(s)  exp(log.hD(s))



S <- function(t) {
  f <- function (s) hazD(s) + hazR(s)
  sapply(t, function (ix) exp(-(integrate(f, lower = 0, upper = ix)$value)) )   #exp(-(integrate(f, lower = 0, upper = ix)$value))
}


CIFr <- function(t) {
  g <- function(s) hazR(s)*S(s)
  sapply(t, function(ix) 
    integrate(g, lower = tt[1], upper = ix, rel.tol = 1e-06,
              subdivisions = 500, stop.on.error = FALSE)$value)
}

CIFd <- function(t) {
  g <- function(s) hazD(s)*S(s)
  sapply(t, function(ix) 
    integrate(g, lower = tt[1], upper = ix, rel.tol = 1e-06,
              subdivisions = 500, stop.on.error = FALSE)$value)
}



