####
# data = long formal data for the continuous longitudinal analysis
# data2 = long formal data for the ordinal longitudinal analysis
# data.id = data for the survival analysis
####

library(JMbayes)
library("rms")
library("arm")


load("data1.Rdata")
load("data.id1.Rdata")


colnames(data) <- c("id", "Time", "y1", "y2")

#################################
# fit the continuous longitudinal outcome using a mixed-effetcs model
lmeObject <- lme(y1 ~ ns(Time,3), data = data,
           na.action = na.exclude,
           random = list(id = pdDiag(form = ~ ns(Time, 3))))
 
#################################
# fit the ordinal longitudinal outcome using a ontinuation ratio mixed-effects model
yOR <- data$y2
u <- crback.setup(yOR)
cohort <- u$cohort
y <- u$y
id <- data$id[u$subs]
Time <- data$Time[u$subs] ### 


data2 <- data.frame(y, cohort, id, Time)

glmObject  <- glmer(y ~ cohort + Time + (Time|id), family = binomial, data= data2)

#################################
### Set
timeVar <- "Time"
param <- "td-value"
survMod <- "spline-PH"
lag <- 0

y1 <- data$y1
y2 <-data2$y

data.id$y1 <- data[tapply(row.names(data), data$id, tail, 1) ,  3]

Time <- data.id$Time ### is the survival time

#################################
# set MCMC details
con <- list(program = "JAGS", n.chains = 1, n.iter = 10000,
            n.burnin = 1000, n.thin = 10, n.adapt = 2000, K = 100,
            C = 5000, working.directory = getwd(), bugs.directory = "C:/Program Files/WinBUGS14/",
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4,
            bugs.seed = 1, quiet = FALSE)

#################################
# for the continuous longitudinal outcome create the design matrices
id <- data$id 

offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)


data.id[[timeVar]] <- pmax(data.id$Time - 0, 0) ### where Time = survival time

mfX.id <- model.frame(TermsX, data = data.id)  
mfZ.id <- model.frame(TermsZ, data = data.id) 
Xtime <- model.matrix(formYx, mfX.id)
Ztime <- model.matrix(formYz, mfZ.id)

#################################
# for the ordinal longitudinal outcome create the design matrices
id2 <- data2$id

offset2 <- as.vector(c(1, 1 + cumsum(tapply(id2, id2, length))))

formYx2 <- y ~ cohort + Time  
TermsX2 <- terms(glmObject)
mfX2 <- model.frame(TermsX2, data = data2)
X2 <- model.matrix(formYx2, mfX2)

formYz2 <- ~ Time 
mfZ2 <- model.frame(terms(formYz2), data = data2)
TermsZ2 <- attr(mfZ2, "terms")
Z2 <- model.matrix(formYz2, mfZ2)


data2.id <- NULL
ttime2 <- tapply(data2$Time,data2$id,max)
id22 <- unique(id2)
data2.id <- data2[data2$id == id22[1]&data2$Time == ttime2[1],]
for (i in 2:dim(data.id)[1]){
  dat <- data2[data2$id == id22[i]&data2$Time == ttime2[i],]
  data2.id <- rbind(data2.id,dat)   
}   

data2.id[[timeVar]] <- rep(Time, tapply(data2.id$Time, data2.id$id,length))

mfX.id2 <- model.frame(TermsX2, data = data2.id)
mfZ.id2 <- model.frame(TermsZ2, data = data2.id)
Xtime2 <- model.matrix(formYx2, mfX.id2)
Ztime2 <- model.matrix(formYz2, mfZ.id2)

###################################
# for the continuous longitudinal outcome - design matrices for the 15-point Gauss-Kronrod quadrature rule approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)

P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))

mfX <- model.frame(TermsX, data = data.id2)  
mfZ <- model.frame(TermsZ, data = data.id2)    
Xs <- model.matrix(formYx, mfX)
Zs <- model.matrix(formYz, mfZ)

#################################
# for the ordinal longitudinal outcome - design matrices for the 15-point Gauss-Kronrod quadrature rule approximation

# first select the appropriate cohort for each of the 15-point Guass-Kronrod quadrature rule
data.id22 <- NULL
kk <- 1

ID <- ID_n <- c(table(data2$Time[data2$id == unique(id2)[1]]))
numbID <- c(ID)
IDnumb <- IDnumb_n <- rep(c(1:length(unique(data2$Time[data2$id == unique(id2)[1]]))), numbID)

dat <- data2[data2$id == unique(id2)[1],]

dat$IDnumb <- IDnumb

vec <- unique(data2$Time[data2$id == unique(id2)[1]])

vec<-c(vec, max(data$Time))
vec[1] <- 0

fInd <- findInterval(c(t(st))[1:15],c(vec))
data.id22 <- dat[dat$IDnumb == fInd[1],]
data.id22$kk <- kk
kk <- kk+1
for(j in 2:length(fInd)){
  dat$kk <- kk
  kk <- kk+1
  data.id22 <- rbind(data.id22,dat[dat$IDnumb == fInd[j],])
}



for (i in 2:dim(data.id)[1]) {
  ID <-  c(table(data2$Time[data2$id == unique(id2)[i]]))
  ID_n <- c(ID_n,ID)
  numbID <- c(ID)
  
  IDnumb <- rep(c(1:length(unique(data2$Time[data2$id == unique(id2)[i]]))), numbID)
  IDnumb_n <- c(IDnumb_n,IDnumb)
  
  
  dat <- data2[data2$id == unique(id2)[i],]
  dat$IDnumb <- IDnumb
  vec <- unique(data2$Time[data2$id == unique(id2)[i]])
  
  vec<-c(vec, max(data$Time))
  vec[1] <- 0
  
  fInd <- findInterval(c(t(st))[((i-1)*15+1):(i*15)],c(vec))
  for(j in 1:length(fInd)){
    dat$kk <- kk
    kk <- kk+1
    data.id22 <- rbind(data.id22,dat[dat$IDnumb == fInd[j],])
  }
}

vv <- c(table(data.id22$kk))

data.id22$Time <- rep(c(t(st)),vv) 


mfX2 <- model.frame(TermsX2, data = data.id22)
mfZ2 <- model.frame(TermsZ2, data = data.id22)
Xs2 <- model.matrix(formYx2, mfX2)
Zs2 <- model.matrix(formYz2, mfZ2)

offset2vv <- as.vector(c(1, 1 + cumsum(vv)))
XX <- c(table(data2.id$id))
offset2XX <- as.vector(c(1, 1 + cumsum(XX)))

#################################
# survival submodel
# design matrices for the survival submodel
WR <- model.matrix(~ -1 + Age, data.id)
WD <- model.matrix(~ -1 + Age, data.id)

eventR <- as.numeric(data.id$Event == 2)
eventD <- as.numeric(data.id$Event == 1)
nT <- length(Time)
zeros <- numeric(nT)

x <- list(X = X, Z = Z, WR = if (survMod == "weibull-PH") {
  if (is.null(WR)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                            WR)
} else {
  if (is.null(WR)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(WR) == 1) cbind(WR, rep(0, nT)) else WR
  }
}, WD = if (survMod == "weibull-PH") {
  if (is.null(WD)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                            WD)
} else {
  if (is.null(WD)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(WD) == 1) cbind(WD, rep(0, nT)) else WD
  }
})

#################################
# design matrices for the baseline hazard (a B-splines baseline hazard function is asssumed)
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- Time
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr
W2R <- splineDesign(rr, Time, ord = con$ordSpline)
W2D <- splineDesign(rr, Time, ord = con$ordSpline)
if (any(colSums(W2R) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")
if (any(colSums(W2D) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")

# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation
W2sR <- splineDesign(rr, c(t(st)), ord = con$ordSpline)
W2sD <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(W2R = W2R, W2D = W2D, W2sR = W2sR, W2sD = W2sD))

#################################
ncX <- ncol(X)
ncZ <- ncol(Z)
ncWR <- ncol(x$WR)
ncWD <- ncol(x$WD)
ncW2R <- ncol(x$W2R)
ncW2D <- ncol(x$W2D)
ncX2 <- ncol(X2)
ncZ2 <- ncol(Z2)
nb <- ncZ + ncZ2


b1 <- data.matrix(ranef(lmeObject))
b2 <- (ranef(glmObject))
b2 <- as.matrix(b2$id)
b <- cbind(b1,b2)
b <- data.frame(b)

nY <- nrow(b)
sigma2 <- lmeObject$sigma^2
nY2 <- nrow(b2)


y.long <- model.response(mfX, "numeric")
y <- list(y = y.long, offset = offset, logT = log(Time),
          eventR = eventR, eventD = eventD, zeros = zeros, lag = lag)
y.long2 <- model.response(mfX2, "numeric")
y2 <- list(y = y.long2, offset = offset2, logT = log(Time),
           eventR = eventR, eventD = eventD, zeros = zeros, lag = lag)

#################################
# priors/hyperpriors
betas <- rep(0, ncX)
var.betas <- rep(con$K, ncX)
betas2 <- rep(0, ncX2)
var.betas2 <- rep(con$K, ncX2)

mu0 <- rep(0,ncZ + ncZ2)

alphasD <- DalphasD <- 0
var.alphasD <- var.DalphasD <- con$K
alphas2D <- Dalphas2D <- 0
var.alphas2D <- var.Dalphas2D <- con$K

alphasR <- DalphasR <- 0
var.alphasR <- var.DalphasR <- con$K
alphas2R <- Dalphas2R <- 0
var.alphas2R <- var.Dalphas2R <- con$K


gammasD <- rep(0,(ncWD))
var.gammasD <- rep(con$K, (ncWD))
gammasR <- rep(0,(ncWR))
var.gammasR <- rep(con$K, (ncWR))

Bs.gammasD <- rep(0, (ncW2D))
var.Bs.gammasD <- rep(con$K/10, (ncW2D))
Bs.gammasR <- rep(0, (ncW2R))
var.Bs.gammasR <- rep(con$K/10, (ncW2R))


priorR.D <- diag(1,ncZ + ncZ2)
priorK.D <- (ncZ + ncZ2)

#################################
# specify parameters of interest
parms <- c("betas", "betas2", "tau", "inv.D", "gammasR", "gammasD","b", "alphasR", "alphas2R", "alphasD", "alphas2D",
           "Bs.gammasR", "Bs.gammasD")



#################################
Data <- list(N = nY, K = K, offset = offset, offset2 = offset2, X = X, X2 = X2, Xtime = Xtime, Xtime2 = Xtime2,
             y = y$y, y2 = y2$y, Xs = Xs, Xs2 = Xs2, Z = Z, Z2 = Z2, Ztime = Ztime, Ztime2 = Ztime2,
             Zs = Zs, Zs2 = Zs2, eventR = eventR, eventD = eventD, zeros = zeros, WR = x$WR, WD = x$WD, ncZ = ncol(Z), 
             ncZ2 = ncol(Z2),
             ncX = ncol(X), ncX2 = ncol(X2), ncWR = ncol(x$WR), ncWD = ncol(x$WD), W2R = W2R, W2D = W2D, W2sR = W2sR, 
             W2sD = W2sD, ncW2R = ncol(x$W2R), ncW2D = ncol(x$W2D), C = C, P = P,
             wk = wk, nb = nb, #nb2 = nb2,
             mu0 = mu0, offset2vv = offset2vv, offset2XX = offset2XX,  #mu02 = mu02,
             priorMean.betas = betas,  priorMean.betas2 = betas2,
             priorTau.betas = diag(1/var.betas), priorTau.betas2 = diag(1/var.betas2),
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, priorMean.gammas = gammasR,
             priorTau.gammas = diag(1/var.gammasR),
             priorMean.alphas = alphasR, priorMean.alphas2 = alphas2R,
             priorTau.alphas = 1/var.alphasR, 
             priorTau.alphas2 = 1/var.alphas2R,
             priorMean.Bs.gammas = Bs.gammasR,
             priorTau.Bs.gammas = diag(1/var.Bs.gammasR),
             priorR.D = priorR.D , priorK.D = priorK.D)

