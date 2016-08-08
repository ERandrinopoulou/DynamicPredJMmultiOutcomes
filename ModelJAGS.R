model <- function ()
{
    for (i in 1:N) {
        for (j in offset[i]:(offset[i + 1] - 1)) {
            muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 
                1:ncZ], Z[j, 1:ncZ])
            y[j] ~ dnorm(muy[j], tau)
        }
        for (j in offset2[i]:(offset2[i + 1] - 1)) {
            muy2[j] <- inprod(betas2[1:ncX2], X2[j, 1:ncX2]) + 
                inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], Z2[j, 1:ncZ2])
            Pr[j] <- max(1.00000E-05, min(0.99999, (exp(muy2[j])/(1 + 
                exp(muy2[j])))))
            y2[j] ~ dbin(Pr[j], 1)
        }
        f.T[i] <- inprod(betas[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 
            1:ncZ], Ztime[i, 1:ncZ])
        for (f in offset2XX[i]:(offset2XX[i + 1] - 1)) {
            ff.T2[i, f] <- inprod(betas2[1:ncX2], Xtime2[f, 1:ncX2]) + 
                inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], Ztime2[f, 
                  1:ncZ2])
        }
        f.T2[i] <- sum(ff.T2[i, offset2XX[i]:(offset2XX[i + 1] - 
            1)])
        etaBaselineR[i] <- inprod(gammasR[1:(ncWR)], WR[i, 1:ncWR])
        etaBaselineD[i] <- inprod(gammasD[1:(ncWD)], WD[i, 1:ncWD])
        log.h0.TR[i] <- inprod(Bs.gammasR[1:(ncW2R)], W2R[i, 
            1:ncW2R])
        log.h0.TD[i] <- inprod(Bs.gammasD[1:(ncW2D)], W2D[i, 
            1:ncW2D])
        log.hazardR[i] <- log.h0.TR[i] + etaBaselineR[i] + alphasR * 
            f.T[i] + alphas2R * f.T2[i]
        log.hazardD[i] <- log.h0.TD[i] + etaBaselineD[i] + alphasD * 
            f.T[i] + alphas2D * f.T2[i]
        for (k in 1:K) {
            log.h0.sR[i, k] <- inprod(Bs.gammasR[1:(ncW2R)], 
                W2sR[K * (i - 1) + k, 1:ncW2R])
            log.h0.sD[i, k] <- inprod(Bs.gammasD[1:(ncW2D)], 
                W2sD[K * (i - 1) + k, 1:ncW2D])
            f.s[i, k] <- inprod(betas[1:ncX], Xs[K * (i - 1) + 
                k, 1:ncX]) + inprod(b[i, 1:ncZ], Zs[K * (i - 
                1) + k, 1:ncZ])
            for (f in offset2vv[k + (15 * (i - 1))]:(offset2vv[(k + 
                1) + (15 * (i - 1))] - 1)) {
                ff.s2[i, f] <- inprod(betas2[1:ncX2], Xs2[f, 
                  1:ncX2]) + inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], 
                  Zs2[f, 1:ncZ2])
            }
            f.s2[i, k] <- sum(ff.s2[i, offset2vv[k + (15 * (i - 
                1))]:(offset2vv[(k + 1) + (15 * (i - 1))] - 1)])
            SurvLongR[i, k] <- wk[k] * exp(log.h0.sR[i, k] + 
                alphasR * f.s[i, k] + alphas2R * f.s2[i, k])
            SurvLongD[i, k] <- wk[k] * exp(log.h0.sD[i, k] + 
                alphasD * f.s[i, k] + alphas2D * f.s2[i, k])
        }
        log.survivalR[i] <- -exp(etaBaselineR[i]) * P[i] * sum(SurvLongR[i, 
            ])
        log.survivalD[i] <- -exp(etaBaselineD[i]) * P[i] * sum(SurvLongD[i, 
            ])
        phi[i] <- C - ((eventR[i] * log.hazardR[i]) + (eventD[i] * 
            log.hazardD[i])) - (log.survivalR[i] + log.survivalD[i])
        zeros[i] ~ dpois(phi[i])
        b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
    }
    betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, 
        ])
    betas2[1:ncX2] ~ dmnorm(priorMean.betas2[], priorTau.betas2[, 
        ])
    tau ~ dgamma(priorA.tau, priorB.tau)
    gammasR[1:(ncWR)] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, 
        ])
    gammasD[1:(ncWD)] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, 
        ])
    alphasR ~ dnorm(priorMean.alphas, priorTau.alphas)
    alphas2R ~ dnorm(priorMean.alphas2, priorTau.alphas2)
    alphasD ~ dnorm(priorMean.alphas, priorTau.alphas)
    alphas2D ~ dnorm(priorMean.alphas2, priorTau.alphas2)
    Bs.gammasR[1:(ncW2R)] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, 
        ])
    Bs.gammasD[1:(ncW2D)] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, 
        ])
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
