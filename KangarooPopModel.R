##########################################################################################################################################
## red kangaroo (Osphranter rufus) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")

##############################
## OSPHRANTER (rufus) (OR)

# mass
OR.mass <- 25 # (Croft & Clacy 2008 Macropus rufus. The Mammals of Australia (eds S. van Dyck & R. Strahan), pp. 352-354. Reed New Holland, Sydney)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
OR.rm.pred <- 10^(0.6914 - (0.2622*log10(OR.mass*1000)))
OR.lm.pred <- exp(OR.rm.pred)

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
OR.D.pred <- (10^(4.196 - (0.74*log10(OR.mass*1000))))/2 # divided by 2 for females only
OR.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
OR.age.max <- round(10^(0.89 + (0.13*log10(OR.mass*1000))), 0)
OR.age.max <- 13 # (Jonzén et al. 2010-J Anim Ecol)

## age vector
OR.age.vec <- 0:OR.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
OR.F.pred <- exp(2.719 - (0.211*log(OR.mass*1000)))/2 # divided by 2 for females
OR.F.pred <- 1.5/2 # 1.5 young/year (Jonzén et al. 2010-J Anim Ecol)

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
OR.alpha <- ceiling(exp(-1.34 + (0.214*log(OR.mass*1000))))
OR.alpha <- 2 # (Jonzén et al. 2010-J Anim Ecol)

## define m function with age
OR.m.vec <- c(rep(0, OR.alpha-1), rep(0.75*OR.F.pred, round(OR.alpha/2,0)), rep(OR.F.pred, (OR.age.max+1-((OR.alpha-1+round(OR.alpha/2,0))))))
OR.m.sd.vec <- 0.05*OR.m.vec
plot(OR.age.vec, OR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# Weibull y = a - b*exp(-c*x^d)
OR.m.dat <- data.frame(OR.age.vec, OR.m.vec)
param.init <- c(0.75, 0.75, 0.69, 5.1)
OR.fit.logp <- nls(OR.m.vec ~ a - b*exp(-c*OR.age.vec^d), 
                   data = OR.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3], d = param.init[4]),
                   trace = TRUE,      
                   nls.control(maxiter = 100000, tol = 1e-05, minFactor = 1/1024))
plot(OR.age.vec, OR.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
OR.age.vec.cont <- seq(0,max(OR.age.vec),1)
OR.pred.p.m <- coef(OR.fit.logp)[1] - coef(OR.fit.logp)[2]*exp(-coef(OR.fit.logp)[3]*OR.age.vec.cont^coef(OR.fit.logp)[4])
OR.pred.p.mm <- ifelse(OR.pred.p.m > 1, 1, OR.pred.p.m)
lines(OR.age.vec.cont, OR.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
OR.s.tran <- ln.a.s + b.s*log(OR.mass*1000) + log(1)
OR.s.ad.yr <- exp(-exp(OR.s.tran))

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.99*OR.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.0 # rate of mortality decline (also known as bt)
a2 <- 1 - 1.01*OR.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.2e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.1 # rate of mortality increase
longev <- OR.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
OR.Sx <- c(0.975*OR.s.ad.yr, 1 - qx)
plot(x, OR.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
OR.s.sd.vec <- 0.05*OR.Sx

## create matrix
OR.popmat <- matrix(data = 0, nrow=OR.age.max+1, ncol=OR.age.max+1)
diag(OR.popmat[2:(OR.age.max+1),]) <- OR.Sx[-(OR.age.max+1)]
OR.popmat[OR.age.max+1,OR.age.max+1] <- OR.Sx[OR.age.max+1]
OR.popmat[1,] <- OR.pred.p.mm
colnames(OR.popmat) <- c(0:OR.age.max)
rownames(OR.popmat) <- c(0:OR.age.max)
OR.popmat.orig <- OR.popmat ## save original matrix

## matrix properties
max.lambda(OR.popmat.orig) ## 1-yr lambda
OR.lm.pred
max.r(OR.popmat.orig) # rate of population change, 1-yr
OR.ssd <- stable.stage.dist(OR.popmat.orig) ## stable stage distribution
plot(OR.age.vec, OR.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(OR.popmat.orig, OR.age.max) # reproductive value
OR.gen.l <- G.val(OR.popmat.orig, OR.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
OR.pop.found <- round(area*OR.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
OR.init.vec <- OR.ssd * OR.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*OR.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

OR.tot.F <- sum(OR.popmat.orig[1,])
OR.popmat <- OR.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
OR.n.mat[,1] <- OR.init.vec

## set up projection loop
for (i in 1:t) {
  OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]
}

OR.n.pred <- colSums(OR.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(OR.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
OR.K.max <- 1*OR.pop.found
OR.K.vec <- c(1, OR.K.max/2, 0.75*OR.K.max, OR.K.max) 
OR.red.vec <- c(1,0.88,0.77,0.63)
plot(OR.K.vec, OR.red.vec,pch=19,type="b")
OR.Kred.dat <- data.frame(OR.K.vec, OR.red.vec)

# logistic power function a/(1+(x/b)^c)
OR.param.init <- c(1, 2*OR.K.max, 2)
OR.fit.lp <- nls(OR.red.vec ~ a/(1+(OR.K.vec/b)^c), 
                 data = OR.Kred.dat,
                 algorithm = "port",
                 start = c(a = OR.param.init[1], b = OR.param.init[2], c = OR.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
OR.fit.lp.summ <- summary(OR.fit.lp)
plot(OR.K.vec, OR.red.vec, pch=19,xlab="N",ylab="reduction factor")
OR.K.vec.cont <- seq(1,2*OR.pop.found,1)
OR.pred.lp.fx <- coef(OR.fit.lp)[1]/(1+(OR.K.vec.cont/coef(OR.fit.lp)[2])^coef(OR.fit.lp)[3])
lines(OR.K.vec.cont, OR.pred.lp.fx, lty=3,lwd=3,col="red")

OR.a.lp <- coef(OR.fit.lp)[1]
OR.b.lp <- coef(OR.fit.lp)[2]
OR.c.lp <- coef(OR.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
OR.n.mat <- matrix(0, nrow=OR.age.max+1, ncol=(t+1))
OR.n.mat[,1] <- OR.init.vec
OR.popmat <- OR.popmat.orig

## set up projection loop
for (i in 1:t) {
  OR.totN.i <- sum(OR.n.mat[,i])
  OR.pred.red <- as.numeric(OR.a.lp/(1+(OR.totN.i/OR.b.lp)^OR.c.lp))
  diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.Sx[-(OR.age.max+1)])*OR.pred.red
  OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.Sx[OR.age.max+1])*OR.pred.red
  OR.popmat[1,] <- OR.pred.p.mm
  OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]
}

OR.n.pred <- colSums(OR.n.mat)
plot(yrs, OR.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=OR.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

OR.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
OR.s.arr <- OR.m.arr <- array(data=NA, dim=c(t+1, OR.age.max+1, iter))

for (e in 1:iter) {
  OR.popmat <- OR.popmat.orig
  
  OR.n.mat <- matrix(0, nrow=OR.age.max+1,ncol=(t+1))
  OR.n.mat[,1] <- OR.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    OR.s.alpha <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$alpha
    OR.s.beta <- estBetaParams(OR.Sx, OR.s.sd.vec^2)$beta
    OR.s.stoch <- rbeta(length(OR.s.alpha), OR.s.alpha, OR.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    OR.fert.stch <- rnorm(length(OR.popmat[,1]), OR.pred.p.mm, OR.m.sd.vec)
    OR.m.arr[i,,e] <- ifelse(OR.fert.stch < 0, 0, OR.fert.stch)
    
    OR.totN.i <- sum(OR.n.mat[,i], na.rm=T)
    OR.pred.red <- OR.a.lp/(1+(OR.totN.i/OR.b.lp)^OR.c.lp)
    
    diag(OR.popmat[2:(OR.age.max+1),]) <- (OR.s.stoch[-(OR.age.max+1)])*OR.pred.red
    OR.popmat[OR.age.max+1,OR.age.max+1] <- (OR.s.stoch[OR.age.max+1])*OR.pred.red
    OR.popmat[1,] <- OR.m.arr[i,,e]
    OR.n.mat[,i+1] <- OR.popmat %*% OR.n.mat[,i]

    OR.s.arr[i,,e] <- OR.s.stoch * OR.pred.red
    
  } # end i loop
  
  OR.n.sums.mat[e,] <- ((as.vector(colSums(OR.n.mat))/OR.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

OR.n.md <- apply(OR.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
OR.n.up <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
OR.n.lo <- apply(OR.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,OR.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(OR.n.lo),1.05*max(OR.n.up)))
lines(yrs,OR.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,OR.n.up,lty=2,col="red",lwd=1.5)

OR.s.add <- OR.m.add  <- rep(0, OR.age.max+1)
for (m in 1:iter) {
  OR.s.add <- rbind(OR.s.add, OR.s.arr[ceiling(OR.gen.l):(t+1),,m])
  OR.m.add <- rbind(OR.m.add, OR.m.arr[ceiling(OR.gen.l):(t+1),,m])
}
OR.s.add <- OR.s.add[-1,]
OR.m.add <- OR.m.add[-1,]

OR.s.md <- apply(OR.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
OR.s.up <- apply(OR.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
OR.s.lo <- apply(OR.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(OR.age.vec,OR.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(OR.s.lo),1.05*max(OR.s.up)))
lines(OR.age.vec,OR.s.lo,lty=2,col="red",lwd=1.5)
lines(OR.age.vec,OR.s.up,lty=2,col="red",lwd=1.5)

OR.m.md <- apply(OR.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
OR.m.up <- apply(OR.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
OR.m.lo <- apply(OR.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(OR.age.vec,OR.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(OR.m.lo),1.05*max(OR.m.up)))
lines(OR.age.vec,OR.m.lo,lty=2,col="red",lwd=1.5)
lines(OR.age.vec,OR.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))
