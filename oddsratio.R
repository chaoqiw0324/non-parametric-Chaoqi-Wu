## Method for estimating odds ratio using doubly robust procedure from
## Tchetgen Tchetgen, Robins, and Rotnitzky (2010) 
## "On doubly robust estimation in a semiparametric odds ratio model"
## and Tan (2019)
## "On doubly robust estimation for logistic partially linear models"

## @author: Daniel Malinsky 
## March 2019

library(BB)

# rm(list=ls())
# cat("\f")
# set.seed(2)
# N <- 100000

## A,Y,L binary
# L.true <- rbinom(N,1,0.5) ## L.true <- cbind(L.true,rep(0,N))
# Z <- 0.5 + 0.5*L.true
# pr <- 1/(1+exp(-Z))
# Y.true <- rbinom(N,1,pr)
# A.true <- rbinom(N,1,pr)
# dat <- data.frame(cbind(Y.true,A.true,L.true)) ## Y and A should be independent given L

## A,Y binary, L empty
# A.true <- rbinom(N,1,0.7)
# Z <- 0.5 + 0.5*A.true
# pr <- 1/(1+exp(-Z))
# Y.true <- rbinom(N,1,pr)
# L.true <- rep(0,N)
# dat <- data.frame(cbind(Y.true,A.true,L.true)) ## Y and A should be dependent given L

## A,Y continuous, L empty
# A.true <- rnorm(N)
# Y.true <- rnorm(N,0.5 + 0.8*A.true)
# L.true <- rep(0,N)
# dat <- data.frame(cbind(Y.true,A.true,L.true)) ## Y and A should be dependent given L

## A,Y,L continuous
# L.true <- rnorm(N)
# A.true <- rnorm(N,0.5 + 0.5*L.true)
# Y.true <- rnorm(N,0.5 + 0.8*L.true)
# dat <- data.frame(cbind(Y.true,A.true,L.true)) ## Y and A should be independent given L
# 
#  colnames(dat) <- c("Y","A")
# 
#  Y <- dat$Y
#  A <- dat$A
#  L <- dat[,-c(1,2)]

## A,Y,L continuous, dim(A)=2
# L1.true <- rnorm(N)
# L2.true <- rnorm(N,4,1)
# A1.true <- rnorm(N,0.5 + 0.5*L1.true + 0.8*L2.true)
# A2.true <- rnorm(N,0.2 + 0.7*L1.true + 0.8*L2.true)
# Y.true <- rnorm(N,0.5 + 0.8*L1.true + 0.8*L2.true)
# data <- data.frame(cbind(Y.true,A1.true,A2.true,L1.true,L2.true)) ## Y and A should be independent given L
# Y <- as.matrix(data[,1])
# A <- as.matrix(data[,c(2,3)])
# #L <- as.matrix(data[,c(4,5)])
# L <- as.matrix(data[,c(4)]) # using this, Y and A should be not independent given L

## A,L continuous, Y binary
# L1.true <- rnorm(N)
# L2.true <- rnorm(N,4,1)
# A1.true <- rnorm(N,0.5 + 0.5*L1.true + 0.8*L2.true)
# A2.true <- rnorm(N,0.2 + 0.7*L1.true + 0.6*L2.true)
# # Z <- 0.5 + 0.8*L1.true + 0.8*L2.true ## Y \indep of A given L
# Z <- 0.5 + 0.4*L1.true + 0.8*L2.true + 0.4*A1.true - 1.5*A2.true
# pr <- 1/(1+exp(-Z))
# Y <- rbinom(N,1,pr)
# data <- data.frame(cbind(Y,A1.true,A2.true,L1.true,L2.true))
# Y <- as.matrix(data[,1])
# A <- as.matrix(data[,c(2,3)])
# L <- as.matrix(data[,c(4,5)])
# # L <- as.matrix(data[,c(4)]) # using this, Y and A should be not independent given L



psi.hat <- function(Y, A, L=c(), subset = NULL, out.bin = TRUE, exp.bin = FALSE, exp.scalar = FALSE){
  ## Function: estimate the odds ratio parameter psi
  ## Input: 1. An outcome nuisance model, onm = f(Y|L,A=0)
  ##        2. An exposure nuisance model, enm = g(A|Y=0,L)
  ## Output: A real number (vector) psi.hat, an estimate of the conditional odds ratio parameter psi
  
  ## Todo: generalize estimating eq to arbitrary dimensions
  
  #temp
  if(!exp.bin){
    A[,1] <- A[,1] - mean(A[,1],na.rm=TRUE)
    A[,2] <- A[,2] - mean(A[,2],na.rm=TRUE)
    L[,1] <- L[,1] - mean(L[,1],na.rm=TRUE)
    if(ncol(L) > 1) L[,2] <- L[,2] - mean(L[,2],na.rm=TRUE)
  }
  
  if(length(L)==0){
    L <- cbind(rep(0,length(Y)),rep(0,length(Y))) ## for cases where conditioning set is empty
    colnames(L) <- c("L1","L2")
  }
  
  if(!is.null(subset)) Y <- Y[subset]
  if(!is.null(subset)) A <- A[subset,]
  if(!is.null(subset)) L <- L[subset,]
  #temp
  
  dat <- data.frame(Y,A,L)
  
  refA <- 0 
  refY <- 0
  
  if(exp.scalar){
  
    if(out.bin && exp.bin){
      h.dag <- 0.25 ## probability f.dag(Y|L) = g.dag(A|L) = 0.5, i.e., Y ~ A ~ Bernoulli(0.5)
    
      outcome <- glm(Y ~ A + L, family = binomial)
      dat1 <- data.frame(cbind(A,L))
      dat1[,1] <- 0 ## setting A=0
      onm <- predict.glm(outcome, newdata=dat1,type="response")
      onm[Y==0] <- 1-onm[Y==0]
    
      exposure <- glm(A ~ Y + L, family = binomial)
      dat2 <- data.frame(cbind(Y,L))
      dat2[,1] <- 0 ## setting Y=0
      enm <- predict.glm(exposure, newdata=dat2,type="response")
      enm[A==0] <- 1-enm[A==0]
    
      d.diff <- (-1)^(Y+A) ## Eric's suggestion, for Y,A binary
    
      U <- function(psi,onm,enm,d.diff){ sum( d.diff*h.dag / (exp(psi*Y*A)*onm*enm) ) }
      est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001, onm = onm, enm = enm, d.diff=d.diff)
      return(est$root)
    }
  
    if(!out.bin){ 
      outcome <- glm(Y ~ A + L, family = gaussian)
    } else outcome <- glm(Y ~ A + L, family = binomial)
    dat1 <- data.frame(cbind(A,L))
    dat1[,1] <- 0 ## setting A=0
    onm <- predict.glm(outcome, newdata=dat1,type="response")
    onm[Y==0] <- 1-onm[Y==0]
  
    if(!exp.bin){
      exposure <- glm(A ~ Y + L, family = gaussian)
    } else exposure <- glm(A ~ Y + L, family = binomial)
    dat2 <- data.frame(cbind(Y,L))
    dat2[,1] <- 0 ## setting Y=0
    enm <- predict.glm(exposure, newdata=dat2,type="response")
    enm[A==0] <- 1-enm[A==0]
  
  
    #U <- function(psi,onm,enm){ sum( (Y - onm)*(A - enm)*exp(-psi*Y*A) )}
    #est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001, onm = onm, enm = enm)
    #return(est$root)
  
    U <- function(psi){ sum( (Y - onm)*(A - enm)*exp(-psi*Y*A) )}
    par.init <- c(1.0)
    sol <- BBsolve(par=par.init, fn=U, quiet=TRUE) 
    est <- sol$par
  
    return(est)
  } ## end if(exposure is scalar)
  
  if(!exp.scalar){
    
    d <- ncol(A)
    exposure <- list()
    enm <- list()
    
    # if(out.bin && exp.bin){
    #   cat("error!") ## do this binary case later
    # }
    
    fmla = as.formula(paste("Y ~", paste0(colnames(A),collapse="+"), "+", paste0(colnames(L),collapse="+")))
    if(!out.bin){ 
      outcome <- glm(fmla, dat, family = gaussian)
    } else outcome <- glm(fmla, dat, family = binomial)  
    dat1 <- dat[,-1] #data.frame(cbind(A,L))
    for(j in 1:d){ dat1[,j] <- refA } ## setting A=0
    onm <- predict.glm(outcome, newdata=dat1, type="response")
    onm[Y==0] <- 1-onm[Y==0]
    
    if(!exp.bin){
      for(j in 1:d){
        fmla_a = as.formula(paste(colnames(A)[j], " ~ Y + ", paste0(colnames(L),collapse="+") ))
        exposure[[j]] <- glm(fmla_a, dat, family = gaussian)
      }
    } else {
      for(j in 1:d){
        fmla_a = as.formula(paste(colnames(A)[j], " ~ Y + ", paste0(colnames(L),collapse="+") ))
        exposure[[j]] <- glm(fmla_a, dat, family = binomial)
      }
    }
    dat2 <- dat[,-c(2:(d+1))] #data.frame(cbind(as.numeric(Y),L))
    colnames(dat2) = c("Y",colnames(L))
    dat2[,1] <- refY ## setting Y=0
    for(j in 1:d){
      enm[[j]] <- predict.glm(exposure[[j]], newdata=dat2, type="response")
      enm[[j]][A[,j]==0] <- 1-enm[[j]][A[,j]==0]
    }
  
    U <- function(psi){
      f <- rep(NA,d)
      for(j in 1:d){
        f[j] <- sum( (Y - onm)*(A[,j] - enm[[j]])*exp( -psi[1]*Y*A[,1] -psi[2]*Y*A[,2] ) ) # need to generalize to arbitrary d!
        #f[j] <- sum( (Y - onm)*(A[,j] - enm[[j]])*exp( ... ) ) ## psi1*Y*A1 + psi2*Y*A2 + ... psid*Y*Ad
      }
      f
    }
    
    par.init <- outcome$coefficients[2:(d+1)]
    sol <- BBsolve(par=par.init, fn=U, control=c(noimp=2000, maxit=5000), quiet=TRUE) 
    est <- sol$par
    
    if(sol$convergence > 0) {
      cat("ERROR: DR odds ratio estimator failed", "\n")
      est <- par.init ## ad hoc fix in case method fails
    }
    
    return(est)
      
  }
  
}

## execute:
#(psi.hat <- psi.hat(Y,A,L))
#(psi.hat <- psi.hat(Y,A,L,exp.scalar=TRUE)) ## if both binary
#(psi.hat <- psi.hat(Y,A,L,out.bin=FALSE,exp.bin=FALSE,exp.scalar = TRUE)) ## if both not binary
#(psi.hat <- psi.hat(Y,A,L,out.bin=FALSE,exp.bin=FALSE)) ## if both continuous
#(psi.hat <- psi.hat(Y,A,L,out.bin=TRUE,exp.bin=FALSE)) ## binary outcome, continuous exposures
#exp(psi.hat) ## the odds ratio
#summary(glm(Y~A+L,family=binomial))
#summary(glm(Y~offset(psi.hat[1]*A[,1] + psi.hat[2]*A[,2])+L,family=binomial))
