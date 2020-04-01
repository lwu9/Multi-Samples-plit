###### contains single/multiple splits, whether to switch #####
###### copy from the krr_parallel.R ######
list.of.packages <- c("foreach", "parallel", "doParallel", "parallel", "listdtr", "ramify",
                      "randomForest", "mvtnorm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(foreach)
library(doParallel)
library(parallel)
library(listdtr)
library(ramify)
library(mvtnorm)
library(randomForest)

aipw <- function(actions, x.out, a.out, y.out, psm) {
  ## psm -- propensity score model: 1--average of a.out; 2--logistic regression
  if (psm==1) {
    gamma.est<-glm(a.out ~ 1, family="binomial")$coef
    phat<-rep(1/(1+exp(-gamma.est )), length(a.out))
    q <- exp(gamma.est) / (1+exp(gamma.est))**2 ## q is needed in var estimation
  } else if (psm==2) {
    gamma.est<-glm(a.out ~ x.out, family="binomial")$coef
    lpshat<-cbind(1,x.out)%*%gamma.est
    phat<-1/(1+exp(-lpshat ))
    q <- exp(lpshat) / (1+exp(lpshat))**2
  } 
  a.out.pred <- phat; a.out.pred[actions==0] <- 1-phat[actions==0]
  C <- (actions==a.out)
  lm_fit.out <- lm(y.out ~ x.out + a.out*x.out)
  betas <- lm_fit.out$coefficients
  Q <- cbind(rep(1, length(a.out)), x.out, actions, actions*x.out)%*% as.matrix(betas, length(betas),1)
  aipw.list <- C*y.out/a.out.pred-(C-a.out.pred)/a.out.pred*Q
  return(list(aipw=aipw.list, Q.hat=lm_fit.out$fitted.values, lm_fit.out=lm_fit.out, ps.hat=phat, 
              gamma.est=gamma.est, q=q))
}

aipw.diff.var.est <- function(X, A, Y, d1.actions, d2.actions, a.out.pred, 
                              v1.out, v2.out, Q.hat) {
  ### this is for the aipw with known propensity score (=0.5) and correctly specified linear model for Q function
  # X: covariate matrix not including intercepts
  # d1.actions: array of actions from policy d1
  # a.out.pred: know propensity score
  # v1.out: array of aipw value estimate from policy d1
  # Q.hat: array of fitted Q values
  # return: the estimated variance of the difference of the value estimtates of two policies d1 and d2
  # From the M-estimation theory to get the asymptotic variance (https://www4.stat.ncsu.edu/~davidian/dtr19/dtr2.pdf, page 79)
  n <- dim(X)[1]
  X <- cbind(rep(1, n), X) # add intercept term
  p <- dim(X)[2]
  An <- mean((v1.out-v2.out-mean(v1.out-v2.out))**2)
  Cn <- matrix(NA, 2*p, 1)
  Cn[1:p,] <- as.matrix(apply(X*((d2.actions==A)-(d1.actions==A))/a.out.pred, 2, mean), p, 1)
  Cn[(p+1):(2*p),] <- as.matrix(apply(X*d2.actions*((d2.actions==A)/a.out.pred-1)-
                                        X*d1.actions*((d1.actions==A)/a.out.pred-1), 2, mean), p, 1)
  S <- cbind(X, A*X)
  Dn <- t(S)%*%S/n
  Dn.inv <- solve(Dn)
  S.prime <- S*(Y-Q.hat)
  Bn <- t(S.prime)%*%S.prime/n
  var.est <- An+t(Cn)%*%Dn.inv%*%Bn%*%Dn.inv%*%Cn
  return(var.est)
}

aipw.diff.var.est.general <- function(X, A, Y, d1.actions, d2.actions, 
                              ps.hat, v1.out, v2.out, q, lm_fit.out) {
  ### this is for a more general aipw with logistic model for propensity score and linear model for Q function
  # X: covariate matrix not including intercepts
  # d1.actions: array of actions from policy d1
  # ps.hat: array of estimated propensity score P(A=1|X) -- \pi(X, A; \hat{\gamma)}
  # v1.out: array of aipw value estimate from policy d1
  # Q.hat: array of fitted Q values
  # q: exp(1+x^t \hat{\gamma}) / ( exp(1+x^t \hat{\gamma}) )^2
  # return: the estimated variance of the difference of the value estimtates of two policies d1 and d2
  # Based on the M-estimation theory to get the asymptotic variance (https://www4.stat.ncsu.edu/~davidian/dtr19/dtr2.pdf, page 79)
  n <- dim(X)[1]
  X <- cbind(rep(1, n), X) # add intercept term
  p <- dim(X)[2]
  p.gamma <- ifelse(length(q)>1, p, 1)
  Q.hat <- lm_fit.out$fitted.values
  betas <- lm_fit.out$coefficients
  Q.d1 <- cbind(X, d1.actions*X) %*% as.matrix(betas, length(betas),1)
  Q.d2 <- cbind(X, d2.actions*X) %*% as.matrix(betas, length(betas),1)
  
  # d1.ps: array of the estimated propensity score based on d1: \pi(X, d1; \hat{\gamma})
  ps.d1 <- ps.d2 <- ps.hat # initialize
  ps.d1[which(d1.actions == 0)] <- 1 - ps.hat[which(d1.actions == 0)]
  ps.d1[which(d2.actions == 0)] <- 1 - ps.hat[which(d2.actions == 0)]
  C.d2 <- (d2.actions==A); C.d1 <- (d1.actions==A)
  M1 <- v1.out-v2.out-mean(v1.out-v2.out)
  
  Fn <- mean((M1)**2)

  Bn.beta <- matrix(NA, 2*p, 1)
  Bn.beta[1:p,] <- as.matrix(apply(X*(C.d2/ps.d2 - C.d1/ps.d1), 2, mean), p, 1)
  Bn.beta[(p+1):(2*p),] <- as.matrix(apply(X*d2.actions*(C.d2/ps.d2-1) - X*d1.actions*(C.d1/ps.d1-1), 2, mean), p, 1)
  Bn.gamma <- matrix(NA, p, 1)
  if (p.gamma==1) {
    Bn.gamma <- mean(X[,p.gamma] * q * c(C.d2*((d2.actions==1)-(d2.actions==0))*(Y-Q.d2) / ps.d2**2 - 
                                                        C.d1*((d1.actions==1)-(d1.actions==0))*(Y-Q.d1) / ps.d1**2))
  } else {
    Bn.gamma <- as.matrix(apply(X * q * c(C.d2*((d2.actions==1)-(d2.actions==0))*(Y-Q.d2) / ps.d2**2 - 
                                                        C.d1*((d1.actions==1)-(d1.actions==0))*(Y-Q.d1) / ps.d1**2), 2, mean), p, 1)
    
  }
  Bn <- t(rbind(Bn.gamma, Bn.beta))
  
  S <- cbind(X, A*X)
  Dn.beta <- t(S) %*% S / n
  Dn.beta.inv <- -solve(Dn.beta)
  Dn.gamma <- t(X[,1:p.gamma]*sqrt(q)) %*% (X[,1:p.gamma]*sqrt(q)) / n
  Dn.gamma.inv <- -solve(Dn.gamma)
  Dn.inv <- rbind(cbind(Dn.gamma.inv, matrix(0, nrow=dim(Dn.gamma.inv)[1], ncol=dim(Dn.beta.inv)[2])), 
                  cbind(matrix(0, nrow=dim(Dn.beta.inv)[1], ncol=dim(Dn.gamma.inv)[2]), Dn.beta.inv))
  
  S.prime <- S*(Y-Q.hat)
  Gn.beta <- t(S.prime) %*% S.prime / n
  X.prime <- X[,1:p.gamma]*(A-ps.hat)
  Gn.gamma1 <- t(X.prime) %*% (X.prime) / n
  Gn.gamma2 <- t(X.prime) %*% (S.prime) / n
  Gn <- rbind(cbind(Gn.gamma1, Gn.gamma2), cbind(t(Gn.gamma2), Gn.beta))
  
  if (p.gamma==1) {
    Hn.gamma <- mean(X.prime * c(M1))
  } else {
    Hn.gamma <- as.matrix(apply(X.prime * c(M1), 2, mean), p.gamma, 1)
  }
  Hn.beta <- as.matrix(apply(S.prime * c(M1), 2, mean), 2*p, 1)
  Hn <- rbind(Hn.gamma, Hn.beta)
  
  term <- Bn %*% Dn.inv %*% Hn
  var.est <- Fn - term - t(term) +  Bn %*% Dn.inv %*% Gn %*% Dn.inv %*% t(Bn)
  if (var.est < 0) browser()
  return(var.est)
}

aipw_krr_Q <- function(actions, x.out, a.out, y.out, a.out.pred) {
  ## using KRR to estimate Q
  C <- (actions==a.out)
  obj <- krr(x.out, y.out, group = a.out)
  pred <- predict(obj, x.out)
  Q <- sapply(1:length(actions), function(i) pred[i, 1+actions[i]])
  fitted.values <-sapply(1:length(y.out), function(i) pred[i, 1+a.out[i]])
  aipw.list <- C*y.out/a.out.pred-(C-a.out.pred)/a.out.pred*Q
  return(list(aipw=aipw.list, Q.hat=fitted.values))
}
  
aipw_krr_Q.diff.var.est <- function(x.out, a.out, y.out, d1.actions, d2.actions, a.out.pred) {
  d1.aipw.list <- aipw_krr_Q (d1.actions, x.out, a.out, y.out, a.out.pred)$aipw
  d2.aipw.list <- aipw_krr_Q (d2.actions, x.out, a.out, y.out, a.out.pred)$aipw
  var.est <- mean((d1.aipw.list - d2.aipw.list)**2)
  return(var.est)
}
  
pValue.oneSplit.krr <- function(ind.in, x, a, y, Switch=F, delta_m=T, true.diff=0, L=10L) {
  n <- length(y)
  ind.out <- setdiff(1:n, ind.in)
  x.in <- x[ind.in,]; a.in <- a[ind.in]; y.in <- y[ind.in]
  x.out <- x[ind.out,]; a.out <- a[ind.out]; y.out <- y[ind.out]
  n.out <- length(y.out)
  delta.m <- log(log10(2*n.out))/(2*n.out)^(1/6)
  stage.x <- rep(1, dim(x)[2]) # consider 1-stage here
  # decision list
  list.fit <- listdtr(y.in, a.in, x.in, stage.x, maxlen = L)
  pred.list <- as.numeric(as.character(predict(list.fit, x.out, 1)))
  # krr
  krr.fit <- krr(x.in, y.in, group = a.in)
  pred <- predict(krr.fit, x.out)
  pred.krr <- argmax(pred) - 1
  
  # estimate Propensity score is more efficient than using the true PS 
  #a.out.pred <- ifelse(a.out==1, mean(a.out), 1-mean(a.out))
  # a.out.pred <- rep(0.5, length(y.out)) # use the true PS
  #table(a.in, a.in.pred)
  if (sum(abs(pred.list-pred.krr))==0) return(1.0) else {
    result <-  aipw(pred.list, x.out, a.out, y.out, psm=1)
    v.out.list <- result$aipw
    Q.hat <- result$Q.hat
    v.out.krr <- aipw(pred.krr, x.out, a.out, y.out, psm=1)$aipw
    diff.v <- mean(v.out.list) - mean(v.out.krr)
    # var.est <- aipw.diff.var.est(x.out, a.out, y.out, pred.list, pred.krr, 
    #                              a.out.pred, v.out.list, v.out.krr, Q.hat)
    var.est <- aipw.diff.var.est.general(x.out, a.out, y.out, pred.list, pred.krr, 
                                         result$ps.hat, v.out.list, v.out.krr, result$q, result$lm_fit.out)
  }
  if (delta_m==T) {
    # use delta_m in the denominator
    test.stat <- (diff.v-true.diff)*sqrt(n.out/max(var.est, delta.m^2))
  } else {
    test.stat <- (diff.v-true.diff)*sqrt(n.out/var.est)
  }
  
#  print(paste(test.stat, diff.v*sqrt(n.out/var.est), var.est, delta_m^2))
  # true values of the estimated policies from the Data
  p.value.raw <- pnorm(test.stat)
  if (Switch) {
    p.value.raw2 <- pValue.oneSplit.krr(ind.out, x, a, y)
    p.value.raw <- 2 * min(c(p.value.raw, p.value.raw2))
  }
  # browser()
  return(list(p.value.raw=p.value.raw, var.est=var.est))
}

Q.gamma <-  function(gamma, pvalues) {
  return(min(1, quantile(pvalues/gamma, gamma)))
}

pValue.multiSplit.krr <- function(x, a, y, nB, Switch, delta_m, true.diff=0, L=10L) {
  # nB: the number of splits
  # true.diff: \in (-inf,0]
  # return: nB p-values
  n <- length(y)
  pValues <- var.ests <- c(); b=0
  while(b < nB) {
    b=b+1
    ind.in <- sample(1:n, size=n/2)
    #pvalue <- min(pValue.oneSplit(ind.in, x, a, y), pValue.oneSplit(setdiff(1:n, ind.in), x, a, y))
    check <- tryCatch({
      pvalue <- pValue.oneSplit.krr(ind.in, x, a, y, Switch=Switch, delta_m, true.diff=true.diff, L=L)
      pvalue.var.est <- pvalue$var.est
      pvalue <- pvalue$p.value.raw
      #print(paste(b, pvalue))
    }, error=function(e) {return("error")} )
    if (typeof(check) == "character") {
      nB = nB+1
      #print(c(b, nB, length(pValues), length(var.ests)))
    } else {
      pValues <- c(pValues, pvalue)
      var.ests <- c(var.ests, pvalue.var.est)
    }
  }
  return(list(pValues=pValues, var.ests=var.ests))
}

pValue.oneSplit.krr_Q <- function(ind.in, x, a, y, Switch=F, L=10L) {
  # Estimate Q usig KRR
  n <- length(y)
  ind.out <- setdiff(1:n, ind.in)
  x.in <- x[ind.in,]; a.in <- a[ind.in]; y.in <- y[ind.in]
  x.out <- x[ind.out,]; a.out <- a[ind.out]; y.out <- y[ind.out]
  n.out <- length(y.out)
  delta_m <- log(log10(2*n.out))/(2*n.out)^(1/6)/10
  stage.x <- rep(1, dim(x)[2]) # consider 1-stage here
  # decision list
  list.fit <- listdtr(y.in, a.in, x.in, stage.x, maxlen = L)
  pred.list <- as.numeric(as.character(predict(list.fit, x.out, 1)))
  # krr
  krr.fit <- krr(x.in, y.in, group = a.in)
  pred <- predict(krr.fit, x.out)
  pred.krr <- argmax(pred) - 1
  
  # estimate Propensity score is more efficient than using the true PS 
  #a.out.pred <- ifelse(a.out==1, mean(a.out), 1-mean(a.out))
  a.out.pred <- rep(0.5, length(y.out)) # use the true PS
  #table(a.in, a.in.pred)
  if (sum(abs(pred.list-pred.krr))==0) return(1.0) else {
    result <- aipw_krr_Q(pred.list, x.out, a.out, y.out, a.out.pred)
    v.out.list <- result$aipw
    Q.hat <- result$Q.hat
    v.out.krr <- aipw_krr_Q(pred.krr, x.out, a.out, y.out, a.out.pred)$aipw
    diff.v <- mean(v.out.list) - mean(v.out.krr)
    var.est <- aipw_krr_Q.diff.var.est(x.out, a.out, y.out, pred.list, pred.krr, 
                                 a.out.pred)
  }
  test.stat <- diff.v*sqrt(n.out/max(var.est, delta_m^2))
  #  print(paste(test.stat, diff.v*sqrt(n.out/var.est), var.est, delta_m^2))
  # true values of the estimated policies from the Data
  p.value.raw <- pnorm(test.stat)
  if (Switch) {
    p.value.raw2 <- pValue.oneSplit.krr(ind.out, x, a, y)
    p.value.raw <- 2 * min(c(p.value.raw, p.value.raw2))
  }
  return(p.value.raw)
}

pValue.multiSplit.krr_Q <- function(x, a, y, nB, Switch=F) {
  # nB: the number of splits
  # return: nB p-values
  n <- length(y)
  pValues <- c(); b=0
  while(b < nB) {
    b=b+1
    ind.in <- sample(1:n, size=n/2)
    #pvalue <- min(pValue.oneSplit(ind.in, x, a, y), pValue.oneSplit(setdiff(1:n, ind.in), x, a, y))
    check <- tryCatch({
      pvalue <- pValue.oneSplit.krr_Q(ind.in, x, a, y, Switch=Switch)
      #print(paste(b, pvalue))
    }, error=function(e) {return("error")} )
    if (typeof(check) == "character") {
      nB = nB+1
      # print(c(b, B, length(pValues)))
    } else
      pValues <- c(pValues, pvalue)
  }
  return(pValues)
}
