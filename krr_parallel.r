###### sims for single and multiple split ######
################################################
#source("/Users/lwu9/Documents/lab_proj/SP/functions2.r")
source("/home/lwu9/sampleSplit/functions2.r") # use this in cluster

p=10
stage.x <- rep(1, p)
cov_matrix <- matrix(0, nrow = p, ncol = p)
for (j1 in 1:p) {
  for (j2 in 1:p) {
    cov_matrix[j1, j2] <- 0.2**abs(j1-j2)*4
  }
}
sigma <- 0.1; #H1=F
# ############ single split ##############
# Replicates = 49
# rs <- 1:Replicates; alpha = 0.05; n <- 300; Switch <- F; 
# registerDoParallel(48) 
# results <- foreach(seed=rs, .combine=c, .errorhandling = 'remove',
#                    .packages=c('listdtr','mvtnorm','randomForest','ramify')) %dopar% {
#                      set.seed(seed)
#                      ind.in <- sample(1:n, size=n/2)
#                      a <- rbinom(n, 1, 0.5)
#                      x <- rmvnorm(n=n, mean=rep(0, p), sigma=cov_matrix)
#                      if (H1==T) {
#                        # under H1
#                        y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * (x[, 1] - x[, 2] + x[, 3] - x[, 4]) + rnorm(n, 0, sigma)
#                      } else {
#                        # under H0
#                        #y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * sin(x[, 1] * pi / 8) + rnorm(n, 0, sigma)
#                        y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * (x[,1] - 0.1) * (x[,2] + 0.1) + rnorm(n, 0, sigma)
#                      }
#                       return(pValue.oneSplit.krr(ind.in, x, a, y, Switch = Switch))
#                    }
# stopImplicitCluster()
# print(paste0("H1=", H1, ", reject rate=", round(mean(results < alpha), 2), 
#              ", replicates=", length(results), ", Switch=", Switch, ", sigma=", sigma))
# opt.d <- rep(0, n)
# opt.d[(x[,1] - 0.1) * (x[,2] + 0.1) > 0] <- 1
# list.fit <- listdtr(y,a,x,rep(1, dim(x)[2]))
# tbl <- table(opt.d, predict(list.fit, x, 1))
# sum(diag(tbl))/sum(tbl)


############ multiple splits ##############
### refer to paper "P-Values for High-Dimensional Regression" ###
contrast.expressions <- list()
contrast.expressions[["H1FALSE"]] <- c("sin(x[, 1] * pi / 8)", "(x[,1] - 0.1) * (x[,2] + 0.1)")
#contrast.expressions[["H1TRUE"]] <- c("(x[, 1] - x[, 2] + x[, 3] - x[, 4])","(x[, 1] + x[, 2] - 1)", 
#                                      "atan(exp(1+x[,1]) - 3*x[,2] - 5)")
contrast.expressions[["H1TRUE"]] <- c("(x[, 1] - x[, 2] + x[, 3] - x[, 4])")
Replicates = 200; nB = 600; alpha = 0.05;
print(paste0("Replicates=", Replicates, "sigma=", sigma))
rs <- 1:Replicates; results.list <- list()
for (n in c(200)) {
  for (H1 in c(F)) {
   for (icont in 1: length(contrast.expressions[[paste0("H1",H1)]])) {
	   sigma <- 2.5
  # if (H1 & (icont>1)) sigma <- 0.1 else
#	   sigma <- 1.0
   # for (icont in c(1)) {
      registerDoParallel(96) 
      pvalues.multi <- foreach(seed=rs, .combine=rbind, #.errorhandling = 'remove',
                               .packages=c('listdtr','mvtnorm','randomForest','ramify')) %dopar% {
                                 set.seed(seed)
                                 a <- rbinom(n, 1, 0.5)
                                 x <- rmvnorm(n=n, mean=rep(0, p), sigma=cov_matrix)
                                 contrast <- eval(parse(text=contrast.expressions[[paste0("H1",H1)]][icont]))
                                 y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * contrast + rnorm(n, 0, sigma)
                                 return(pValue.multiSplit.krr(x, a, y, nB=nB, Switch=F, delta_m = T, L=3L))
                               }
      stopImplicitCluster()
      print(paste0("n",n,"_H1",H1,"_expr",icont))
      results.list[[paste0("n",n,"_H1",H1,"_expr",icont)]] <- pvalues.multi
    }
  }
}
save(results.list, file=paste0("/home/lwu9/sampleSplit/pvalues_multi_list_var_est_general_sigma1_2dot5_H1," H1,"_n", n, "_icont,", icont, "_nB", nB, ".RData"))

# registerDoParallel(48) 
# pvalues.multi <- foreach(seed=rs, .combine=rbind, .errorhandling = 'remove',
#                    .packages=c('listdtr','mvtnorm','randomForest','ramify')) %dopar% {
#                      set.seed(seed)
#                      a <- rbinom(n, 1, 0.5)
#                      x <- rmvnorm(n=n, mean=rep(0, p), sigma=cov_matrix)
#                      if (H1==T) {
#                        # under H1
#                        # y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * (x[, 1] + x[, 2] - 1) + rnorm(n, 0, sigma)
#                        y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * (x[, 1] - x[, 2] + x[, 3] - x[, 4]) + rnorm(n, 0, sigma)
#                        #y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * atan(exp(1+x[,1]) - 3*x[,2] - 5) + rnorm(n, 0, sigma)
#                      } else {
#                        # under H0
#                        #y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * sin(x[, 1] * pi / 8) + rnorm(n, 0, sigma)
#                        y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * (x[,1] - 0.1) * (x[,2] + 0.1) + rnorm(n, 0, sigma)
#                      }
#                      return(pValue.multiSplit.krr(x, a, y, nB=nB, Switch=F))
#                    }
# stopImplicitCluster()

print(dim(pvalues.multi))
nRep = dim(pvalues.multi)[1]
pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=nRep, ncol=nB)
for (gamma.min in seq(0.01,0.02,0.005)) {
  reject.rate <- mean(pvalues[,1] < alpha)
  for (b in 2:nB) {
    p.values <- c()
    for (r in 1:nRep) {
      pValues <- pvalues[r, 1:b]
      Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
      p.value <- min(1, (1-log(gamma.min))*Q_inf) ## adaptive way to find gamma
      # p.value <- Q.gamma(0.5, pValues)
      p.values <-c(p.values, p.value)
    }
    reject.rate <- c(reject.rate, mean(p.values < alpha))
  }
  cat(gamma.min, reject.rate,"\n")
  plot(reject.rate, main=paste0("H1=", H1,", sigma=", sigma))
}

