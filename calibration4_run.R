# ##############################################################################################
source("/home/lwu9/sampleSplit/functions2.r")
source("/home/lwu9/sampleSplit/calibration4.R")

p=10
stage.x <- rep(1, p)
cov_matrix <- matrix(0, nrow = p, ncol = p)
for (j1 in 1:p) {
  for (j2 in 1:p) {
    cov_matrix[j1, j2] <- 0.2**abs(j1-j2)*4
  }
}
# Contrast function expressions
contrast.expressions <- list()
contrast.expressions[["H1FALSE"]] <- c("sin(x[, 1] * pi / 8)", "(x[,1] - 0.1) * (x[,2] + 0.1)")
contrast.expressions[["H1TRUE"]] <- c("(x[, 1] - x[, 2] + x[, 3] - x[, 4])","(x[, 1] + x[, 2] - 1)",
                                      "atan(exp(1+x[,1]) - 3*x[,2] - 5)")

Replicates <- 100; nB <- 100; gamma.min <- 0.01; alpha_target <- 0.05
tune <- TRUE; nfold <- 10; Switch <- F; maxL <- 10L
print(paste0("Replicates=", Replicates, ", tune=", tune, ", nB=", nB, ", gamma.min=adapt"))
rs <- 1:Replicates; results.list <- list()
for (n in c(200, 400)) {
  for (H1 in c(F)) {
    for (icont in 2: length(contrast.expressions[[paste0("H1",H1)]])) {
      # sigma <- 0.1
      if (H1 & (icont>1)) sigma <- 0.1 else
        sigma <- 2
      registerDoParallel(Replicates)
      pvalues.multi <- foreach(seed=rs, .combine=rbind, #.errorhandling = 'remove',
                               .packages=c('listdtr','mvtnorm','ramify')) %dopar% {
                                 set.seed(seed)
                                 a <- rbinom(n, 1, 0.5)
                                 x <- mvtnorm::rmvnorm(n=n, mean=rep(0, p), sigma=cov_matrix)
                                 contrast <- eval(parse(text=contrast.expressions[[paste0("H1",H1)]][icont]))
                                 y <- 2 + apply(x[, c(1,3,5,7)], 1, sum) + a * contrast + rnorm(n, 0, sigma)
                                 dlist <- listdtr::listdtr(y, a, x, stage.x = rep(1, dim(x)[2]), maxlen=maxL)
                                 d <- as.numeric(as.character(predict(dlist, x, 1)))

                                 w <- tune_parameters_find_coef_vectors(x, a, y, d, type="rbf", tune=tune, nfold=nfold)
                                 y_hat <- w$y_hat
                                 # y_hat <-  2 + apply(x[, c(1,3,5,7)], 1, sum) + a * contrast
                                 residual <- y - y_hat
                                 cat("Mean and sd of residuals: ", mean(residual), sd(residual), "\n")
                                 # paras <- find_alpha_tilde(residual, y_hat, 1:100, nB=nB, x, a,
                                                           # alpha_target, Switch=Switch, L=maxL, gamma.min=gamma.min)
                                 paras <- find_alpha_tilde3(residual, y_hat, 1:200, nB=nB, x, a,
                                                            alpha_target, Switch=Switch, L=maxL)
                                 alpha.tilde <- paras$alpha.tilde
                                 gamma.mins <- paras$gamma.min
                                 # gamma.mins <- gamma.min
                                 diff_min <- paras$diff_min
                                 pValues <- pValue.multiSplit.krr(x, a, y, nB=nB, Switch=Switch, delta_m=T, L=maxL)[[1]]
                                 return(list(pValues=pValues, alpha.tilde=alpha.tilde, residual=residual, 
                                             gamma.mins=gamma.mins, diff_min=diff_min))
                               }
      stopImplicitCluster()
      print(paste0("n",n,"_H1",H1,"_expr",icont))
      if (nB > 1) {
        p_adapt <- compute_adative_pvalues(pvalues.multi)
      } else {
        p_adapt <- unlist(pvalues.multi[,1])
      }
      results.list[[paste0("n",n,"_H1",H1,"_expr",icont)]] <- p_adapt
      print("reject rate fixed alpha 0.05:")
      print(apply(p_adapt < alpha_target, 2, mean))
      print("reject rate fixed alpha tilde calibrated:")
      print(apply(p_adapt < unlist(pvalues.multi[,2]), 2, mean))
      cat('alpha.tilde', unlist(pvalues.multi[,2]), "\n")
      cat('gamma.min', unlist(pvalues.multi[,4]), "\n")
      cat('diff_min', unlist(pvalues.multi[,5]), "\n")
      # print(apply(pvalues.multi, 2, mean))
      # print(mean(pvalues.multi[,1] <= pvalues.multi[,2]))
    }
  }
}
result <- list(pvalues.multi=pvalues.multi, results.list=results.list)
save(result,
     file=paste0("/home/lwu9/sampleSplit/pvalues_multi_list_var_est_general_sigma_",sigma,"_H1",H1, "_n", n, "_icont",
                 icont, "_B", nB, "_rep", Replicates,"_tune", tune, "_nfold", nfold, "_Switch", Switch, "_maxL", maxL,
                 "_gamMinTune", "_calib4.RData"))














