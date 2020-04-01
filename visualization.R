setwd("/Users/lwu9/Documents/lab_proj/SP")
source("/Users/lwu9/Documents/lab_proj/SP/functions2.r")
Results <- list()
load("pvalues_multi_list_var_est_general_sigma1_dot1.RData")
Results <- c(Results, results.list)
load("pvalues_multi_list_var_est_general_sigma1_dot1_H1T_n1200.RData")
Results <- c(Results, results.list)
load("pvalues_multi_list_var_est_general_sigma1_dot1_H1T_n100.RData")
Results <- c(Results, results.list)
results.list <- Results; rm(Results)
(NAMES <- names(results.list))
alpha <- 0.05
ns <- H1s <- exprs <- reject.rates <- c()
for (k in 1:length(results.list)) {
  pvalues.multi <- results.list[[k]]
  name <- NAMES[k]
  splitstr <-  strsplit(name,split = '_')[[1]]
  n <- gsub('n', 'n=',splitstr[1])
  H1 <- splitstr[2]
  expr <- splitstr[3]
  reject.rate <- rej.rate(pvalues.multi)
  reject.rates <- c(reject.rates, reject.rate)
  nB <- length(reject.rate)
  ns <- c(ns, rep(n, nB))
  H1s <- c(H1s, rep(H1, nB))
  exprs <- c(exprs, rep(expr, nB))
  print(k)
}


ns <- factor(ns, levels=c("n=100","n=200","n=400","n=800","n=1200"))
result <- data.frame(reject.rates, ns, H1s, exprs)
result$B <- rep(1:nB, length(NAMES))

ind <- (result$exprs=="expr3")&(result$H1s=="H1TRUE")&(result$ns=="n=100")
result0 <- result[!ind,]
ind <- (result0$exprs=="expr2")&(result0$H1s=="H1TRUE")&(result0$ns=="n=100")
result0 <- result0[!ind,]
ind <- (result0$exprs=="expr1")&(result0$H1s=="H1TRUE")&(result0$ns=="n=1200")
result0 <- result0[!ind,]
library(ggplot2)
for (expr in paste0("expr",1:3)) {
  for (H1 in paste0("H1",c(T,F))) {
    ind <- (result0$exprs==expr)&(result0$H1s==H1)
    if (sum(ind) > 0) {
      result1 <- result0[ind, ]
      gp <- ggplot(data=result1, aes(x=B, y=reject.rates)) +
        geom_line()+ labs(x = "The number of splits")+labs(y="reject rate")+
        geom_point(size = 1)+facet_wrap(~ ns, ncol=1)
      print(gp)
      if ((H1=="H1TRUE") & (expr!="expr1")) file_name <- paste0(H1, expr, "sigma0dot1.eps") else
        file_name <- paste0(H1, expr, "sigma1.eps")
      ggsave(gp, file=file_name)
    }
  }
}


rej.rate <- function(pvalues.multi, gamma.min=0.01) {
  nRep = dim(pvalues.multi)[1]
  nB <- length(pvalues.multi[,1][[1]])
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=nRep, ncol=nB)
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
  # cat(gamma.min, reject.rate,"\n")
  # plot(reject.rate, main=names(results.list)[k])
  # plot(reject.rate, main=paste0("H1=", H1,", sigma=", sigma))
  return(reject.rate)
}

load("pvalues_multi_list_var_est_general_sigma1_dot1_H1T_n200_nB600.RData")
names(results.list)
alpha <- 0.05
for (k in 1:length(results.list)) {
  pvalues.multi <- results.list[[k]]
  nRep = dim(pvalues.multi)[1]
  nB <- length(pvalues.multi[,1][[1]])
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=nRep, ncol=nB)
  for (gamma.min in c(0.01)) {
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
    plot(reject.rate, main=names(results.list)[k])
    # plot(reject.rate, main=paste0("H1=", H1,", sigma=", sigma))
  }
}


