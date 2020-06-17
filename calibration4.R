compute_kernel_and_feature_matrix <- function(x, sigmas, type, return_Phi=FALSE) {
  # Computes the kernel matrix (and feature) matrix.
  # 
  # Args:
  #   x: The covariate matrix with sample size n by the number of covaraites p.
  #   sigmas: The vector of the tuning parameters in the kernels.
  #   type: The string of the kernel type. If "rbf", rbf kernels; if "rec", rectangular kernels.
  #   return_Phi: If "TRUE", calculate the feature matrix; if "FALSE", not calculate the feature matrix.
  # 
  # Returns:
  #   The kernel matrix (and feature) matrix.
  n <- dim(x)[1]
  if (type == "rec") {
    K <- diag(nrow=n, ncol=n)/2
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        K[i,j] <- exp(- max(abs(x[i, ] - x[j, ])/sigmas))
      }
    }
    K <- K + t(K)
  } else if (type == "rbf") {
    sigx <- x %*% diag(sqrt(sigmas))  # to multiply matrix columns with vector elements in 'sigmas'
    xx <- apply(sigx^2, 1, sum)
    kxx <- matrix(xx, ncol=1) %*% matrix(1, nrow=1, ncol=n)
    kxy <- sigx %*% t(sigx)
    K <- kxx + t(kxx) - 2*kxy
    K <- exp(-K)
  } else {
    stop("No such choice of kernel types")
  }
  if (return_Phi) {
    # Eigenvalue decomposition to get feature vector
    decom <- eigen(K, symmetric=TRUE)
    Phi <- sqrt(decom$values) * t(decom$vectors) 
    return(list(K=K, Phi=Phi))
  } 
  list(K=K)
}

# Rcplex
solve_mix_program <- function(a, y, lambda_1, lambda_0, d, K) {
  # Solve the mix integer programming using the R package Rcplex.
  # 
  # Args:
  #   a: The vector of the actions in the observations.
  #   y: The vector of outcomes in the sample.
  #   lambda_1: The coefficient in fron of the pernalty term ||w_1||^2.
  #   lambda_0: The coefficient in fron of the pernalty term ||w_0||^2.
  #   d: The vector of actions assigned by the decision list.
  #   K: The kernel matrix
  #
  # Returns:
  #   The coefficient vectors w_1, w_0.
  n <- length(a)
  yK <- y * K
  a_ind1 <- (a==1)
  a_ind0 <- (a==0)
  a_n1 <- sum(a)
  a_n0 <- sum(1-a)
  cvec <- -c(colSums(yK[a_ind1, a_ind1]), colSums(yK[a_ind0, a_ind0]))
  Amat1 <- (1-2*d) * K[, a_ind1]
  Amat2 <- (2*d-1) * K[, a_ind0]
  Amat <- cbind(Amat1, Amat2)
  bvec <- rep(0, n)
  phi_phi_transpose1 <- matrix(0, a_n1, a_n1)
  phi_phi_transpose0 <- matrix(0, a_n0, a_n0)
  for (i in 1:n) {
    if (a[i]==1) {
      phi_phi_transpose1 <- phi_phi_transpose1 + matrix(K[i, a_ind1], ncol=1) %*% matrix(K[i, a_ind1], nrow=1)
    } else {
      phi_phi_transpose0 <- phi_phi_transpose0 + matrix(K[i, a_ind0], ncol=1) %*% matrix(K[i, a_ind0], nrow=1)
    }
  }
  Qmat <- matrix(0, n, n)
  Qmat[1:a_n1, 1:a_n1] <- lambda_1 * diag(a_n1) + phi_phi_transpose1
  Qmat[(a_n1+1):n, (a_n1+1):n] <- lambda_0 * diag(a_n0) + phi_phi_transpose0
  res <- Rcplex(cvec, Amat, bvec, Qmat=Qmat, sense=c("L"), objsense="min", control=list(trace=0))
  list(w1=res$xopt[1:a_n1], w0=res$xopt[(a_n1+1):n])
}

cv_nfold <- function(parameters, other_argu) {
  # Campute the mean square error from n-fold cross-validation.
  # 
  # Args:
  #   parameters: The vector of tuning parameters sigmas and lambda.
  #   a: The vector of the actions in the observations.
  #   y: The vector of outcomes in the sample.
  #   w1: The cofficient vector corresponsing action 1.
  #   w0: The cofficient vector corresponsing action 0.
  #   Phi: The matrix with size n by n, where n is the sample size, where the jth column vector is feature vector \phi(x_j).
  #   lambda_1: The coefficient in fron of the pernalty term ||w_1||^2.
  #   lambda_0: The coefficient in fron of the pernalty term ||w_0||^2.
  #   gam: The coefficent in front of the last penalty term.
  #   d: The vector of actions assigned by the decision list.
  #   maxT: The maximum number of iterations in the subgradient descent.  #   type: The string of the kernel type. If "rbf", rbf kernels; if "rec", rectangular kernels.
  #
  #   parallel: If "TRUE", use parallel in the loocv; if "FALSE", not use.
  #   njob: The number of the cores used in the parallel.
  #
  # Returns:
  #   The mean square error from leave-one-out cross-validation.
  x <- other_argu$x
  a <- other_argu$a
  y <- other_argu$y
  d <- other_argu$d
  type <- other_argu$type
  nfold <- other_argu$nfold
  n <- length(a)
  p <- dim(x)[2]
  sigmas <- parameters[1:p]
  lambda_1 <- lambda_0 <- parameters[p+1]
  K <- compute_kernel_and_feature_matrix(x, sigmas, type)$K
  ind <- 1:n
  size_each_fold <- n/nfold
  square_error <- 0
  a_ind1 <- which(a==1)
  a_ind0 <- which(a==0)
  for (i in 1:nfold) {
    ind_each_fold <- sample(ind, size_each_fold)
    # w <- execute_subgradient_descent(a[-ind_each_fold], y[-ind_each_fold], w1, w0, Phi[,-ind_each_fold], 
    #                                  lambda_1, lambda_0, gam, d[-ind_each_fold], maxT) 
    w <- solve_mix_program(a[-ind_each_fold], y[-ind_each_fold], lambda_1, lambda_0, 
                           d[-ind_each_fold], K[-ind_each_fold, -ind_each_fold]) 
    w1 <- matrix(w$w1, ncol=1)
    w0 <- matrix(w$w0, ncol=1)
    # browser()
    fitted_test1 <- K[ind_each_fold, setdiff(a_ind1, ind_each_fold)] %*% w1
    fitted_test0 <- K[ind_each_fold, setdiff(a_ind0, ind_each_fold)] %*% w0
    fitted_test <- fitted_test1 
    fitted_test[a[ind_each_fold] == 0] <- fitted_test0[a[ind_each_fold] == 0]
    square_error_test <- mean((fitted_test - y[ind_each_fold])**2)
    ind <- setdiff(ind, ind_each_fold)
    square_error <- square_error + square_error_test
  }
  square_error/nfold
}

tune_parameters_find_coef_vectors <- function(x, a, y, d, type="rbf", tune=TRUE, nfold=10) {
  # Implement the Algorithm 2: tuning and finding the solutions to the constrained kernel ridge regression where 
  # the solutions are consistent with the decsion list.
  # 
  # Args:
  #   x: The covariate matrix with sample size n by the number of covaraites p.
  #   a: The vector of the actions in the observations.
  #   y: The vector of outcomes in the sample.
  #   d: The vector of actions assigned by the decision list.
  #   maxT: The maximum number of iterations in the subgradient descent.
  #   type: The string of the kernel type. If "rbf", rbf kernels; if "rec", rectangular kernels.
  #   eps: The tolenrance error used for the last penalty term.
  #   eps.gam: The step size factor multiplied to the base gamma.
  #   parallel: If "TRUE", use parallel in the loocv; if "FALSE", not use.
  #   njob: The number of the cores used in the parallel.
  #   tune: If "TRUE", execute the tuning procedure; if "FALSE", use fixed tuning parameters.
  #
  # Returns:
  #   The square error for the ith observation and the other data is used for training.
  p <- dim(x)[2]
  n <- dim(x)[1]
  other_argu <- list()
  other_argu$x <- x
  other_argu$a <- a
  other_argu$y <- y 
  other_argu$d <- d
  other_argu$type <- type
  other_argu$nfold <- nfold
  sigmas.upper <- 20/sum((x[1,]-x[2,])^2)
  if (tune) {
    # op <- spg(c(rep(sigmas.upper/10, p), 1/(mean(y)+sd(y)*3)/2/n), cv_nfold, other_argu=other_argu, lower=c(rep(0.1/p, p), 0.00001),
    #       upper=c(rep(sigmas.upper, p), 1000), quiet = T)
    
    op <- optim(c(rep(sigmas.upper/5, p), 1/(mean(y)+sd(y)*3)/2/n), cv_nfold, gr=NULL, method="L-BFGS-B",
                lower=c(rep(0.1/p, p), 0.00001), upper=c(rep(sigmas.upper, p), 1000), other_argu=other_argu)
    tune.para <- op$par
    sigmas <- tune.para[1:p]
    lambda_1 <- lambda_0 <- tune.para[p+1]
  } else {
    sigmas <- rep(0.1, p)
    lambda_1 <- lambda_0 <- 1
  }
  K <- compute_kernel_and_feature_matrix(x, sigmas, type)$K
  w <- solve_mix_program(a, y, lambda_1, lambda_0, d, K)
  w1 <- w$w1
  w0 <- w$w0
  a_ind1 <- (a==1)
  a_ind0 <- (a==0)
  y_hat <- rep(0, length(a))
  y_hat[a_ind1] <- K[a_ind1, a_ind1] %*% w1
  y_hat[a_ind0] <- K[a_ind0, a_ind0] %*% w0
  # To test whether constrains are satisfied
  Amat1 <- (1-2*d) * K[, a_ind1] %*% w1
  Amat2 <- (2*d-1) * K[, a_ind0] %*% w0
  condition.value <- mean(Amat1 + Amat2 <= 1e-6)  # should be 1 if satisfying constrains
  cat(", tuned paras:", c(sigmas,lambda_1),", constrains value:", condition.value, "\n")
  list(w1=w1, w0=w0, y_hat=y_hat)
}

find_alpha_tilde <- function(residual, y_hat, rs, nB, x, a, alpha_target, Switch, L, gamma.min) {
  # Compute the adaptive multiple-sample-splitting p-values from the new data sets constructed by boostrapping the residuals. 
  # And find the calibrated type I error level, alpha tilde, which satisfies the reject rates is closest alpha_target
  # 
  # Args:
  #   residual: The residuals got from the Algorithm 2.
  #   y_hat: The fitting function esitmated by the constrained Kernel Ridge Regression.
  #   rs: The vector of deciding the number of the bootstraps.
  #   nB: The number of the sample splits.
  #   x: The covariate matrix with sample size n by the number of covaraites p.
  #   a: The vector of the actions in the observations.
  #   alpha_target: The nominal type I error level alpha.
  #
  # Returns:
  #   The calibrated type I error level.
  pvalues.multi <- c()
  n <- length(y_hat)
  for (i in rs) {
    resi.boot <- sample(residual, size=n, replace=T)
    y <- y_hat + resi.boot
    pvalues.multi <- rbind(pvalues.multi, pValue.multiSplit.krr(x, a, y, nB=nB, Switch=Switch, delta_m=T, L=L))
  }
  opt <- optimize(compute_reject_rate, c(0,1), tol=0.00001, pvalues.multi=pvalues.multi, nB=nB, 
                  alpha_target=alpha_target, gamma.min=gamma.min)
  list(alpha.tilde=opt$minimum, diff_min=opt$objective)
}

compute_adative_pvalues <- function(pvalues.multi) {
  # Campute the adaptive p-values using multiple sample splits.
  # 
  # Args:
  #   pvalues.multi: The matrix containing the p-values for each split in each replicate, 
  #                  with size the number of replicates "n_rep" by the number of sample splits "nB".
  #   gamma.min: The value used for adative p-value finding.
  #
  # Returns:
  #   The matrix with size n_rep by nB containing the adaptive p-values for different 
  #   number of splits from 1 to nB in each replciate.
  n_rep <- dim(pvalues.multi)[1]  # the number of experiment replicates
  nB <- length(pvalues.multi[, 1][[1]])
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=n_rep, ncol=nB)
  gamma.mins <- unlist(pvalues.multi[, 4])
  p_adapt <- pvalues
  for (b in 2:nB) {
    p_values <- c()
    for (r in 1:n_rep) {
      pValues <- pvalues[r, 1:b]
      Q_inf <- optimize(Q.gamma, c(gamma.mins[r], 1), tol = 0.00001, pvalues = pValues)$objective
      p_value <- min(1, (1-log(gamma.mins[r]))*Q_inf)  # adaptive way to find gamma
      # p_value <- Q.gamma(gamma.mins[r], pValues)
      p_values <-c(p_values, p_value)
    }
    p_adapt[, b] <- p_values
  }
  p_adapt
}

compute_reject_rate <- function(alpha, pvalues.multi, nB, alpha_target, gamma.min) {
  # Compute the calibrated alpha using boostrap.
  # 
  # Args:
  #   alpha: The type I error level needed to calibrate.
  #   pvalues.multi: The matrix with p-values for each split in each replicate.
  #   nB: The number of the sample splits.
  #   alpha_target: The nominal type I error level alpha.
  #   gamma.min: The value used for adative p-value finding.
  #
  # Returns:
  #   The absolute difference between the calibrated reject rate and the 
  #   nominal type-I-error level alpha_target.
  n_rep <- dim(pvalues.multi)[1]
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=n_rep, ncol=nB)
  p_values <- c()
  if (nB == 1) {
    p_values <- pvalues
  } else {
    for (r in 1:n_rep) {
      pValues <- pvalues[r, ]
      # Adaptive way to find gamma:
      Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
      p_value <- min(1, (1-log(gamma.min))*Q_inf)
      # Fixed gamma to calculate the p-value:
      # p_value <- Q.gamma(gamma.min, pValues)
      p_values <- c(p_values, p_value)
    }
  }
  reject_rate <- mean(p_values < alpha)
  abs(reject_rate-alpha_target)
}

find_alpha_tilde2 <- function(residual, y_hat, rs, nB, x, a, alpha_target, Switch, L) {
  # Compute the adaptive multiple-sample-splitting p-values from the new data sets constructed by boostrapping the residuals. 
  # And find the calibrated type I error level, alpha tilde, which satisfies the reject rates is closest alpha_target
  # 
  # Args:
  #   residual: The residuals got from the Algorithm 2.
  #   y_hat: The fitting function esitmated by the constrained Kernel Ridge Regression.
  #   rs: The vector of deciding the number of the bootstraps.
  #   nB: The number of the sample splits.
  #   x: The covariate matrix with sample size n by the number of covaraites p.
  #   a: The vector of the actions in the observations.
  #   alpha_target: The nominal type I error level alpha.
  #
  # Returns:
  #   The calibrated type I error level.
  pvalues.multi <- c()
  n <- length(y_hat)
  for (i in rs) {
    resi.boot <- sample(residual, size=n, replace=T)
    y <- y_hat + resi.boot
    pvalues.multi <- rbind(pvalues.multi, pValue.multiSplit.krr(x, a, y, nB=nB, Switch=Switch, delta_m=T, L=L))
  }
  n_rep <- dim(pvalues.multi)[1]
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=n_rep, ncol=nB)
  diff_min <- Inf
  for (gamma.min in seq(0.02, 0.98, 0.02)) {
    p_values <- c()
    if (nB == 1) {
      p_values <- pvalues
    } else {
      for (r in 1:n_rep) {
        pValues <- pvalues[r, ]
        # Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
        # p_value <- min(1, (1-log(gamma.min))*Q_inf)  # adaptive way to find gamma
        p_value <- Q.gamma(gamma.min, pValues)
        p_values <- c(p_values, p_value)
      }
    }
    for (alpha in seq(0.01, 0.99, 0.005)) {
      reject_rate <- mean(p_values < alpha)
      diff_temp <- abs(reject_rate-alpha_target)
      if (diff_min > diff_temp) {
        diff_min <- diff_temp
        alpha.tilde <- alpha
        gamma.min_chosen <- gamma.min
        # print(list(gamma.min=gamma.min_chosen, alpha.tilde=alpha.tilde, diff_min=diff_min))
      }
    }
  }
  list(gamma.min=gamma.min_chosen, alpha.tilde=alpha.tilde, diff_min=diff_min)
  # opt <- optim(c(0.05,0.05), compute_reject_rate2, method = "L-BFGS-B", pvalues.multi=pvalues.multi, nB=nB, 
  #                 alpha_target=alpha_target)
  # alpha.tilde <- opt$minimum
  # p_values <- c()
  # if (nB == 1) {
  #   p_values <- pvalues
  # } else {
  #   for (r in 1:n_rep) {
  #     pValues <- pvalues[r, ]
  #     Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
  #     p_value <- min(1, (1-log(gamma.min))*Q_inf)  # adaptive way to find gamma
  #     # p_value <- Q.gamma(0.1, pValues)
  #     p_values <- c(p_values, p_value)
  #   }
  # }
  # alpha.tilde2 <- sort(p_values)[ceil(n_rep*alpha_target)]
  # list(alpha.tilde=alpha.tilde, alpha.tilde2=alpha.tilde2)
}

find_alpha_tilde3 <- function(residual, y_hat, rs, nB, x, a, alpha_target, Switch, L) {
  # Compute the adaptive multiple-sample-splitting p-values from the new data sets constructed by boostrapping the residuals. 
  # And find the calibrated type I error level, alpha tilde, which satisfies the reject rates is closest alpha_target
  # 
  # Args:
  #   residual: The residuals got from the Algorithm 2.
  #   y_hat: The fitting function esitmated by the constrained Kernel Ridge Regression.
  #   rs: The vector of deciding the number of the bootstraps.
  #   nB: The number of the sample splits.
  #   x: The covariate matrix with sample size n by the number of covaraites p.
  #   a: The vector of the actions in the observations.
  #   alpha_target: The nominal type I error level alpha.
  #
  # Returns:
  #   The calibrated type I error level.
  pvalues.multi <- c()
  n <- length(y_hat)
  for (i in rs) {
    resi.boot <- sample(residual, size=n, replace=T)
    y <- y_hat + resi.boot
    pvalues.multi <- rbind(pvalues.multi, pValue.multiSplit.krr(x, a, y, nB=nB, Switch=Switch, delta_m=T, L=L))
  }
  n_rep <- dim(pvalues.multi)[1]
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=n_rep, ncol=nB)
  diff_min <- Inf
  alpha <- alpha.tilde <- 0.05
  for (gamma.min in seq(0.01, 0.99, 0.01)) {
    p_values <- c()
    if (nB == 1) {
      p_values <- pvalues
    } else {
      for (r in 1:n_rep) {
        pValues <- pvalues[r, ]
        Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
        p_value <- min(1, (1-log(gamma.min))*Q_inf)  # adaptive way to find gamma
        # p_value <- Q.gamma(gamma.min, pValues)
        p_values <- c(p_values, p_value)
      }
    }
    reject_rate <- mean(p_values < alpha)
    diff_temp <- abs(reject_rate-alpha_target)
    if (diff_min > diff_temp) {
      diff_min <- diff_temp
      gamma.min_chosen <- gamma.min
      # print(list(gamma.min=gamma.min_chosen, alpha.tilde=alpha.tilde, diff_min=diff_min))
    }
  }
  list(gamma.min=gamma.min_chosen, alpha.tilde=alpha.tilde, diff_min=diff_min)
  # opt <- optim(c(0.05,0.05), compute_reject_rate2, method = "L-BFGS-B", pvalues.multi=pvalues.multi, nB=nB, 
  #                 alpha_target=alpha_target)
  # alpha.tilde <- opt$minimum
  # p_values <- c()
  # if (nB == 1) {
  #   p_values <- pvalues
  # } else {
  #   for (r in 1:n_rep) {
  #     pValues <- pvalues[r, ]
  #     Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
  #     p_value <- min(1, (1-log(gamma.min))*Q_inf)  # adaptive way to find gamma
  #     # p_value <- Q.gamma(0.1, pValues)
  #     p_values <- c(p_values, p_value)
  #   }
  # }
  # alpha.tilde2 <- sort(p_values)[ceil(n_rep*alpha_target)]
  # list(alpha.tilde=alpha.tilde, alpha.tilde2=alpha.tilde2)
}


compute_reject_rate2 <- function(alpha_gamma.min, pvalues.multi, nB, alpha_target) {
  # Compute the calibrated alpha using boostrap.
  # 
  # Args:
  #   alpha: The type I error level needed to calibrate.
  #   pvalues.multi: The matrix with p-values for each split in each replicate.
  #   nB: The number of the sample splits.
  #   alpha_target: The nominal type I error level alpha.
  #   gamma.min: The value used for adative p-value finding.
  #
  # Returns:
  #   The absolute difference between the calibrated reject rate and the 
  #   nominal type-I-error level alpha_target.
  alpha <- alpha_gamma.min[1]
  gamma.min <- alpha_gamma.min[2]
  n_rep = dim(pvalues.multi)[1]
  pvalues <- matrix(Reduce(rbind, pvalues.multi[,1]), nrow=n_rep, ncol=nB)
  p_values <- c()
  if (nB == 1) {
    p_values <- pvalues
  } else {
    for (r in 1:n_rep) {
      pValues <- pvalues[r, ]
      Q_inf <- optimize(Q.gamma, c(gamma.min, 1), tol = 0.00001, pvalues = pValues)$objective
      p_value <- min(1, (1-log(gamma.min))*Q_inf)  # adaptive way to find gamma
      # p_value <- Q.gamma(0.1, pValues)
      p_values <- c(p_values, p_value)
    }
  }
  reject_rate <- mean(p_values < alpha)
  abs(reject_rate-alpha_target)
}

