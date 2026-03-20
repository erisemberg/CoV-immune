buildZ <- function(flow_mat, X_mat){
  sex_int_terms <- flow_mat*X_mat[,"sex"]
  colnames(sex_int_terms) <- paste0(colnames(flow_mat), "_sex")
  SARS_int_terms <- flow_mat*X_mat[,"infectionSARSCoV"]
  colnames(SARS_int_terms) <- paste0(colnames(flow_mat), "_SARS")
  SARS2_int_terms <- flow_mat*X_mat[,"infectionSARS2CoV"]
  colnames(SARS2_int_terms) <- paste0(colnames(flow_mat), "_SARS2")
  flow_wInt <- cbind(flow_mat, sex_int_terms, SARS_int_terms, SARS2_int_terms)
  return(flow_wInt)
}

jaccard <- function(set1, set2) {
  intersection <- sum(set1 & set2)
  union <- sum(set1 | set2)
  if (union == 0) return(NA)  # avoid division by zero
  return(intersection / union)
}

extract_draws <- function(res){
  gammas <- do.call(rbind, lapply(res, function(x) x$gamma))
  betas <- do.call(rbind, lapply(res, function(x) x$beta))
  deltas <- do.call(rbind, lapply(res, function(x) x$delta))
  sigma2s <- as.numeric(do.call(rbind, lapply(res, function(x) x$sigma2)))
  flow_miss <- do.call(rbind, lapply(res, function(x) x$flow_test))
  return(list(gammas = gammas,
              betas = betas,
              deltas = deltas,
              sigma2s = sigma2s,
              flow_miss = flow_miss))
}

run_sampler_once <- function(y, X, flow, X_test, flow_test, 
                             a_psi = 2, b_psi = 7, tau2 = 2.5e-05, c0 = 300, 
                             a0 = 3, b0 = 0.2,S = 1000, burn_in = 100, n_chains = 4, 
                             seed = 123){
  set.seed(seed)
  omega2 <- tau2*c0
  q <- ncol(flow) # number of effects subject to selection 
  p <- ncol(X) # number of fixed effects
  n <- length(y)
  n_test <- nrow(flow_test)
  
  # imputation priors 
  # prior on Bimp ~ MN(M0, Lambda0, Sigma)
  M0 <- matrix(0, p, q) # initialize beta_imp and set M0 (prior mean)  
  Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
  Lambda0inv <- chol2inv(chol(Lambda0)) 
  # prior on Sigma ~ IW(S0, nu0), centers the prior at E(Sigma) = I_q 
  nu0 <- q + 2
  S0 <- diag(q) * (nu0 - q - 1)
  Sigma_init <- S0
  # initialize missing flow data to means 
  flow_full_init <- flow # leave original dataset unchanged 
  O <- 1*(!is.na(flow))
  # initialize all missing values in flow_full to phenotype average
  for (j in 1:q){
    flow_full_init[is.na(flow_full_init[,j]), j] <- colMeans(flow_full_init, na.rm = TRUE)[j]
  }
  
  # initialize missing flow data to means 
  flow_full_test_init <- flow_test # leave original dataset unchanged 
  O_test <- 1*(!is.na(flow_test))
  miss_ix_test <- which(O_test == 0)
  # initialize all missing values in flow_full to phenotype average
  for (j in 1:q){
    flow_full_test_init[is.na(flow_full_test_init[,j]), j] <- colMeans(flow_full_test_init, na.rm = TRUE)[j]
  }
  
  # regression priors 
  muB0 <- rep(0, p) # prior mean
  V <- diag(p) # identity matrix 
  phi2 <- var(y)*10 # prior variance on betas 
  b0 <- var(y)*(a0-1)
  sigma2_init <- 1 # initialize sigma  
  # initialize delta 
  delta_ols <- matrix(lm(y ~ flow)$coefficients[-1], nrow = q, ncol = 1) # initialize delta to OLS estimates
  int_deltas <- matrix(0, nrow = 3*q, 1)  
  delta_init <- rbind(delta_ols, int_deltas)  
  # warm start for gamma 
  marg_corr <- abs(cor(flow_full_init, y))
  top_k <- order(marg_corr, decreasing = TRUE)[1:10]
  gamma_init <- rep(0, (q*4))  
  gamma_init[top_k] <- 1
  
  res <- vector("list", n_chains)
  for (chain in 1:n_chains){
    flow_full <- flow_full_init
    flow_full_test <- flow_full_test_init
    Sigma <- Sigma_init
    sigma2 <- sigma2_init
    delta <- delta_init
    gamma <- gamma_init
    
    nsave <- (S - burn_in)
    flow_test_samples <- matrix(NA, nrow = nsave, ncol = length(miss_ix_test)) # missing values
    beta_samples <- matrix(NA, nrow = nsave, ncol = p)  # regression effects (fixed)
    gamma_samples <- matrix(NA, nrow = nsave, ncol = (q*4)) # inclusion indicator 
    delta_samples <- matrix(NA, nrow = nsave, ncol = (q*4)) # regression effects (shrinkage)
    sigma2_samples <- numeric(nsave) # residual variance 
    loglik <- matrix(NA, nrow = nsave, ncol = n)
    
    for (s in 1:S){
      # ----------------------update imputation params------------------------ #
      # update Beta: regression estimates for predicting immune traits
      Lambda_n <- chol2inv(chol(Lambda0inv + t(X) %*% X)) # p x p
      Mn <- Lambda_n %*% (Lambda0inv %*% M0 + t(X) %*% flow_full) # p x q 
      Bimp <- matrix(mvtnorm::rmvnorm(1, mean = c(Mn), sigma = kronecker(Sigma, Lambda_n)), nrow = p, ncol = q) # p x q
      
      # update Sigma: residual covariance matrix in likelihood for Z, and column covariance matrix in B prior
      # Sn <- S0 + Y'Y + M0' Lambda0^{-1} M0 - Mn' Lambda_n^{-1} Mn
      Sn <- S0 +
        t(flow_full) %*% flow_full +
        crossprod(backsolve(chol(Lambda0), M0, transpose = TRUE)) -
        crossprod(backsolve(chol(Lambda_n), Mn, transpose = TRUE))
      if (!isSymmetric(Sn)) { Sn <- (Sn + t(Sn))/2 } # symmetrize
      if (inherits(try(chol(Sn), silent = TRUE), "try-error")){ Sn <- as.matrix(Matrix::nearPD(Sn)$mat) } # make PD
      Sn_inv <- chol2inv(chol(Sn))
      if (!isSymmetric(Sn_inv)) { Sn_inv <- (Sn_inv + t(Sn_inv))/2 } # symmetrize
      if (inherits(try(chol(Sn_inv), silent = TRUE), "try-error")){ Sn_inv <- as.matrix(Matrix::nearPD(Sn_inv)$mat) } # make PD
      Omega <- rWishart(1, nu0 + n, Sn_inv)[,,1] # draw precision from Wishart
      if (inherits(try(chol(Omega), silent = TRUE), "try-error")){ Omega <- as.matrix(Matrix::nearPD(Omega)$mat) } # make PD
      Sigma <- chol2inv(chol(Omega)) # recover Sigma from precision using Cholesky
      
      # ------------------------------imputation------------------------------ #
      for (i in 1:n) {
        b <- (O[i, ] == 0) # missing/masked phenotypes
        if (!any(b)) next
        a <- !b # observed phenotypes
        mu_i <- X[i,] %*% Bimp # expected missing values including fixed effects
        Sa <- Sigma[a, a, drop = FALSE] # Nobs x Nobs
        ch <- chol(Sa)
        beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch) # Nmis x Nobs
        W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
        Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W # Nmis x Nmis
        if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){ Sigma_j <- make_spd(Sigma_j) }
        theta_j <- mu_i[b] + beta_j %*% t(flow_full[i, a, drop = FALSE] - mu_i[a])
        flow_full[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
      }
      
      # ------------------------imputation of test data------------------------ #
      for (i in 1:n_test) {
        b <- (O_test[i, ] == 0) # missing/masked phenotypes
        if (!any(b)) next
        a <- !b # observed phenotypes
        mu_i <- X_test[i,] %*% Bimp # expected missing values including fixed effects
        Sa <- Sigma[a, a, drop = FALSE] # Nobs x Nobs
        ch <- chol(Sa)
        beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch) # Nmis x Nobs
        W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
        Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W # Nmis x Nmis
        if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){ Sigma_j <- make_spd(Sigma_j) }
        theta_j <- mu_i[b] + beta_j %*% t(flow_full_test[i, a, drop = FALSE] - mu_i[a])
        flow_full_test[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
      }
      
      # -------------------------------regression------------------------------- #
      # create interaction terms from imputed data
      Z <- buildZ(flow_full, X)
      
      # update fixed effects, beta
      A <- (t(X) %*% X)/sigma2 + diag(p)/phi2
      B <- (t(X) %*% (y - Z %*% delta))/sigma2 + (solve(V) %*% muB0)/phi2 
      Ainv <- chol2inv(chol(A))
      beta <- t(rmvnorm(1, Ainv %*% B, Ainv))
      
      # update psi (probability that a variable belongs in the important state)
      if (s <= 50){ # fix psi in early iterations, giving model stable inclusion probability during early learning
        psi <- a_psi/(a_psi + b_psi) # prior mean of psi ~ Beta(a_psi, b_psi)
      } else {
        psi <- rbeta(1, a_psi + sum(gamma), b_psi + (q*4) - sum(gamma))
      }
      
      # update gamma using normal densities of the current delta_j under spike vs slab 
      fj_norm2 <- colSums(Z^2)
      for (j in 1:(4*q)){
        if (s <= 50 & j %in% top_k){ # warm start; keep gamma[j] = 1 for first 50 iterations
          gamma[j] <- 1
          next 
        }
        
        # calculate residuals with the current j removed (so we assess j's marginal contribution)
        fj <- Z[, j, drop = FALSE]
        rj <- y - X %*% beta - Z %*% delta + fj*delta[j] 
        
        # posterior variance/mean of delta_j under spike prior
        s2_spike <- 1 / (fj_norm2[j]/sigma2 + 1/tau2)
        m_spike <- s2_spike * as.numeric(crossprod(fj, rj)) / sigma2
        # ... and under the SLAB prior
        s2_slab <- 1 / (fj_norm2[j]/sigma2 + 1/omega2)
        m_slab <- s2_slab * as.numeric(crossprod(fj, rj)) / sigma2
        # conditional Bayes factor (integrating out delta_j)
        logBF_j <- 0.5 * ((log(s2_slab) - log(omega2)) -
                            (log(s2_spike) - log(tau2)) +
                            (m_slab^2)/s2_slab - (m_spike^2)/s2_spike)
        
        # posterior inclusion probability
        prob_imp <- plogis(log(psi) - log1p(-psi) + logBF_j)
        prob_imp <- max(min(prob_imp, 0.99), 0.01) # optional stabilizer
        gamma[j] <- rbinom(1, 1, p = prob_imp)
      }
      
      # update immune effects (effects subject to shrinkage), delta 
      D_gamma <- ifelse(gamma == 1, omega2, tau2) # var is omega2 if included, tau2 if not
      prior_prec <- 1/D_gamma 
      A <- (t(Z) %*% Z)/sigma2 + diag(prior_prec)  
      B <- (t(Z) %*% (y - X%*%beta))/sigma2
      Ainv <- chol2inv(chol(A))
      delta <- t(rmvnorm(1, Ainv %*% B, Ainv))
      
      # update residual variance, sigma2
      sigma2 <- 1/rgamma(1, 
                         a0 + n/2, 
                         b0 + 0.5*crossprod(y - X%*%beta - Z%*%delta)) # SSE/2
      
      # save log-likelihood if past burn in 
      if (s > burn_in) { 
        s_ix <- (s - burn_in)
        flow_test_samples[s_ix,] <- flow_full_test[miss_ix_test]
        beta_samples[s_ix,] <- beta
        gamma_samples[s_ix,] <- gamma
        delta_samples[s_ix,] <- delta
        sigma2_samples[s_ix] <- sigma2 
        mu_s <- as.numeric(X %*% beta + Z %*% delta)
        loglik[s_ix,] <- dnorm(y, mean = mu_s, sd = sqrt(sigma2), log = TRUE)
      }
    }
    
    # save results for chain 
    res[[chain]] <- list(
      flow_test = flow_test_samples,
      miss_ix_test = miss_ix_test,
      beta = beta_samples,
      gamma = gamma_samples,
      delta = delta_samples,
      sigma2 = sigma2_samples,
      loglik = loglik
    )
  }
  return(res)
}