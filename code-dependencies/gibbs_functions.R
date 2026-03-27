library(readr)
library(mvtnorm)

# -------------------------------misc functions------------------------------- #
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# load_cross_as_df: loads Rqtl cross object as a regular data.frame 
# Input: 
#     file_name: csv file, formatted for Rqtl 
#     n_geno_start: the column corresponding to the first genotype column in csv file 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
load_cross_as_df <- function(file_name, n_geno_start){
  cross_data <- read_csv(file_name, show_col_types = FALSE)
  cross_data <- cross_data[-c(1:2),] # remove rows with chr # and marker position 
  
  # convert categorical variables (infection, batch) to factors 
  cross_data$infection <- as.factor(cross_data$infection) 
  cross_data$batch <- as.factor(cross_data$batch)
  cross_data$flow_batch <- as.factor(cross_data$flow_batch)
  
  # code genotypes as c(0,1,2) - additive model 
  for (c in n_geno_start:ncol(cross_data)){
    cross_data[c] <- as.numeric(mapvalues(cross_data[[c]], from = c("AA", "AB", "BB"), to = c(0, 1, 2)))
  }
  
  return(cross_data)
}


# --------------------------------CV functions-------------------------------- #
# function for creating cell-wise folds over observed cells (i.e., exclude originally 
# missing cells where `O_base == 0`. Don't mask entire rows, because imputation model 
# (conditional multivariate normal) relies on at least some observed data per row. 
make_cellwise_folds <- function(O_base, K = 10, seed = 123, keep_one_per_row = TRUE) {
  stopifnot(is.matrix(O_base))
  set.seed(seed)
  n <- nrow(O_base)
  q <- ncol(O_base)
  
  # 0 = never masked; 1..K = fold id for masking
  folds <- matrix(0L, n, q)
  
  for (i in 1:n) {
    obs_cols <- which(O_base[i, ] == 1L)
    if (!length(obs_cols)) next
    
    # optionally reserve one "anchor" observed cell that is never masked
    if (keep_one_per_row && length(obs_cols) >= 1) {
      anchor <- sample(obs_cols, 1)
      folds[i, anchor] <- 0L
      obs_cols <- setdiff(obs_cols, anchor)
    }
    
    if (length(obs_cols) > 0) {
      base <- rep(1:K, length.out = length(obs_cols))
      folds[i, obs_cols] <- sample(base)  # permute across folds
    }
  }
  return(folds)
}

# function to run $K$-fold cross-validation
impute_cv <- function(Y, K = 10, seed = 123, 
                      RF = FALSE, X = NULL, geno_imp = NULL,
                      chain_fn_name = NULL, chains = 4, ncores = chains, S = 1000, 
                      mu0 = NULL, L0 = NULL, S0 = NULL, nu0 = NULL, M0 = NULL, 
                      Lambda0 = NULL, G = NULL, Sg0 = NULL, nug0 = NULL,
                      v_fixed = NULL, p0 = NULL, tau2 = NULL, lambda2 = NULL, 
                      Psi0 = NULL, nu_aux = NULL, m = NULL, xi_aux = NULL,
                      nu_tau0 = NULL, S_tau0 = NULL, batch = NULL){
  
  if (RF == FALSE){ # Bayesian imputation 
    chain_fn <- match.fun(chain_fn_name)
    
    # build export list dynamically based on non-NULL values 
    exported_vars <- c("make_spd", "seed", "Y", "O", "S")
    optional_vars <- list(X = X, mu0 = mu0, L0 = L0, S0 = S0, nu0 = nu0, B0 = B0, 
                          V0 = V0, G = G, Sg0 = Sg0, nug0 = nug0, v_fixed = v_fixed,
                          p0 = p0, tau2 = tau2, lambda2 = lambda2, Psi0 = Psi0,
                          nu_aux = nu_aux, m = m, xi_aux = xi_aux,
                          nu_tau0 = nu_tau0, S_tau0 = S_tau0, batch = batch)
    for (var_name in names(optional_vars)){
      if (!is.null(optional_vars[[var_name]])){
        exported_vars <- c(exported_vars, var_name)
      }
    } 
  }
  
  set.seed(seed)
  n <- nrow(Y)
  q <- ncol(Y)
  
  # Y may already contain missing data. Only cells with O_base == 1 are eligible for masking. 
  # we evaluate predictive performance only on the cells masked by CV - not on those originally 
  # missing (not that that would be possible anyway bc no observed value to compare to). 
  O_base <- matrix(1L, n, q) # O indicates whether data is observed or missing
  og_missing <- which(is.na(Y)) 
  O_base[is.na(Y)] <- 0L 
  folds <- make_cellwise_folds(O_base, K = K, seed = seed, keep_one_per_row = FALSE) # define folds
  
  fold_metrics <- data.frame(fold = 1:K, 
                             n_holdout = NA_integer_,
                             RMSE = NA_real_,
                             R2 = NA_real_)
  
  # collect all holdout errors for overall metrics
  all_true <- numeric(0)
  all_pred <- numeric(0)
  
  # collect all predicted values for originally missing data 
  #og_miss_pred <- matrix(NA, nrow = K, ncol = length(og_missing))
  pred_data <- matrix(NA, nrow = K, ncol = length(Y))
  
  for (k in 1:K) {
    O <- O_base
    to_mask <- (folds == k) # subset of currently observed cells
    O[to_mask] <- 0L # add masking on top of existing missingness
    
    # safety check: no row should be fully unobserved 
    if (any(rowSums(O) == 0)) {
      stop(sprintf("Fold %d would fully mask a row. Try keep_one_per_row=TRUE or reduce K.", k))
    }
    set.seed(seed + k*1000)
    
    if (RF == TRUE){
      Yk <- Y
      Yk[O==0] <- NA # create masked df
      # all data to pass into ranger 
      datk <- as.data.frame(cbind(Yk, X, geno_imp)) # passing in imputed genotype data makes this much faster 
      imp <- missRanger(data = datk, pmm.k = 0, num.trees = 500, maxiter = 5, verbose = 0)
      Y_hat <- as.matrix(imp[,1:q]) # imputed flow data 
    } else { # Bayesian imputation 
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      chain_res <- foreach(ch = 1:chains, 
                           .export = exported_vars, 
                           .packages = c("mvtnorm")) %dopar% chain_fn(ch)
      stopCluster(cl)
      
      # aggregate posterior means over S and chains
      miss_ix <- chain_res[[1]]$miss_ix
      # average within chain, then across chains
      chain_means <- lapply(chain_res, function(cr) colMeans(cr$flow_miss))
      mean_vec <- Reduce("+", chain_means) / length(chain_means)
      Y_hat <- Y # predicted matrix
      Y_hat[miss_ix] <- mean_vec # missing/masked filled by posterior mean 
    }
    
    # calculate RMSE and R-squared only on held-out cells
    y_true <- Y[to_mask]
    y_pred <- Y_hat[to_mask]
    all_true <- c(all_true, y_true)
    all_pred <- c(all_pred, y_pred)
    
    # metrics
    sse <- sum((y_true - y_pred)^2)
    rmse <- sqrt(mean((y_true - y_pred)^2))
    sst <- sum((y_true - mean(y_true))^2)
    r2 <- if (sst > 0) 1 - sse / sst else NA_real_
    
    fold_metrics$n_holdout[k] <- length(y_true)
    fold_metrics$RMSE[k] <- rmse
    fold_metrics$R2[k] <- r2
    
    message(sprintf("Fold %d/%d: RMSE=%.4f, R2=%.4f", k, K, rmse, r2))
    
    # fill pred_data with predicted data (originally missing + masked)
    pred_data[k, og_missing] <- Y_hat[og_missing]
    pred_data[k, to_mask] <- Y_hat[to_mask] 
  }
  
  # overall metrics across all folds
  overall_rmse <- sqrt(mean((all_true - all_pred)^2))
  overall_r2 <- {
    sse <- sum((all_true - all_pred)^2)
    sst <- sum((all_true - mean(all_true))^2)
    if (sst > 0) 1 - sse / sst else NA_real_
  }
  
  return(list(per_fold = fold_metrics,
              overall = list(RMSE = overall_rmse, R2 = overall_r2),
              pred_data = pred_data)) 
}


# ------------------------------Gibbs functions------------------------------- #
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

# a tiny nudge used as a practical guard against floating point noise; shouldn't be
# the mechanism that makes the sampler run - if it frequently ramps up noise, chain
# will be biased
make_spd <- function(M, eps0 = 1e-12, max_eps = 1e-2) {
  # Symmetrize + add diagonal jitter until Cholesky succeeds; return SPD matrix
  k <- nrow(M)
  M <- (M + t(M)) / 2
  eps <- 0
  repeat {
    test <- try(chol(M + diag(eps, k)), silent = TRUE)
    if (!inherits(test, "try-error")) return(M + diag(eps, k))
    eps <- if (eps == 0) eps0 else eps * 10
    if (eps > max_eps) stop("make_spd: could not make matrix SPD")
  }
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
  # initialize all missing values in Yhat to phenotype average
  for (j in 1:q){
    flow_full_init[is.na(flow_full_init[,j]), j] <- colMeans(flow_full_init, na.rm = TRUE)[j]
  }
  
  # initialize missing flow data to means 
  flow_full_test_init <- flow_test # leave original dataset unchanged 
  O_test <- 1*(!is.na(flow_test))
  miss_ix_test <- which(O_test == 0)
  # initialize all missing values in Yhat to phenotype average
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


# -------------------------Gibbs imputation functions------------------------- #
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# run_chain0: function to run one chain for imputation 
# based only on Sigma (variance-covariance matrix between 
# phenotypes) 
# Input: 
#     chain_id: chain identifier
#+++++++++++++++++++++++++++++++++++++++++++++++++++
run_chain0 <- function(chain_id){
  
  set.seed(seed + chain_id)
  n <- nrow(Y)
  q <- ncol(Y)
  
  # initialize variables 
  Yhat <- Y
  for (j in 1:q){
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j] # missing values with column means 
  }
  Lambda <- solve(S0) # precision 
  Lambda0 <- chol2inv(chol(L0)) 
  
  # matrices to store samples
  nsave <- S - burn_in
  theta_samples <- matrix(NA, nsave, q) # immune phenotype means
  Sigma_samples <- matrix(NA, nsave, q*q) # immune phenotype variance-covariance matrix
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = nsave, ncol = length(miss_ix))
  
  for (s in 1:S){
    # update theta (predictor means) 
    Ln_inv <- Lambda0 + n * Lambda
    Ln_inv <- (Ln_inv + t(Ln_inv))/2 # symmetrize
    if (inherits(try(chol(Ln_inv), silent = TRUE), "try-error")){ Ln_inv <- as.matrix(Matrix::nearPD(Ln_inv)$mat) }
    Ln <- chol2inv(chol(Ln_inv)) 
    Ybar <- apply(Yhat, 2, mean) # get current phenotype averages 
    mu_n <- Ln %*% (Lambda0 %*% mu0 + n * Lambda %*% Ybar)
    theta <- rmvnorm(1, mu_n, Ln)  
    
    # update Sigma (variance-covariance matrix for predictors)
    Sn <- S0 + (t(Yhat) - c(theta)) %*% t(t(Yhat) - c(theta)) 
    Sn <- (Sn + t(Sn))/2 # symmetrize  
    if (inherits(try(chol(Sn), silent = TRUE), "try-error")){ Sn <- as.matrix(Matrix::nearPD(Sn)$mat) }
    Sn_inv <- chol2inv(chol(Sn))
    Sn_inv <- (Sn_inv + t(Sn_inv))/2 # symmetrize 
    if (inherits(try(chol(Sn_inv), silent = TRUE), "try-error")){ Sn_inv <- as.matrix(Matrix::nearPD(Sn_inv)$mat) }
    Lambda <- rWishart(1, nu0 + n, Sn_inv)[,,1] # draw precision from Wishart 
    if (inherits(try(chol(Lambda), silent = TRUE), "try-error")){ Lambda <- as.matrix(Matrix::nearPD(Lambda)$mat) }
    Sigma <- chol2inv(chol(Lambda)) # recover Sigma from precision using Cholesky
    
    # imputation 
    for (i in 1:n) {
      b <- (O[i, ] == 0) # missing/masked phenotypes 
      if (!any(b)) next
      a <- !b # observed phenotypes 
      Sa <- Sigma[a, a, drop = FALSE] 
      ch <- chol(Sa)
      beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch)
      W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
      Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W
      theta_j <- theta[b] + beta_j %*% t(Yhat[i, a, drop = FALSE] - theta[a])
      if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){ Sigma_j <- make_spd(Sigma_j) }
      Yhat[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    
    if (s > burn_in){ # save draws
      s_ix <- (s - burn_in)
      theta_samples[s,] <- theta 
      Sigma_samples[s_ix,] <- c(Sigma)  
      flow_miss_samples[s_ix,] <- Yhat[miss_ix] 
    }
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(theta = theta_samples,
                             Sigma = Sigma_samples)))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# run_chain1: function to run one chain for imputation 
# based on Sigma (variance-covariance matrix between 
# phenotypes) and fixed effects for sex and infection
# Input: 
#     chain_id: chain identifier
#+++++++++++++++++++++++++++++++++++++++++++++++++++
run_chain1 <- function(chain_id){
  
  set.seed(seed + chain_id)
  n <- nrow(Y) # number of individuals 
  q <- ncol(Y) # number of outcomes (variables with missing data to impute)
  p <- ncol(X) # number of covariates in model predicting imputed values 
  
  # initialize variables
  Yhat <- Y
  for (j in 1:q){
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j] # missing values with column means
  }
  Sigma <- S0
  Lambda0inv <- chol2inv(chol(Lambda0)) 
  
  # matrices to store samples 
  nsave <- S - burn_in
  Bimp_samples <- array(NA, dim = c(nsave, p, q)) # fixed effects in imputation model 
  Sigma_samples <- matrix(NA, nsave, q*q) # immune phenotype variance-covariance matrix 
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = nsave, ncol = length(miss_ix))
  
  for (s in 1:S){
    # update Beta (conjugate matrix-normal update)
    Lambda_n <- chol2inv(chol(Lambda0inv + t(X) %*% X)) # p x p
    Mn <- Lambda_n %*% (Lambda0inv %*% M0 + t(X) %*% Yhat) # p x q
    Bimp <- matrix(mvtnorm::rmvnorm(1, mean = c(Mn), sigma = kronecker(Sigma, Lambda_n)), nrow = p, ncol = q)
    
    # update Sigma (conjugate normal-Wishart update)
    # Sn <- S0 + Y'Y + M0' Lambda0^{-1} M0 - Mn' Lambda_n^{-1} Mn 
    Sn <- S0 + 
      t(Yhat) %*% Yhat + 
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
    
    # imputation 
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
      theta_j <- mu_i[b] + beta_j %*% t(Yhat[i, a, drop = FALSE] - mu_i[a]) 
      Yhat[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    
    if (s > burn_in){ # save draws
      s_ix <- (s - burn_in)
      Bimp_samples[s_ix,,] <- Bimp
      Sigma_samples[s_ix,] <- c(Sigma)  
      flow_miss_samples[s_ix,] <- Yhat[miss_ix] 
    }
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(Beta_imp = Bimp_samples,
                             Sigma = Sigma_samples)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# run_chain2: function to run one chain for imputation 
# based on Sigma (variance-covariance matrix between 
# phenotypes), fixed effects for sex and infection, and
# a random polygenic effect
# Input: 
#     chain_id: chain identifier
#+++++++++++++++++++++++++++++++++++++++++++++++++++
run_chain2 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  
  Yhat <- Y
  for (j in 1:q){ # initialize missing values with column means
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j]
  }
  Z <- diag(n) # individual-level random effects 
  
  Sigma_samples <- matrix(NA, S, q*q) # immune phenotype variance-covariance matrix 
  Bimp_samples <- array(NA, dim = c(S, p, q)) # fixed effects in imputation model 
  Sigma_g_samples <- matrix(NA, S, q*q) # genetic variance-covariance matrix for immune phenotypes 
  U_samples <- array(NA, dim = c(S, n, q)) # each U matrix is n x q 
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = S, ncol = length(miss_ix))
  
  # initialize parameters
  # also initialize theta <- colMeans(Yhat) if including intercept 
  Sigma <- S0
  Sigma_g <- Sg0 # q x q genetic variance-covariance matrix 
  U <- matrix(0, n, q) # n x q matrix of random effects 
  
  Gsym <- (G + t(G))/2 # symmetrize G 
  eg <- eigen(Gsym, symmetric = TRUE) # compute eigenvectors of grm  
  Q <- eg$vectors # assign eigenvectors to Q 
  lam <- pmax(eg$values, 1e-8) # lam is eigenvalues, unless they are smaller than 1e-8 (to avoid issues with inverting)
  invlam <- 1/lam # invert 
  
  V0inv <- chol2inv(chol(V0)) 
  Ginv <- chol2inv(chol(G)) 
  
  for (s in 1:S){
    # update Beta 
    Vn <- chol2inv(chol(V0inv + t(X) %*% X)) # p x p
    flow_residR <- Yhat - Z%*%U # residuals after random effects 
    Bn <- Vn %*% (V0inv %*% B0 + t(X) %*% flow_residR) # p x q  
    Bimp <- t(matrix(mvtnorm::rmvnorm(1, c(t(Bn)), kronecker(Sigma, Vn)), q, p)) # transpose to get p x q
    Bimp_samples[s,,] <- Bimp
    
    # update Sigma 
    # Sn <- S0 + (Y-U)'(Y-U) + B0' V0^{-1} B0 - Bn' Vn^{-1} Bn 
    Sn <- S0 + 
      t(flow_residR) %*% flow_residR + 
      crossprod(backsolve(chol(V0), B0, transpose = TRUE)) - 
      crossprod(backsolve(chol(Vn), Bn, transpose = TRUE))
    if (!isSymmetric(Sn)) { Sn <- (Sn + t(Sn))/2 } # symmetrize  
    if (inherits(try(chol(Sn), silent = TRUE), "try-error")){ Sn <- as.matrix(Matrix::nearPD(Sn)$mat) } # make PD
    Sn_inv <- chol2inv(chol(Sn))
    if (!isSymmetric(Sn_inv)) { Sn_inv <- (Sn_inv + t(Sn_inv))/2 } # symmetrize 
    if (inherits(try(chol(Sn_inv), silent = TRUE), "try-error")){ Sn_inv <- as.matrix(Matrix::nearPD(Sn_inv)$mat) } # make PD
    Lambda <- rWishart(1, nu0 + n, Sn_inv)[,,1] # draw precision from Wishart 
    if (inherits(try(chol(Lambda), silent = TRUE), "try-error")){ Lambda <- as.matrix(Matrix::nearPD(Lambda)$mat) } # make PD
    Sigma <- chol2inv(chol(Lambda)) # recover Sigma from precision using Cholesky
    Sigma_samples[s,] <- c(Sigma) 
    
    # update U 
    # vec(U) ~ N(0, kronecker(Sigma_g, G)) 
    ### CONFIRM 
    flow_residF <- Yhat - X %*% Bimp # residuals after fixed effects
    flow_residF_star <- t(Q) %*% flow_residF # rotate rows
    Sigma_inv <- chol2inv(chol(Sigma))
    Sigma_g_inv <- chol2inv(chol(Sigma_g))
    U_star <- matrix(0, n, q)
    for (i in 1:n){ # update per-row 
      Sigma_i <- chol2inv(chol(Sigma_g_inv * invlam[i] + Sigma_inv)) # q x q 
      mu_i <- Sigma_i %*% (Sigma_inv %*% flow_residF_star[i,])
      U_star[i,] <- as.numeric(mvtnorm::rmvnorm(1, mean = mu_i, sigma = Sigma_i))
    }
    U <- Q %*% U_star
    U_samples[s,,] <- U
    
    # update Sigma_g
    # Su_n <- S0_u + t(U) %*% Ginv %*% U 
    U_star_scaled <- U_star * invlam 
    Sg_n <- Sg0  + t(U_star) %*% U_star_scaled
    if (!isSymmetric(Sg_n)) { Sg_n <- (Sg_n + t(Sg_n))/2 }
    if (inherits(try(chol(Sg_n), silent = TRUE), "try-error")) { Sg_n <- as.matrix(Matrix::nearPD(Sg_n)$mat) }
    Sg_n_inv <- chol2inv(chol(Sg_n))
    Lambda_g <- rWishart(1, nug0 + n, Sg_n_inv)[,,1]
    if (inherits(try(chol(Lambda_g), silent = TRUE), "try-error")) { Lambda_g <- as.matrix(Matrix::nearPD(Lambda_g)$mat) }
    Sigma_g <- chol2inv(chol(Lambda_g))
    Sigma_g_samples[s,] <- c(Sigma_g)
    
    # imputation 
    for (i in 1:n) {
      b <- (O[i, ] == 0) # missing/masked phenotypes 
      if (!any(b)) next
      a <- !b # observed phenotypes 
      # incl theta if intercept separate: theta + as.vector(t(Bimp) %*% X[i,])
      mu_i <- X[i,] %*% Bimp + U[i,] # expected missing values including fixed effects
      Sa <- Sigma[a, a, drop = FALSE] # Nobs x Nobs
      ch <- chol(Sa)
      beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch) # Nmis x Nobs
      W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
      Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W # Nmis x Nmis
      if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){ Sigma_j <- make_spd(Sigma_j) }
      theta_j <- mu_i[b] + beta_j %*% t(Yhat[i, a, drop = FALSE] - mu_i[a]) 
      Yhat[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    flow_miss_samples[s,] <- Yhat[miss_ix] # save results
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(Beta_imp = Bimp_samples,
                             Sigma = Sigma_samples, 
                             Sigma_g = Sigma_g_samples,
                             U = U_samples)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# run_chain3: function to run one chain for imputation 
# based on Sigma (variance-covariance matrix between 
# phenotypes), fixed effects for sex and infection, and
# a random batch effect
# Input: 
#     chain_id: chain identifier
#+++++++++++++++++++++++++++++++++++++++++++++++++++
run_chain3 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  G <- max(batch)
  
  Yhat <- Y
  # initialize missing values with column means 
  for (j in 1:q){
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j]
  }
  
  Sigma_samples <- matrix(NA, S, q*q) # immune phenotype variance-covariance matrix 
  Bimp_samples <- array(NA, dim = c(S, p, q)) # p x q fixed effects in imputation model 
  Tau_samples <- matrix(NA, S, q*q) # batch effect covariance matrix 
  U_samples <- array(NA, dim = c(S, G, q)) # G x q random batch effects
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = S, ncol = length(miss_ix))
  
  # priors / initial values 
  Sigma <- S0
  Tau <- S_tau0 / (nu_tau0 - q - 1)
  Tau_prec <- chol2inv(chol(Tau))
  U <- matrix(0, G, q)
  
  V0inv <- chol2inv(chol(V0)) 
  
  # pre-compute index lists ### CONFIRM 
  idx_by_batch <- split(seq_len(n), batch)
  n_by_batch <- vapply(idx_by_batch, length, 1L)
  
  for (s in 1:S){
    # update Beta (conjugate matrix-normal update)
    flow_resid <- Yhat - U[batch,,drop = FALSE] # n x q ### CONFIRM U[batch,]
    Vn <- chol2inv(chol(V0inv + t(X) %*% X)) # p x p
    Bn <- Vn %*% (V0inv %*% B0 + t(X) %*% flow_resid) # p x q ### CONFIRM: replaced Yhat w flow_resid  
    Bimp <- t(matrix(mvtnorm::rmvnorm(1, c(t(Bn)), kronecker(Sigma, Vn)), q, p)) # transpose to get p x q
    Bimp_samples[s,,] <- Bimp
    
    # update Sigma (conjugate normal-Wishart update)
    flow_resid <- Yhat - X %*% Bimp - U[batch, , drop = FALSE] # n x q
    Sn <- S0 + 
      t(flow_resid) %*% flow_resid + 
      crossprod(backsolve(chol(V0), B0, transpose = TRUE)) - 
      crossprod(backsolve(chol(Vn), Bn, transpose = TRUE))
    if (!isSymmetric(Sn)) { Sn <- (Sn + t(Sn))/2 }   
    if (inherits(try(chol(Sn), silent = TRUE), "try-error")){ Sn <- as.matrix(Matrix::nearPD(Sn)$mat) } 
    Sn_inv <- chol2inv(chol(Sn))
    if (!isSymmetric(Sn_inv)) { Sn_inv <- (Sn_inv + t(Sn_inv))/2 } 
    if (inherits(try(chol(Sn_inv), silent = TRUE), "try-error")){ Sn_inv <- as.matrix(Matrix::nearPD(Sn_inv)$mat) } 
    Lambda <- rWishart(1, nu0 + n, Sn_inv)[,,1] # draw precision from Wishart 
    if (!isSymmetric(Lambda)) { Lambda <- (Lambda + t(Lambda))/2 }
    if (inherits(try(chol(Lambda), silent = TRUE), "try-error")){ Lambda <- as.matrix(Matrix::nearPD(Lambda)$mat) } 
    Sigma <- chol2inv(chol(Lambda)) # recover Sigma from precision using Cholesky
    Sigma_samples[s,] <- c(Sigma) 
    
    # update Ug for g = 1,.., G
    for (g in seq_len(G)){
      ids <- idx_by_batch[[g]]
      ng <- n_by_batch[g]
      # sufficient statistic: sum of residuals without Ug
      Rg <- Yhat[ids,,drop = FALSE] - X[ids,,drop = FALSE] %*% Bimp 
      rg_sum <- colSums(Rg)
      mat <- Tau_prec + ng*Lambda
      if (!isSymmetric(mat)) { mat <- (mat + t(mat))/2 }
      if (inherits(try(chol(mat), silent = TRUE), "try-error")){ mat <- as.matrix(Matrix::nearPD(mat)$mat) } 
      Vg <- chol2inv(chol(mat)) # q x q 
      mg <- Vg %*% (Lambda %*% rg_sum)
      Ug <- mvtnorm::rmvnorm(1, mean = c(mg), sigma = Vg)
      U[g,] <- Ug
    }
    U_samples[s,,] <- U
    
    # update Tau (batch-effect covariance) 
    S_tau <- S_tau0 + crossprod(U)
    if (!isSymmetric(S_tau)) { S_tau <- (S_tau + t(S_tau))/2 }
    if (inherits(try(chol(S_tau), silent = TRUE), "try-error")){ S_tau <- as.matrix(Matrix::nearPD(S_tau)$mat) }
    S_tau_inv <- chol2inv(chol(S_tau))
    if (!isSymmetric(S_tau_inv)) { S_tau_inv <- (S_tau_inv + t(S_tau_inv))/2 }
    if (inherits(try(chol(S_tau_inv), silent = TRUE), "try-error")){ S_tau_inv <- as.matrix(Matrix::nearPD(S_tau_inv)$mat) }
    Tau_prec <- rWishart(1, nu_tau0 + G, S_tau_inv)[,,1]
    if (!isSymmetric(Tau_prec)) { Tau_prec <- (Tau_prec + t(Tau_prec))/2 }
    if (inherits(try(chol(Tau_prec), silent = TRUE), "try-error")){ Tau_prec <- as.matrix(Matrix::nearPD(Tau_prec)$mat) }
    Tau <- chol2inv(chol(Tau_prec))
    Tau_samples[s,] <- c(Tau)
    
    # imputation 
    for (i in 1:n) {
      b <- (O[i, ] == 0) # missing/masked phenotypes 
      if (!any(b)) next
      a <- !b # observed phenotypes 
      # incl theta if intercept separate: theta + as.vector(t(Bimp) %*% X[i,])
      mu_i <- X[i,] %*% Bimp + U[batch[i], ] # expected missing values including fixed effects
      Sa <- Sigma[a, a, drop = FALSE] # Nobs x Nobs
      ch <- chol(Sa)
      beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch) # Nmis x Nobs
      W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
      Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W # Nmis x Nmis
      if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){ Sigma_j <- make_spd(Sigma_j) }
      theta_j <- mu_i[b] + beta_j %*% t(Yhat[i, a, drop = FALSE] - mu_i[a]) 
      Yhat[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    flow_miss_samples[s,] <- Yhat[miss_ix] # save results
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(Beta_imp = Bimp_samples,
                             Sigma = Sigma_samples,
                             U = U_samples,
                             Tau = Tau_samples)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# run_chain3: function to run one chain for imputation 
# based on Sigma (variance-covariance matrix between 
# phenotypes), fixed effects for sex and infection, and
# horseshoe-shrunk effects for all genotypes 
# Input: 
#     chain_id: chain identifier
### CODE NOT WORKING YET 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
run_chain4 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  Yhat <- Y
  for (j in 1:q){ # initialize missing values with column means
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j]
  }
  
  #Sigma_samples <- matrix(NA, S, q*q) # immune phenotype variance-covariance matrix 
  #Sigma_mean <- matrix(0, q, q)
  #Bimp_samples <- array(NA, dim = c(S, p, q)) # fixed effects in imputation model 
  #Bimp_mean <- matrix(0, p, q)
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = (S/10), ncol = length(miss_ix))
  #flow_miss_mean <- numeric(length(miss_ix))
  
  # inverse-gamma sampler 
  #rinv_gamma <- function(shape, rate) 1 / rgamma(length(shape), shape = shape, rate = rate)
  rinv_gamma <- function(shape, rate, eps = 1e-12, cap = 1e12) {
    rate <- pmax(rate, eps)
    n <- length(rate)
    out <- 1 / rgamma(n, shape = shape, rate = rate)
    out[!is.finite(out)] <- cap
    pmin(out, cap)
  }
  
  for (s in 1:S){
    print(s)
    # build V (row variances)
    v_min <- 1e-8; v_max <- 1e8
    v <- c(rep(v_fixed, p0), tau2*lambda2)
    v <- pmin(pmax(v, v_min), v_max)
    
    # compute B without pxp inverse (Woodbury)
    # Bn = B0 + V X' A^{-1} (Y - X B0),  A = I + X V X'
    R  <- Yhat - X %*% B0                # n x q
    Xv <- sweep(X, 2, v, `*`)                 # n x p  (X * v_j by columns)
    A  <- diag(n) + Xv %*% t(X)               # n x n
    # solve A %*% Z = R  via Cholesky
    if (inherits(try(chol(A), silent = TRUE), "try-error")){ A <- as.matrix(Matrix::nearPD(A)$mat) }
    cholA <- chol(A)
    Z  <- backsolve(cholA, forwardsolve(t(cholA), R, upper.tri = TRUE, transpose = TRUE))
    Bn <- B0 + t(X) %*% Z                     # p x q, because X' A^{-1} R, then multiply each row by v:
    Bn <- sweep(Bn, 1, v, `*`)                   # (V X' A^{-1} R) + B0
    
    # update Sigma | Y, V  ~ IW(Psi_n, nu0+n) 
    # Use (★): Psi_n = Psi0 + (Y - X Bn)'(Y - X Bn) + (Bn - B0)' Λ0 (Bn - B0).
    # Λ0 = diag(1/v). Avoid forming Λ0 explicitly.
    E <- Yhat - X %*% Bn                    # n x q
    term1 <- crossprod(E)                        # q x q
    Db  <- Bn - B0                               # p x q
    # (Bn - B0)' Λ0 (Bn - B0) = sum_j (1/v_j) * (Db_j·' Db_j·)
    # compute as weighted row crossproducts:
    term2 <- t(Db) %*% (Db / v)                  # q x q  (divides each row j by v_j before crossprod)
    Psi_n <- Psi0 + term1 + term2
    # symmetrize / PD guard
    if (!isSymmetric(Psi_n)) Psi_n <- 0.5*(Psi_n + t(Psi_n))
    if (inherits(try(chol(Psi_n), silent=TRUE), "try-error")) { Psi_n <- as.matrix(Matrix::nearPD(Psi_n)$mat) }
    # draw Sigma
    Sigma <- chol2inv(chol(rWishart(1, nu0 + n, chol2inv(chol(Psi_n)))[,,1]))
    #Sigma_samples[s,] <- c(Sigma)
    #Sigma_mean <- Sigma_mean + Sigma 
    LS <- chol(Sigma)                            # q x q (upper)
    
    # update B | Sigma, Y, V  ~ MN(Bn, Vn, Sigma) 
    # Sample using low-rank representation: Δ = Zp - V X' A^{-1}(X Zp - Un),
    # with Zp ~ MN(0, V, Σ), Un ~ MN(0, I_n, Σ). Then B = Bn + Δ.
    # Draw Zp:
    Zstd_p <- matrix(rnorm(p * q), p, q)
    Zp <- sweep(Zstd_p %*% t(LS), 1, sqrt(v), `*`)  # diag(sqrt(v)) * N * LS'
    # Draw Un:
    Ustd_n <- matrix(rnorm(n * q), n, q)
    Un <- Ustd_n %*% t(LS)
    # Compute A^{-1}(X Zp - Un) using same chol(A)
    XZp <- X %*% Zp
    W   <- XZp - Un
    Slin <- backsolve(cholA, forwardsolve(t(cholA), W, upper.tri = TRUE, transpose = TRUE))
    # Δ = Zp - V X' Slin
    Delta <- Zp - sweep(t(X) %*% Slin, 1, v, `*`)
    Bimp  <- Bn + Delta
    #Bimp_samples[s,,] <- Bimp
    #Bimp_mean <- Bimp_mean + Bimp
    
    # update horseshoe scales | Bimp, Sigma
    Siginv <- chol2inv(LS)
    Bg <- Bimp[(p0+1):p, , drop = FALSE]         # m x q
    # κ_j = row-wise quadratic form Bg Siginv Bg'
    Kappa <- rowSums((Bg %*% Siginv) * Bg)       # length m
    
    # local scales
    shape_lam <- (q + 1)/2
    rate_lam  <- (1/nu_aux) + 0.5 * Kappa / pmax(tau2, 1e-12)
    lambda2   <- rinv_gamma(shape_lam, rate_lam)
    lambda2   <- pmin(pmax(lambda2, 1e-6), 1e6)
    # aux for locals
    nu_aux    <- rinv_gamma(0.5, 1 + 1/pmax(lambda2, 1e-12))
    
    # global scale
    shape_tau <- (m * q + 1)/2
    rate_tau  <- (1/xi_aux) + 0.5 * sum(Kappa / pmax(lambda2, 1e-12))
    tau2      <- rinv_gamma(shape_tau, rate_tau)
    tau2      <- min(max(tau2, 1e-6), 1e6)
    # aux for global
    xi_aux    <- rinv_gamma(0.5, 1 + 1/pmax(tau2, 1e-12))
    
    # imputation 
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
      theta_j <- mu_i[b] + beta_j %*% t(Yhat[i, a, drop = FALSE] - mu_i[a]) 
      Yhat[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
      if (any(!is.finite(Yhat[i,b])) || any(abs(Yhat[i,b]) > 1e8)){
        Yhat[i,b] <- pmax(pmin(Yhat[i,b], 1e8), -1e8)
      }
    }
    
    if (s%%10 == 0){ flow_miss_samples[s/10,] <- Yhat[miss_ix] } # save results every 10
    #flow_miss_mean <- flow_miss_mean + Yhat[miss_ix]
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples))
}