library(readr)
library(readxl)
library(mvtnorm)
library(AGHmatrix)
source("/Users/ellenrisemberg/Documents/ValdarFerris/scripts/qtl_functions.R")
setwd("/Users/ellenrisemberg/Documents/ValdarFerris/Coronavirus/cov-immune")

# ---------------------------------load data---------------------------------- #
# load phenotypes and covariates 
cross_data <- read_csv("derived_data/cross_data.csv", show_col_types = FALSE)

# define flow cols 
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
flow_cols <- pheno_names$flow_col_name[-1] # remove titer 
q <- length(flow_cols)
# predictors 
mice_w_any_flow <- rowSums(is.na(cross_data[,flow_cols])) < q
flow <- as.matrix(cross_data[mice_w_any_flow,flow_cols]) 
n <- nrow(flow) 
# design matrix 
options(na.action = "na.pass")
#X <- model.matrix(~ sex + infection + y, data = cross_data) # for cases where downstream analyses are genetic 
X <- model.matrix(~ sex + infection, data = cross_data) # for cases where downstream analyses involve y 
options(na.action = "na.omit")
X <- X[mice_w_any_flow,]
p <- ncol(X)

# outcome 
y <- scale(cross_data$weight_aac, center = TRUE, scale = TRUE)
y <- y[mice_w_any_flow]

# genotype data
geno <- cross_data[mice_w_any_flow, 
                   c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]

# imputed genotype data 
cross_data_imp <- load_cross_as_df("derived_data/cross_imputed.csv", n_geno_start = 122)
geno_imp <- cross_data_imp[mice_w_any_flow, 
                           c(which(colnames(cross_data_imp) == "gUNC145909"):which(colnames(cross_data_imp) == "mJAX00187167"))]

# ------------------------------basic imputation------------------------------ #  
mu0 <- colMeans(flow, na.rm = TRUE) # average of each immune phenotype 
sd0 <- mu0/2 # chosen to keep most of prior mass on values above zero 
L0 <- matrix(0.1, q, q) # lambda0 = prior var-covar matrix on means of immune phenotypes  
diag(L0) <- 1
L0 <- L0*outer(sd0, sd0)
S0 <- L0
nu0 <- ncol(flow) + 2
# ^ should prior scale matrix and prior covariance matrix also be sample covariance matrix? 
#S0 <- cov(flow, use = "pairwise.complete.obs")

# default weak priors
#L0 <- diag(q) * 1e3
#S0 <- diag(q)

# starting values 
#Sigma <- cov(flow, use = "pairwise.complete.obs") # Sigma initalized to sample covariance matrix - should S0 also be sample covariance matrix? 
Sigma <- S0

flow_full <- flow # leave original dataset unchanged 
O <- 1*(!is.na(flow))
# initialize all missing values in flow_full to phenotype average
for (j in 1:q){
  flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
}

# code for imputation only, with none of the numerical instability solutions / cholesky 
# stuff. Compare to results from cross-validation code which has a bunch of extra stuff 
# for avoiding numerically unstable Sigma. 
S <- 1000
set.seed(123)
cl <- makeCluster(4) # 4 chains, run in parallel
registerDoParallel(cl)

imp_res0 <- foreach(ch = 1:4, .packages = c("mvtnorm")) %dopar% {
  # matrices to store samples in 
  flow_miss_samples <- matrix(NA, nrow = S, ncol = length(flow_full[O == 0])) # missing values 
  theta_samples <- matrix(NA, nrow = S, ncol = q) # means of immune phenotypes
  Sigma_samples <- matrix(NA, nrow = S, ncol = q*q) # covariance of immune phenotypes
  
  for (s in 1:S){
    # update theta (predictor means) 
    flow_bar <- apply(flow_full, 2, mean) # get current phenotype averages 
    Ln <- solve(solve(L0) + n*solve(Sigma)) 
    mu_n <- Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% flow_bar) 
    theta <- rmvnorm(1, mu_n, Ln)  
    theta_samples[s,] <- theta # save results 
    
    # update Sigma (variance-covariance matrix for predictors)
    Sn <- S0 + (t(flow_full) - c(theta)) %*% t(t(flow_full) - c(theta))  
    Sigma <- solve(rWishart(1, nu0 + n, solve(Sn))[,,1]) # inverse-Wishart  
    Sigma_samples[s,] <- c(Sigma) # save results 
    
    # update missing data 
    for (i in 1:n){
      b <- (O[i,] == 0) # phenotypes for which data is missing 
      a <- (O[i,] == 1) # phenotypes for which data is observed 
      if(sum(b) == 0){ next } # no missing data
      iSa <- solve(Sigma[a,a]) # 99 x 99 (when 6 missing, 99 observed variables for mouse i)
      beta_j <- Sigma[b,a] %*% iSa # 6 x 99 
      Sigma_j <- Sigma[b,b] - Sigma[b,a] %*% iSa %*% Sigma[a,b]
      theta_j <- theta[b] + beta_j %*% t(t(flow_full[i,a]) - theta[a])
      flow_full[i,b] <- rmvnorm(1, theta_j, Sigma_j) # update missing phenotypes for mouse i 
    }
    flow_miss_samples[s,] <- flow_full[O == 0] # save results
  }
  
  # save results for chain 
  chain_res <- list(
    flow_miss = flow_miss_samples,
    thetaD = theta_samples, 
    SigmaD = Sigma_samples
  )
  return(chain_res)
}

stopCluster(cl)



# ----------------------------CV helper functions----------------------------- # 
# a tiny nudge used as a practical guard against floating point noise; shouldn't be
# the mechanism that makes your sampler run - if it frequently ramps up noise, chain
# will be biases. 
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

# Function for creating cell-wise folds over observed cells (i.e., exclude originally 
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
  
  fold_metrics <- data.frame(
    fold = 1:K,
    n_holdout = NA_integer_,
    RMSE = NA_real_,
    R2 = NA_real_
  )
  
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


# -------------------------imputation based on Sigma-------------------------- #  
# function to run one chain for imputation based on variance-covariance matrix 
# between phenotypes (Sigma) only
run_chain0 <- function(chain_id){
  set.seed(seed + chain_id)
  flow_full <- Y
  n <- nrow(Y)
  q <- ncol(Y)
  # initialize missing with column means 
  for (j in 1:q){
    flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
  }
  
  theta_samples <- matrix(NA, S, q) # immune phenotype means
  Sigma_samples <- matrix(NA, S, q*q) # immune phenotype variance-covariance matrix
  
  # initialize 
  theta <- colMeans(flow_full)
  Sigma <- S0
  Lambda <- solve(S0) # precision 
  
  # store missing draws each iteration 
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = S, ncol = length(miss_ix))
  
  Lambda0 <- chol2inv(chol(L0)) # Lambda0 = L0^{-1}
  
  for (s in 1:S){
    # update theta (predictor means) 
    flow_bar <- apply(flow_full, 2, mean) # get current phenotype averages 
    Ln_inv <- Lambda0 + n * Lambda
    Ln_inv <- (Ln_inv + t(Ln_inv))/2 # symmetrize
    if (inherits(try(chol(Ln_inv), silent = TRUE), "try-error")){
      Ln_inv <- as.matrix(Matrix::nearPD(Ln_inv)$mat)
    }
    Ln <- chol2inv(chol(Ln_inv)) 
    mu_n <- Ln %*% (Lambda0 %*% mu0 + n * Lambda %*% flow_bar)
    theta <- rmvnorm(1, mu_n, Ln)  
    theta_samples[s,] <- theta # save results 
    
    # update Sigma (variance-covariance matrix for predictors)
    Sn <- S0 + (t(flow_full) - c(theta)) %*% t(t(flow_full) - c(theta)) 
    Sn <- (Sn + t(Sn))/2 # symmetrize  
    if (inherits(try(chol(Sn), silent = TRUE), "try-error")){
      Sn <- as.matrix(Matrix::nearPD(Sn)$mat)
    }
    Sn_inv <- chol2inv(chol(Sn))
    Sn_inv <- (Sn_inv + t(Sn_inv))/2 # symmetrize 
    if (inherits(try(chol(Sn_inv), silent = TRUE), "try-error")){
      Sn_inv <- as.matrix(Matrix::nearPD(Sn_inv)$mat)
    }
    Lambda <- rWishart(1, nu0 + n, Sn_inv)[,,1] # draw precision from Wishart 
    if (inherits(try(chol(Lambda), silent = TRUE), "try-error")){
      Lambda <- as.matrix(Matrix::nearPD(Lambda)$mat)
    }
    Sigma <- chol2inv(chol(Lambda)) # recover Sigma from precision using Cholesky
    Sigma_samples[s,] <- c(Sigma) # save results 
    
    for (i in 1:n) {
      b <- (O[i, ] == 0) # missing/masked phenotypes 
      if (!any(b)) next
      a <- !b # observed phenotypes 
      Sa <- Sigma[a, a, drop = FALSE] 
      ch <- chol(Sa)
      beta_j <- Sigma[b, a, drop = FALSE] %*% chol2inv(ch)
      W <- backsolve(ch, Sigma[a, b, drop = FALSE], transpose = TRUE)
      Sigma_j <- Sigma[b, b, drop = FALSE] - t(W) %*% W
      theta_j <- theta[b] + beta_j %*% t(flow_full[i, a, drop = FALSE] - theta[a])
      if (inherits(try(chol(Sigma_j), silent = TRUE), "try-error")){
        Sigma_j <- make_spd(Sigma_j)
      }
      flow_full[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    
    flow_miss_samples[s,] <- flow_full[miss_ix] # save results
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(theta = theta_samples,
                             Sigma = Sigma_samples)))
}

# run 10-fold cross-validation
cv_res0 <- impute_cv(chain_fn_name = "run_chain0", Y = flow, K = 10, S = 1000, 
                     chains = 4, mu0 = mu0, L0 = L0, S0 = L0, nu0 = nu0,
                     seed = 123, ncores = 4)

# compare predicted values from cross-validation code to those from original code 
imp_res0_flowmiss <- do.call(rbind, lapply(imp_res0, function(x) x$flow_miss))
plot(colMeans(imp_res0_flowmiss), colMeans(cv_res0$pred_data[,is.na(flow)])) 
abline(a = 0, b = 1)



# -----------------------imputation with fixed effects------------------------ #
# imputation based on var-covar matrix Sigma and fixed effects for sex and infection 
run_chain1 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  Yhat <- Y
  # initialize missing values with column means 
  for (j in 1:q){
    Yhat[is.na(Yhat[,j]), j] <- colMeans(Yhat, na.rm = TRUE)[j]
  }
  Sigma <- S0
  Lambda0inv <- chol2inv(chol(Lambda0)) 
  
  nsave <- S - burn_in
  Sigma_samples <- matrix(NA, nsave, q*q) # immune phenotype variance-covariance matrix 
  Bimp_samples <- array(NA, dim = c(nsave, p, q)) # fixed effects in imputation model 
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

# define priors and hyperparameters 
S0 <- diag(q) * (nu0 - q - 1)
nu0 <- q + 2
M0 <- matrix(0, p, q) # initialize beta_imp and set M0 (prior mean) - ((p-1) x q) matrix 
Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)

# 10-fold cross_validation
cv_res1 <- impute_cv(chain_fn_name = "run_chain1", Y = flow, K = 10, S = 1000, X = X,
                     chains = 4, seed = 123, ncores = 4, S0 = S0, nu0 = nu0, M0 = M0, Lambda0 = Lambda0)
cv_res1$per_fold

# ---------------------------create imputed phenos---------------------------- #
# create standard imputed version of phenotypes from the imputation model above 
# (run_chain1), based on fixed effects for sex and infection and var-covar matrix Sigma 

# **Note**: Because these imputed variables are being used for QTL mapping *only*, 
# I include weight loss as a predictor in the design matrix X. For the combined 
# imputation / variable selection Gibbs sampler, weight loss will *not* be included
# in the design matrix, as that would induce data leakage / circular reasoning 
# (that is, using weight loss to impute missing data, then using imputed data to 
# determine how important each variable is in the prediction of weight loss).
S <- 3000
burn_in <- 300
Y <- flow
seed = 123
ncores <- chains <- 4
O <- matrix(1L, n, q) # O indicates whether data is observed or missing
O[is.na(Y)] <- 0L 
# prior hyperparameters 
M0 <- matrix(0, p, q) # initialize beta_imp and set B0 (prior mean) - ((p-1) x q) matrix 
Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
S0 <- diag(q) * (nu0 - q - 1)
nu0 <- q + 2

cl <- makeCluster(ncores)
registerDoParallel(cl)
imp_res1 <- foreach(ch = 1:chains, .packages = c("mvtnorm")) %dopar% run_chain1(ch)
stopCluster(cl)

# Create imputed phenotypes for comparison with SynSurr
flow_miss1 <- do.call(rbind, lapply(imp_res1, function(x) x$flow_miss))
flow_miss1_means <- colMeans(flow_miss1)
flow_imp <- flow
flow_imp[is.na(flow_imp)] <- flow_miss1_means
colnames(flow_imp) <- paste0(colnames(flow), "_imp")

# create data frame for QTL mapping 
cross_dataI <- cross_data[mice_w_any_flow,]
cross_dataI <- cbind(cross_dataI, flow_imp)
write_csv(cross_dataI, "derived_data/cross_data_flow_imp.csv")


# ------------------imputation with random polygenic effect------------------- #
# function to run one chain of Gibbs sampler based on fixed effects for sex and 
# infection, a random polygenic effect, and var-covar matrix Sigma 
run_chain2 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  
  flow_full <- Y
  for (j in 1:q){ # initialize missing values with column means
    flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
  }
  Z <- diag(n) # individual-level random effects 
  
  Sigma_samples <- matrix(NA, S, q*q) # immune phenotype variance-covariance matrix 
  Bimp_samples <- array(NA, dim = c(S, p, q)) # fixed effects in imputation model 
  Sigma_g_samples <- matrix(NA, S, q*q) # genetic variance-covariance matrix for immune phenotypes 
  U_samples <- array(NA, dim = c(S, n, q)) # each U matrix is n x q 
  miss_ix <- which(O == 0)
  flow_miss_samples <- matrix(NA, nrow = S, ncol = length(miss_ix))
  
  # initialize parameters
  # also initialize theta <- colMeans(flow_full) if including intercept 
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
    #----------------(1) Update B | Sigma, Y, U (matrix-normal)----------------#
    Vn <- chol2inv(chol(V0inv + t(X) %*% X)) # p x p
    flow_residR <- flow_full - Z%*%U # residuals after random effects 
    Bn <- Vn %*% (V0inv %*% B0 + t(X) %*% flow_residR) # p x q  
    Bimp <- t(matrix(mvtnorm::rmvnorm(1, c(t(Bn)), kronecker(Sigma, Vn)), q, p)) # transpose to get p x q
    Bimp_samples[s,,] <- Bimp
    
    #--------------(2) Update Sigma | Y, U (factorize via Sigma)---------------#
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
    
    #-----------(3) Update U | B, Sigma, Y, Sigma_g (spectral trick)-----------#
    # vec(U) ~ N(0, kronecker(Sigma_g, G)) 
    ### CONFIRM 
    flow_residF <- flow_full - X %*% Bimp # residuals after fixed effects
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
    
    #-----------------------(4) Update Sigma_g | U ----------------------------#
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
      theta_j <- mu_i[b] + beta_j %*% t(flow_full[i, a, drop = FALSE] - mu_i[a]) 
      flow_full[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    flow_miss_samples[s,] <- flow_full[miss_ix] # save results
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(Beta_imp = Bimp_samples,
                             Sigma = Sigma_samples, 
                             Sigma_g = Sigma_g_samples,
                             U = U_samples)))
}

# genetic relatedness matrix
G <- Gmatrix(as.matrix(geno), method = "VanRaden")
all_mouseIDs <- cross_data$mouse_ID[mice_w_any_flow]
colnames(G) <- rownames(G) <- all_mouseIDs

# specify priors 
Sg0 <- diag(q) # prior for genetic variance 
# set nug0 to something else? prior df for genetic variance 

# pass in variables that are passed into run_chain1 as well as ones above 
cv_res2 <- impute_cv(chain_fn_name = "run_chain2", Y = flow, K = 10, S = 1000, 
                     G = G, chains = 4, ncores = 4, seed = 123, S0 = S0, 
                     nu0 = nu0, B0 = B0, V0 = V0, X = X, Sg0 = Sg0, nug0 = nu0)

# --------------------imputation with random batch effects-------------------- #
# imputation based on fixed effects for sex and infection, random effect for 
# batch, and var-covar matrix Sigma 
run_chain3 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  G <- max(batch)
  
  flow_full <- Y
  # initialize missing values with column means 
  for (j in 1:q){
    flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
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
    flow_resid <- flow_full - U[batch,,drop = FALSE] # n x q ### CONFIRM U[batch,]
    Vn <- chol2inv(chol(V0inv + t(X) %*% X)) # p x p
    Bn <- Vn %*% (V0inv %*% B0 + t(X) %*% flow_resid) # p x q ### CONFIRM: replaced flow_full w flow_resid  
    Bimp <- t(matrix(mvtnorm::rmvnorm(1, c(t(Bn)), kronecker(Sigma, Vn)), q, p)) # transpose to get p x q
    Bimp_samples[s,,] <- Bimp
    
    # update Sigma (conjugate normal-Wishart update)
    flow_resid <- flow_full - X %*% Bimp - U[batch, , drop = FALSE] # n x q
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
      Rg <- flow_full[ids,,drop = FALSE] - X[ids,,drop = FALSE] %*% Bimp 
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
      theta_j <- mu_i[b] + beta_j %*% t(flow_full[i, a, drop = FALSE] - mu_i[a]) 
      flow_full[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
    }
    flow_miss_samples[s,] <- flow_full[miss_ix] # save results
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples,
              samples = list(Beta_imp = Bimp_samples,
                             Sigma = Sigma_samples,
                             U = U_samples,
                             Tau = Tau_samples)))
}

# define priors / hyperparameters
#Xr <- X[,-1]
p <- ncol(Xr)

S0 <- diag(q) * (nu0 - q - 1)
B0 <- matrix(0, p, q) # initialize beta_imp and set B0 (prior mean) - ((p-1) x q) matrix 
V0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)

nu_tau0 <- q + 3
c <- 0.3 # prior batch SD per trait
S_tau0 <- (nu_tau0 - q - 1) * (c^2) * diag(q)

batch <- as.factor(cross_data$flow_batch[mice_w_any_flow]) # length n (int))
recode_df <- data.frame(original = levels(as.factor(batch)),
                        new = seq(1, length(levels(as.factor(batch)))))
recoded_batch <- as.data.frame(batch) %>%
  mutate(batch = as.character(batch)) %>%                    
  left_join(recode_df %>% mutate(original = as.character(original)),
            by = c("batch" = "original")) %>%
  mutate(batch_int = as.integer(new)) %>%                    
  select(-new)
batch <- recoded_batch$batch_int

# run 10-fold cross-validation 
cv_res3 <- impute_cv(chain_fn_name = "run_chain3", Y = flow, K = 10, S = 1000,
                     chains = 4, seed = 123, ncores = 4, 
                     S0 = S0, nu0 = nu0, B0 = B0, V0 = V0, X = Xr,
                     nu_tau0 = nu_tau0, S_tau0 = S_tau0, batch = batch)
cv_res3$per_fold



# --------------imputation with horseshoe-shrunk genetic effects-------------- #
# Imputation with fixed effects for sex, infection, horseshoe-shrunk effects 
# for all genotypes, and var-covar matrix Sigma

### CODE NOT VALIDATED 
run_chain4 <- function(chain_id){
  set.seed(seed + chain_id)
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  flow_full <- Y
  for (j in 1:q){ # initialize missing values with column means
    flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
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
    #-------------------------build V (row variances)--------------------------#
    v_min <- 1e-8; v_max <- 1e8
    v <- c(rep(v_fixed, p0), tau2*lambda2)
    v <- pmin(pmax(v, v_min), v_max)
    
    #-----------------compute B without pxp inverse (Woodbury)-----------------#
    # Bn = B0 + V X' A^{-1} (Y - X B0),  A = I + X V X'
    R  <- flow_full - X %*% B0                # n x q
    Xv <- sweep(X, 2, v, `*`)                 # n x p  (X * v_j by columns)
    A  <- diag(n) + Xv %*% t(X)               # n x n
    # solve A %*% Z = R  via Cholesky
    if (inherits(try(chol(A), silent = TRUE), "try-error")){ A <- as.matrix(Matrix::nearPD(A)$mat) }
    cholA <- chol(A)
    Z  <- backsolve(cholA, forwardsolve(t(cholA), R, upper.tri = TRUE, transpose = TRUE))
    Bn <- B0 + t(X) %*% Z                     # p x q, because X' A^{-1} R, then multiply each row by v:
    Bn <- sweep(Bn, 1, v, `*`)                   # (V X' A^{-1} R) + B0
    
    # ---------- Sigma | Y, V  ~ IW(Psi_n, nu0+n) ------------------------
    # Use (★): Psi_n = Psi0 + (Y - X Bn)'(Y - X Bn) + (Bn - B0)' Λ0 (Bn - B0).
    # Λ0 = diag(1/v). Avoid forming Λ0 explicitly.
    E <- flow_full - X %*% Bn                    # n x q
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
    
    # ---------- B | Sigma, Y, V  ~ MN(Bn, Vn, Sigma) --------------------
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
    
    # ---------- Horseshoe scales given (Bimp, Sigma) --------------------
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
      theta_j <- mu_i[b] + beta_j %*% t(flow_full[i, a, drop = FALSE] - mu_i[a]) 
      flow_full[i, b] <- mvtnorm::rmvnorm(1, mean = c(theta_j), sigma = Sigma_j)
      if (any(!is.finite(flow_full[i,b])) || any(abs(flow_full[i,b]) > 1e8)){
        flow_full[i,b] <- pmax(pmin(flow_full[i,b], 1e8), -1e8)
      }
    }
    
    if (s%%10 == 0){ flow_miss_samples[s/10,] <- flow_full[miss_ix] } # save results every 10
    #flow_miss_mean <- flow_miss_mean + flow_full[miss_ix]
  }
  
  return(list(miss_ix = miss_ix, 
              flow_miss = flow_miss_samples))
}

# define priors / hyperparameters
Psi0 <- S0
v_fixed <- 10 # variance for fixed effects not subject to shrinkage 
# design matrices 
geno_imp_ctr <- scale(geno_imp)
Xall <- as.matrix(cbind(X, geno_imp_ctr))
p0 <- ncol(X)
m <- ncol(geno_imp)
p <- ncol(Xall)
B0 <- matrix(0, p, q)

lambda2 <- rep(1,m) # local scales 
nu_aux <- rep(1,m) # auxiliaries for lambda2
tau2 <- 1 # global scale  
xi_aux <- 1 # auxiliary for tau2 

# run 10-fold CV 
cv_res3 <- impute_cv(chain_fn_name = "run_chain4", Y = flow, K = 10, seed = 123, 
                     S = 1000, chains = 4, ncores = 4, X = Xall, B0 = B0, 
                     v_fixed = v_fixed, p0 = p0, tau2 = tau2, lambda2 = lambda2,
                     Psi0 = Psi0, nu0 = nu0, nu_aux = nu_aux, m = m, xi_aux = xi_aux)


# -----------------------------compare CV results----------------------------- #
# extended design matrix for RF imputation 
options(na.action = "na.pass")
Xrf <- model.matrix(~ sex + infection + HS + Titer + weight_aac + as.factor(flow_batch) + 
                      as.factor(batch) + pgm + F1_dam_geno + Gater, data = cross_data)[,-1]
options(na.action = "na.omit")
Xrf[,"weight_aac"] <- scale(Xrf[,"weight_aac"])
Xrf[,"HS"] <- scale(Xrf[,"HS"])
Xrf[,"Titer"] <- scale(Xrf[,"Titer"])
Xrf <- Xrf[mice_w_any_flow,]

# CV with random forest 
cv_resRF <- impute_cv(Y = flow, K = 10, seed = 123, RF = TRUE, X = Xrf, geno_imp = geno_imp)
cv_resRFnoG <- impute_cv(Y = flow, K = 10, seed = 123, RF = TRUE, X = Xrf)

# compare R2 and RMSE 
message(sprintf("Overall R-squared for random forest with genotypes: %.4f", cv_resRF$overall$R2))
message(sprintf("Overall R-squared for random forest without genotypes: %.4f", cv_resRFnoG$overall$R2))
message(sprintf("Overall R-squared for model based on Sigma: %.4f", cv_res0$overall$R2))
message(sprintf("Overall R-squared for model based on Sigma + sex + infection + weight loss: %.4f", cv_res1$overall$R2))
message(sprintf("Overall R-squared for model based on Sigma + sex + infection + weight loss + polygenic effect: %.4f", cv_res2$overall$R2))

# plot results 
r2_df <- data.frame(Fold = seq(1,10),
                    RandomForest = cv_resRF$per_fold$R2,
                    RandomForestnoG = cv_resRFnoG$per_fold$R2,
                    Sigma = cv_res0$per_fold$R2,
                    SigmaSexInf = cv_res1$per_fold$R2,
                    SigmaSexInfRandG = cv_res2$per_fold$R2)
r2_df_long <- pivot_longer(r2_df, 
                           cols = c("RandomForest", "RandomForestnoG", "Sigma", "SigmaSexInf", "SigmaSexInfRandG"), 
                           names_to = "Model", 
                           values_to = "R-squared")

labs <- c("Random Forest", "Random Forest w/o genotypes", "Bayesian: Sigma", "Bayesian: Sigma + Sex/Infection", "Bayesian: Sigma + Sex/Infection + GRM")
png("flow/imputation/model_comparison.png", width = 700)
ggplot(data = r2_df_long, mapping = aes(x = Fold, y = `R-squared`, color = Model)) + 
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3) +
  scale_x_continuous(breaks = seq(1,10)) +
  scale_color_discrete(labels = function(x) str_wrap(labs, width = 20)) +
  ylim(0,1) + 
  big_theme
dev.off()


### Imputation results 
Y <- flow 
q <- ncol(Y); n <- nrow(Y); p <- ncol(X)
O <- matrix(1L, n, q) 
O[which(is.na(Y))] <- 0L
# X already defined - intercept, sex, infection 
S <- 1000
chains <- ncores <- 4
S0 <- diag(q) * (nu0 - q - 1)
B0 <- matrix(0, p, q) # initialize beta_imp and set B0 (prior mean) - ((p-1) x q) matrix 
V0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
nu0 <- ncol(flow) + 2

exported_vars <- c("make_spd", "seed", "Y", "O", "S", "X", "S0", "B0", "V0", "nu0")

cl <- makeCluster(ncores)
registerDoParallel(cl)
chain_res <- foreach(ch = 1:chains, 
                     .export = exported_vars, 
                     .packages = c("mvtnorm")) %dopar% run_chain1(ch)
stopCluster(cl)


# traceplot 
plot(1:1000, chain_res[[1]]$samples$Beta_imp[,2,1], type = "l", col = "green", main = "Sex effect on L_pctCD45ofLive")
lines(1:1000, chain_res[[2]]$samples$Beta_imp[,2,1], type = "l", col = "orange")
lines(1:1000, chain_res[[3]]$samples$Beta_imp[,2,1], type = "l", col = "lightblue")
lines(1:1000, chain_res[[4]]$samples$Beta_imp[,2,1], type = "l", col = "purple")


