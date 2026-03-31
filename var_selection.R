library(readr)
library(readxl)
library(mvtnorm)
library(doParallel)
library(foreach)
library(coda)
library(nimble)
library(statmod)
library(AGHmatrix)
library(missRanger)
library(qtl)
library(MASS)
library(heatmaply)
library(ComplexHeatmap)
library(ggrepel)
library(corrplot)
library(cowplot)
library(factoextra)
library(patchwork)
library(loo)
library(matrixStats)
library(ggridges)
#source("code-dependencies/cov_qtl_functions.R")
source("code-dependencies/gibbs_functions.R")

# ---------------------------------Load data---------------------------------- #
cross_data <- read_csv("derived_data/cross_data.csv", show_col_types = FALSE)

# define flow cols 
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
flow_cols <- pheno_names$flow_col_name[-1] # remove titer 

# create flow matrix / remove mice with no flow data 
q <- length(flow_cols)
J <- q*4 # total number of predictors, incl. interaction terms 
mice_w_any_flow <- rowSums(is.na(cross_data[,flow_cols])) < q
flow <- as.matrix(cross_data[mice_w_any_flow,flow_cols]) 
n <- nrow(flow) # 567 mice with any flow observations 

all_colnames <- c(colnames(flow), 
                  paste0(colnames(flow), "_sex"), 
                  paste0(colnames(flow), "_SARS"),
                  paste0(colnames(flow), "_SARS2"))

# outcome 
y <- scale(cross_data$weight_aac, center = TRUE, scale = TRUE)

# design matrix 
options(na.action = "na.pass")
#X <- model.matrix(~ sex + infection + y, data = cross_data) # for cases where downstream analyses are genetic 
X <- model.matrix(~ sex + infection, data = cross_data) # for cases where downstream analyses involve y 
options(na.action = "na.omit")
X <- X[mice_w_any_flow,]
p <- ncol(X)
y <- y[mice_w_any_flow]

# genotype data
geno <- cross_data[mice_w_any_flow, 
                   c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]

# imputed genotype data 
cross_data_imp <- load_cross_as_df("derived_data/cross_imputed.csv", n_geno_start = 122)
geno_imp <- cross_data_imp[mice_w_any_flow, 
                           c(which(colnames(cross_data_imp) == "gUNC145909"):which(colnames(cross_data_imp) == "mJAX00187167"))]

cross_dataI <- read_csv("derived_data/cross_data_flow_imp.csv")
flow_imp <- cross_dataI[,paste0(flow_cols, "_imp")]

# ---------------------------priors/initial values---------------------------- #
n <- nrow(flow) # number of mice 
q <- ncol(flow) # number of effects subject to selection 
p <- ncol(X) # number of fixed effects

# imputation priors 
# prior on Bimp ~ MN(M0, Lambda0, Sigma)
M0 <- matrix(0, p, q) # initialize beta_imp and set M0 (prior mean)  
Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
Lambda0inv <- chol2inv(chol(Lambda0)) 
# prior on Sigma ~ IW(S0, nu0)
# centers the prior at E(Sigma) = I_q 
nu0 <- ncol(flow) + 2
S0 <- diag(q) * (nu0 - q - 1)
Sigma <- S0


muB0 <- rep(0, p) # prior mean
V <- diag(p) # identity matrix 
phi2 <- var(y)*10 # prior variance on betas 

flow_full <- flow # leave original dataset unchanged 
O <- 1*(!is.na(flow))
# initialize all missing values in flow_full to phenotype average
for (j in 1:q){
  flow_full[is.na(flow_full[,j]), j] <- colMeans(flow_full, na.rm = TRUE)[j]
}

a0 <- 3
#b0 <- 0.2
b0 <- var(y)*(a0-1)
sigma2 <- 1 # initialize sigma 

# beta prior on psi (inclusion probability) 
a_psi <- 2
b_psi <- 7
# a_psi <- 1
# b_psi <- 13

# initialize delta to OLS estimates for main effects, 0 for interaction effects
delta_ols <- matrix(lm(y ~ flow)$coefficients[-1], nrow = q, ncol = 1) 
int_deltas <- matrix(0, nrow = 3*q, 1)  
delta <- rbind(delta_ols, int_deltas)  

# warm start for gamma 
marg_corr <- abs(cor(flow_full, y))
top_k <- order(marg_corr, decreasing = TRUE)[1:10]
gamma <- rep(0, (q*4))  
gamma[top_k] <- 1

# SSVS 
tau2 <- 2.5e-05 # spike variance 
c0 <- 300
omega2 <- c0*tau2 # slab variance 

# -------------------------------Gibbs sampler-------------------------------- #
S <- 5000
burn_in <- 1000
set.seed(123)
nchains <- ncores <- 4
dnm_fname <- "results/var_selection/res_dnm.RDS"

if (file.exists(dnm_fname)){
  message(paste("Variable selection results already exist in", dnm_fname))
  res_dnm <- readRDS("results/var_selection/res_dnm.RDS")
} else {
  message(paste("Running variable selection for", S, "iterations over", nchains, "chains"))
  cl <- makeCluster(ncores) 
  registerDoParallel(cl)
  res_dnm <- foreach(ch = 1:nchains, .packages = c("mvtnorm")) %dopar% {
    # matrices to store samples in 
    nsave <- (S - burn_in)
    miss_ix <- which(O == 0)
    flow_miss_samples <- matrix(NA, nrow = nsave, ncol = length(flow_full[O == 0])) # missing values 
    beta_samples <- matrix(NA, nrow = nsave, ncol = p)  # regression effects (fixed)
    psi_samples <- numeric(nsave) # inclusion probability 
    gamma_samples <- matrix(NA, nrow = nsave, ncol = (q*4)) # inclusion indicator 
    delta_samples <- matrix(NA, nrow = nsave, ncol = (q*4)) # regression effects (shrinkage)
    sigma2_samples <- numeric(nsave) # residual variance 
    loglik <- matrix(NA, nrow = nsave, ncol = n)
    
    for (s in 1:S){
      # -------------------------------imputation------------------------------- #
      # update Beta: regression estimates for predicting immune traits
      Lambda_n <- chol2inv(chol(Lambda0inv + t(X) %*% X)) 
      Mn <- Lambda_n %*% (Lambda0inv %*% M0 + t(X) %*% flow_full) 
      Bimp <- matrix(mvtnorm::rmvnorm(1, mean = c(Mn), sigma = kronecker(Sigma, Lambda_n)), nrow = p, ncol = q) 
      
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
      
      # imputation
      for (i in 1:n) {
        b <- (O[i, ] == 0) # missing/masked phenotypes
        if (!any(b)) next
        a <- !b # observed phenotypes
        # incl theta if intercept separate: theta + as.vector(t(Bimp) %*% X[i,])
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
      
      # -------------------------------regression------------------------------- #
      Z <- buildZ(flow_full, X) # create predictor matrix, incl. interaction terms
      
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
        Zj <- Z[, j, drop = FALSE]
        rj <- y - X %*% beta - Z %*% delta + Zj*delta[j] 
        
        # posterior variance/mean of delta_j under spike prior
        s2_spike <- 1 / (fj_norm2[j]/sigma2 + 1/tau2)
        m_spike <- s2_spike * as.numeric(crossprod(Zj, rj)) / sigma2
        # ... and under the SLAB prior
        s2_slab <- 1 / (fj_norm2[j]/sigma2 + 1/omega2)
        m_slab <- s2_slab * as.numeric(crossprod(Zj, rj)) / sigma2
        # conditional Bayes factor (integrating out delta_j)
        logBF_j <- 0.5 * ((log(s2_slab) - log(omega2)) -
                            (log(s2_spike) - log(tau2)) +
                            (m_slab^2)/s2_slab - (m_spike^2)/s2_spike)
        
        # posterior inclusion probability
        prob_imp <- plogis(log(psi) - log1p(-psi) + logBF_j)
        prob_imp <- max(min(prob_imp, 0.99), 0.01) # stabilizer
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
      
      # save samples if past burn in 
      if (s > burn_in) { 
        s_ix <- (s - burn_in)
        flow_miss_samples[s_ix,] <- flow_full[miss_ix]
        beta_samples[s_ix,] <- beta
        psi_samples[s_ix] <- psi
        gamma_samples[s_ix,] <- gamma
        delta_samples[s_ix,] <- delta
        sigma2_samples[s_ix] <- sigma2 
        loglik[s_ix,] <- dnorm(y, mean = as.numeric(X %*% beta + Z %*% delta), sd = sqrt(sigma2), log = TRUE)
      }
    }
    
    # save results for chain 
    chain_res <- list(
      miss_ix = miss_ix, 
      flow_miss = flow_miss_samples,
      beta = beta_samples, 
      psi = psi_samples, 
      gamma = gamma_samples,
      delta = delta_samples, 
      sigma2 = sigma2_samples,
      loglik = loglik
    )
    return(chain_res)
  }
  stopCluster(cl)
  
  saveRDS(res_dnm, file = "results/var_selection/res_dnm.RDS")
}

# ----------------------------sampler performance----------------------------- #
params <- c("beta", "psi", "gamma", "delta", "sigma2")
for (param in params){
  mcmc_list <- mcmc.list(lapply(res_dnm, function(ch) mcmc(ch[[param]])))
  #traceplot(mcmc_list)
  #autocorr.plot(mcmc_list)
  print(effectiveSize(mcmc_list))
  print(gelman.diag(mcmc_list, autoburnin = FALSE))
}

# mixing quality 
rhat_all <- numeric(0)
lab_all  <- character(0)
for (param in params) {
  arr <- to_arr(res_dnm, param)
  rh  <- as.numeric(rhat(arr))
  
  K <- length(rh)
  lab <- if (K == 1) param else paste0(param, "[", seq_len(K), "]")
  
  rhat_all <- c(rhat_all, rh)
  lab_all  <- c(lab_all, lab)
}
df <- data.frame(param = lab_all, rhat = rhat_all)
ecdf <- ggplot(df, aes(x = rhat)) +
  stat_ecdf(size = 0.8) +
  geom_vline(xintercept = 1.01, linetype = 2) +
  geom_vline(xintercept = 1.05, linetype = 3) +
  labs(
    x = expression(hat(R)),
    y = "ECDF",
    title = expression(hat(R)~"distribution across all"~beta~", "~sigma^2~", "~psi~", "~gamma~", "~delta))

# descriptive fit diagnostic (prediction after selection): in-sample Bayesian R-squared 
bayesR2 <- bayes_R2(res_dnm, y, X, flow)
bayesR2_df <- data.frame(r2 = bayesR2$r2)
R2_hist <- ggplot(bayesR2_df, aes(x = r2)) +
  geom_histogram(binwidth = 0.01) + 
  labs(x = expression(Bayesian~R^2), title = expression("In-sample"~Bayesian~R^2)) +
  xlim(c(0,1)) +
  theme_minimal()

# ------------------------------cross-validation------------------------------ #
PIP_thresh <- 0.23
R <- 5 # number of repeats of K-fold CV
K <- 10 # number of folds 
seed1 <- 123

cv_res_fname <- "results/var_selection/diagnostics/cv_res_list.rds"
if (file.exists(cv_res_fname)){
  message(paste("CV results already exist in", cv_res_fname))
  cv_res_list <- readRDS(cv_res_fname)
} else {
  cv_dat <- data.frame(mouse_id = cross_data$mouse_ID[mice_w_any_flow],
                       sex_inf = factor(paste0(cross_data$sex[mice_w_any_flow], 
                                               cross_data$infection[mice_w_any_flow])))
  cl <- makeCluster(R)
  registerDoParallel(cl)
  cv_res_list <- foreach(r = 1:R, .packages = c("mvtnorm", "matrixStats")) %dopar% {
    # K-fold split, stratified by sex/infection
    set.seed(seed1 + r)
    fold_id <- vector("integer", length = nrow(cv_dat))
    for (grp in levels(cv_dat$sex_inf)) {   
      idx <- which(cv_dat$sex_inf == grp)
      fold_id[idx] <- sample(rep(1:K, length.out = length(idx)))
    }
    
    # initialize storage 
    elpd_R <- numeric(K)
    num_sel_R <- numeric(K)
    keep_R <- matrix(FALSE, nrow = K, ncol = J)
    pip_R <- matrix(NA, nrow = K, ncol = J)
    bayesR2_R <- vector("list", K) # holds matrices of varying size for bayesian R2 
    
    for (k in 1:K){
      message(paste0("Repeat", r, ": Fold", k))
      # split data into test/train                     
      test_idx  <- which(fold_id == k)
      train_idx <- setdiff(seq_along(fold_id), test_idx)
      
      # fold-specific y scaling
      y_train <- y[train_idx]
      y_mean <- mean(y_train)
      y_sd   <- sd(y_train)
      y_train <- (y_train - y_mean) / y_sd 
      y_test  <- (y[test_idx] - y_mean) / y_sd # center/scale based on mean/sd from y_train
      
      X_train <- X[train_idx, ,drop=FALSE]
      X_test  <- X[test_idx, ,drop=FALSE]
      
      # immune datasets (with missing data)
      flow_train <- flow[train_idx, , drop = FALSE]  
      flow_test <- flow[test_idx, , drop = FALSE]
      
      # run sampler on training only (1 chain, 1000 iters, 100 burn-in)
      res <- run_sampler_once(y = y_train, X = X_train, flow = flow_train, 
                              X_test = X_test, flow_test = flow_test, 
                              S = 1000, burn_in = 100, n_chains = 1, seed = seed1 + k)
      
      draws <- extract_draws(res)  
      pips <- colMeans(draws$gammas)
      pip_R[k,] <- pips
      
      # get imputed flow_test data, build Z_test 
      miss_ix_test <- res[[1]]$miss_ix_test
      flow_test[miss_ix_test] <- colMeans(draws$flow_miss)
      Z_test <- buildZ(flow_test, X_test)
      
      nsave <- length(draws$sigma2)
      logS <- log(nsave) 
      
      keep <- (pips >= PIP_thresh)
      num_sel_R[k] <- sum(keep)
      keep_R[k,] <- keep
      
      lpd_i <- numeric(length(test_idx))
      y_pred_k <- matrix(NA, nrow = length(y_test), ncol = nsave)
      bayesR2_R[[k]] <- numeric(nsave)
  
      for (i in seq_along(test_idx)) {
        mu <- as.numeric(draws$betas %*% X_test[i,]) # nsave-length vector of predictive means for mouse i
        if (any(keep)) {
          mu <- mu + as.numeric(draws$deltas[,keep, drop = FALSE] %*% as.numeric(Z_test[i,keep]))
        }
        y_pred_k[i,] <- mu # for calculation of Bayesian R-squared
        # compute log predictive density for each test obs i
        ll_s <- dnorm(y_test[i], mean = mu, sd = sqrt(draws$sigma2s), log = TRUE) # log-likelihood
        # subtract log(n) bc we are computing the log of the average predictive density across posterior draws 
        lpd_i[i] <- matrixStats::logSumExp(ll_s) - logS
      }
      elpd_R[k] <- sum(lpd_i)
      
      for (iter in 1:nsave){
        mu <- y_pred_k[,iter] # n(tst)-length vector of predictions for sampler iteration iter 
        var_res <- var(y_test - mu)
        var_fit <- var(mu)
        # each row of bayesR2_R[[k]][t,] is the distribution of bayes R2 for that fold and threshold 
        bayesR2_R[[k]][iter] <- var_fit / (var_fit + var_res)  
      }
    }
    return(list(pip = pip_R,
                keep_array = keep_R,
                num_sel = num_sel_R,
                elpd = elpd_R,
                bayesR2 = bayesR2_R))
  }
  stopCluster(cl)
  
  ensure_directory("results/var_selection/diagnostics")
  saveRDS(cv_res_list, "results/var_selection/diagnostics/cv_res_list.rds")
}

# predictive accuracy: out-of-sample Bayesian R-squared
all_R2 <- NULL
for (r in 1:R){
  for (k in 1:K){
    all_R2 <- c(all_R2, cv_res_list[[r]]$bayesR2[[k]])
  }
}
oos_bayesR2_df <- data.frame(r2 = all_R2)
oosR2_hist <- ggplot(oos_bayesR2_df, aes(x = r2)) +
  geom_histogram(binwidth = 0.01) + 
  labs(x = expression(Bayesian~R^2), 
       title = expression("Out-of-sample"~Bayesian~R^2)) +
  xlim(c(0,1)) +
  theme_minimal()

# variable stability
J <- ncol(cv_res_list[[1]]$pip) # number of variables
stability_mat <- matrix(NA, nrow = R*K, ncol = J)
pip_mat <- matrix(NA, nrow = R*K, ncol = J)
row_idx <- 0
for (r in 1:R) {
  for (k in 1:K) {
    row_idx <- row_idx + 1
    stability_mat[row_idx,] <- cv_res_list[[r]]$keep_R[k,] # was variable j selected in repeat r, fold k?
    pip_mat[row_idx,] <- cv_res_list[[r]]$pip[k,] # PIPs for this fold
  }
}
variable_stability <- data.frame(variable = 1:J,
                                 mean_pip = colMeans(pip_mat, na.rm = TRUE),
                                 sd_pip = apply(pip_mat, 2, sd, na.rm = TRUE),
                                 stability = colMeans(stability_mat), # proportion selected across all folds
                                 n_selected = colSums(stability_mat), # raw count
                                 total_folds = R * K)
variable_stability$var_name <- all_colnames
variable_stability <- variable_stability[order(-variable_stability$stability),]
variable_stability$highly_stable <- variable_stability$stability >= 0.70
stable_vars <- all_colnames[variable_stability$variable[variable_stability$highly_stable]]
stable_phenos <- unique(sub("(_sex|_SARS2|_SARS)$", "", stable_vars))

# stability vs mean PIP 
stability_vs_pip <- ggplot(variable_stability, aes(x = mean_pip, y = stability)) +
  geom_point(aes(color = highly_stable), alpha = 0.6, size = 3) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  geom_vline(xintercept = PIP_thresh, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "darkgreen"),
                     name = "Stable\n(≥70%)") +
  labs(x = "Mean PIP (across all folds)", 
       y = paste0("Selection Stability at threshold = ", PIP_thresh),
       title = "Variable Selection Stability") +
  theme_bw() +
  theme(legend.position = "right")

# calculate Jaccard similarity across all fold pairs (within and across repeats) 
n_folds_total <- R * K
n_pairs <- choose(n_folds_total, 2) # total number of fold pairs
pairs <- comb(1:n_folds_total, 2)
Jvals <- apply(pairs, 2, function(idx) {
  jaccard(stability_mat[idx[1], ], stability_mat[idx[2], ])
})
# calculate overlap sizes for each pair
overlaps <- apply(pairs, 2, function(idx) {
  sum(stability_mat[idx[1], ] & stability_mat[idx[2], ])
})
# histogram of jaccard similarity
jacc_hist <- ggplot(data.frame(jacc = Jvals), aes(x = jacc)) +
  geom_histogram(binwidth = 0.01) + 
  labs(x = "Jaccard similarity", 
       title = "Pairwise Jaccard similarity") +
  xlim(c(0,1)) +
  theme_minimal()

# -------------------------------sampler results------------------------------ #
gammas_dnm <- do.call(rbind, lapply(res_dnm, function(x) x$gamma))

PIPs <- colMeans(gammas_dnm)  # q-vector of inclusion probabilities
ranked_vars <- order(PIPs, decreasing = TRUE)
PIPdf <- data.frame(phenotype = all_colnames[ranked_vars],
                    PIP = PIPs[ranked_vars])
important_vars <- PIPdf$phenotype[PIPdf$PIP > PIP_thresh] 
important_phenos <- unique(sub("(_sex|_SARS2|_SARS)$", "", important_vars))
important_vars <- important_vars[important_vars %in% stable_vars]
important_phenos <- important_phenos[important_phenos %in% stable_phenos]

imp_phenos_df <- data.frame(phenotype = important_phenos, 
                            stable = important_phenos %in% stable_phenos)
write_csv(imp_phenos_df, "results/var_selection/important_phenos.csv")

# calculate joint PIPs
base_name <- gsub("(_sex|_SARS|_SARS2)$", "", all_colnames)
idx_by_base <- split(seq_along(all_colnames), base_name)
joint_pip_vec <- vapply(idx_by_base, function(idx) mean(rowSums(gammas_dnm[, idx, drop = FALSE]) > 0), numeric(1))
jointPIPdf <- data.frame(base = names(joint_pip_vec), 
                         jointPIP = as.numeric(joint_pip_vec), 
                         row.names = NULL)
jointPIPdf <- merge(jointPIPdf, multiR2_df, by.x = "base", by.y = "phenotype") 
jointPIPdf$VIF <- 1 / (1 - jointPIPdf$multiR2)
jointPIPdf <- jointPIPdf[order(jointPIPdf$jointPIP, decreasing = TRUE),]

# calculate multiple correlation coefficient 
multiR2s <- sapply(1:q, function(j) summary(lm(flow[,j] ~ flow[,-j]))$r.squared)
multiR2_df <- data.frame(phenotype = colnames(flow),
                         multiR2 = multiR2s)
# is there a relationship between multiple R-squared and joint PIP?
summary(lm(jointPIP ~ multiR2, data = jointPIPdf))
plot(jointPIPdf$multiR2, jointPIPdf$jointPIP)

# wide table with PIPs and joint PIPs
PIP_wide <- PIPdf %>% 
  mutate(base = sub("(_sex|_SARS2|_SARS)$", "", phenotype),
         term = case_when(str_ends(phenotype, "_sex") ~ "sex",
                          str_ends(phenotype, "_SARS") ~ "SARS",
                          str_ends(phenotype, "_SARS2") ~ "SARS2",
                          TRUE ~ "main")) %>%
  select(base, term, PIP) %>% 
  pivot_wider(names_from = term, values_from = PIP) %>% 
  select(base, main, sex, SARS, SARS2)

PIP_wide <- merge(PIP_wide, jointPIPdf[,1:3], by = "base")

# summarize deltas (effects)
deltas_dnm <- do.call(rbind, lapply(res_dnm, function(x) x$delta))
colnames(deltas_dnm) <- all_colnames

col_mean <- colMeans(deltas_dnm)
# col_lwr <- apply(deltas_dnm, 2, quantile, probs = 0.025)
# col_upr <- apply(deltas_dnm, 2, quantile, probs = 0.975)
col_pip <- colMeans(gammas_dnm) 
# conditional mean (given inclusion) 
conditional_mean <- map_dbl(seq_len(ncol(deltas_dnm)), function(j) {
  incl <- gammas_dnm[,j] > 0
  mean(deltas_dnm[incl, j])
})
# lower bound of 95% CI, conditional on inclusion 
cond_lwr <- vapply(seq_len(ncol(deltas_dnm)), function(j){
  incl <- gammas_dnm[,j] > 0
  as.numeric(quantile(deltas_dnm[incl, j], probs = 0.025, names = FALSE))
}, numeric(1))
# lower bound of 50% CI, conditional on inclusion 
cond_q25 <- vapply(seq_len(ncol(deltas_dnm)), function(j){
  incl <- gammas_dnm[,j] > 0
  as.numeric(quantile(deltas_dnm[incl, j], probs = 0.25, names = FALSE))
}, numeric(1))
# upper bound of 50% CI, conditional on inclusion 
cond_q75 <- vapply(seq_len(ncol(deltas_dnm)), function(j){
  incl <- gammas_dnm[,j] > 0
  as.numeric(quantile(deltas_dnm[incl, j], probs = 0.75, names = FALSE))
}, numeric(1))
# upper bound of 95% CI, conditional on inclusion 
cond_upr <- vapply(seq_len(ncol(deltas_dnm)), function(j){
  incl <- gammas_dnm[,j] > 0
  as.numeric(quantile(deltas_dnm[incl, j], probs = 0.975, names = FALSE))
}, numeric(1))

base <- gsub("(_sex|_SARS|_SARS2)$", "", all_colnames)
term <- case_when(str_ends(all_colnames, "_sex") ~ "Sex",
                  str_ends(all_colnames, "_SARS") ~ "SARS-CoV",
                  str_ends(all_colnames, "_SARS2") ~ "SARS-CoV-2",
                  TRUE ~ "Main")

effect_sum <- data.frame(phenotype = all_colnames,
                         base = base,
                         term = factor(term, levels = c("Main", "SARS-CoV", "SARS-CoV-2", "Sex")),
                         mean = as.numeric(col_mean),
                         lower = as.numeric(cond_lwr),
                         upper = as.numeric(cond_upr),
                         PIP = as.numeric(col_pip),
                         conditional_mean = conditional_mean,
                         cond_q25 = cond_q25,
                         cond_q75 = cond_q75)

# filter to important phenotypes 
effect_sum <- effect_sum %>% 
  filter(base %in% important_phenos) %>% 
  left_join(select(PIP_wide, base, jointPIP), by = "base") %>% 
  mutate(base = reorder(base, jointPIP))

# are SARS2-interaction CIs larger than the others?
ci_wide <- effect_sum %>%
  filter(term %in% c("Main", "SARS-CoV", "SARS-CoV-2", "Sex")) %>%
  mutate(ci_width = upper - lower) %>%
  select(base, term, ci_width) %>%
  distinct() %>%
  pivot_wider(names_from = term, values_from = ci_width)
wilcox.test(ci_wide$`SARS-CoV-2`, ci_wide$Main, paired = TRUE, alternative = "greater")
wilcox.test(ci_wide$`SARS-CoV-2`, ci_wide$`SARS-CoV`, paired = TRUE, alternative = "greater")
wilcox.test(ci_wide$`SARS-CoV-2`, ci_wide$Sex, paired = TRUE, alternative = "greater")

ci_wide2 <- ci_wide %>%
  mutate(other_mean = (Main + `SARS-CoV` + Sex) / 3)
wilcox.test(ci_wide2$`SARS-CoV-2`, ci_wide2$other_mean, paired = TRUE, alternative = "greater")

# ------------------(virus-specific) posterior mean effects------------------- #
# map columns to base + term
col_map <- tibble(
  phenotype = colnames(deltas_dnm),
  base = gsub("(_sex|_SARS|_SARS2)$", "", colnames(deltas_dnm)),
  term = case_when(
    str_ends(colnames(deltas_dnm), "_SARS")  ~ "SARS",
    str_ends(colnames(deltas_dnm), "_SARS2") ~ "SARS2",
    str_ends(colnames(deltas_dnm), "_sex")   ~ "Sex",
    TRUE ~ "Main"))

# helper: get column index for a given base + term
get_col <- function(b, t) {
  ph <- col_map$phenotype[col_map$base == b & col_map$term == t]
  if (length(ph) == 0) return(NA_integer_)
  match(ph, colnames(deltas_dnm))
}

#bases_keep <- important_phenos
bases_keep <- flow_phenos

virus_conditional_means <- map_dfr(bases_keep, function(b) {
  
  j_main  <- get_col(b, "Main")
  j_sars  <- get_col(b, "SARS")
  j_sars2 <- get_col(b, "SARS2")
  
  if (anyNA(c(j_main, j_sars, j_sars2))) {
    return(tibble(
      base = b,
      beta_SARS_cond_mean = NA_real_,
      beta_SARS2_cond_mean = NA_real_,
      n_sars_used = 0,
      n_sars2_used = 0
    ))
  }
  
  tot_sars  <- deltas_dnm[, j_main] + deltas_dnm[, j_sars]
  tot_sars2 <- deltas_dnm[, j_main] + deltas_dnm[, j_sars2]
  
  # "total active" gating
  active_sars  <- (gammas_dnm[, j_main] > 0) | (gammas_dnm[, j_sars] > 0)
  active_sars2 <- (gammas_dnm[, j_main] > 0) | (gammas_dnm[, j_sars2] > 0)
  
  tibble(
    base = b,
    
    beta_SARS_cond_mean  = mean(tot_sars[active_sars]),
    beta_SARS2_cond_mean = mean(tot_sars2[active_sars2]),
    
    n_sars_used  = sum(active_sars),
    n_sars2_used = sum(active_sars2)
  )
})

# ----------------------consistency of effect direction----------------------- #
# identify columns for Main / SARS / SARS2
base <- gsub("(_sex|_SARS|_SARS2)$", "", colnames(deltas_dnm))
term <- case_when(
  str_ends(colnames(deltas_dnm), "_SARS") ~ "SARS",
  str_ends(colnames(deltas_dnm), "_SARS2") ~ "SARS2",
  str_ends(colnames(deltas_dnm), "_sex") ~ "Sex",
  TRUE ~ "Main"
)

# restrict to important phenos (optional)
# bases_keep <- important_phenos 
# col_map <- tibble(phenotype = colnames(deltas_dnm),
#                   base = base,
#                   term = term) %>%
#   filter(base %in% bases_keep) %>%
#   filter(term %in% c("Main", "SARS", "SARS2"))

col_map <- tibble(phenotype = colnames(deltas_dnm),
                  base = base,
                  term = term) %>%
  filter(term %in% c("Main", "SARS", "SARS2"))

bases <- sort(unique(col_map$base))

effect_direction <- function(deltas_dnm, gammas_dnm, bases, get_col, eps = 0.01) {
  
  purrr::map_dfr(bases, function(b) {
    
    j_main  <- get_col(b, "Main")
    j_sars  <- get_col(b, "SARS")
    j_sars2 <- get_col(b, "SARS2")
    
    if (anyNA(c(j_main, j_sars, j_sars2))) {
      return(tibble::tibble(
        base = b,
        prop_same_sign = NA_real_,
        n_draws_used = 0
      ))
    }
    
    tot_sars  <- deltas_dnm[, j_main] + deltas_dnm[, j_sars]
    tot_sars2 <- deltas_dnm[, j_main] + deltas_dnm[, j_sars2]
    
    # total is "active" if Main OR interaction is included
    active_sars  <- (gammas_dnm[, j_main]  > 0) | (gammas_dnm[, j_sars]  > 0)
    active_sars2 <- (gammas_dnm[, j_main]  > 0) | (gammas_dnm[, j_sars2] > 0)
    
    keep <- active_sars & active_sars2 &
      (abs(tot_sars)  > eps) &
      (abs(tot_sars2) > eps)
    
    if (sum(keep) == 0) {
      return(tibble::tibble(
        base = b,
        prop_same_sign = NA_real_,
        n_draws_used = 0
      ))
    }
    
    tibble::tibble(
      base = b,
      prop_same_sign = mean(sign(tot_sars[keep]) == sign(tot_sars2[keep])),
      n_draws_used = sum(keep)
    )
  })
}

dir_summary <- effect_direction(
  deltas_dnm = deltas_dnm,
  gammas_dnm = gammas_dnm,
  bases = bases,
  get_col = get_col,
  eps = 0.01
)

dir_summary
sum(dir_summary$prop_same_sign > 0.8)
dir_summary[which.min(dir_summary$prop_same_sign),]

# -------------------------------Supp. Table 1-------------------------------- #
varcomp_data <- read_csv("results/var_comp_res.csv", show_col_types = FALSE)
clust_ord <- read_csv("results/qtl_mapping/hclust_pheno_order.csv")
stbl1 <- varcomp_data %>% 
  left_join(pheno_names, by = c("pheno" = "flow_col_name")) %>%
  left_join(jointPIPdf, by = c("pheno" = "base")) %>%
  left_join(virus_conditional_means, by = c("pheno" = "base")) %>%
  left_join(dir_summary, by = c("pheno" = "base")) %>% 
  left_join(clust_ord, by = "pheno") %>%
  select(flow_display_name, multiR2, VIF, h2, h2_lwr, h2_upr, jointPIP, beta_SARS_cond_mean, beta_SARS2_cond_mean, prop_same_sign, n_draws_used, order) %>%
  mutate(across(starts_with("h2"), ~ .x*100)) %>%
  mutate(across(where(is.numeric) & !all_of("order"), ~ round(.x, 2))) %>%
  arrange(flow_display_name)
writexl::write_xlsx(stbl1, "results/var_selection/STable1.xlsx")


# ---------------------------------------------------------------------------- #
# ----------------------------------figures----------------------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------DNM prior---------------------------------- #
dnm_d1 <- density(rnorm(100000, mean = 0, sd = sqrt(omega2)))
dnm_d2 <- density(rnorm(100000, mean = 0, sd = sqrt(tau2)))

df <- rbind(data.frame(x = dnm_d1$x, y = dnm_d1$y, component = "omega2"),
            data.frame(x = dnm_d2$x, y = dnm_d2$y, component = "tau2"))
df$component <- factor(df$component, levels = c("tau2", "omega2"))

dnm_prior <- ggplot(df, aes(x = x, y = y, color = component)) +
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Density", title = "Discrete normal mixture prior", color = NULL) +
  scale_color_discrete(labels = c("omega2" = expression("slab (var = "*omega^2*")"), "tau2" = expression("spike (var = "*tau^2*")"))) + 
  theme_minimal() + 
  theme(legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13))

# --------------------------------Supp. Fig. 1-------------------------------- #
plot_grid(
  plot_grid(dnm_prior, ecdf, R2_hist, ncol = 3, labels = c("A", "B", "C")),
  ggdraw() + draw_label("Cross-validation results", fontface = "bold", x = 0.01, hjust = 0),
  plot_grid(oosR2_hist, stability_vs_pip, jacc_hist, ncol = 3, 
            rel_widths = c(1,1.2,1), labels = c("D", "E", "F")),
  nrow = 3, rel_heights = c(1,0.15,1)
)
ggsave("figures/supplemental/gibbs_performance.png", width = 11, height = 7, bg = "white")

# -------------------------post dist for missing val-------------------------- #
miss1 <- unlist(lapply(res_dnm, function(x) x$flow_miss[,1]))
df <- data.frame(value = miss1)
ggplot(df, aes(x = value)) + 
  geom_density(fill = "gray70", alpha = 0.5, linewidth = 1) + 
  geom_vline(xintercept = median(miss1), lty = 2) + 
  labs(x = NULL, y = "Density", title = "Posterior distribution of first missing value") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 15))
ggsave("figures/supplemental//miss1.png")

# ---------------------------post effects for j=33---------------------------- #
j = 33
main_deltas <- unlist(lapply(res_dnm, function(x) x$delta[,j]))
sex_deltas <- unlist(lapply(res_dnm, function(x) x$delta[,(j+105)]))
sars1_deltas <- unlist(lapply(res_dnm, function(x) x$delta[,(j+210)]))
sars2_deltas <- unlist(lapply(res_dnm, function(x) x$delta[,(j+315)]))

df <- rbind(data.frame(value = main_deltas, effect = "main"),
            data.frame(value = sars1_deltas, effect = "sars_int"),
            data.frame(value = sars2_deltas, effect = "sars2_int"),
            data.frame(value = sex_deltas, effect = "sex_int"))
df$effect <- factor(df$effect, levels = rev(c("main", "sars_int", "sars2_int", "sex_int")))
set2_n5 <- brewer.pal(n = 5, "Set2")
effect_pal <- c("main" = set2_n5[5], "sars_int" = set2_n5[2], "sars2_int" = set2_n5[3], "sex_int" = set2_n5[4])

ggplot(df, aes(x = value, y = effect, fill = effect)) + 
  geom_density_ridges(scale = 0.9, alpha = 0.8, color = "black", linewidth = 0.25, rel_min_height = 0.01) +
  labs(x = "Value", y = NULL, fill = "Effect", title = "Posterior distribution of effects\nfor DN T cells (% CM)") + 
  theme_classic() +
  scale_y_discrete(labels = c("main" = "Main", "sars_int" = "SARS-CoV\ninteraction",
                              "sars2_int" = "SARS-CoV-2\ninteraction", "sex_int" = "Sex\ninteraction"),
                   expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(-0.051, 0.351)) + # breaks = seq(-0.05, 0.35, by = 0.1)
  scale_fill_manual(values = effect_pal) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 13))
ggsave("figures/supplemental/dnt_post.png", width = 5.5, height = 6)

# ------------------------------------PCA------------------------------------- #
pc <- prcomp(flow_imp)
hab <- ifelse(cross_dataI$infection == "PBS", "Control", 
              ifelse(cross_dataI$infection == "SARSCoV", "SARS-CoV", "SARS-CoV-2"))
pca_plot <- fviz_pca_ind(pc, geom = "point", habillage = hab, addEllipses = TRUE) +
  scale_color_brewer(name = "Infection", palette = "Set2") +
  scale_shape_discrete(name = "Infection") + 
  scale_fill_brewer(name = "Infection", palette = "Set2") +
  guides(color = guide_legend(title = "Infection"),
         shape = guide_legend(title = "Infection"),
         fill = guide_legend(title = "Infection")) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 15))

# ----------------------------------heatmap----------------------------------- #
corr <- cor(flow, use = "pairwise.complete.obs")
dcols <- dist(1-corr, method = "euclidian")
hc_cols <- hclust(dcols, method = "complete")
hm <- Heatmap(corr, show_row_dend = FALSE, show_row_names = FALSE, 
              show_column_names = FALSE, col = heatmaply::cool_warm(256), 
              column_dend_height = unit(0.7,"inches"), column_title = "Immune trait correlation heatmap",
              cluster_columns = as.dendrogram(hc_cols), cluster_rows = as.dendrogram(hc_cols),
              heatmap_legend_param = list(title = "",
                                          title_gp = gpar(fontsize = 14, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 12),
                                          grid_width = unit(0.6, "cm"),
                                          grid_height = unit(0.6, "cm"),
                                          legend_height = unit(5, "cm"),
                                          at = seq(-1, 1, by = 0.5),
                                          direction = "vertical")) 
hmgrob <- grid.grabExpr(draw(hm))

# Labeled heatmap for supplement 
pdf("figures/supplemental/heatmap_labeled.pdf", width = 16, height = 16)
labeled_hm <- Heatmap(corr, show_row_dend = FALSE, 
                      row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                      col = heatmaply::cool_warm(256), column_dend_height = unit(1,"inches"), 
                      cluster_columns = as.dendrogram(hc_cols), cluster_rows = as.dendrogram(hc_cols),
                      heatmap_legend_param = list(title = "",
                                                  title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                  labels_gp = gpar(fontsize = 8),
                                                  grid_width = unit(0.6, "cm"),
                                                  grid_height = unit(0.6, "cm"),
                                                  legend_height = unit(5, "cm"),
                                                  at = seq(-1, 1, by = 0.5),
                                                  direction = "vertical"))
draw(labeled_hm)
dev.off()

# -------------------------------PIP histograms------------------------------- #
hist_df <- data.frame(PIP = c(PIPdf$PIP, PIP_wide$jointPIP),
                      "Type" = c(rep("Individual", nrow(PIPdf)), rep("Joint", nrow(PIP_wide))))
pip_hist <- ggplot(hist_df, aes(x = PIP, fill = Type)) + 
  geom_histogram(alpha = 0.6, position = "identity") + 
  geom_vline(xintercept = c(PIP_thresh), linetype = "dashed") +
  theme_minimal() + 
  theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill = "white", color = NA))
# 0.469 is the minimum joint PIP of all "important" variables 
# (as determined by having an individual term with PIP > PIP_thresh)

pip_hist_ind <- ggplot(hist_df[hist_df$Type == "Individual",], aes(x = PIP)) +
  geom_histogram() + 
  geom_vline(xintercept = c(PIP_thresh), linetype = "dashed") + 
  xlim(0,1) + 
  labs(x = "Posterior inclusion probability (PIP)", title = "Individual PIPs") +
  theme_minimal() 

pip_hist_joint <- ggplot(hist_df[hist_df$Type == "Joint",], aes(x = PIP)) +
  geom_histogram() +
  xlim(0,1) +
  labs(x = "Posterior inclusion probability (PIP)", title = "Joint PIPs") +
  theme_minimal() 

both_hists <- plot_grid(pip_hist_ind, pip_hist_joint, nrow = 2, labels = c("C", "D"))

# -----------------------------joint PIP barplot------------------------------ #
pheno_names <- pheno_names %>%
  mutate(flow_display_name2 = flow_display_name %>% 
           str_replace_all("⁺", "+") %>%
           str_replace_all("⁻", "-"))

PIP_wide_filt <- PIP_wide[PIP_wide$base %in% important_phenos,]
topN <- length(important_phenos)
plot_df <- PIP_wide_filt %>% 
  left_join(pheno_names, by = c("base" = "flow_col_name")) %>%
  mutate(base = flow_display_name2,
         total_comp = main + sex + SARS + SARS2) %>%
  select(-flow_display_name, flow_display_name2) %>%
  mutate(w_main = main / total_comp,
         w_sex = sex / total_comp,
         w_SARS = SARS/ total_comp,
         w_SARS2 = SARS2 / total_comp) %>%
  pivot_longer(cols = starts_with("w_"), names_to = "term", values_to = "weight") %>%
  mutate(term = case_when(term == "w_main" ~ "Main",
                          term == "w_sex" ~ "Sex",
                          term == "w_SARS" ~ "SARS-CoV",
                          term == "w_SARS2" ~ "SARS-CoV-2"),
         segment = jointPIP*weight) %>% 
  group_by(base) %>% 
  ungroup() %>% 
  slice_max(order_by = jointPIP, n = topN*4, with_ties = TRUE) %>%
  mutate(base = reorder(base, jointPIP),
         term = factor(term, levels = c("Sex", "SARS-CoV-2", "SARS-CoV", "Main")))

set2_colors <- c(brewer.pal(5, "Set2")[2:5])
names(set2_colors) <- c("SARS-CoV", "SARS-CoV-2", "Sex", "Main")

jpip_barplot <- ggplot(plot_df, aes(x = base, y = segment, fill = term)) +
  geom_col(width = 0.8, color = "grey30", linewidth = 0.2) +
  coord_flip() + 
  scale_y_continuous(name = "Joint inclusion probability",
                     limits = c(0, max(plot_df$jointPIP) * 1.08),
                     breaks = pretty_breaks(6),
                     labels = number_format(accuracy = 0.01)) +
  scale_fill_manual(values = set2_colors, 
                    breaks = c("Main", "SARS-CoV", "SARS-CoV-2", "Sex"), 
                    labels = str_wrap(c("Main effect", "SARS-CoV interaction", "SARS-CoV-2 interaction", "Sex interaction"),
                                      width = 12),
                    name = "Component") +
  labs(x = "Phenotype") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = c(0.8, 0.3),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.height = unit(1.4, "lines"))

# ------------------------effect size caterpillar plot------------------------ #
lvl <- c("Main", "SARS-CoV", "SARS-CoV-2", "Sex")
eff_df <- effect_sum %>%
  mutate(term = factor(as.character(term), levels = lvl),
         facet = factor(as.character(term), levels = c("Joint inclusion", lvl)))

term_labels <- c("Main" = "Main effect",
                 "SARS-CoV" = "SARS-CoV interaction",
                 "SARS-CoV-2" = "SARS-CoV-2 interaction",
                 "Sex" = "Sex interaction")

effect_plot <- ggplot(eff_df, aes(y = base, color = term)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, size = 0.4, alpha = 0.9) +
  geom_errorbar(aes(xmin = cond_q25, xmax = cond_q75), width = 0, size = 1) +
  xlim(-0.3,0.3) +
  geom_point(aes(x = conditional_mean), size = 1.6) +
  facet_wrap(~ term, ncol = 4, scales = "free_x",
             labeller = as_labeller(term_labels)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
  scale_color_manual(values = set2_colors) +
  labs(x = "Posterior effects given inclusion", y = NULL) +
  theme_minimal(base_size = 12) + 
  theme(legend.position = "none",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

# ----------------------------------Figure 2---------------------------------- #
f2_top <- plot_grid(hmgrob, pca_plot, both_hists, ncol = 3, labels = c("A", "B", ""), rel_widths = c(0.9, 1, 0.7))
f2_bottom <- (jpip_barplot | effect_plot) + 
  plot_layout(widths = c(1.2, 3)) +
  plot_annotation(tag_levels = list(c("E", "F"))) &
  theme(plot.tag = element_text(face = "bold"))
f2 <- plot_grid(f2_top, f2_bottom, nrow = 2, rel_heights = c(1,1.5))
ggsave("figures/Figure2.png", bg = "white", width = 14, height = 11)
