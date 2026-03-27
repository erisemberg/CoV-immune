library(readr)
library(readxl)
library(mvtnorm)
library(AGHmatrix)
#source("code-dependencies/cov_qtl_functions.R")
source("code-dependencies/cmdline.R")
source("code-dependencies/gibbs_functions.R")

# ---------------------------------load data---------------------------------- #
run_cv <- cmdline.logical("runCV", default = FALSE)
comp_cv <- cmdline.logical("compareCV", default = FALSE)

# load phenotypes and covariates 
cross_data <- read_csv("derived_data/cross_data.csv", show_col_types = FALSE)

# define flow cols 
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
flow_cols <- pheno_names$flow_col_name[-1] # remove titer 
# filter to mice with any flow data
mice_w_any_flow <- rowSums(is.na(cross_data[,flow_cols])) < length(flow_cols)
flow <- as.matrix(cross_data[mice_w_any_flow, flow_cols]) 

# covariate design matrix 
options(na.action = "na.pass")
# for cases where downstream analyses are genetic only, can include y as a predictor 
# X <- model.matrix(~ sex + infection + y, data = cross_data)
# for cases where downstream analyses involve y, e.g. variable selection, including 
# y as a predictor would cause data leakage 
X <- model.matrix(~ sex + infection, data = cross_data) 
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

# ---------------------------create imputed phenos---------------------------- #
# create standard imputed version of phenotypes from the imputation model based 
# on fixed effects for sex and infection and var-covar matrix Sigma (run_chain1)
#
# **Note**: Because these imputed variables are being used for QTL mapping *only*, 
# I include weight loss as a predictor in the design matrix X. For the combined 
# imputation / variable selection Gibbs sampler, weight loss will *not* be included
# in the design matrix, as that would induce data leakage / circular reasoning 
# (that is, using weight loss to impute missing data, then using imputed data to 
# determine how important each variable is in the prediction of weight loss).

q <- ncol(flow)
n <- nrow(flow)
p <- ncol(X)

# prior hyperparameters 
S0 <- diag(q) * (nu0 - q - 1)
nu0 <- q + 2
M0 <- matrix(0, p, q) # initialize beta_imp and set B0 (prior mean) - ((p-1) x q) matrix 
Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)

S <- 3000
burn_in <- 300
seed = 123
ncores <- chains <- 4

Y <- flow
O <- matrix(1L, n, q) # O indicates whether data is observed or missing
O[is.na(Y)] <- 0L 

# run 4 chains 
cl <- makeCluster(ncores)
registerDoParallel(cl)
imp_res1 <- foreach(ch = 1:chains, .packages = c("mvtnorm")) %dopar% run_chain1(ch)
stopCluster(cl)

# take posterior mean for imputed dataset 
flow_miss1 <- do.call(rbind, lapply(imp_res1, function(x) x$flow_miss))
flow_miss1_means <- colMeans(flow_miss1)
flow_imp <- flow
flow_imp[is.na(flow_imp)] <- flow_miss1_means
colnames(flow_imp) <- paste0(colnames(flow), "_imp")

# create data frame for QTL mapping 
cross_dataI <- cross_data[mice_w_any_flow,]
cross_dataI <- cbind(cross_dataI, flow_imp)
write_csv(cross_dataI, "derived_data/cross_data_flow_imp.csv")


# ------------------------------cross-validation------------------------------ #  
S <- 1000 
set.seed(seed)

# run cross-validation for chosen imputation method (i.e., imputation based on 
# variance-covariance matrix and fixed effect covariates, in run_chain1)
if (run_cv | comp_cv){
  # ----------------------imputation with fixed effects----------------------- #
  # define priors and hyperparameters 
  S0 <- diag(q) * (nu0 - q - 1)
  nu0 <- q + 2
  M0 <- matrix(0, p, q) # initialize beta_imp and set M0 (prior mean) - ((p-1) x q) matrix 
  Lambda0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
  
  # 10-fold cross_validation
  cv_res1 <- impute_cv(chain_fn_name = "run_chain1", Y = flow, K = 10, S = 1000, X = X,
                       chains = 4, seed = 123, ncores = 4, 
                       S0 = S0, nu0 = nu0, M0 = M0, Lambda0 = Lambda0)
  cv_res1$per_fold  
}

# also run cross-validation for other imputation methods:
#   1. imputation based on variance-covariance matrix only 
#   2. imputation based on Sigma, fixed effects and random polygenic effect 
#   3. imputation based on Sigma, fixed effects and random batch effects 
#   4. imputation based on Sigma, fixed effects and horseshoe-shrunk genetic effects 
if (comp_cv){ 
  # -------------------------imputation based on Sigma-------------------------- #  
  # define additional hyperparameters 
  mu0 <- colMeans(flow, na.rm = TRUE) # average of each immune phenotype 
  sd0 <- mu0/2 # chosen to keep most of prior mass on values above zero 
  L0 <- matrix(0.1, q, q) # lambda0 = prior var-covar matrix on means of immune phenotypes  
  diag(L0) <- 1
  L0 <- L0*outer(sd0, sd0) # alt: diag(q) * 1e3
  S0 <- L0 # alt: diag(q)
  
  # run 10-fold cross-validation
  cv_res0 <- impute_cv(chain_fn_name = "run_chain0", Y = flow, K = 10, S = 1000, 
                       chains = 4, mu0 = mu0, L0 = L0, S0 = L0, nu0 = nu0,
                       seed = 123, ncores = 4)
  
  # compare predicted values from cross-validation code to those from original code 
  imp_res0_flowmiss <- do.call(rbind, lapply(imp_res0, function(x) x$flow_miss))
  plot(colMeans(imp_res0_flowmiss), colMeans(cv_res0$pred_data[,is.na(flow)])) 
  abline(a = 0, b = 1)
  
  # -----------------imputation with random polygenic effect------------------ #
  # genetic relatedness matrix
  G <- Gmatrix(as.matrix(geno), method = "VanRaden")
  all_mouseIDs <- cross_data$mouse_ID[mice_w_any_flow]
  colnames(G) <- rownames(G) <- all_mouseIDs
  
  # specify priors 
  Sg0 <- diag(q) # prior for genetic variance 
  B0 <- matrix(0, p, q) # initialize beta_imp and set B0 (prior mean) - ((p-1) x q) matrix 
  V0 <- diag(p)*n # large variance / tune as desired ### try diag(1,k)
  # set nug0 to something else? prior df for genetic variance 
  
  # pass in variables that are passed into run_chain1 as well as ones above 
  cv_res2 <- impute_cv(chain_fn_name = "run_chain2", Y = flow, K = 10, S = 1000, 
                       G = G, chains = 4, ncores = 4, seed = 123, S0 = S0, 
                       nu0 = nu0, B0 = B0, V0 = V0, X = X, Sg0 = Sg0, nug0 = nu0)
  
  # --------------------imputation with random batch effects-------------------- #
  # # define priors / hyperparameters
  # Xr <- X[,-1]
  # p <- ncol(Xr)
  # 
  # nu_tau0 <- q + 3
  # c <- 0.3 # prior batch SD per trait
  # S_tau0 <- (nu_tau0 - q - 1) * (c^2) * diag(q)
  # 
  # batch <- as.factor(cross_data$flow_batch[mice_w_any_flow]) # length n (int))
  # recode_df <- data.frame(original = levels(as.factor(batch)),
  #                         new = seq(1, length(levels(as.factor(batch)))))
  # recoded_batch <- as.data.frame(batch) %>%
  #   mutate(batch = as.character(batch)) %>%                    
  #   left_join(recode_df %>% mutate(original = as.character(original)),
  #             by = c("batch" = "original")) %>%
  #   mutate(batch_int = as.integer(new)) %>%                    
  #   select(-new)
  # batch <- recoded_batch$batch_int
  # 
  # # run 10-fold cross-validation 
  # cv_res3 <- impute_cv(chain_fn_name = "run_chain3", Y = flow, K = 10, S = 1000,
  #                      chains = 4, seed = 123, ncores = 4, 
  #                      S0 = S0, nu0 = nu0, B0 = B0, V0 = V0, X = Xr,
  #                      nu_tau0 = nu_tau0, S_tau0 = S_tau0, batch = batch)
  # cv_res3$per_fold
  
  # -------------imputation with horseshoe-shrunk genetic effects------------- #
  # # define priors / hyperparameters
  # Psi0 <- S0
  # v_fixed <- 10 # variance for fixed effects not subject to shrinkage 
  # # design matrices 
  # geno_imp_ctr <- scale(geno_imp)
  # Xall <- as.matrix(cbind(X, geno_imp_ctr))
  # p0 <- ncol(X)
  # m <- ncol(geno_imp)
  # p <- ncol(Xall)
  # lambda2 <- rep(1,m) # local scales 
  # nu_aux <- rep(1,m) # auxiliaries for lambda2
  # tau2 <- 1 # global scale  
  # xi_aux <- 1 # auxiliary for tau2 
  # 
  # # run 10-fold CV 
  # cv_res4 <- impute_cv(chain_fn_name = "run_chain4", Y = flow, K = 10, seed = 123, 
  #                      S = 1000, chains = 4, ncores = 4, X = Xall, B0 = B0, 
  #                      v_fixed = v_fixed, p0 = p0, tau2 = tau2, lambda2 = lambda2,
  #                      Psi0 = Psi0, nu0 = nu0, nu_aux = nu_aux, m = m, xi_aux = xi_aux)
  
  # -------------------------------RF imputation------------------------------ #
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
  
  # ----------------------------compare CV results---------------------------- #
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
}


