#' Convert preset options to log prior case probabilities function for GxT analysis
#'
#' This function takes the log prior case probabilities, and if a preset is provided, converts it into the formal log prior case
#' probability.
#'
#' @param ln_prior_c Log prior case probabilities. If posterior_summary() is being used for a non-default posterior odds
#' summary, the log prior case probabilities used with bmediatR() are stored in its output.
#' @return \code{return_ln_prior_c_from_presets_gxt} returns a vector of length 64 consisting of
#' \code{0}s for models to include and \code{-Inf}s for models to exclude. The order of models
#' is the same as c1-c64 described in \code{model_info()}.
#' @export
#' @examples return_ln_prior_c_from_presets_gxt()
return_ln_prior_c_from_presets_gxt <- function(ln_prior_c) {
  if (is.character(ln_prior_c) && ln_prior_c[1] == "complete_gxt") {
    ln_prior_c <- rep(0,64) # uniform prior: log(1/64) after normalization 
  } else {
    ln_prior_c <- ln_prior_c # user-provided custom prior 
  }
  ln_prior_c
}

#' Summaries of posterior model probabilities function
#'
#' This function takes the log posterior probability of the data (posterior likelihood) for the various cases, the log prior case probabilities, and
#' returns log posterior odds.
#'
#' @param ln_prob_data Log posterior likelihoods under the various models, returned by bmediatR().
#' @param ln_prior_c Log prior case probabilities. If posterior_summary_gxt() is being used for a non-default posterior odds
#' summary, the log prior case probabilities used with bmediatR() are stored in its output. NOTE: not implemented yet for gxt mediation 
#' @param c_numerator The index of cases to be summed in the numerator of the posterior odds. Cases, their order, and likelihoods
#' are provided in model_info().
#' @return \code{posterior_summary_gxt} returns a list containing the following components:
#'  \item{ln_post_c}{a matrix with posterior probabilities of each causal model for each candidate mediator.}
#'  \item{ln_post_odds}{a matrix with posterior odds of individual models or combinations of models for each candidate mediator.}
#'  \item{ln_prior_odds}{a single row matrix with prior odds of individual models or combinations of models.}
#'  \item{ln_ml}{the natural log of the marginal likelihood.}
#' @export
#' @examples posterior_summary_gxt()
posterior_summary_gxt <- function(ln_prob_data, ln_prior_c, c_numerator) {
  
  # function to compute log odds from log probabilities
  ln_odds <- function(ln_p, numerator){
    ln_odds_numerator <- apply(ln_p[,numerator, drop=F], 1, matrixStats::logSumExp)
    ln_odds_denominator <- apply(ln_p[,-numerator, drop=F], 1, matrixStats::logSumExp)
    ln_odds <- ln_odds_numerator -ln_odds_denominator
  }
  
  # ensure c_numerator is a list
  if (!is.list(c_numerator)) {
    c_numerator <- list(c_numerator)
  }
  
  # presets for ln_prior_c
  ln_prior_c <- return_ln_prior_c_from_presets_gxt(ln_prior_c = ln_prior_c)
  
  # ensure ln_prior_c sum to 1 on probability scale and that it is a matrix
  if (is.matrix(ln_prior_c)){
    ln_prior_c <- t(apply(ln_prior_c, 1, function(x){x - matrixStats::logSumExp(x)}))
  } else {
    ln_prior_c <- ln_prior_c - matrixStats::logSumExp(ln_prior_c)
    ln_prior_c <- matrix(ln_prior_c, nrow(ln_prob_data), length(ln_prior_c), byrow=T)
  }
  
  # compute posterior probabilities for all cases
  # ln_prob_data already contains joint likelihoods, so just add prior for each model 
  ln_post_c <- ln_prob_data + ln_prior_c 
  
  # compute marginal likelihood and normalize posterior 
  ln_ml <- apply(ln_post_c, 1, matrixStats::logSumExp)
  ln_post_c <- ln_post_c - ln_ml
  rownames(ln_post_c) <- rownames(ln_prob_data) # might not be necessary 
  
  # compute prior odds for each combination of cases
  ln_prior_odds <- sapply(c_numerator, ln_odds, ln_p = ln_prior_c)
  ln_prior_odds <- matrix(ln_prior_odds, ncol = length(c_numerator))
  rownames(ln_prior_odds) <- rownames(ln_post_c)
  
  # compute posterior odds for each combination of cases
  ln_post_odds <- sapply(c_numerator, ln_odds, ln_p = ln_post_c)
  ln_post_odds <- matrix(ln_post_odds, ncol = length(c_numerator))
  rownames(ln_post_odds) <- rownames(ln_post_c)
  
  if (is.null(c_numerator)) {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- c_numerator
  } else {
    colnames(ln_post_odds) <- colnames(ln_prior_odds) <- names(c_numerator)
  }
  
  # return results
  list(ln_post_c = ln_post_c, 
       ln_prior_c = ln_prior_c, 
       ln_post_odds = ln_post_odds, 
       ln_prior_odds = ln_prior_odds,
       ln_prob_data = ln_prob_data,
       ln_ml=ln_ml)
}

#' Column indices for commonly used posterior odds
#'
#' This helper function returns the columns of the log posterior case probabilities to be summed for
#' commonly desired log posterior odds summaries.
#'
#' @param odds_type The desired posterior odds.
#' @return \code{return_preset_odds_index_gxt} returns a list with indices for individual models and combinations of models.
#' @export
#' @examples return_preset_odds_index_gxt()
return_preset_odds_index_gxt <- function(odds_type = c("mediation_any", "mediation_grp1", "mediation_grp2",
                                                       "partial_any", "partial_grp1", "partial_grp2",
                                                       "complete_any", "complete_grp1", "complete_grp2",
                                                       "colocal_any", "colocal_grp1", "colocal_grp2")) {
  
  presets <- list(
    "mediation_any" = c(18, 20, 22, 24, 26, 28, 30, 32, 50, 54, 58, 62, # group 1
                        35, 36, 39, 40, 43, 44, 47, 48, 51, 55, 59, 63, # group 2
                        52, 56, 60, 64), # both groups
    "mediation_grp1" = c(18, 20, 22, 24, 26, 28, 30, 32, 50, 54, 58, 62, # group 1
                         52, 56, 60, 64),  # both groups
    "mediation_grp2" = c(35, 36, 39, 40, 43, 44, 47, 48, 51, 55, 59, 63,  # group 2 
                         52, 56, 60, 64),  # both groups 
    "partial_any" = c(22, 24, 30, 32, 54, 56, 62, # group 1 partial
                      43, 44, 47, 48, 59, 60, 63, # group 2 partial
                      64),  # both groups partial
    "partial_grp1" = c(22, 24, 30, 32, 54, 56, 62, # group 1 partial 
                       64),  # both groups partial
    "partial_grp2" = c(43, 44, 47, 48, 59, 60, 63, # group 2 partial
                       64),  # both groups partial
    "complete_any" = c(18, 20, 26, 28, 50, 58, 60, # group 1 complete
                       35, 36, 39, 40, 51, 55, 56, # group 2 complete 
                       52),  # both groups complete
    "complete_grp1" = c(18, 20, 26, 28, 50, 58, 60, # group 1 complete
                        52),  # both groups complete
    "complete_grp2" = c(35, 36, 39, 40, 51, 55, 56, # group 2 complete
                        52),  # both groups complete
    "colocal_any" = c(21, 23, 29, 31, 53, 55, 63, # group 1 colocal 
                      41, 42, 45, 46, 57, 58, 62, # group 2 colocal 
                      61),  # both groups colocal
    "colocal_grp1" = c(21, 23, 29, 31, 53, 55, 63, # group 1 colocal
                       61),  # both groups colocal
    "colocal_grp2" = c(41, 42, 45, 46, 57, 58, 62,  # group 2 colocal
                       61)  # both groups colocal
  )
  
  index_list <- presets[odds_type]
  index_list
}

## Function to process data and optionally align them
process_data <- function(y, M, X, group = NULL,
                         Z = NULL, Z_y = NULL, Z_M = NULL,
                         w = NULL, w_y = NULL, w_M = NULL,
                         align_data = TRUE,
                         verbose = TRUE) {
  
  # Ensure y is a vector
  if (is.matrix(y)) { y <- y[,1] }
  
  # Ensure X, M, Z, Z_y, and Z_M are matrices
  X <- as.matrix(X)
  M <- as.matrix(M)
  if (!is.null(Z)) { Z <- as.matrix(Z) }
  if (!is.null(Z_y)) { Z_y <- as.matrix(Z_y) }
  if (!is.null(Z_M)) { Z_M <- as.matrix(Z_M) }
  
  # ensure group is binary and numeric
  if (!is.null(group)) {
    unique_groups <- sort(unique(group[!is.na(group)]))
    if (length(unique_groups) != 2) {
      stop("group must contain exactly 2 unique values")
    }
    if (!all(unique_groups %in% c(0,1))) {
      stop("group must be coded as 0/1")
    }
    group <- as.numeric(group) 
  }
  
  # Process covariate matrices
  if (is.null(Z_y)) { Z_y <- matrix(NA, length(y), 0); rownames(Z_y) <- names(y) }
  if (is.null(Z_M)) { Z_M <- matrix(NA, nrow(M), 0); rownames(Z_M) <- rownames(M) }
  
  if (!is.null(Z)) {
    if (align_data) {
      Z <- Z[(rownames(Z) %in% names(y)) & (rownames(Z) %in% rownames(M)), , drop = FALSE]
      Z_y <- cbind(Z, Z_y[rownames(Z),])
      Z_M <- cbind(Z, Z_M[rownames(Z),])
    } else {
      Z_y <- cbind(Z, Z_y)
      Z_M <- cbind(Z, Z_M)
    }
  }
  
  # Process weight vectors
  if (is.null(w)) {
    if (is.null(w_y)) { w_y <- rep(1, length(y)); names(w_y) <- names(y) }
    if (is.null(w_M)) { w_M <- rep(1, nrow(M)); names(w_M) <- rownames(M) }
  } else {
    w_y <- w_M <- w
  }
  
  if (align_data) {
    # M and Z_M can have NAs
    overlapping_samples <- Reduce(f = intersect, x = list(names(y),
                                                          rownames(X),
                                                          rownames(Z_y),
                                                          names(w_y),
                                                          names(group)))
    
    if (length(overlapping_samples) == 0 | !any(overlapping_samples %in% unique(c(rownames(M), rownames(Z_M), names(w_M), names(group))))) {
      stop("No samples overlap. Check rownames of M, X, Z (or Z_y and Z_M) and names of y and w (or w_y and w_M) and group.", call. = FALSE)
    } else if (verbose) {
      writeLines(text = c("Number of overlapping samples:", length(overlapping_samples)))
    }
    
    # Ordering
    y <- y[overlapping_samples]
    M <- M[overlapping_samples,, drop = FALSE]
    X <- X[overlapping_samples,, drop = FALSE]
    Z_y <- Z_y[overlapping_samples,, drop = FALSE]
    Z_M <- Z_M[overlapping_samples,, drop = FALSE]
    w_y <- w_y[overlapping_samples]
    w_M <- w_M[overlapping_samples]
    group <- group[overlapping_samples]
  }
  
  # Drop observations with missing y or X and update n
  complete_y <- !is.na(y)
  complete_X <- !apply(is.na(X), 1, any)
  
  y <- y[complete_y & complete_X]
  M <- M[complete_y & complete_X,, drop = FALSE]
  X <- X[complete_y & complete_X,, drop = FALSE]
  Z_y <- Z_y[complete_y & complete_X,, drop = FALSE]
  Z_M <- Z_M[complete_y & complete_X,, drop = FALSE]
  w_y <- w_y[complete_y & complete_X]
  w_M <- w_M[complete_y & complete_X]
  group <- group[complete_y & complete_X]
  
  # Drop columns of Z_y and Z_M that are invariant
  Z_y_drop <- which(apply(Z_y, 2, function(x) var(x)) == 0)
  Z_M_drop <- which(apply(Z_M, 2, function(x) var(x)) == 0)
  if (length(Z_y_drop) > 0) {
    if (verbose) {
      writeLines(paste("Dropping invariants columns from Z_y:", colnames(Z_y)[Z_y_drop]))
    }
    Z_y <- Z_y[,-Z_y_drop, drop = FALSE]
  }
  if (length(Z_M_drop) > 0) {
    if (verbose) {
      writeLines(paste("Dropping invariants columns from Z_M:", colnames(Z_M)[Z_M_drop]))
    }
    Z_M <- Z_M[,-Z_M_drop, drop = FALSE]
  }
  
  # Return processed data
  list(y = y,
       M = M,
       X = X,
       group = group,
       Z_y = Z_y, Z_M = Z_M,
       w_y = w_y, w_M = w_M)
}

#' Bayesian model selection for mediation analysis function
#'
#' This function takes an outcome (y), candidate mediators (M), and a driver as a design matrix (X) to perform a
#' Bayesian model selection analysis for mediation.
#'
#' @param y Vector or single column matrix of an outcome variable. Single outcome variable expected.
#' Names or rownames must match across M, X, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param M Vector or matrix of mediator variables. Multiple mediator variables are supported.
#' Names or rownames must match across y, X, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs.
#' @param X Design matrix of the driver. Names or rownames must match across y, M, Z, Z_y, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE.
#' If align_data = FALSE, dimensions and order must match across inputs. One common application is for X to represent genetic information at a QTL,
#' either as founder strain haplotypes or variant genotypes, though X is generalizable to other types of variables.
#' @param Z DEFAULT: NULL. Design matrix of covariates that influence the outcome and mediator variables.
#' Names or rownames must match to those of y, M, X, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data=FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_y DEFAULT: NULL. Design matrix of covariates that influence the outcome variable.
#' Names or rownames must match to those of y, M, X, Z_M, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param Z_M DEFAULT: NULL. Design matrix of covariates that influence the mediator variables.
#' Names or rownames must match across y, M, X, Z_y, w, w_y, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. If Z is provided, it supercedes Z_y and Z_M.
#' @param w DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis that applies to both
#' y and M. Names must match across y, M, X, Z, Z_y, and Z_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where w
#' is a vector of the number of individuals per strain. If no w, w_y, or w_M is given, observations are equally weighted as 1s for y and M.
#' If w is provided, it supercedes w_y and w_M.
#' @param w_y DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of y. Names must match across y, M, X, Z, Z_y, Z_M, and w_M (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_y is a vector of the number of individuals per strain used to
#' measure y. If no w_y (or w) is given, observations are equally weighted as 1s for y.
#' @param w_M DEFAULT: NULL. Vector or single column matrix of weights for individuals in analysis, specific to the measurement
#' of M. Names must match across y, M, X, Z, Z_y, Z_M, and w_y (if provided) when align_data = TRUE. If align_data = FALSE,
#' dimensions and order must match across inputs. A common use would be for an analysis of strain means, where y and M
#' are summarized from a different number of individuals per strain. w_M is a vector of the number of individuals per strain use to
#' measure M. If no w_M (or w) is given, observations are equally weighted as 1s for M.
#' @param kappa DEFAULT: c(0.001, 0.001). Shape hyperparameters for the precisions of m and y. The DEFAULT represents uninformative priors on the precisions,
#' in combination with the default prior for lambda.
#' @param lambda DEFAULT: c(0.001, 0.001). Rate hyperparameters for the precisions of m and y. The DEFAULT represents uninformative priors on the precisions,
#' in combination with the default prior for kappa.
#' @param tau_sq_mu DEFAULT: c(1000, 1000). Variance components for the intercepts mu_m and mu_y. The DEFAULT represents diffuse priors, analogous to
#' fixed effect terms.
#' @param tau_sq_Z DEFAULT: c(1000, 1000). Variance components for the covariates encoded in Z_m and Z_y. The DEFAULT represents diffuse priors, analogous
#' to fixed effect terms.
#' @param phi_sq DEFAULT: c(1, 1, 1). Each element of (a, b, c) represents one of the relationships being evaluated for mediation,
#' specifically the ratio of signal to noise. a is the effect of X on M, b is the effect of M on y, and c is the effect of X on y.
#' The DEFAULT represents relationships that explain 50\% of the variation in the outcome variable.
#' @param ln_prior_c DEFAULT: "complete". The prior log case probabilities. See model_info() for description of likelihoods and their
#' combinations into cases. Simplified pre-set options are available, including "complete", "partial", and "reactive".
#' @param options_X DEFAULT: list(sum_to_zero = TRUE, center = FALSE, scale = FALSE). Optional transformations for the X design matrix. Sum_to_zero imposes
#' a sum-to-zero contrast on the columns of X. Center sets the mean of each column to 0. Scale sets the variance of each column to 1.
#' @param align_data DEFAULT: TRUE. If TRUE, expect vector and matrix inputs to have names and rownames, respectively. The overlapping data
#' will then be aligned, allowing the user to not have to reduce data to overlapping samples and order them.
#' @return \code{bmediatR} returns a list containing the following components:
#' \describe{
#'  \item{\code{ln_prob_data}}{a matrix with likelihoods of each hypothesis H1-H8 for each candidate mediator.}
#'  \item{\code{ln_post_c}}{a matrix with posterior probabilities of each causal model for each candidate mediator.}
#'  \item{\code{ln_post_odds}}{a matrix with posterior odds of individual models or combinations of models for each candidate mediator.}
#'  \item{\code{ln_prior_c}}{a single row matrix with posterior probabilities of each causal model.}
#'  \item{\code{ln_prior_odds}}{a single row matrix with prior odds of individual models or combinations of models.}
#'  \item{\code{ln_ml}}{the natural log of the marginal likelihood.}
#' }
#' @note See examples in vignettes with \code{vignette("use_bmediatR", "bmediatR")} or \code{vignette("bmediatR_in_DO", "bmediatR")}.
#' @export
#' @examples bmediatR()
bmediatR_GxT <- function(y, M, X, group,
                         Z = NULL, Z_y = NULL, Z_M = NULL,
                         w = NULL, w_y = NULL, w_M = NULL,
                         kappa = c(0.001, 0.001),
                         lambda = c(0.001, 0.001),
                         tau2_mu = c(1000, 1000),
                         tau2_Z = c(1000, 1000),
                         phi2 = c(1, 1, 1, 1, 1, 1),
                         ln_prior_c = "complete_gxt",
                         options_X = list(sum_to_zero = TRUE, center = FALSE, scale = FALSE),
                         align_data = TRUE,
                         verbose = TRUE) {
  
  # process and optionally align data
  processed_data <- process_data(y = y, M = M, X = X, 
                                 group = group, Z = Z,
                                 Z_y = Z_y, Z_M = Z_M,
                                 w_y = w_y, w_M = w_M,
                                 align_data = align_data,
                                 verbose = verbose)
  y <- processed_data$y
  M <- processed_data$M
  X <- processed_data$X
  group <- processed_data$group
  Z_y <- processed_data$Z_y
  Z_M <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_M <- processed_data$w_M
  
  # dimension of y
  n <- length(y)
  
  # dimension of Z's - add 1 for group
  p_y <- ncol(Z_y) + 1
  p_M <- ncol(Z_M) + 1
  
  # scale y, M, and Z
  y <- c(scale(y))
  M <- apply(M, 2, scale)
  if (p_y > 0) { Z_y <- apply(Z_y, 2, scale) }
  if (p_M > 0) { Z_M <- apply(Z_M, 2, scale) }
  
  # optionally use sum-to-zero contrast for X
  # recommended when X is a matrix of factors, with a column for every factor level
  if (options_X$sum_to_zero == TRUE) {
    if (ncol(X)<2){
      stop("Cannot use sum-to-zero contrast when X has one column")
    } else {
      C <- sumtozero_contrast(ncol(X))
      X <- X%*%C
    }
  }
  
  # optionally center and scale X
  X <- apply(X, 2, scale, center = options_X$center, scale = options_X$scale)
  
  # dimension of X
  d <- ncol(X)
  
  # column design matrix for mu
  ones <- matrix(1, n)
  colnames(ones) <- "intercept"
  
  # begin Bayesian calculations
  if (verbose) { print("Initializing", quote = FALSE) }
  
  # 4 M models and 16 y models = 64 unique combinations 
  y_models <- expand.grid(b1=0:1, b2=0:1, c1=0:1, c2=0:1)
  M_models <- expand.grid(a1 = 0:1, a2 = 0:1)
  
  # likelihood models for all hypothesis
  # hypotheses encoded by presence (1) or absence (0) of 'x->m, m->y, x->y' edges on the DAG
  # H1:  '-,-,0,0,0,0' / y does not depend on X or m 
  # H2:  '-,-,1,0,0,0' / y depends on m in group 1 
  # H3:  '-,-,0,1,0,0' / y depends on m in group 2
  # H4:  '-,-,1,1,0,0' / y depends on m in groups 1 and 2 
  # H5:  '-,-,0,0,1,0' / y depends on X in group 1 
  # H6:  '-,-,1,0,1,0' / y depends on m and X in group 1
  # H7:  '-,-,0,1,1,0' / y depends on m in group 2 and X in group 1
  # H8:  '-,-,1,1,1,0' / y depends on m in groups 1 and 2 and X in group 1
  # H9:  '-,-,0,0,0,1' / y depends on X in group 2
  # H10: '-,-,1,0,0,1' / y depends on m in group 1 and X in group 2
  # H11: '-,-,0,1,0,1' / y depends on m and X in group 2
  # H12: '-,-,1,1,0,1' / y depends on m in groups 1 and 2 and X in group 2
  # H13: '-,-,0,0,1,1' / y depends on X in groups 1 and 2
  # H14: '-,-,1,0,1,1' / y depends on m in group 1 and X in groups 1 and 2
  # H15: '-,-,0,1,1,1' / y depends on m in group 2 and X in groups 1 and 2
  # H16: '-,-,1,1,1,1' / y depends on m and X in groups 1 and 2
  # H17: '0,0,-,-,-,-' / m does not depend on X  
  # H18: '1,0,-,-,-,-' / m depends on X in group 1
  # H19: '0,1,-,-,-,-' / m depends on X in group 2
  # H20: '1,1,-,-,-,-' / m depends on X in both groups 
  # all include covariates Z
  
  # iterate over batches of M with same pattern of missing values
  if (verbose) { print("Iterating", quote = FALSE) }
  counter <- 0
  
  # object to store likelihoods
  ln_prob_data <- matrix(NA, nrow = ncol(M), ncol = nrow(M_models)*nrow(y_models)) 
  rownames(ln_prob_data) <- colnames(M)
  y_labels <- paste0("H", 1:16)
  M_labels <- paste0("H", 17:20)
  colnames(ln_prob_data) <- paste(rep(M_labels, each = 16), "x", rep(y_labels, times = 4))
  
  # for each mediator
  for (i in 1:ncol(M)) {
    counter <- counter + 1
    if (counter%%1000==0 & verbose) { print(paste(counter, "of", ncol(M)), quote=F) }
    
    m <- M[,i,drop = FALSE]
    
    # compute likelihoods for 4 mediator models 
    ln_prob_M <- numeric(nrow(M_models))
    for (m_idx in 1:nrow(M_models)) {
      X_m <- cbind(ones, group, Z_M) # base design matrix 
      v_m <- c(tau2_mu[2], rep(tau2_Z[2], p_M)) # base prior variance matrix (diagonal)
      # build X_m and v_m based on M_models[m_idx,]
      if (M_models$a1[m_idx] == 1) { # add X effect for samples in group 1
        X_grp1 <- X*(group == 0)
        colnames(X_grp1) <- "X_grp1"
        X_m <- cbind(X_m, X_grp1)  
        v_m <- c(v_m, phi2[1]) # phi2[a1]
      }
      if (M_models$a2[m_idx] == 1) { # add X effect for samples in group 2 
        X_grp2 <- X*(group == 1)
        colnames(X_grp2) <- "X_grp2"
        X_m <- cbind(X_m, X_grp2) 
        v_m <- c(v_m, phi2[2]) # phi2[a2]
      }
      
      # build scale matrix 
      sigma <- crossprod(sqrt(lambda[2]*v_m) * t(X_m)) 
      diag(sigma) <- diag(sigma) + lambda[2]/w_M
      sigma_chol <- chol(sigma)
      
      # compute likelihood
      ln_prob_M[m_idx] <- bmediatR:::dmvt_chol(m, sigma_chol = sigma_chol, df = kappa[2])
    }
    
    # compute likelihoods for 16 outcome models 
    ln_prob_y <- numeric(nrow(y_models))
    for (y_idx in 1:nrow(y_models)) {
      X_y <- cbind(ones, group, Z_y) # base design matrix 
      v_y <- c(tau2_mu[1], rep(tau2_Z[1], p_y)) # base prior variance matrix (diagonal)
      # build X_y and v_y based on y_models[y_idx,]
      if (y_models$b1[y_idx] == 1) { # add m effect for samples in group 1
        m_grp1 <- m*(group == 0)
        colnames(m_grp1) <- "m_grp1"
        X_y <- cbind(X_y, m_grp1)  
        v_y <- c(v_y, phi2[3]) # phi2[3] = phi2_b1
      }
      if (y_models$b2[y_idx] == 1) { # add m effect for samples in group 2 
        m_grp2 <- m*(group == 1)
        colnames(m_grp2) <- "m_grp2"
        X_y <- cbind(X_y, m_grp2) 
        v_y <- c(v_y, phi2[4]) # phi2[4] = phi2_b2
      }
      if (y_models$c1[y_idx] == 1) { # add X effect for samples in group 1
        X_grp1 <- X*(group == 0)
        colnames(X_grp1) <- "X_grp1"
        X_y <- cbind(X_y, X_grp1)  
        v_y <- c(v_y, phi2[5]) # phi2[5] = phi2_c1
      }
      if (y_models$c2[y_idx] == 1) { # add X effect for samples in group 2 
        X_grp2 <- X*(group == 1)
        colnames(X_grp2) <- "X_grp2"
        X_y <- cbind(X_y, X_grp2) 
        v_y <- c(v_y, phi2[6]) # phi2[6] = phi2_c2
      }
      
      # build scale matrix 
      sigma <- crossprod(sqrt(lambda[1]*v_y) * t(X_y)) 
      diag(sigma) <- diag(sigma) + lambda[1]/w_y
      sigma_chol <- chol(sigma)
      
      # compute likelihood
      ln_prob_y[y_idx] <- bmediatR:::dmvt_chol(y, sigma_chol = sigma_chol, df = kappa[1])
    }
    
    # combine into 64 joint likelihoods (4 × 16 = 64)
    model_ix <- 0
    for (m_idx in 1:nrow(M_models)){
      for (y_idx in 1:nrow(y_models)){
        model_ix <- model_ix + 1
        ln_prob_data[i,model_ix] <- ln_prob_M[m_idx] + ln_prob_y[y_idx]
      }
    }
  }
  
  
  # compute posterior probabilities for all cases
  # compute posterior odds for specified combinations of cases
  # cases encoded by presence (1) or absence (0) of 'X->m, X->m, X->y' edges on the DAG
  # c1:  '0,0,0,0,0,0' = H17 x H1
  # c2:  '0,0,1,0,0,0' = H17 x H2 
  # c3:  '0,0,0,1,0,0' = H17 x H3
  # c4:  '0,0,1,1,0,0' = H17 x H4 
  # c5:  '0,0,0,0,1,0' = H17 x H5
  # c6:  '0,0,1,0,1,0' = H17 x H6 
  # c7:  '0,0,0,1,1,0' = H17 x H7
  # c8:  '0,0,1,1,1,0' = H17 x H8
  # c9:  '0,0,0,0,0,1' = H17 x H9
  # c10: '0,0,1,0,0,1' = H17 x H10 
  # c11: '0,0,0,1,0,1' = H17 x H11
  # c12: '0,0,1,1,0,1' = H17 x H12
  # c13: '0,0,0,0,1,1' = H17 x H13
  # c14: '0,0,1,0,1,1' = H17 x H14
  # c15: '0,0,0,1,1,1' = H17 x H15
  # c16: '0,0,1,1,1,1' = H17 x H16
  # c17: '1,0,0,0,0,0' = H18 x H1
  # c18: '1,0,1,0,0,0' = H18 x H2 - complete mediation, group 1
  # c19: '1,0,0,1,0,0' = H18 x H3
  # c20: '1,0,1,1,0,0' = H18 x H4 - complete mediation, group 1
  # c21: '1,0,0,0,1,0' = H18 x H5 - colocalization, group 1
  # c22: '1,0,1,0,1,0' = H18 x H6 - partial mediation, group 1
  # c23: '1,0,0,1,1,0' = H18 x H7 - colocalization, group 1 
  # c24: '1,0,1,1,1,0' = H18 x H8 - partial mediation, group 1
  # c25: '1,0,0,0,0,1' = H18 x H9 
  # c26: '1,0,1,0,0,1' = H18 x H10 - complete mediation, group 1
  # c27: '1,0,0,1,0,1' = H18 x H11 
  # c28: '1,0,1,1,0,1' = H18 x H12 - complete mediation, group 1
  # c29: '1,0,0,0,1,1' = H18 x H13 - colocalization, group 1
  # c30: '1,0,1,0,1,1' = H18 x H14 - partial mediation, group 1 
  # c31: '1,0,0,1,1,1' = H18 x H15 - colocalization, group 1
  # c32: '1,0,1,1,1,1' = H18 x H16 - partial mediation, group 1
  # c33: '0,1,0,0,0,0' = H19 x H1
  # c34: '0,1,1,0,0,0' = H19 x H2  
  # c35: '0,1,0,1,0,0' = H19 x H3 - complete mediation, group 2
  # c36: '0,1,1,1,0,0' = H19 x H4 - complete mediation, group 2
  # c37: '0,1,0,0,1,0' = H19 x H5 
  # c38: '0,1,1,0,1,0' = H19 x H6 
  # c39: '0,1,0,1,1,0' = H19 x H7 - complete mediation, group 2
  # c40: '0,1,1,1,1,0' = H19 x H8 - complete mediation, group 2
  # c41: '0,1,0,0,0,1' = H19 x H9 - colocalization, group 2
  # c42: '0,1,1,0,0,1' = H19 x H10 - colocalization, group 2
  # c43: '0,1,0,1,0,1' = H19 x H11 - partial mediation, group 2
  # c44: '0,1,1,1,0,1' = H19 x H12 - partial mediation, group 2
  # c45: '0,1,0,0,1,1' = H19 x H13 - colocalization, group 2
  # c46: '0,1,1,0,1,1' = H19 x H14 - colocalization, group 2
  # c47: '0,1,0,1,1,1' = H19 x H15 - partial mediation, group 2
  # c48: '0,1,1,1,1,1' = H19 x H16 - partial mediation, group 2
  # c49: '1,1,0,0,0,0' = H20 x H1
  # c50: '1,1,1,0,0,0' = H20 x H2 - complete mediation, group 1
  # c51: '1,1,0,1,0,0' = H20 x H3 - complete mediation, group 2
  # c52: '1,1,1,1,0,0' = H20 x H4 - complete mediation, both groups 
  # c53: '1,1,0,0,1,0' = H20 x H5 - colocalization, group 1
  # c54: '1,1,1,0,1,0' = H20 x H6 - partial mediation, group 1
  # c55: '1,1,0,1,1,0' = H20 x H7 - colocalization, group 1; complete mediation, group 2
  # c56: '1,1,1,1,1,0' = H20 x H8 - partial mediation, group 1; complete mediation, group 2
  # c57: '1,1,0,0,0,1' = H20 x H9 - colocalization, group 2
  # c58: '1,1,1,0,0,1' = H20 x H10 - complete mediation, group 1; colocalization, group 2
  # c59: '1,1,0,1,0,1' = H20 x H11 - partial mediation, group 2
  # c60: '1,1,1,1,0,1' = H20 x H12 - complete mediation, group 1; partial mediation, group 2
  # c61: '1,1,0,0,1,1' = H20 x H13 - colocalization, both groups 
  # c62: '1,1,1,0,1,1' = H20 x H14 - partial mediation, group 1; colocalization, group 2 
  # c63: '1,1,0,1,1,1' = H20 x H15 - colocalization, group 1; partial mediation, group 2
  # c64: '1,1,1,1,1,1' = H20 x H16 - partial mediation, both groups 
  preset_odds_index <- return_preset_odds_index_gxt()
  output <- posterior_summary_gxt(ln_prob_data, ln_prior_c, preset_odds_index)
  colnames(output$ln_post_odds) <- colnames(output$ln_prior_odds) <- colnames(output$ln_post_odds)
  colnames(output$ln_prior_c) <- colnames(output$ln_post_c)
  rownames(output$ln_prior_c) <- rownames(output$ln_post_c)
  
  if (verbose) { print("Done", quote = FALSE) }
  return(output)
}