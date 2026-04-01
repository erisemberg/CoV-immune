library(readr)
library(readxl)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)
library(qtl)
library(nlme)
library(lme4qtl)
source("code-dependencies/cov_qtl_functions.R")
source("code-dependencies/cmdline.R")

run_mode <- cmdline.option("mode", default = "local", allowed.values = c("local", "slurm"))
# if running via slurm, set working directory to HPC directory
if (run_mode == "slurm"){
  hpc_dir <- cmdline.string("hpc_dir", default = "/work/users/e/r/erisembe/Coronavirus/cov-immune")
  setwd(hpc_dir)
}

# ---------------------------------functions---------------------------------- #
# defining locally so that this script can be standalone / without dependencies 
# in case being run on the cluster 

# define row weights (precision weights): w_i = 1 / s2_g(i) 
make_row_weights <- function(covar, var_by_group, mask_rows) {
  g <- as.character(covar$trt[mask_rows])
  w <- 1 / var_by_group[g]
  return(w)
}

# permutation test for 
univar_perm_gxt <- function(cross, target, covar, num_perms = 1000, var_by_group){
  n <- length(target) 
  max_logPs <- rep(NA, num_perms)
  
  for (i in 1:num_perms){
    samp <- sample(x = 1:n, size = n)
    # define weights
    covar_i <- covar[samp,]
    y_i <- target[samp]
    has_y <- !is.na(y_i)
    wts <- make_row_weights(covar_i, var_by_group, has_y)
    
    res <- try(univar_scan(cross, target = y_i, covar = covar_i, gxt = TRUE, wts = wts), silent = TRUE)
    if (inherits(res, "try-error")) { message("Skipping permutation ", i); next }
    
    # get max log P per marker (between p_overall, p_G, and p_GxT)
    max_logPs_permarker <- apply(res[,3:5], 1, max, na.rm = TRUE)
    # max of max log Ps
    max_logPs[i] <- max(max_logPs_permarker, na.rm = TRUE) 
  }
  
  max_logPs <- max_logPs[!is.na(max_logPs)] # remove NAs
  # fit GEV distribution 
  fitgev <- fevd(max_logPs, type = "GEV")
  fitgev <- add_attribute(fitgev, "class", c("permgev", "fevd"))
  return(fitgev)
}

# ---------------------------------load data---------------------------------- #
# raw phenotype data, covariates, genotype data
cross <- read.cross(format="csv", 
                    file="derived_data/RqtlCC006xCC044_ctrlAndSARS.csv",
                    na.strings=c("-","NA", "na", "<NA>", "not tested"), 
                    genotypes=c("AA","AB","BB"),
                    alleles=c("A","B")) # A = CC006, B = CC044
cross <- jittermap(cross)

# processed phenotype data: logit-transformed, de-BLUPd batch effects (based on
# fully specified model with polygenic term), and centered/scaled 
cross_data <- read_csv("derived_data/cross_data.csv", show_col_types = FALSE)
geno <- cross_data[,c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]

# covar <- data.frame(trt = cross_data$infection,
#                     sex = cross_data$sex,
#                     pgm = cross_data$pgm)
covar <- data.frame(sex = as.numeric(pull.pheno(cross, "sex") == "M"),
                    trt = pull.pheno(cross, "infection"),
                    pgm = pull.pheno(cross, "pgm"))

# define flow cols 
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
flow_cols <- pheno_names$flow_col_name[-1] # remove titer 
q <- length(flow_cols)

groups <- c("PBS", "SARSCoV", "SARS2CoV", "GxT")
group_dirs <- c("PBS" = "PBS", "SARSCoV" = "SARS", "SARS2CoV" = "SARS2", "GxT" = "GxT")

mod_rds_dir <- "/results/qtl_mapping/modRDS/"
perm_rds_dir <- "/results/qtl_mapping/permRDS/"
scan_dir <- "/results/qtl_mapping/scans/"

num_perms <- 1000

# ----------------------------immune QTL mapping------------------------------ #
if (run_mode == "slurm"){
  ncores <- q
  cl <- parallel::makeForkCluster(ncores)
} else {
  ncores <- parallel::detectCores() - 4 # if local, use available cores 
  cl <- makeCluster(ncores)
}
doParallel::registerDoParallel(cl)

foreach(i = 1:q) %dopar% {
  pheno_name <- flow_cols[i]
  for (group in groups){
    group_dir <- group_dirs[[group]]
    if (group %in% c("PBS", "SARSCoV", "SARS2CoV")){
      pheno <- cross_data[[pheno_name]][cross_data$infection == group]
      cross_g <- subset(cross, ind = (cross$pheno$infection == group))
      covar_g <- covar[covar$trt == group,]
      gxt <- FALSE
      wts <- rep(1, nrow(covar_g))
    } else {
      pheno <- cross_data[[pheno_name]]
      cross_g <- cross
      covar_g <- covar 
      gxt <- TRUE
      
      # define weights 
      has_pheno <- !is.na(pheno)
      covar_i <- covar[has_pheno, , drop = FALSE]
      covar_i$pheno <- pheno[has_pheno]
      
      fit_w <- lm(pheno ~ sex + trt, data = covar_i)
      res <- resid(fit_w)
      s2 <- tapply(res, covar_i$trt, var)
      wts <- 1 / s2[as.character(covar_i$trt)]
    }
    
    mod <- univar_scan(cross_g, pheno, covar_g, gxt = gxt, wts = wts) # QTL scan (univariate) 
    saveRDS(mod, paste0(mod_rds_dir, group_dir, "/", pheno_name, ".rds")) # save mod to RDS 
    
    # permutation test 
    if (group == "GxT"){
      perm <- univar_perm_gxt(cross_g, pheno, covar_g, num_perms = num_perms, var_by_group = s2)
    } else {
      perm <- univar_perm(cross_g, pheno, covar_g, num_perms) # permutation test (univariate)
    }
    saveRDS(perm, paste0(perm_rds_dir, group_dir, "/", pheno_name, ".rds")) # save perm to RDS
    
    # save univariate genome scan to png  
    lodcols <- if(group == "GxT") { c(1:3) } else { 1 }
    png(paste0(scan_dir, group_dir, "/", pheno_name, ".png"), width = 750)
    plot(mod, ylab = "", xlab = "", lodcolumn = lodcols,
         main = paste0(pheno_name, " (", group_dir, "), v2"), 
         alternate.chrid = T, 
         cex.main = 2, cex.axis = 2)
    title(ylab = "-log(P)", line = 2.5, cex.lab = 2)
    title(xlab = "Chromosome", cex.lab = 2.3, line = 3)
    abline(h = c(summary(perm)[1], summary(perm)[2]), lty=1:2)
    dev.off()
  }
}

parallel::stopCluster(cl)

# -----------------------------titer QTL mapping------------------------------ #
pheno_name <- "Titer"
cross$pheno$Titer <- log(cross$pheno$Titer+1)
groups <- c("SARSCoV", "SARS2CoV")
for (group in groups){
  group_dir <- group_dirs[[group]]
  cross_g <- subset(cross, ind = (cross$pheno$infection == group))

  mod <- scanone(cross_g, pheno.col = "Titer", model = "2part")
  saveRDS(mod, paste0(mod_rds_dir, group_dir, "/Titer.rds"))
  perm <- scanone(cross_g, pheno.col = "Titer", model = "2part", n.perm = 1000, n.cluster = 4)
  saveRDS(perm, paste0(perm_rds_dir, group_dir, "/Titer.rds"))
  
  png(paste0(scan_dir, group_dir, "/", pheno_name, ".png"), width = 750)
  plot(mod, lodcolumn = 1:3, main = "Viral titer (SARS-CoV)", alternate.chrid = TRUE)
  abline(h = summary(perm)[1,], lty = 2, col = c("black", "blue", "red"))
  dev.off()
}

