library(qtl)
library(readr)
library(readxl)
library(AGHmatrix)
library(lme4)
library(lme4qtl)
library(r2glmm)
library(parallel)
library(doParallel)
library(foreach)
source("code-dependencies/cmdline.R")

run_mode <- cmdline.option("mode", default = "local", allowed.values = c("local", "slurm"))

# ---------------------------------Load data---------------------------------- #
# raw data
cross <- read.cross(format = "csv", 
                    file = "derived_data/RqtlCC006xCC044_ctrlAndSARS.csv",
                    na.strings=c("-","NA", "na", "<NA>", "not tested"), 
                    genotypes=c("AA","AB","BB"),
                    alleles=c("A","B")) # A = CC006, B = CC044
cross <- jittermap(cross)

# define flow cols 
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
all_phenos <- pheno_names$flow_col_name
p <- length(all_phenos)
flow_cols <- all_phenos[-1] # remove titer 
q <- length(flow_cols)

# normalize / transform, but don't deBLUP batch effects (for variance components analysis)
cross$pheno$Titer <- log(cross$pheno$Titer+1)
cross$pheno$batch <- paste0("batch", cross$pheno$batch)
cross$pheno$sex <- as.numeric(cross$pheno$sex == "M")
for (flow_col in flow_cols){
  cross$pheno[[flow_col]] <- car::logit(cross$pheno[[flow_col]])
}

# load genotypes 
cross_data <- read_csv("derived_data/cross_data.csv")
geno <- cross_data[,c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]

# define genomic relationship matrix 
grm <- Gmatrix(as.matrix(geno), method = "VanRaden")
colnames(grm) <- rownames(grm) <- cross_data$mouse_ID

set.seed(123)

if (run_mode == "slurm"){
  ncores <- q
  cl <- parallel::makeForkCluster(ncores)
} else {
  ncores <- parallel::detectCores() - 4 # if local, use available cores 
  cl <- makeCluster(ncores)
}
doParallel::registerDoParallel(cl)

message(paste("estimating variance components for immune phenotypes in parallel using", ncores, "cores."))
varcomp_data <- foreach(i = 2:p, .combine = rbind, .packages = c("lme4qtl", "r2glmm")) %dopar% {
  pheno_name <- all_phenos[i]
  form <- as.formula(paste0(pheno_name, " ~ sex + infection + sex:infection + (1|flow_batch) + (1|mouse_ID)"))
  mod <- relmatLmer(form, data = cross$pheno, relmat = list(mouse_ID = grm))
  
  # -------------------estimate random effect components-------------------- #
  vp <- VarProp(mod)
  h2 <- vp[vp$grp == "mouse_ID", "prop"]
  batch_vp <- vp[vp$grp == "flow_batch", "prop"] 
  
  # --------------est confidence intervals for ranef varcomps--------------- #
  # try using varpropProf(profile()), which is more theoretically correct and uses
  # proper likelihood-based profiling for heritability, but fails when h2 = 0 or 1 
  prof_prop <- try(varpropProf(profile(mod)), silent = TRUE)
  if (inherits(prof_prop, "try-error")){ # varpropProf won't work, use bootstrapping 
    print(paste("using bootstrap for", pheno_name))
    # estimate confidence intervals 
    vc_fun <- function(fit) {
      vc <- as.data.frame(VarCorr(fit))
      sigma2  <- sigma(fit)^2
      sigma2G <- vc$vcov[vc$grp == "mouse_ID"]
      sigma2F <- vc$vcov[vc$grp == "flow_batch"]
      h2 <- sigma2G / (sigma2G + sigma2F + sigma2)
      batch_vp <- sigma2F / (sigma2G + sigma2F + sigma2)
      return(c(h2 = h2, batch_vp = batch_vp))
    }
    boot <- bootMer(mod, FUN = vc_fun, nsim = 1000, use.u = FALSE, type = "parametric", parallel = "multicore", ncpus = 5)
    ci <- quantile(boot$t[,"h2"], probs = c(0.025, 0.975))
    h2_data <- c(h2, ci[1], ci[2])
    ci <- quantile(boot$t[,"batch_vp"], probs = c(0.025, 0.975))
    batch_vp_data <- c(batch_vp, ci[1], ci[2])
  } else {
    vc <- as.data.frame(VarCorr(mod))
    ci <- confint(prof_prop)
    ix <- which(vc$grp == "mouse_ID")
    h2_data <- c(h2, ci[ix,1], ci[ix,2])
    ix <- which(vc$grp == "flow_batch")
    batch_vp_data <- c(batch_vp, ci[ix,1], ci[ix,2])
  }
  
  # --------------------estimate fixed effect components-------------------- #
  # use Nakagawa and Schielzeth method to estimate variance explained by fixed effects
  fixed_r2 <- as.data.frame(r2beta(mod, method = "nsj"))
  sex_vp_data <- fixed_r2[fixed_r2$Effect == "sex", c("Rsq", "lower.CL", "upper.CL")]
  inf_vp_data <- fixed_r2[fixed_r2$Effect == "infection", c("Rsq", "lower.CL", "upper.CL")]
  int_vp_data <- fixed_r2[fixed_r2$Effect == "sex:infection", c("Rsq", "lower.CL", "upper.CL")]
  
  return(c(pheno_name, sex_vp_data, inf_vp_data, int_vp_data, h2_data, batch_vp_data))
}
stopCluster(cl)

colnames(varcomp_data) <- c("pheno", 
                            "sex_vp", "sex_vp_lwr", "sex_vp_upr",
                            "inf_vp", "inf_vp_lwr", "inf_vp_upr",
                            "int_vp", "int_vp_lwr", "int_vp_upr",
                            "h2", "h2_lwr", "h2_upr",
                            "batch_vp", "batch_vp_lwr", "batch_vp_upr")

message(paste("estimating variance components for titer"))
mod <- relmatLmer(Titer ~ sex + infection + sex:infection + (1|batch) + (1|mouse_ID), 
                  data = cross$pheno, relmat = list(mouse_ID = grm))
vp <- VarProp(mod)
h2 <- vp[vp$grp == "mouse_ID", "prop"]
batch_vp <- vp[vp$grp == "batch", "prop"]
prof_prop <- varpropProf(profile(mod))
vc <- as.data.frame(VarCorr(mod))
ci <- confint(prof_prop)
ix <- which(vc$grp == "mouse_ID")
h2_data <- c(h2, ci[ix,1], ci[ix,2])
ix <- which(vc$grp == "batch")
batch_vp_data <- c(batch_vp, ci[ix,1], ci[ix,2])
fixed_r2 <- as.data.frame(r2beta(mod, method = "nsj"))
sex_vp_data <- fixed_r2[fixed_r2$Effect == "sex", c("Rsq", "lower.CL", "upper.CL")]
inf_vp_data <- fixed_r2[fixed_r2$Effect == "infection", c("Rsq", "lower.CL", "upper.CL")]
int_vp_data <- fixed_r2[fixed_r2$Effect == "sex:infection", c("Rsq", "lower.CL", "upper.CL")]
varcomp_data <- rbind(varcomp_data, c("Titer", sex_vp_data, inf_vp_data, int_vp_data, h2_data, batch_vp_data))

write_csv(as.data.frame(varcomp_data), file = "results/var_comp_res.csv")


