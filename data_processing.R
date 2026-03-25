library(readxl)
library(qtl)
library(car)
library(lme4qtl)
library(AGHmatrix)
source("code-dependencies/cov_qtl_functions.R")

# -------------------------process phenotype data----------------------------- #
# load original cross_data in as regular df
cross_data <- load_cross_as_df("derived_data/RqtlCC006xCC044_ctrlAndSARS.csv", n_geno_start = 122)
# get genotypes for GRM estimation 
geno <- cross_data[,c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]

# define flow columns
pheno_names <- read_xlsx("source_data/pheno_names.xlsx")
flow_cols <- pheno_names$flow_col_name[-1] # remove titer, will transform separately 

# save raw data 
flow_untr <- cross_data[,flow_cols]
write_csv(flow_untr, "derived_data/data_processing/raw_immune_data.csv")

# add weight area above the curve phenotype
cross_data$pd0 <- rep(100, nrow(cross_data))
cross_data$pd1 <- cross_data$d1/cross_data$d0*100
cross_data$pd2 <- cross_data$d2/cross_data$d0*100
cross_data$pd3 <- cross_data$d3/cross_data$d0*100
cross_data$pd4 <- cross_data$d4/cross_data$d0*100
cross_data <- calc_auc(cross_data,
                       steps = c(0,1,2,3,4),
                       col.name = "weight_aac",
                       phenos = c("pd0", "pd1", "pd2", "pd3", "pd4"),
                       is_df = TRUE)
# mice_w_NAweights <- is.na(cross_data$weight_aac)
# cross_data <- cross_data[!mice_w_NAweights,] # four mice

cross_data$sex <- as.numeric(cross_data$sex == "M")
cross_data$Titer <- log(cross_data$Titer + 1) # log(x+1) transform Titer data

# create GRM
grm <- Gmatrix(as.matrix(geno), method = "VanRaden")
colnames(grm) <- rownames(grm) <- cross_data$mouse_ID

# logit-transform, deBLUP batch effects, and center/scale flow data
for (flow_col in flow_cols){
  tr_pheno <- car::logit(cross_data[[flow_col]])

  # fit <- lmer(tr_pheno ~ sex + infection + sex:infection + (1|flow_batch), data = cross_data)
  fit <- relmatLmer(tr_pheno ~ sex + infection + sex:infection + (1|flow_batch) + (1|mouse_ID),
                    data = cross_data, relmat = list(mouse_ID = grm)) # LMM including polygenic term
  blup <- ranef(fit)$flow_batch[as.character(cross_data$flow_batch), "(Intercept)"]
  deblupd <- tr_pheno - blup
  
  cross_data[,flow_col] <- scale(deblupd, center = TRUE, scale = TRUE)
}

write_csv(cross_data, "derived_data/cross_data.csv") # save to csv - use this for Gibbs, SynSurr, etc. 


# --------------------------impute genotype data------------------------------ #
f2cross <- read.cross(format="csv",
                      file="derived_data/RqtlCC006xCC044_ctrlAndSARS.csv",
                      na.strings=c("-","NA", "na", "<NA>", "not tested"), 
                      genotypes=c("AA","AB","BB"),
                      alleles=c("A","B")) # A = CC006, B = CC044
f2cross <- fill.geno(f2cross, method = "argmax") # impute genotypes 
#f2cross <- subset(f2cross, ind = !mice_w_NAweights)
write.cross(f2cross, format = "csv", filestem = "derived_data/cross_imputed")


