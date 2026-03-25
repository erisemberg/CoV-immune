library(tidyverse)
library(readxl) 
library(lme4)

# load all lymphoid datasets (see Summary tab in spreadsheet for details) 
lymph20 <- read_xlsx('source_data/lymphoid-combined.xlsx', sheet = '082520-LYMPHOID-KJ-forR', 
                     na = c('not measured', 'not assayed', 'not tested', "")) 
lymph22 <- read_xlsx('source_data/lymphoid-combined.xlsx', sheet = '010322-LYMPHOID-KJ-forR', 
                     na = c('not measured', 'not assayed', 'not tested', ""))
lymphkk <- read_xlsx('source_data/lymphoid-combined.xlsx', sheet = 'KK-combined')

# function to split up DATE_SAMPLE column into date and sample number columns 
split_date_sample <- function(df, date_sample_col){
  date_sample <- df[,date_sample_col]
  df$DATE <- str_split_fixed(date_sample[[1]], "_", 2)[,1]
  df$DATE <- as.factor(df$DATE)
  df$SAMPLE_ID <-str_split_fixed(date_sample[[1]], "_", 2)[,2]
  return(df)
}

lymph20 <- split_date_sample(lymph20, "DATE_SAMPLE")
lymph22 <- split_date_sample(lymph22, "DATE_SAMPLE")
lymphkk <- split_date_sample(lymphkk, "DATE_SAMPLE")

# define lymphoid phenotypes 
phenotypes <- read_delim(file = "source_data/phenotype_names.txt", delim = "\n", col_names = FALSE)[[1]]
phenotypes <- phenotypes[startsWith(phenotypes, "L_")]

# add empty columns that don't exist in each spreadsheet so we have a common set 
# of columns between datasets
for (i in 1:length(phenotypes)){
  pheno <- phenotypes[i]
  if (!(pheno %in% colnames(lymph20))){
    lymph20[,pheno] <- rep(NA, nrow(lymph20))
  }
  if (!(pheno %in% colnames(lymph22))){
    lymph22[,pheno] <- rep(NA, nrow(lymph22))
  }
  if (!(pheno %in% colnames(lymphkk))){
    lymphkk[,pheno] <- rep(NA, nrow(lymphkk))
  }
}

# Select relevant columns, all in same order, and add `Gater` column (KJ1 = Kara's 
# first dataset, KJ2 = Kara's second dataset, KK = Kalika's data) and `Gater_date` 
# column (for random effect modeling). 
lymph20 <- lymph20 %>% dplyr::select(c('ORIGINAL_FILE_NAME', 'DATE_SAMPLE', 'DATE', 
                                       'SAMPLE_ID', 'NOTES', any_of(phenotypes)))
lymph20$Gater <- 'KJ1'
lymph20$Gater_date <- paste0(lymph20$Gater, lymph20$DATE)

lymph22 <- lymph22 %>% dplyr::select(c('ORIGINAL_FILE_NAME', 'DATE_SAMPLE', 'DATE', 
                                       'SAMPLE_ID', 'NOTES', any_of(phenotypes)))
lymph22$Gater <- 'KJ2'
lymph22$Gater_date <- paste0(lymph22$Gater, lymph22$DATE)

lymphkk <- lymphkk %>% dplyr::select(c('ORIGINAL_FILE_NAME', 'DATE_SAMPLE', 'DATE', 
                                       'SAMPLE_ID', 'NOTES', any_of(phenotypes)))
lymphkk$Gater <- 'KK'
lymphkk$Gater_date <- paste0(lymphkk$Gater, lymphkk$DATE)

# now that all datasets have the same columns, rbind 
all_lymph_og <- rbind(lymph20, lymph22, lymphkk) # un-corrected data
all_lymph_og$Gater <- as.factor(all_lymph_og$Gater)
all_lymph <- all_lymph_og # to store residuals

# for EDA, correct phenotypes for gater and date
gater_date <- all_lymph_og$Gater_date
for (pheno in phenotypes){
  y <- all_lymph_og[,pheno][[1]]
  
  mod <- lmer(y ~ (1|gater_date), na.action = na.exclude)
  mu <- coef(summary(mod))['(Intercept)','Estimate']
  y_hat <- mu + resid(mod)
  
  all_lymph[,pheno] <- y_hat
}

# protocol: starting with raw data (all_lymph_og):
#   1. For measurements with data from KJ1 and KJ2, correct for slope between the 

# plan: correct for slope (KJ1 vs KJ2 or KJ1 vs KK) in one dataset, and not correct in the other. 
# In the corrected dataset, I will correct the KJ2 data for the KJ1-KJ2 slope, and correct the KK data for the KJ1-KK slope, 
# effectively projecting all of the data onto the KJ1 space.

#In both datasets, I will average the phenotypes that have >1 measurements. 

# NOTE: Do not use these residuals when analyzing SARS-CoV-2 infection data. 
# Since date is completely confounded with SARS-CoV-2 infection 

# Harmonization 
flowIDs <- unique(all_lymph$DATE_SAMPLE)

pheno3 <- data.frame(matrix(NA, nrow = length(flowIDs), ncol = length(phenotypes)))
rownames(pheno3) <- flowIDs
names(pheno3) <- phenotypes
pheno3$Gater <- rep(NA, nrow = length(flowIDs))

# loop over each phenotype 
for (phenotype in phenotypes){
  df <- all_lymph_og[,c('DATE_SAMPLE', phenotype, 'Gater')] # filter relevant data
  widedf <- pivot_wider(df, id_cols = 'DATE_SAMPLE', 'names_from' = 'Gater', values_from = phenotype) # pivot wider 
  overlaps <- rowSums(!is.na(widedf[,2:3])) # count rows with data in both KJ2 and KJ2 
  
  if(sum(overlaps > 1) > 0){ # there is data overlapping between KJ1 and KJ2, calculate intercept and slope
    overlap_df <- widedf[which(overlaps > 1), c('KJ1', 'KJ2')] # subset relevant data
    mod <- lm(KJ2 ~ KJ1, data = overlap_df) # Linear model KJ1 vs KJ2
    mu <- summary(mod)$coefficients['(Intercept)', 'Estimate'] # intercept
    b <- summary(mod)$coefficients['KJ1', 'Estimate'] # slope
  }
  
  # Loop over samples (DATE_SAMPLE) and determine strategy based on presence of overlapping data 
  for (flowID in flowIDs){
    if (overlaps[widedf$DATE_SAMPLE == flowID] > 1){
      # if there is overlapping data between KJ1 and KJ2, correct KJ2 data for KJ1-KJ2 slope and calculate average
      KJ1val <- as.numeric(widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ1'])
      KJ2val <- as.numeric(widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ2'])
      correctedKJ2val <- mu + b*KJ2val # calculate corrected KJ2 value 
      avg <- mean(KJ1val, correctedKJ2val)
      pheno3[flowID, phenotype] <- avg 
      pheno3[flowID, 'Gater'] <- 'KJ1avg'
      next 
    } else if (!is.na(widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ1'])){
      # no overlapping data, but we have KJ1 data, so use KJ1 data
      # this will include all data from 11/13/19-1/24/20, except for central mem, effector mem and naive from Kalika (eventually)
      pheno3[flowID, phenotype] <- widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ1']
      pheno3[flowID, 'Gater'] <- 'KJ1'
      next 
    } else if (!is.na(widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ2'])){
      # no overlapping data, no KJ1 data, but we have KJ2 data, so use KJ2 data 
      # this will be data from 06/20/19-11/12/19 that Kara re-gated with naive, central mem, effector mem 
      pheno3[flowID, phenotype] <- widedf[which(widedf$DATE_SAMPLE == flowID), 'KJ2']
      pheno3[flowID, 'Gater'] <- 'KJ2'
      next
    } else if (!is.na(widedf[which(widedf$DATE_SAMPLE == flowID), 'KK'])){
      # no KJ1 or KJ2 data, but we have KK data, use KK data 
      pheno3[flowID, phenotype] <- widedf[which(widedf$DATE_SAMPLE == flowID), 'KK']
      pheno3[flowID, 'Gater'] <- 'KK'
      next 
    } else {
      # no data 
      pheno3[flowID, phenotype] <- NA
    }
  }
}

write.csv(pheno3, 'source_data/lymph-harmonization-3.csv')






