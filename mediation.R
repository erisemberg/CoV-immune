library(readr)
library(readxl)
library(qtl)
library(ggplot2)
library(cowplot)
source("code-dependencies/cov_qtl_functions.R")
source("code-dependencies/bmediatR-master/R/bmediatR.R")

# ---------------------------------load data---------------------------------- #
cross_data <- read_csv("derived_data/cross_data.csv")

# load cross for QTL mapping 
cross <- read.cross(format = "csv", 
                    file = "derived_data/RqtlCC006xCC044_ctrlAndSARS.csv",
                    na.strings=c("-","NA", "na", "<NA>", "not tested"), 
                    genotypes=c("AA","AB","BB"),
                    alleles=c("A","B")) # A = CC006, B = CC044
cross <- jittermap(cross)

pheno_names <- read_excel("source_data/pheno_names.xlsx")
all_phenos <- pheno_names$flow_col_name
q <- length(all_phenos)
flow_phenos <- all_phenos[-1]

disease_phenos <- c("weight_aac", "HS")

groups <- c("PBS", "SARSCoV", "SARS2CoV", "GxT")
strat_groups <- c("PBS", "SARSCoV", "SARS2CoV")
group_disp_names <- c("PBS" = "Control", "SARSCoV" = "SARS-CoV", "SARS2CoV" = "SARS-CoV-2")
group_dirs <- list("PBS" = "PBS", "SARSCoV" = "SARS", "SARS2CoV" = "SARS2", "GxT" = "GxT")
group_lookup <- c("C" = "PBS", "1" = "SARSCoV", "2" = "SARS2CoV", "A" = "GxT")

geno_map <- list(A = "CC006", B = "CC044")

qtl_map <- read_csv("results/qtl_map.csv")
dis_qtl_map <- read_csv("results/dis_qtl_map.csv")

condensed_qtl <- read_csv("results/qtl_mapping/condensed_qtl.csv")

imp_phenos_df <- read_csv("results/var_selection/important_phenos.csv")
important_phenos <- imp_phenos_df$phenotype[imp_phenos_df$stable]

# ---------------------------load Leist (2024) QTL---------------------------- #
# load QTL from Leist (2024)
disease_qtl <- read_csv("results/qtl_mapping/leist_qtl.csv")
wt_loss_phenos <- c("Weight loss (2 dpi)", "Weight loss (3 dpi)", "Weight loss (4 dpi)", "Weight loss (AAC)")

# condense across infections 
cond_dis_qtl <- disease_qtl %>%
  mutate(chr = sub("Chr ([0-9]+):.*", "\\1", chr_pos),
         pos = as.numeric(sub(".*: ([0-9.]+) \\(.*", "\\1", chr_pos)),
         CI_lwr = as.numeric(sub(".*\\(([0-9.]+)-.*", "\\1", chr_pos)),
         CI_upr = as.numeric(sub(".*-([0-9.]+)\\)", "\\1", chr_pos))) %>%
  group_by(qtl_id, phenotype) %>%
  summarise(infection = paste(infection, collapse = ","),
            chr = chr[which.min(adjP)],
            pos = pos[which.min(adjP)],               
            CI_lwr_conservative = min(CI_lwr, na.rm = TRUE),       
            CI_upr_conservative = max(CI_upr, na.rm = TRUE),       
            CI_lwr_liberal = max(CI_lwr, na.rm = TRUE),
            CI_upr_liberal = min(CI_upr, na.rm = TRUE),
            var_expl = var_expl[which.min(adjP)],
            adjP = min(adjP),
            .groups =) %>%
  ungroup() %>%
  select(qtl_id, disease_pheno = phenotype, infection, chr, pos, adjP) %>%
  mutate(disease_pheno = case_when(disease_pheno %in% wt_loss_phenos ~ "weight_aac",
                                   disease_pheno == "Congestion score" ~ "HS")) %>%
  group_by(qtl_id, disease_pheno) %>%
  summarise(infection = paste(sort(unique(trimws(unlist(strsplit(infection, ","))))), collapse = ","),
            chr = chr[which.min(adjP)],
            pos = pos[which.min(adjP)],
            adjP = min(adjP),
            .groups = "drop")

# -------------------------------GxT mediation-------------------------------- #
# assume titer = 0 for control mice so we can incl. as a covariate in 3-group mediation 
cross_data$Titer[cross_data$infection == "PBS"] <- 0 

med_fname <- "results/mediation/med_gxt_all.csv"
if(file.exists(med_fname)){
  message(paste("Mediation results already exist in", med_fname))
  med_gxt_all <- read_csv("results/mediation/med_gxt_all.csv")
} else {
  # -------------------------mediation on immune QTL-------------------------- #
  message("Running mediation analysis on immune QTL...")
  med_gxt_all <- data.frame(
    chr = character(), marker = character(), 
    immune_pheno = character(), disease_pheno = character(), analysis = character(),
    BFmed_any = numeric(), BFmed_grp1 = numeric(), BFmed_grp2 = numeric(), BFmed_grp3 = numeric(),
    BFcomp_any = numeric(), BFcomp_grp1 = numeric(), BFcomp_grp2 = numeric(), BFcomp_grp3 = numeric(),
    BFpart_any = numeric(), BFpart_grp1 = numeric(), BFpart_grp2 = numeric(), BFpart_grp3 = numeric())
  ix = 0
  for (i in 1:nrow(condensed_qtl)){ # for all immune QTL
    chr <- condensed_qtl$chr[i]
    marker <- find.marker(cross, chr, as.numeric(condensed_qtl$pos[i]))
    immune_pheno <- condensed_qtl$phenotype[i]
    
    for (disease_pheno in disease_phenos){
      ix = ix + 1
      df <- data.frame(X = pull.geno(cross, chr)[,marker],
                       M = cross_data[[immune_pheno]],
                       y = cross_data[[disease_pheno]],
                       infection = cross_data$infection,
                       sex = cross_data$sex,
                       Titer = cross_data$Titer)
      df <- na.omit(df) # columns don't change w/in loop, so we can do this outside of loop 
      
      if (immune_pheno == "Titer" | disease_pheno == "HS"){ # 2-way GxT mediation 
        analysis <- "2group"
        df <- df[df$infection %in% c("SARSCoV", "SARS2CoV"),]
        
        if (immune_pheno == "Titer"){
          Z <- model.matrix(~ sex, data = df)[,-1, drop = FALSE] 
        } else { # incl. titer as a covariate if titer isn't the mediator 
          Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE] 
        }
        M <- as.matrix(df$M) 
        colnames(M) <- immune_pheno
        X <- model.matrix(~ X, data = df)[,-1, drop = FALSE] 
        group_med <- as.numeric(df$infection == "SARS2CoV") 
        
        # define weights for y and M to handle variance heterogeneity
        df$y <- scale(df$y)
        fit_y <- lm(y ~ sex + infection, data = df)
        res <- resid(fit_y)
        s2_y1 <- var(res[group_med == 0])
        s2_y2 <- var(res[group_med == 1])
        wts_y <- ifelse(group_med == 0, 1/s2_y1, 1/s2_y2)
        
        fit_m <- lm(M ~ sex + infection, data = df)
        res <- resid(fit_m)
        s2_m1 <- var(res[group_med == 0])
        s2_m2 <- var(res[group_med == 1])
        wts_M <- ifelse(group_med == 0, 1/s2_m1, 1/s2_m2)
        
        # run bmediatR 
        med_gxt <- bmediatR_GxT(y = df$y, M = M, X = X, Z = Z, 
                                group = group_med, w_y = wts_y, w_M = wts_M,
                                align_data = FALSE,
                                options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
      
      } else { # 3-way GxT mediation 
        analysis <- "3group"
        Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE] # covariate matrix for both M and y 
        M <- as.matrix(df$M) # mediator matrix 
        colnames(M) <- immune_pheno
        X <- model.matrix(~ X, data = df)[,-1, drop = FALSE] # exposure matrix 
        group_med <- ifelse(df$infection == "SARS2CoV", 2, 
                            ifelse(df$infection == "SARSCoV", 1, 0))
        
        # define weights for y and M to handle variance heterogeneity
        df$y <- scale(df$y, center = TRUE, scale = TRUE)
        fit_y <- lm(y ~ sex + infection, data = df)
        res <- resid(fit_y)
        s2_y0 <- var(res[group_med == 0])
        s2_y1 <- var(res[group_med == 1])
        s2_y2 <- var(res[group_med == 2])
        wts_y <- ifelse(group_med == 2, 1/s2_y2, 
                        ifelse(group_med == 1, 1/s2_y1, 1/s2_y0))
        
        fit_m <- lm(M ~ sex + infection, data = df)
        res <- resid(fit_m)
        s2_m0 <- var(res[group_med == 0])
        s2_m1 <- var(res[group_med == 1])
        s2_m2 <- var(res[group_med == 2])
        wts_M <- ifelse(group_med == 2, 1/s2_m2, 
                        ifelse(group_med == 1, 1/s2_m1, 1/s2_m0))
        
        # run bmediatR 
        med_gxt <- bmediatR_GxT3(y = df$y, M = M, X = X, Z = Z, 
                                 group = group_med, w_y = wts_y, w_M = wts_M,
                                 align_data = FALSE,
                                 options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
      }
      
      BFmed <- signif(log10(exp(med_gxt$ln_post_odds) / exp(med_gxt$ln_prior_odds)), 2)
      med_gxt_all[ix,] <- c(chr, marker, immune_pheno, disease_pheno, analysis,
                            BFmed[1,"mediation_any"], BFmed[1,"mediation_grp1"], BFmed[1,"mediation_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"mediation_grp3"], NA),
                            BFmed[1,"complete_any"], BFmed[1,"complete_grp1"], BFmed[1,"complete_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"complete_grp3"], NA),
                            BFmed[1,"partial_any"], BFmed[1,"partial_grp1"], BFmed[1,"partial_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"partial_grp3"], NA))
      
      # save posterior bar plot to png 
      save_fname <- paste0("figures/med_wTiterCov/", immune_pheno, "_", disease_pheno, "_chr", chr, ".png")
      png(save_fname, width = 450, height = 300)
      if (analysis == "2group"){
        print(plot_posterior_bar_gxt(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                     grp_labels = c(grp1 = group_disp_names[[2]], 
                                                    grp2 = group_disp_names[[3]])))
      } else {
        print(plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                      grp_labels = c(grp1 = group_disp_names[[1]], 
                                                     grp2 = group_disp_names[[2]],
                                                     grp3 = group_disp_names[[3]])))
      }
      dev.off()
    }
  }
  
  # ------------------------mediation on disease QTL-------------------------- #
  for (i in 1:nrow(cond_dis_qtl)){ # for all disease QTL from Leist (2024)
    chr <- cond_dis_qtl$chr[i]
    marker <- find.marker(cross, chr, as.numeric(cond_dis_qtl$pos[i]))
    qtl_id <- cond_dis_qtl$qtl_id[i]
    disease_pheno <- cond_dis_qtl$disease_pheno[i] 
    
    for (immune_pheno in all_phenos){ # check for mediation by all immune traits
      df <- data.frame(X = pull.geno(cross, chr)[,marker],
                       M = cross_data[[immune_pheno]],
                       y = cross_data[[disease_pheno]],
                       infection = cross_data$infection,
                       sex = cross_data$sex,
                       Titer = cross_data$Titer,
                       d0 = cross_data$d0)
      df <- na.omit(df) 
      
      # only test for mediation if X-M association is significant? ### try without titer
      XMp <- summary(lm(M ~ infection + Titer + sex + X, data = df))$coefficients["X","Pr(>|t|)"]
      if (XMp > 0.05){ next }
      ix = ix + 1
      
      if (immune_pheno == "Titer" | disease_pheno == "HS"){ # 2-group analysis 
        analysis <- "2group"
        df <- df[df$infection %in% c("SARSCoV", "SARS2CoV"),]
        if (immune_pheno == "Titer"){
          Z <- model.matrix(~ sex, data = df)[,-1, drop = FALSE]
        } else { # incl. titer as a covariate if titer is not the mediator
          Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE]
        }
        M <- as.matrix(df$M)  
        colnames(M) <- immune_pheno
        X <- model.matrix(~ X, data = df)[,-1, drop = FALSE] 
        group_med <- as.numeric(df$infection == "SARS2CoV") 
        
        # define weights for y and M to handle variance heterogeneity
        df$y <- scale(df$y)
        fit_y <- lm(y ~ sex + infection, data = df)
        res <- resid(fit_y)
        s2_y1 <- var(res[group_med == 0])
        s2_y2 <- var(res[group_med == 1])
        wts_y <- ifelse(group_med == 0, 1/s2_y1, 1/s2_y2)
        
        fit_m <- lm(M ~ sex + infection, data = df)
        res <- resid(fit_m)
        s2_m1 <- var(res[group_med == 0])
        s2_m2 <- var(res[group_med == 1])
        wts_M <- ifelse(group_med == 0, 1/s2_m1, 1/s2_m2)
        
        med_gxt <- bmediatR_GxT(y = df$y, M = M, X = X, Z = Z, 
                                group = group_med, w_y = wts_y, w_M = wts_M,
                                align_data = FALSE,
                                options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
      } else { # 3 way GxT mediation
        analysis <- "3group"
        Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE] 
        M <- as.matrix(df$M) 
        colnames(M) <- immune_pheno
        X <- model.matrix(~ X, data = df)[,-1, drop = FALSE] 
        group_med <- ifelse(df$infection == "SARS2CoV", 2, 
                            ifelse(df$infection == "SARSCoV", 1, 0))
        
        # define weights for y and M to handle variance heterogeneity
        df$y <- scale(df$y, center = TRUE, scale = TRUE)
        fit_y <- lm(y ~ sex + infection, data = df)
        res <- resid(fit_y)
        s2_y0 <- var(res[group_med == 0])
        s2_y1 <- var(res[group_med == 1])
        s2_y2 <- var(res[group_med == 2])
        wts_y <- ifelse(group_med == 2, 1/s2_y2, 
                        ifelse(group_med == 1, 1/s2_y1, 1/s2_y0))
        
        fit_m <- lm(M ~ sex + infection, data = df)
        res <- resid(fit_m)
        s2_m0 <- var(res[group_med == 0])
        s2_m1 <- var(res[group_med == 1])
        s2_m2 <- var(res[group_med == 2])
        wts_M <- ifelse(group_med == 2, 1/s2_m2, 
                        ifelse(group_med == 1, 1/s2_m1, 1/s2_m0))
        
        med_gxt <- bmediatR_GxT3(y = df$y, M = M, X = X, Z = Z, 
                                 group = group_med, w_y = wts_y, w_M = wts_M,
                                 align_data = FALSE,
                                 options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
      }
      
      BFmed <- signif(log10(exp(med_gxt$ln_post_odds) / exp(med_gxt$ln_prior_odds)), 2)
      med_gxt_all[ix,] <- c(chr, marker, immune_pheno, disease_pheno, analysis,
                            BFmed[1,"mediation_any"], BFmed[1,"mediation_grp1"], BFmed[1,"mediation_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"mediation_grp3"], NA),
                            BFmed[1,"complete_any"], BFmed[1,"complete_grp1"], BFmed[1,"complete_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"complete_grp3"], NA),
                            BFmed[1,"partial_any"], BFmed[1,"partial_grp1"], BFmed[1,"partial_grp2"],
                            ifelse(analysis == "3group", BFmed[1,"partial_grp3"], NA))
      
      # save posterior bar plot to png 
      save_fname <- paste0("figures/med_wTiterCov/", cond_dis_qtl$qtl_id[i], "_", immune_pheno, "_", disease_pheno, ".png")
      png(save_fname, width = 450, height = 300)
      if (analysis == "2group"){
        print(plot_posterior_bar_gxt(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                     grp_labels = c(grp1 = group_disp_names[[2]], 
                                                    grp2 = group_disp_names[[3]])))
      } else {
        print(plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                      grp_labels = c(grp1 = group_disp_names[[1]], 
                                                     grp2 = group_disp_names[[2]],
                                                     grp3 = group_disp_names[[3]])))
      }
      dev.off()
    }
  }
  write_csv(med_gxt_all, med_fname)
}


# ----------------------------------Table 2----------------------------------- #
est_XM_effect <- function(cross_data, cross, marker, immune_pheno, 
                          group_var = "infection", group_value, covars = c("sex")) {
  d <- cross_data[cross_data[[group_var]] == group_value, , drop = FALSE]
  
  if (which_chr(cross, marker) != "X"){ # autosome
    f <- as.formula(paste(immune_pheno, "~", paste(c(covars, marker), collapse = " + ")))
    
    fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(NA_character_)
    beta <- summary(fit)$coefficients[marker, "Estimate"]
    sig <- summary(fit)$coefficients[marker, "Pr(>|t|)"]
    
    return(paste0(signif(beta, 2), ifelse(sig < 0.05, "*", "")))
  } else { # X-chromosome marker 
    covars <- setdiff(covars, "sex")
    f <- as.formula(paste(immune_pheno, "~", paste(covars, marker), collapse = " + "))
    
    # females
    d_f <- d[d$sex == 0,]
    fit <- tryCatch(lm(f, data = d_f), error = function(e) NULL)
    if (is.null(fit)) return(NA_character_)
    betaF <- summary(fit)$coefficients[marker, "Estimate"]
    sigF <- summary(fit)$coefficients[marker, "Pr(>|t|)"]
    
    # males - recode 0/2 -> 0/1 so on same additive scale as females 
    d_m <- d[d$sex == 1,]
    d_m[[marker]] <- d_m[[marker]] / 2 # 0 stays 0, 2 becomes 1 
    fit <- tryCatch(lm(f, data = d_m), error = function(e) NULL)
    if (is.null(fit)) return(NA_character_)
    betaM <- summary(fit)$coefficients[marker, "Estimate"]
    sigM <- summary(fit)$coefficients[marker, "Pr(>|t|)"]
    
    return(paste0("female: ", signif(betaF, 2), ifelse(sigF < 0.05, "*", ""), "; male: ", signif(betaM, 2), ifelse(sigM < 0.05, "*", "")))
  }
}

est_MY_effect <- function(cross_data, immune_pheno, disease_pheno, 
                          group_var = "infection", group_value, covars = c("sex")) {
  d <- cross_data[cross_data[[group_var]] == group_value, , drop = FALSE]
  f <- as.formula(paste(disease_pheno, "~", paste(c(covars, immune_pheno), collapse = " + ")))
  fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA)
  beta <- summary(fit)$coefficients[immune_pheno, "Estimate"]
  sig <- summary(fit)$coefficients[immune_pheno, "Pr(>|t|)"]
  return(paste0(signif(beta, 2), ifelse(sig < 0.05, "*", "")))
}

# for mapping disease QTL to IDs
pos_tbl <- find.markerpos(cross, med_gxt_all$marker)
pos_tbl$marker <- rownames(pos_tbl)

dis_qtl_map <- dis_qtl_map %>% 
  mutate(pos = signif(pos, 3))

med_gxt_all$source <- c(rep("imm_QTL", 136), rep("dis_QTL", 103))
  
table2 <- med_gxt_all %>% 
  as_tibble() %>%
  mutate(across(starts_with("BF"), ~ suppressWarnings(as.numeric(.x)))) %>%
  left_join(pos_tbl, by = c("chr", "marker")) %>%
  mutate(pos = signif(pos, 3)) %>%
  mutate(grp_control = "PBS",
         grp_sars    = "SARSCoV",
         grp_sars2   = "SARS2CoV") %>%
  mutate(BFmed_control  = if_else(analysis == "3group", BFmed_grp1, NA_real_),
         BFcomp_control = if_else(analysis == "3group", BFcomp_grp1, NA_real_),
         BFpart_control = if_else(analysis == "3group", BFpart_grp1, NA_real_),
         BFmed_sars     = if_else(analysis == "3group", BFmed_grp2, BFmed_grp1),
         BFcomp_sars    = if_else(analysis == "3group", BFcomp_grp2, BFcomp_grp1),
         BFpart_sars    = if_else(analysis == "3group", BFpart_grp2, BFpart_grp1),
         BFmed_sars2    = if_else(analysis == "3group", BFmed_grp3, BFmed_grp2),
         BFcomp_sars2   = if_else(analysis == "3group", BFcomp_grp3, BFcomp_grp2),
         BFpart_sars2   = if_else(analysis == "3group", BFpart_grp3, BFpart_grp2)) %>%
  rowwise() %>%
  mutate(
    XM_control = if (analysis == "3group") est_XM_effect(cross_data, cross, marker, immune_pheno, group_value = grp_control) else NA,
    XM_sars    = est_XM_effect(cross_data, cross, marker, immune_pheno, group_value = grp_sars),
    XM_sars2   = est_XM_effect(cross_data, cross, marker, immune_pheno, group_value = grp_sars2),
    MY_control = if (analysis == "3group") est_MY_effect(cross_data, immune_pheno, disease_pheno, group_value = grp_control) else NA,
    MY_sars    = est_MY_effect(cross_data, immune_pheno, disease_pheno, group_value = grp_sars),
    MY_sars2   = est_MY_effect(cross_data, immune_pheno, disease_pheno, group_value = grp_sars2)
  ) %>%
  ungroup() %>%
  left_join(pheno_names, by = c("immune_pheno" = "flow_col_name")) %>%
  left_join(qtl_map %>% transmute(chr, immune_pheno = phenotype, qtl_id_imm = qtl_id), by = c("chr", "immune_pheno")) %>%
  left_join(dis_qtl_map %>% transmute(chr = as.character(chr), pos, qtl_id_dis = qtl_id), by = c("chr", "pos")) %>%
  mutate(qtl_id = case_when(source == "imm_QTL" ~ qtl_id_imm,
                            source == "dis_QTL" ~ qtl_id_dis,
                            TRUE ~ NA)) %>%
  select(chr, pos, qtl_id, source, marker, immune_pheno, flow_display_name, disease_pheno, analysis,
         grp_control, XM_control, MY_control, BFmed_control, BFcomp_control, BFpart_control,
         grp_sars, XM_sars, MY_sars, BFmed_sars, BFcomp_sars, BFpart_sars,
         grp_sars2, XM_sars2, MY_sars2, BFmed_sars2, BFcomp_sars2, BFpart_sars2) %>%
  filter(if_any(starts_with("BF"), ~ !is.na(.x) & .x > bf_thresh)) %>%
  arrange(as.numeric(chr))

# table for publication 
table2_long <- bind_rows(
  table2 %>%
    transmute(
      chr, qtl_id, immune_pheno, flow_display_name, disease_pheno,
      group = grp_control,
      XM_effect = XM_control,
      MY_effect = MY_control,
      BFmed = BFmed_control,
      BFcomp = BFcomp_control,
      BFpart = BFpart_control),
  table2 %>%
    transmute(
      chr, qtl_id, immune_pheno, flow_display_name, disease_pheno,
      group = grp_sars,
      XM_effect = XM_sars,
      MY_effect = MY_sars,
      BFmed = BFmed_sars,
      BFcomp = BFcomp_sars,
      BFpart = BFpart_sars),
  table2 %>%
    transmute(
      chr, qtl_id, immune_pheno, flow_display_name, disease_pheno,
      group = grp_sars2,
      XM_effect = XM_sars2,
      MY_effect = MY_sars2,
      BFmed = BFmed_sars2,
      BFcomp = BFcomp_sars2,
      BFpart = BFpart_sars2)) %>%
  filter(!is.na(BFmed) | !is.na(BFcomp) | !is.na(BFpart)) %>% # drop NA BFs, which are control rows for 2group analysis 
  filter(BFmed > bf_thresh | BFcomp > bf_thresh | BFpart > bf_thresh) %>% # Keep only groups with evidence in at least one context
  mutate(group = case_when(group == "PBS" ~ "Control",
                           group == "SARSCoV" ~ "SARS-CoV",
                           group == "SARS2CoV" ~ "SARS-CoV-2",
                           TRUE ~ as.character(group)),
         group = factor(group, levels = c("Control", "SARS-CoV", "SARS-CoV-2")),
         disease_pheno = case_when(disease_pheno == "weight_aac" ~ "Weight loss",
                                   disease_pheno == "HS" ~ "Lung CS"),
         disease_pheno = factor(disease_pheno, levels = c("Weight loss", "Lung CS"))) %>%
  arrange(as.numeric(chr), qtl_id, immune_pheno, disease_pheno) %>%
  select(-immune_pheno, -chr)

writexl::write_xlsx(table2_long, "results/mediation/SuppTable3-mediation.xlsx")


# --------------------------summarize shared vs not--------------------------- #
table(table2$immune_pheno)

flag_med <- function(bfmed, bfcomp, bfpart, thr = 0.5) {
  all_missing <- is.na(bfmed) & is.na(bfcomp) & is.na(bfpart)
  out <- pmax(bfmed, bfcomp, bfpart, na.rm = TRUE) > thr
  ifelse(all_missing, NA, out)
}

table2_flagged <- table2 %>%
  mutate(med_sars = flag_med(BFmed_sars,    BFcomp_sars,    BFpart_sars),
         med_sars2 = flag_med(BFmed_sars2,   BFcomp_sars2,   BFpart_sars2))

trait_level <- table2_flagged %>%
  group_by(qtl_id, immune_pheno, flow_display_name) %>%
  summarise(sars_any = any(med_sars    %in% TRUE, na.rm = TRUE),
            sars2_any = any(med_sars2   %in% TRUE, na.rm = TRUE),
            .groups = "drop") %>%              
  mutate(shared = sars_any & sars2_any)

# how many are shared vs virus-specific?
trait_summary <- trait_level %>%
  filter(sars_any | sars2_any) %>%
  summarise(n_traits = n(),
            n_shared = sum(shared),
            prop_shared = n_shared / n_traits,
            n_sars_only  = sum(sars_any & !sars2_any),
            n_sars2_only = sum(!sars_any & sars2_any),
            n_neither = sum(!sars_any & !sars2_any))
trait_summary

# how many unique QTL and immune phenotypes?
length(unique(trait_level$immune_pheno))
length(unique(trait_level$qtl_id))

# which are shared?
trait_level %>%
  filter(shared) %>%
  arrange(qtl_id, immune_pheno)

# QTL level sharing 
qtl_share <- table2 %>%
  mutate(med_sars_row = pmax(BFmed_sars,    BFcomp_sars,    BFpart_sars,    na.rm = TRUE) > 0.5,
         med_sars2_row = pmax(BFmed_sars2,   BFcomp_sars2,   BFpart_sars2,   na.rm = TRUE) > 0.5) %>%
  group_by(qtl_id) %>%
  summarise(sars_med_any = any(med_sars_row,    na.rm = TRUE),
            sars2_med_any = any(med_sars2_row,   na.rm = TRUE),
            .groups = "drop") %>%
  mutate(shared_infections = sars_med_any & sars2_med_any)

qtl_share_summary <- qtl_share %>%
  summarise(n_qtl = n(),
            n_shared = sum(shared_infections, na.rm = TRUE),
            prop_shared = n_shared / n_qtl,
            n_sars_only = sum(sars_med_any & !sars2_med_any, na.rm = TRUE),
            n_sars2_only = sum(!sars_med_any & sars2_med_any, na.rm = TRUE),
            n_neither = sum(!sars_med_any & !sars2_med_any, na.rm = TRUE))
qtl_share_summary

qtl_share %>%
  filter(shared_infections) %>%
  pull(qtl_id)

trait_share <- table2 %>%
  mutate(
    med_sars_row    = pmax(BFmed_sars,    BFcomp_sars,    BFpart_sars,    na.rm = TRUE) > 0.5,
    med_sars2_row   = pmax(BFmed_sars2,   BFcomp_sars2,   BFpart_sars2,   na.rm = TRUE) > 0.5
  ) %>%
  group_by(immune_pheno) %>%
  summarise(
    sars_med_any    = any(med_sars_row,    na.rm = TRUE),
    sars2_med_any   = any(med_sars2_row,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(shared_infections = sars_med_any & sars2_med_any)
trait_share_summary <- trait_share %>%
  summarise(
    n_qtl = n(),
    n_shared = sum(shared_infections, na.rm = TRUE),
    prop_shared = n_shared / n_qtl,
    n_sars_only = sum(sars_med_any & !sars2_med_any, na.rm = TRUE),
    n_sars2_only = sum(!sars_med_any & sars2_med_any, na.rm = TRUE),
    n_neither = sum(!sars_med_any & !sars2_med_any, na.rm = TRUE)
  )
trait_share_summary



# ----------------------------------Figures----------------------------------- #
inf_pal <- list("PBS" = "#66C2A5", "SARSCoV" = "#FC8D62", "SARS2CoV" = "#8DA0CB")
gxt_pal <- c("black", rev(brewer.pal(n = 3, "Set1")[1:2]))
grouplabels = list("PBS" = "Control", "SARSCoV" = "SARS-CoV", "SARS2CoV" = "SARS-CoV-2")

# ---------------------------------------------------------------------------- #
# -----------------------------------Fig 4------------------------------------ #
# ---------------------------------------------------------------------------- #



# ----------------------------Naive helper T cells---------------------------- #
immune_pheno <- "L_HelpTCells_pctNaive"
marker <- "S6T025109239"
disease_pheno <- "weight_aac"
df <- data.frame(mouseID = cross_data$mouse_ID,
                 geno = pull.geno(cross, which_chr(cross, marker))[,marker],
                 M = cross$pheno[[immune_pheno]], # untransformed data for plotting 
                 Mtr = cross_data[[immune_pheno]], # transformed data for regression 
                 y = cross_data[[disease_pheno]], 
                 infection = cross_data$infection,
                 sex = cross_data$sex)

helpN_scan <- make_qtl_scan(flow_pheno = immune_pheno, palette = inf_pal, legend.x = 1700, legend.y = 6, ylim = c(0,10.1))
plot_grid(helpN_scan)

# --------------------------------X->M effect--------------------------------- #
# significance 
geno_terms <- get_add_dom(cross, marker)
df$add <- geno_terms$add
df$dom <- geno_terms$dom
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "PBS",]), lm(Mtr ~ sex, data = df[df$infection == "PBS",])) 
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARSCoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARSCoV",])) 
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARS2CoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARS2CoV",])) 
sig_df <- data.frame(group = strat_groups,
                     label = c("p = 3.4e-9", "p = 2.9e-6", "p = 2.9e-3"),
                     vjust = c(1.2, 2.6, 4.0))
# plot 
helpN_pxg <- pxg_group(cross, pheno = df$Mtr, marker = marker, 
                  geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("Irq4")~"genotype"), 
                  ylab = "Helper T cells (% naive) (transformed)", type = "boxplot",
                  groupcolumn = "infection", groupcolors = inf_pal, grouplabels = grouplabels) + 
  geom_text(data = sig_df, 
            aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) + 
  scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)])
helpN_pxg
# box_legend <- get_legend(AMp2)
# AMp2 <- AMp2 +
#   theme(legend.position = "none")

# --------------------------------M->y effect--------------------------------- #
# significance 
summary(lm(y ~ sex + Mtr, data = df[df$infection == "PBS",])) # p = 0.95
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARSCoV",])) # p = 0.0017
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARS2CoV",])) # p = 0.02
sig_df <- data.frame(group = strat_groups,
                     label = c("ns", "p = 0.0017", "p = 0.02"),
                     vjust = c(1.2, 2.6, 4.0))
helpN_My <- ggplot(data = df, aes(x = Mtr, y = y, color = infection)) + 
  geom_point() + 
  scale_color_manual(values = inf_pal, breaks = names(inf_pal), name = "", labels = grouplabels[names(inf_pal)]) + 
  geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
  labs(x = "Helper T cells (% naive) (transformed)", y = "Weight loss (AAC)") + 
  geom_text(data = sig_df, 
            mapping = aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  theme_minimal()
helpN_My
# scat_legend <- get_legend(AMp3)
# AMp3 <- AMp3 + 
#   theme(legend.position = "none")
# AMp3

# ---------------------------------mediation---------------------------------- # 
df <- na.omit(df)
Z <- model.matrix(~ sex, data = df)[,-1, drop = FALSE] 
M <- as.matrix(df$Mtr) # mediator matrix 
colnames(M) <- immune_pheno
X <- model.matrix(~ geno, data = df)[,-1, drop = FALSE] # exposure matrix 
group_med <- ifelse(df$infection == "SARS2CoV", 2, 
                    ifelse(df$infection == "SARSCoV", 1, 0))
# define weights for y and M to handle variance heterogeneity
df$y <- scale(df$y, center = TRUE, scale = TRUE)
fit_y <- lm(y ~ sex + infection, data = df)
res <- resid(fit_y)
s2_y0 <- var(res[group_med == 0])
s2_y1 <- var(res[group_med == 1])
s2_y2 <- var(res[group_med == 2])
wts_y <- ifelse(group_med == 2, 1/s2_y2, 
                ifelse(group_med == 1, 1/s2_y1, 1/s2_y0))
fit_m <- lm(Mtr ~ sex + infection, data = df)
res <- resid(fit_m)
s2_m0 <- var(res[group_med == 0])
s2_m1 <- var(res[group_med == 1])
s2_m2 <- var(res[group_med == 2])
wts_M <- ifelse(group_med == 2, 1/s2_m2, 
                ifelse(group_med == 1, 1/s2_m1, 1/s2_m0))

med_gxt <- bmediatR_GxT3(y = df$y, M = M, X = X, Z = Z, 
                         group = group_med, w_y = wts_y, w_M = wts_M,
                         align_data = FALSE,
                         options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
helpN_med <- plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                      grp_labels = c(grp1 = "Control", grp2 = "SARS-CoV", grp3 = "SARS-CoV-2"))
helpN_med
# ccr2Lp_med <- ccr2Lp_med +
#   theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
#         legend.key.height = unit(0.5, 'cm'), #change legend key height
#         legend.key.width = unit(0.5, 'cm'), #change legend key width
#         legend.title = element_text(size=10), #change legend title font size
#         legend.text = element_text(size=9))


# ---------------------------------------------------------------------------- #
# ---------------------------CD86+ Ly6C+ mono/macs---------------------------- #
# ---------------------------------------------------------------------------- #
immune_pheno <- "M_Ly6Cp_MonoMac_pctCD86"
marker <- "S3N094839317" 
disease_pheno <- "weight_aac"
df <- data.frame(mouseID = cross_data$mouse_ID,
                 geno = pull.geno(cross, which_chr(cross, marker))[,marker],
                 M = cross$pheno[[immune_pheno]], # untransformed data for plotting 
                 Mtr = cross_data[[immune_pheno]], # transformed for regression 
                 y = cross_data[[disease_pheno]],
                 infection = cross_data$infection,
                 sex = cross_data$sex,
                 Titer = cross_data$Titer)

# --------------------------------genome scan--------------------------------- #
# no genome scan / no significant QTL - associated with HrS43

# --------------------------------X->M effect--------------------------------- #
geno_terms <- get_add_dom(cross, marker)
df$add <- geno_terms$add
df$dom <- geno_terms$dom
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "PBS",]), lm(Mtr ~ sex, data = df[df$infection == "PBS",])) # p = 0.49
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARSCoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARSCoV",])) # p = 0.03
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARS2CoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARS2CoV",])) # p = 0.12

sig_df <- data.frame(group = strat_groups,
                     label = c("p = 0.5", "p = 0.03", "p = 0.1"),
                     vjust = c(1.2, 2.6, 4.0)) # ,hjust = c(1.1, 1, 1.1)
cd86Lp_pxg <- pxg_group(cross, pheno = df$Mtr, marker = marker, 
                        geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("HrS43")~"genotype"), 
                        ylab = expression("Ly6C"^"+"~"monocytes (% CD86"^"+"*") (transformed)"), groupcolumn = "infection", 
                        groupcolors = inf_pal, grouplabels = grouplabels, type = "boxplot") + 
  geom_text(data = sig_df, 
            aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) + 
  scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)])
cd86Lp_pxg

# for defense slides
# pxg_group(subset(cross, ind = (cross$pheno$infection == "SARS2CoV")), pheno = df$M[df$infection == "SARS2CoV"], marker = marker, 
#           geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("HrS43")~"genotype"), 
#           ylab = expression("Inflammatory monocytes (% CD86"^"+"~")"), groupcolumn = "infection", 
#           groupcolors = inf_pal, grouplabels = grouplabels, type = "boxplot") + 
#   scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)]) + 
#   labs(y = NULL, title = expression("Inflammatory monocytes (% CD86"^"+"~")")) + 
#   theme(legend.position = "none",
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 14))
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter4/HrS43_cd86Lp.png", width = 5, height = 4)

# --------------------------------X->y effect--------------------------------- #
# not necessary - this QTL was already reported in Leist (2024)
# anova(lm(y ~ sex + add + dom, data = df[df$infection == "PBS",]), lm(y ~ sex, data = df[df$infection == "PBS",])) # p = 0.076
# anova(lm(y ~ sex + add + dom, data = df[df$infection == "SARSCoV",]), lm(y ~ sex, data = df[df$infection == "SARSCoV",])) # 2.6e-06
# anova(lm(y ~ sex + add + dom, data = df[df$infection == "SARS2CoV",]), lm(y ~ sex, data = df[df$infection == "SARS2CoV",])) # p = 8.0e-05
# sig_df <- data.frame(group = strat_groups,
#                      label = c("p = 0.08", "p = 2.6e-06", "p = 8.0e-05"),
#                      y = c(Inf, Inf, Inf),
#                      vjust = c(1.2, 2.6, 4.0))
# pxg_group(cross, pheno = df$y, marker = marker,
#                         geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("HrS43")~"genotype"),
#                         ylab = "Weight loss (AAC)", type = "boxplot",
#                         groupcolumn = "infection", groupcolors = inf_pal, grouplabels = grouplabels) +
#   geom_text(data = sig_df, 
#             aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
#             hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) + 
#   scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)]) +
#   theme_minimal() + 
#   theme(legend.position = "none")

# --------------------------------M->y effect--------------------------------- #
# significance
summary(lm(y ~ sex + Mtr, data = df[df$infection == "PBS",])) # p = 0.6
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARSCoV",])) # p = 0.08
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARS2CoV",])) # p = 0.02
sig_df <- data.frame(group = strat_groups,
                     label = c("p = 0.6", "p = 0.08", "p = 0.02"),
                     y = c(Inf, Inf, Inf),
                     vjust = c(1.2, 2.6, 4.0))
cd86Lp_wt <- ggplot(data = df, aes(x = Mtr, y = y, color = infection)) +
  geom_point() +
  scale_color_manual(values = inf_pal) +
  geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
  labs(x = expression("Ly6C"^"+"~"monocytes (% CD86"^"+"*") (transformed)"), y = "Weight loss (AAC)") +
  geom_text(data = sig_df, 
            mapping = aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust),
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  theme_minimal() 
cd86Lp_wt

# for defense slides 
# ggplot(data = df[df$infection == "SARS2CoV",], aes(x = M, y = y, color = infection)) +
#   geom_point() +
#   scale_color_manual(values = inf_pal, name = NULL) +
#   geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
#   labs(x = expression("Inflammatory monocytes (% CD86"^"+"~")"), y = "Weight loss (AAC)") +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 14),
#         legend.position = "none")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter4/cd86Lp_wt.png", width = 4, height = 4)

# ---------------------------------mediation---------------------------------- # 
df <- na.omit(df)
Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE] 
M <- as.matrix(df$Mtr) # mediator matrix 
colnames(M) <- immune_pheno
X <- model.matrix(~ geno, data = df)[,-1, drop = FALSE] # exposure matrix 
group_med <- ifelse(df$infection == "SARS2CoV", 2, 
                    ifelse(df$infection == "SARSCoV", 1, 0))
# define weights for y and M to handle variance heterogeneity
df$y <- scale(df$y, center = TRUE, scale = TRUE)
fit_y <- lm(y ~ sex + infection, data = df)
res <- resid(fit_y)
s2_y0 <- var(res[group_med == 0])
s2_y1 <- var(res[group_med == 1])
s2_y2 <- var(res[group_med == 2])
wts_y <- ifelse(group_med == 2, 1/s2_y2, 
                ifelse(group_med == 1, 1/s2_y1, 1/s2_y0))
fit_m <- lm(Mtr ~ sex + infection, data = df)
res <- resid(fit_m)
s2_m0 <- var(res[group_med == 0])
s2_m1 <- var(res[group_med == 1])
s2_m2 <- var(res[group_med == 2])
wts_M <- ifelse(group_med == 2, 1/s2_m2, 
                ifelse(group_med == 1, 1/s2_m1, 1/s2_m0))

med_gxt <- bmediatR_GxT3(y = df$y, M = M, X = X, Z = Z, 
                         group = group_med, w_y = wts_y, w_M = wts_M,
                         align_data = FALSE,
                         options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))
cd86Lp_med <- plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                      grp_labels = c(grp1 = "Control", grp2 = "SARS-CoV", grp3 = "SARS-CoV-2"))
cd86Lp_med

cd86Lp_med <- cd86Lp_med +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=9))
cd86Lp_med

# for defense 
# plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
#                         grp_labels = c(grp1 = "Control", grp2 = "SARS-CoV", grp3 = "SARS-CoV-2")) + 
#   guides(alpha = guide_legend(override.aes = list(fill = "seagreen4"))) +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14))
# ggsave("~/Documents/ValdarFerris/defense/figures/mediation/GxT_cd86_Ly6Cp.png")


# ---------------------------------------------------------------------------- #
# -------------------------------CM DN T cells-------------------------------- #
# ---------------------------------------------------------------------------- #
immune_pheno <- "L_DNTCells_pctCentMem"
marker <- "S3N094839317"
disease_pheno <- "weight_aac"
df <- data.frame(mouseID = cross_data$mouse_ID,
                 geno = pull.geno(cross, which_chr(cross, marker))[,marker],
                 M = cross$pheno[[immune_pheno]], # untransformed data for plotting 
                 Mtr = cross_data[[immune_pheno]], # transformed for regression 
                 y = cross_data[[disease_pheno]],
                 infection = cross_data$infection,
                 sex = cross_data$sex,
                 Titer = cross_data$Titer)

# --------------------------------X->M effect--------------------------------- #
# significance 
geno_terms <- get_add_dom(cross, marker)
df$add <- geno_terms$add
df$dom <- geno_terms$dom
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "PBS",]), lm(Mtr ~ sex, data = df[df$infection == "PBS",])) # p = 0.62
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARSCoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARSCoV",])) # p = 0.02
anova(lm(Mtr ~ sex + add + dom, data = df[df$infection == "SARS2CoV",]), lm(Mtr ~ sex, data = df[df$infection == "SARS2CoV",])) # p = 0.58
sig_df <- data.frame(group = strat_groups,
                     label = c("p = 0.6", "p = 0.02", "p = 0.6"),
                     vjust = c(1.2, 2.6, 4.0))
cmDNT_pxg <- pxg_group(cross, pheno = df$Mtr, marker = marker, 
                       geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("HrS43")~"genotype"), 
                       ylab = "DN T cells (% CM) (transformed)", groupcolumn = "infection", 
                       groupcolors = inf_pal, grouplabels = grouplabels, type = "boxplot") + 
  geom_text(data = sig_df, 
            aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) + 
  scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)])
cmDNT_pxg

box_legend <- get_legend(cmDNT_pxg)
cmDNT_pxg <- cmDNT_pxg + 
  theme(legend.position = "none")

# for defense slides
# pxg_group(subset(cross, ind = (cross$pheno$infection == "SARSCoV")), pheno = df$M[df$infection == "SARSCoV"], marker = marker, 
#           geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("HrS43")~"genotype"), 
#           ylab = "DN T cells (% CM)", groupcolumn = "infection", 
#           groupcolors = inf_pal, grouplabels = grouplabels, type = "boxplot") + 
#   scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)]) + 
#   labs(y = NULL, title = "DN T cells (% CM)") + 
#   theme(legend.position = "none",
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 14))
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter4/HrS43_cmDNT.png", width = 5, height = 4)

# --------------------------------X->y effect--------------------------------- #
# significance 
anova(lm(y ~ sex + add + dom, data = df[df$infection == "PBS",]), lm(y ~ sex, data = df[df$infection == "PBS",])) # p = 0.07
anova(lm(y ~ sex + add + dom, data = df[df$infection == "SARSCoV",]), lm(y ~ sex, data = df[df$infection == "SARSCoV",])) # p = 2.6e-06
anova(lm(y ~ sex + add + dom, data = df[df$infection == "SARS2CoV",]), lm(y ~ sex, data = df[df$infection == "SARS2CoV",])) # p = 8.0e-05
# does the SARS relationship attenuate at all when accounting for CM DNTs? yes 
anova(lm(y ~ sex + Mtr + add + dom, data = df[df$infection == "SARSCoV",]), lm(y ~ sex + Mtr, data = df[df$infection == "SARSCoV",])) # p = 0.48
anova(lm(y ~ sex + Mtr + add + dom, data = df[df$infection == "SARS2CoV",]), lm(y ~ sex + Mtr, data = df[df$infection == "SARS2CoV",]))

summary(lm(y ~ sex + geno, data = df[df$infection == "SARSCoV",]))
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARSCoV",]))
summary(lm(y ~ sex + geno + Mtr, data = df[df$infection == "SARSCoV",]))

summary(lm(y ~ sex + geno, data = df[df$infection == "SARS2CoV",]))
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARS2CoV",]))
summary(lm(y ~ sex + geno + Mtr, data = df[df$infection == "SARS2CoV",]))

sig_df <- data.frame(group = strat_groups,
                     label = c("p = 0.07", "p = 2.6e-06", "p = 8.0e-05"),
                     y = c(Inf, Inf, Inf),
                     vjust = c(1.2, 2.6, 4.0),
                     hjust = c(1.1, 1.1, 1))
cmDNT_yxg <- pxg_group(cross, pheno = df$y, marker = marker,
                       geno.map = list(A = "CC006", B = "CC044"), xlab = expression(italic("Irq24")~"genotype"),
                       ylab = "Weight loss (AAC)", type = "boxplot",
                       groupcolumn = "infection", groupcolors = inf_pal, grouplabels = grouplabels) +
  geom_text(data = sig_df, 
            aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust, hjust = hjust), 
            size = 4, inherit.aes = FALSE, show.legend = FALSE) + 
  scale_color_manual(values = inf_pal, breaks = names(inf_pal), labels = grouplabels[names(inf_pal)]) +
  theme_minimal() + 
  theme(legend.position = "none")
cmDNT_yxg

# --------------------------------M->y effect--------------------------------- #
# significance
summary(lm(y ~ sex + Mtr, data = df[df$infection == "PBS",])) # p = 0.5
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARSCoV",])) # p = 1.2e-05
summary(lm(y ~ sex + Mtr, data = df[df$infection == "SARS2CoV",])) # p = 0.4
sig_df <- data.frame(group = strat_groups,
                     label = c("p = 0.5", "p = 1.2e-5", "p = 0.4"),
                     vjust = c(1.2, 2.6, 4.0))
cmDNT_wt <- ggplot(data = df, aes(x = Mtr, y = y, color = infection)) +
  geom_point() +
  scale_color_manual(values = inf_pal) +
  geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
  labs(x = expression("DN T cells (% CM) (transformed)"), y = "Weight loss (AAC)", color = NULL) +
  geom_text(data = sig_df, 
            mapping = aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust),
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  theme_minimal()
cmDNT_wt

# M effect on y in SARS1 mice only 
# cmDNT_wt <-  ggplot(data = df[df$infection == "SARSCoV",], aes(x = M, y = y, color = infection)) +
#   geom_point() +
#   scale_color_manual(values = inf_pal, name = NULL) +
#   geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
#   labs(x = expression("DN T cells (% CM)"), y = "Weight loss (AAC)") +
#   theme_minimal() +
#   annotate("text", x = Inf, y = Inf, label = "p = 0.000012", hjust = 1.05, vjust = 1.5, size = 4)
# cmDNT_wt

scat_legend <- get_legend(cmDNT_wt)
cmDNT_wt <- cmDNT_wt + 
  theme(legend.position = "none")

# for defense slides 
# ggplot(data = df[df$infection == "SARSCoV",], aes(x = M, y = y, color = infection)) +
#   geom_point() +
#   scale_color_manual(values = inf_pal, name = NULL) +
#   geom_smooth(method = "lm", se = FALSE) + # fullrange = TRUE
#   labs(x = expression("DN T cells (% CM)"), y = "Weight loss (AAC)") +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 14),
#         legend.position = "none")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter4/cmDNT_wt.png", width = 4, height = 4)

# ---------------------------------mediation---------------------------------- #
df <- na.omit(df)
Z <- model.matrix(~ sex + Titer, data = df)[,-1, drop = FALSE] 
M <- as.matrix(df$Mtr) # mediator matrix 
colnames(M) <- immune_pheno
X <- model.matrix(~ geno, data = df)[,-1, drop = FALSE] # exposure matrix 
group_med <- ifelse(df$infection == "SARS2CoV", 2, 
                    ifelse(df$infection == "SARSCoV", 1, 0))
# define weights for y and M to handle variance heterogeneity
df$y <- scale(df$y, center = TRUE, scale = TRUE)
fit_y <- lm(y ~ sex + infection, data = df)
res <- resid(fit_y)
s2_y0 <- var(res[group_med == 0])
s2_y1 <- var(res[group_med == 1])
s2_y2 <- var(res[group_med == 2])
wts_y <- ifelse(group_med == 2, 1/s2_y2, 
                ifelse(group_med == 1, 1/s2_y1, 1/s2_y0))
fit_m <- lm(Mtr ~ sex + infection, data = df)
res <- resid(fit_m)
s2_m0 <- var(res[group_med == 0])
s2_m1 <- var(res[group_med == 1])
s2_m2 <- var(res[group_med == 2])
wts_M <- ifelse(group_med == 2, 1/s2_m2, 
                ifelse(group_med == 1, 1/s2_m1, 1/s2_m0))

med_gxt <- bmediatR_GxT3(y = df$y, M = M, X = X, Z = Z, 
                         group = group_med, w_y = wts_y, w_M = wts_M,
                         align_data = FALSE,
                         options_X = list(sum_to_zero = FALSE, center = FALSE, scale = FALSE))

cmDNT_med <- plot_posterior_bar_gxt3(med_gxt, mediator_id = immune_pheno, grp_name = "Infection",
                                     grp_labels = c(grp1 = "Control", grp2 = "SARS-CoV", grp3 = "SARS-CoV-2"))
cmDNT_med

cmDNT_med <- cmDNT_med +
  guides(fill = "none") + 
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size=9), 
        legend.text = element_text(size=8),
        legend.position = c(0.5, 0.7))
cmDNT_med


# ---------------------------------Figure 4----------------------------------- #
# box_legend <- get_legend(cmDNT_pxg)
# scat_legend <- get_legend(cmDNT_wt)
# med_legend <- get_legend(cmDNT_med)
# 
# cmDNT_pxg <- cmDNT_pxg + theme(legend.position = "none")
# cmDNT_wt <- cmDNT_wt + theme(legend.position = "none")
# cmDNT_med <- cmDNT_med + theme(legend.position = "none")
cd86Lp_pxg <- cd86Lp_pxg + theme(legend.position = "none")
cd86Lp_wt <- cd86Lp_wt + theme(legend.position = "none")
# cd86Lp_med
helpN_pxg <- helpN_pxg + theme(legend.position = "none")
helpN_My <- helpN_My + theme(legend.position = "none")
helpN_med <- helpN_med + theme(legend.position = "none")

newFig4 <- plot_grid(
  ggdraw() + draw_label("Central memory double negative T cells (SARS-CoV)", fontface = "bold", x = 0.01, hjust = 0),
  plot_grid(NULL, cmDNT_pxg, 
            plot_grid(NULL, box_legend, scat_legend, NULL, nrow = 4, rel_heights = c(0.5,1,1,0.5), align = "v"), 
            cmDNT_wt, cmDNT_med, 
            ncol = 5, rel_widths = c(1, 1, 0.5, 1, 1), labels = c("A", "B", "", "C", "D")),
  ggdraw() + draw_label("Activated inflammatory monocytes (SARS-CoV-2)", fontface = "bold", x = 0.01, hjust = 0),
  plot_grid(NULL, cd86Lp_pxg, cd86Lp_wt, cd86Lp_med, ncol = 4, 
            rel_widths = c(1, 1, 1, 1.5), labels = c("E", "F", "G", "H")),
  ggdraw() + draw_label("Naive helper T cells", fontface = "bold", x = 0.01, hjust = 0),
  plot_grid(helpN_scan, helpN_pxg, helpN_My, helpN_med, ncol = 4, rel_widths = c(1.5, 1, 1, 1), labels = c("I", "J", "K", "L")),
  nrow = 6, rel_heights = c(0.15, 1, 0.15, 1, 0.15, 1)
)
ggsave("figures/Figure4.png", bg = "white", width = 16, height = 12)
