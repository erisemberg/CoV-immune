### Functions for QTL analysis (SARS and peanut allergy)
library(lme4) # for modeling random effects
library(MESS) # for AUC calculation 
library(plyr)
library(tidyverse)
#library(dplyr)
library(readr)
library(shades)
library(ranger)
library(scales)
library(extRemes)
library(RColorBrewer)
library(parallel)
library(snow)
library(car)
library(pzfx)
library(ggh4x)
source("code-dependencies/lmmultiresponse.R")
source("code-dependencies/plot_modX.R")

#-----------------------------to be categorized--------------------------------#
make_interpolate_fn <- function(df) {
  ord <- order(df$pval)
  p_sorted <- df$pval[ord]
  q_sorted <- df$qval[ord]
  
  approxfun(x = p_sorted, y = q_sorted, method = "linear")
}

# parse chr_pos column, example: "Chr 1: 168.50 (167.68-175.00)"
parse_chr_pos <- function(x) {
  # extract chromosome number
  chr <- sub("Chr ([0-9XY]+):.*", "\\1", x)
  # extract position
  pos <- as.numeric(sub(".*: ([0-9.]+) \\(.*", "\\1", x)) 
  # extract lower and upper bounds
  lwr <- as.numeric(sub(".*\\(([0-9.]+)-.*", "\\1", x))
  upr <- as.numeric(sub(".*-([0-9.]+)\\).*", "\\1", x))
  
  return(list(chr = chr, pos = pos, lwr = lwr, upr = upr))
}

make_lab_expr <- function(label) {
  m <- regexpr("[\u207A\u207B]", label)
  if (m[1] == -1) return(label)
  sup_char <- substr(label, m, m)
  sup_sym  <- if (sup_char == "\u207A") "+" else "-"
  prefix <- substr(label, 1, m-1)
  suffix <- substr(label, m+1, nchar(label))
  expr_txt <- sprintf('"%s"^"%s"~"%s"', prefix, sup_sym, suffix)
  return(parse(text = expr_txt))
}

test_nominal_sig <- function(chr, pos, pheno_name, sig_groups, cross_data, sig_thresh = 0.05) {
  dat <- data.frame(
    pheno = cross_data[[pheno_name]],
    sex = as.numeric(cross_data$sex),
    pgm = cross_data$pgm, 
    infection = cross_data$infection
  )
  
  # get genotype terms 
  marker <- find.marker(cross, chr = chr, pos = pos)
  if (chr == "X"){
    geno_terms <- get_Xchr_geno_terms(cross, marker)
    dat$ABf <- geno_terms$ABf
    dat$BB <- geno_terms$BB
    dat$BY <- geno_terms$BY
  } else {
    geno_terms <- get_add_dom(cross, marker)
    dat$add <- geno_terms$add
    dat$dom <- geno_terms$dom
  }
  
  # test in each group that's not already genome-wide significant
  group_lookup <- c("C" = "PBS", "1" = "SARSCoV", "2" = "SARS2CoV", "A" = "GxT")
  nominal_groups <- c()
  test_groups <- c("C", "1", "2", "A")
  already_sig <- gsub("[*]+", "", unlist(strsplit(sig_groups, ", ")))
  test_these <- setdiff(test_groups, already_sig)
  test_these <- group_lookup[test_these]
  
  for (group in test_these) {
    if (group == "PBS" & pheno_name == "Titer") { next }
    if (group == "GxT") { # combined analysis, test interaction model
      if (chr == "X"){
        fit_null <- lm(pheno ~ sex + infection + pgm, data = dat)
        fit_full <- lm(pheno ~ sex + infection + pgm + ABf + BB + BY + ABf:infection + BB:infection + BY:infection, data = dat)
      } else {
        fit_null <- lm(pheno ~ sex + infection, data = dat)
        fit_full <- lm(pheno ~ sex + infection + add + dom + add:infection + dom:infection, data = dat)
      }
    } else { # stratified analysis 
      dat_g <- dat[cross_data$infection == group,]
      if (chr == "X"){
        fit_null <- lm(pheno ~ sex, data = dat_g)
        fit_full <- lm(pheno ~ sex + ABf + BB + BY, data = dat_g)
      } else {
        fit_null <- lm(pheno ~ sex, data = dat_g)
        fit_full <- lm(pheno ~ sex + add + dom, data = dat_g)
      }
    }
    p_val <- anova(fit_full, fit_null)$`Pr(>F)`[2]
    message(paste(pheno_name, "in", group, ":", p_val))
    
    # Add to nominal list if p < sig_thresh
    if (!is.na(p_val) && p_val < sig_thresh) {
      group_code <- case_when(
        group == "PBS" ~ "C",
        group == "SARSCoV" ~ "1",
        group == "SARS2CoV" ~ "2",
        group == "GxT" ~ "A"
      )
      nominal_groups <- c(nominal_groups, group_code)
    }
  }
  
  return(paste(nominal_groups, collapse = ", "))
}

#-------------------------------Miscellaneous----------------------------------#
which_chr <- function(cross, marker) {
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o] 
  return(chr)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# has_xchr_perms: indicates whether permutation object has x-chr-specific permutations  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
has_xchr_perms <- function(perm){
  return('xchr' %in% names(attributes(perm)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# add_attributes: add an attribute to an object  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_attribute <- function(obj, attr_name, attr_value){
  attr(obj, attr_name) <- attr_value
  return(obj)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# logit transformation (for proportions/percentages)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
logit <- function(p){log(p/(1-p))}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# RINT transformation from miQTL
#+++++++++++++++++++++++++++++++++++++++++++++++++++
rint <- function(phenotype, prop=0.5){
  rint_phenotype <- qnorm((rank(phenotype, na.last="keep")-prop)/sum(!is.na(phenotype)))
  return(rint_phenotype)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Will's truncated-RINT transformation 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
trint <- function(y, theta=0.01){
  p <- theta + (rank(y)-1)*(1-2*theta)/(length(y)-1)
  qnorm(p)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to ensure directory exists before creating files in it 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create logger function 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_logger <- function(filename, sep="\n"){
  if(file.exists(filename)){
    file.remove(filename);
  }
  function(...){
    text <- sprintf(...);
    cat(text, file=filename, sep=sep, append=T);
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to add phenotype data to R/qtl data table 
# Input:
#   pheno.name = name of phenotype column in phenotype spreadsheet 
#   new.col = empty column 
#   col.pos = position of new column in R/qtl table 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_pheno <- function(pheno.name, new.col, col.pos){
  rqtl <- rqtl %>% 
    add_column(placeholder.name = new.col, .before = col.pos)
  names(rqtl)[names(rqtl) == "placeholder.name"] <- pheno.name
  for (i in 3:nrow(rqtl)){
    rqtl[[pheno.name]][i] <- ifelse(rqtl$mouse_ID[i] %in% pheno$Geno_ID, 
                                    pheno[[pheno.name]][which(pheno$Geno_ID == rqtl$mouse_ID[i])], NA)
  }
  return(rqtl)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# readGRMBin: R script to read the GRM binary file
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
readGRMBin <- function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N, grm=grm))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# fill_in_res:  when a phenotype is missing values (e.g. flow data in SARS 
#               experiments), residuals(fit) will be shorter than the original 
#               phenotype vector/column. This function expands the residuals 
#               vector to the original length, filling in the missing spaces 
#               with NAs. 
# Input:
#     pheno: vector of phenotype values (original length)
#     fit: fit with residuals to use in new data vector 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
fill_in_res <- function(pheno, fit){
  missing <- which(is.na(pheno))
  has_data <- !(seq_along(pheno) %in% missing)
  
  r <- residuals(fit)
  
  tmp <- rep(NA, length(pheno))
  tmp[has_data] <- r
  
  return(tmp)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_summary_for_print: function for extracting p-value from MANOVA 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_summary_for_print <- car:::print.Anova.mlm
body(get_summary_for_print) <- 
  local({
    tmp <- body(get_summary_for_print) 
    tmp[-(length(tmp)-(0:1))]
  })


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# load_themes: loads ggplot themes 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
load_themes <- function(){
  # To find default colors
  #show_col(hue_pal()(2))
  
  poster_theme <<- theme(axis.title.x = element_text(size = 20), 
                   axis.text.x = element_text(size = 18), 
                   axis.title.y = element_text(size = 20),
                   axis.text.y = element_text(size = 18),
                   legend.title = element_text(size = 20),
                   legend.text = element_text(size = 18),
                   plot.title = element_text(size = 22))
  
  poster_theme2 <<- theme(axis.title.x = element_text(size = 22), 
                    axis.text.x = element_text(size = 20), 
                    axis.title.y = element_text(size = 22),
                    axis.text.y = element_text(size = 20),
                    legend.title = element_text(size = 22),
                    legend.text = element_text(size = 20),
                    plot.title = element_text(size = 24))
  
  # For things needing big text  
  big_theme <<- theme(axis.title = element_text(size = 24), 
                      axis.text = element_text(size = 22),
                      plot.title = element_text(size = 24),
                      legend.title = element_text(size = 22),
                      legend.text = element_text(size = 18))
  
  # Same but white background
  bw_biggish_theme <<- theme(axis.title = element_text(size = 20), 
                      axis.text = element_text(size = 18),
                      plot.title = element_text(size = 22),
                      legend.title = element_text(size = 20),
                      legend.text = element_text(size = 18),
                      legend.key = element_rect(fill = "white"),
                      panel.background = element_rect(fill = "white",
                                                      colour = "white"),
                      panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                      panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                      panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
  # Slightly bigger text
  bw_big_theme <<- theme(axis.title = element_text(size = 22), 
                         axis.text = element_text(size = 20),
                         plot.title = element_text(size = 24),
                         legend.title = element_text(size = 20),
                         legend.text = element_text(size = 18),
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_rect(fill = "white",
                                                         colour = "white"),
                         panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                         panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                         panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
  # For things needing big text but smaller x-axis labels   
  # (like HS and titer plots with x-axis = strain)
  big_theme2 <<- theme(axis.title = element_text(size = 24), 
                      axis.text.y = element_text(size = 22),
                      axis.text.x = element_text(size = 18),
                      plot.title = element_text(size = 24),
                      legend.title = element_text(size = 22),
                      legend.text = element_text(size = 20))
  
  # same but white bg
  bw_big_theme2 <<- theme(axis.title = element_text(size = 24), 
                       axis.text.y = element_text(size = 22),
                       axis.text.x = element_text(size = 18),
                       plot.title = element_text(size = 24),
                       legend.title = element_text(size = 22),
                       legend.text = element_text(size = 20),
                       legend.key = element_rect(fill = "white"),
                       panel.background = element_rect(fill = "white",
                                                       colour = "white"),
                       panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                       panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                       panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
  # side-by-side PxG theme 
  sbs_pxg_theme <<- theme(plot.title = element_text(size = 26),
                          axis.text.y = element_text(size = 18),
                          axis.text.x = element_text(size = 17),
                          axis.title = element_text(size = 24))
  
  # For PxG plots needing big text  
  big_pxg_theme <<- theme(axis.title = element_text(size = 28), #26
                      axis.text = element_text(size = 26),
                      legend.title = element_text(size = 26),
                      legend.text = element_text(size = 24))
  
  # themes for Rmarkdown
  rmd_theme <<- theme(axis.title.x = element_text(size = 16), 
                    axis.text.x = element_text(size = 14), 
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14),
                    plot.title = element_text(size = 18))
  
  ppt_theme <<- theme(axis.title.x = element_text(size = 14), 
                      axis.text.x = element_text(size = 12), 
                      axis.title.y = element_text(size = 14),
                      axis.text.y = element_text(size = 12),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12),
                      plot.title = element_text(size = 16))
  
  # theme for publication (multiple plots per figure)
  pub_theme <<- theme(legend.key.size = unit(0.03, 'npc'), # legend key is 0.03 of the plot size (basically)
                      axis.title = element_text(size = 8), 
                      axis.text = element_text(size = 6), 
                      legend.title = element_text(size = 8),
                      legend.text = element_text(size = 6),
                      legend.box.spacing = unit(0, 'pt')) # no space between plot and legend 
  
  # theme for publication (one plot)
  pub_theme2 <<- theme(axis.title = element_text(size = 22), 
                     axis.text = element_text(size = 20), 
                     plot.title = element_text(size = 24))
  
  # theme for publication PxG plots 
  pub_pxg_theme <<- theme(axis.title = element_text(size = 22), 
                          axis.text = element_text(size = 18))
  # add axis.title.y = element_text(hjust = 0) or something to move y-axis away from plot
  
}



#-----------------------------Phenotype plotting-------------------------------#
load_themes()

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# trajectory_plot: plots trajectory data from a r/QTL cross object
# Input:
#     dat: dataframe with data to plot 
#     steps: x-axis data 
#     phenos: y-axis data - can either be the name of a column in dat, 
#             or an integer which will be repeated as the phenotype. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
cov_trajectory_plot <- function(dat, phenos, title = NULL, ylab = NULL, 
                                xlab = "Days post-infection",
                                incl.parents = TRUE, parent.lty = 1, ylim = NULL,...){
  # load parent data 
  if (incl.parents == TRUE){
    CC_weight <- read_pzfx("source_data/screen/CC_mice.pzfx", "CC_weight")
    CC006_weight <- CC_weight[,c(10:13)] #16750
    CC044_weight <- CC_weight[,c(22:25)] #4410
    CC006dat <- data.frame(mouse_ID = 'CC006', 
                           dpi = 0:4, 
                           pct_weight = rowMeans(CC006_weight, na.rm=T))
    CC044dat <- data.frame(mouse_ID = 'CC044', 
                           dpi = 0:4, 
                           pct_weight = rowMeans(CC044_weight, na.rm=T))
  }
  
  # data setup for ggplot 
  df <- dat$mouse_ID
  for (i in 1:length(phenos)){
    df <- cbind(df, dat[phenos[i]])
  }
  df <- as.data.frame(df)
  colnames(df) <- c('mouse_ID', '0', '1', '2', '3', '4')
  df <- pivot_longer(df, cols = c(2:6), names_to = "dpi", values_to = "pct_weight") %>%
    mutate(dpi = as.integer(dpi)) %>%
    group_by(mouse_ID) %>%
    filter(n_distinct(dpi) == length(0:4), all(!is.na(pct_weight))) %>% 
    ungroup()
  
  if (incl.parents == TRUE){
    wtloss <- ggplot(df, aes(x = dpi, y = pct_weight)) +
      geom_line(aes(group = mouse_ID, color = "F2", linetype = "F2", size = "F2")) +
      geom_line(data = CC006dat, aes(x = dpi, y = pct_weight, group = mouse_ID, color = "CC006", linetype = "CC006", size = "CC006")) + 
      geom_line(data = CC044dat, aes(x = dpi, y = pct_weight, group = mouse_ID, color = "CC044", linetype = "CC044", size = "CC044")) +
      scale_linetype_manual(values = c("F2" = 1, "CC006" = parent.lty, "CC044" = parent.lty)) +
      scale_size_manual(values = c("F2" = 0.2, "CC006" = 1, "CC044" = 1)) + # previously 1.3
      scale_color_manual(values = c("F2" = "gray63", "CC006" = "#E41A1C", "CC044" = "#377EB8")) +
      {if(!is.null(ylim)) ylim(ylim) } +
      labs(x = xlab, y = ylab, title = title, color = "Legend", linetype = "Legend", size = "Legend") +
      scale_x_continuous(breaks = 0:4, expand = c(0.03,0.03))
  } else { # don't include parent trajectories 
    wtloss <- ggplot(df, aes(x = dpi, y = pct_weight)) +
      geom_line(aes(group = mouse_ID), color = "gray63", size = 0.2) +
      {if(!is.null(ylim)) ylim(ylim) } +
      labs(x = xlab, y = ylab, title = title) +
      scale_x_continuous(breaks = 0:4, expand = c(0.03,0.03))
  }

  return(wtloss)
}


#--------------------------Calculate derived measures--------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# calc_auc: this function calculates an "area under the curve" metric, which is 
# really the area above the curve and below the horizontal line at the first data 
# point (updates cross object)
# Input:
#     cross: r/qtl cross object 
#     col.name: column name to be used for the new aggregate phenotype 
#     steps: x-axis data 
#     phenos: y-axis data - phenotypes to be used to calculate the new aggregate 
#             phenotype
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc_auc <- function(dat, steps, col.name, phenos, is_df = FALSE){
  if (is_df){
    df <- dat
  } else { # dat is a cross object, extract phenotype data 
    df <- dat$pheno
  }
  
  df[,col.name] <- rep(NA,nrow(df))
  n <- nrow(df)
  
  for (i in 1:n){
    line <- df[i,phenos]
    auc <- auc(x = steps, y = line)
    aac <- line[1]*tail(steps, n=1) - auc
    df[i,col.name] <- aac 
  }
  
  if (is_df){
    dat <- df
  } else { # dat is a cross object, insert updated phenotype data 
    dat$pheno <- df 
  }
  
  return(dat)
}



#----------------------------Working with QTL objects--------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_models: this function creates all single-QTL models specified in the 
#                models table
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes to create single-QTL 
#             models for 
#     covar: dataframe containing covariates 
#     method: whether to use existing R/qtl functions ("rqtl", default);
#             linear mixed modeling code ("lmm"); or
#             gene-by-treatment effect modeling code ("gxt")
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_models <- function(cross.obj, models, covar = NULL, use = "rqtl", 
                          trt = NULL, mod.dir){
  ensure_directory(mod.dir)
  
  for (i in 1:nrow(models)){
    mod_name <- as.character(models[i,'obj'])
    filepath = paste(c(mod.dir, mod_name, ".Rdata"), collapse = "")
    if (file.exists(filepath)){ # model exists already, load .Rdata object 
      print(paste("Loading scanone object from", filepath))
      load(filepath, envir = .GlobalEnv)
    } else { # file doesn't exist, create model 
      print(paste("Creating scanone object for", models[i,'name'], "using", models[i,'type'], "model."))
      
      pheno.col <- models$colname[i]
      mod.type <- models$type[i]
      
      if ((use == "lmm") & (mod.type == "normal")){ # use custom haley-knott function that handles random effect covariates
        assign(mod_name, hk(cross.obj, pheno.col = pheno.col, envir = .GlobalEnv))
        
      } else if ((use == "gxt") & (mod.type == "normal")){ # use custom gene-by-treatment qtl mapping function
        # define covariates, if specified
        covlist <- str_split(models$cov[i], ",")[[1]]
        addcovar <- NULL
        if(!is.na(models$cov[i])){
          addcovar <- covar[,covlist]
        }
        
        assign(mod_name, gxt(cross.obj, pheno.col = pheno.col, addcovar = addcovar, trt = trt), envir = .GlobalEnv)
      } else { # use rqtl::scanone
        # define covariates, if specified
        covlist <- str_split(models$cov[i], ",")[[1]]
        addcovar <- NULL
        if(!is.na(models$cov[i])){
          addcovar <- covar[,covlist]
        }
        
        # create model 
        assign(mod_name, 
               scanone(cross.obj, pheno.col = pheno.col, model = mod.type, method = "hk", addcovar = addcovar), 
               envir = .GlobalEnv) 
      }
      
      print(paste("Saving object to", filepath))
      save(list = mod_name, file = filepath)
    }
  }
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# summary.permgev: summary() function for permutation object
#+++++++++++++++++++++++++++++++++++++++++++++++++++
summary.permgev <- function(object, ...){
  # Calculate 5% and 10% significance thresholds 
  thresholds <- matrix(qevd(p = c(0.95, 0.90), 
                            loc = object$results$par[1], 
                            scale = object$results$par[2], 
                            shape = object$results$par[3]))
  colnames(thresholds) <- c('lod')
  rownames(thresholds) <- c('5%', '10%')
  
  return(thresholds)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_perms: this function creates permutation objects for all single-QTL models 
#               specified in the models table 
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes with single-QTL models 
#     covar: covariates 
#     perm.dir: directory with existing perm objects / to save perm objects to 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_perms <- function(cross.obj, models, perm.dir, n.perm = 1000, gxt = FALSE,
                         perm.strata = NULL, covar = NULL, perm.Xsp = FALSE, ...){

  ensure_directory(perm.dir)
  
  #--------------------Create or load permutation object-----------------------#
  for (i in 1:nrow(models)){
    perm_name <- models$perm.obj[i]
    perm.fname <- gsub("\\.", "", perm_name) # remove period from permutation object name, if present
    permfile <- paste(perm.dir, perm.fname, ".Rdata", sep = "")
    
    # create or load permutation object
    if(file.exists(permfile)){ # permutation object already exists, load from Rdata file 
      load(permfile, envir = .GlobalEnv)
      print(paste("Loaded", perm_name, "from file."))
    } else { # create permutation object 
      print(paste("Creating", perm_name, "..."))
      
      covlist <- str_split(models$cov[i], ",")[[1]]
      addcovar <- NULL
      if(!is.na(models$cov[i])){
        addcovar <- covar[,covlist]
      }
      
      if (gxt == TRUE & models$type[i] == 'normal'){
        assign(perm_name,
               gxt_perm(cross.obj, pheno.col = models$colname[i], covar = addcovar, 
                        n.perm = n.perm, trt = 'infection'),
               envir = .GlobalEnv)
      } else {
        assign(perm_name, 
               scanone(cross.obj, pheno.col = models$colname[i], method = "hk",
                       model = models$type[i], addcovar = addcovar, n.perm = n.perm,
                       perm.strata = perm.strata, perm.Xsp = perm.Xsp, n.cluster = 4), 
               envir = .GlobalEnv)
      }

      save(list = perm_name, file = permfile)
    }

    
    #------------------------Fit GEV distribution(s)---------------------------#  
    if (gxt == TRUE & models$type[i] == 'normal'){
      gevname <- paste0(perm_name, '.gev')
      print(paste("Creating GEV fit", gevname, "..."))
      fitgev <- fevd(as.numeric(get(perm_name)), type = "GEV")
      assign(gevname, fitgev, envir = .GlobalEnv)
      # add custom 'permgev' class so that we can use custom summary function to get GEV-defined thresholds  
      assign(gevname,
             add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
             envir = .GlobalEnv)
    } else if (!has_xchr_perms(get(perm_name))){ # no X-chr-specific thresholds 
      # 1 GEV distribution
      gevname <- paste0(perm_name, '.gev')
      print(paste("Creating GEV fit", gevname, "..."))
      fitgev <- fevd(as.numeric(get(perm_name)), type = "GEV")
      assign(gevname, fitgev, envir = .GlobalEnv)
      # add custom 'permgev' class so that we can use custom summary function to get GEV-defined thresholds  
      assign(gevname,
             add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
             envir = .GlobalEnv)
    } else if ((has_xchr_perms(get(perm_name))) & (ncol(get(perm_name)$A) == 1)) { # X-chr thresholds, 1 lod column
      # 2 GEV distributions - 1 for autosomes, 1 for X-chr
      for (chr in c('A', 'X')){
        gevname <- paste0(perm_name, chr, '.gev')
        print(paste("Creating GEV fit", gevname, "..."))
        fitgev <- fevd(as.data.frame(get(perm_name)[[chr]])$lod, type = "GEV")
        assign(gevname, fitgev, envir = .GlobalEnv)
        # add custom 'permgev' class for use of custom summary function to get GEV-defined thresholds  
        assign(gevname,
               add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
               envir = .GlobalEnv)
      }
    } else if ((has_xchr_perms(get(perm_name))) & (ncol(get(perm_name)$A) > 1)){ # X-chr specific thresholds, >1 LOD column (e.g. 2part model)
      # 6 GEV distributions - 3 for autosomes (1 for each LOD score); 3 for X-chr
      for (chr in c('A', 'X')){
        for (lod.type in colnames(get(perm_name)[[chr]])){
          gevname <- paste0(perm_name, chr, '.', lod.type, '.gev')
          print(paste("Creating GEV fit", gevname, "..."))
          fitgev <- fevd(as.data.frame(get(perm_name)[[chr]])[[lod.type]], type = "GEV")
          assign(gevname, fitgev, envir = .GlobalEnv)
          # add custom 'permgev' class for use of custom summary function to get GEV-defined thresholds  
          assign(gevname,
                 add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
                 envir = .GlobalEnv)
        }
      }
    }
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# doc_peaks: this function finds significant LOD peaks, calculates 95% Bayes credible 
#           intervals for those peaks, and creates a summary table with data for 
#           each one, including marker, chromosome, position, LOD and positions of 
#           Bayes credible intervals.
# Input: 
#     models: table containing information about phenotypes with single-QTL models 
#     sig.level: significance level (can be 0.05 or 0.10)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
doc_peaks <- function(models, gxt = FALSE, sig.level = 0.10){
  
  #--------------------Make sure all permutations exist------------------------#
  perms.exist = T
  for (m in 1:nrow(models)){
    if (!exists(models$perm.obj[m])){
      print(paste(models$perm.obj[m], "does not exist. All permutation objects must exist before attempting to identify significant peaks."))
      perms.exist = F
    }
  }
  if (perms.exist == F){stop("Missing permutation objects, execution halted.")}
  
  
  #-------------------------Count significant peaks----------------------------# 
  num.peaks = 0
  num.peaks.per.mod <- rep(NA, nrow(models))
  gevnames <- vector(mode = 'list', length = nrow(models))
  thresholds <- vector(mode = 'list', length = nrow(models))
  sig.ix <- ifelse(sig.level == 0.10, 2, ifelse(sig.level == 0.05, 1, NA)) ### Handle more than just 0.05 and 0.10?
  for (m in 1:nrow(models)){
    
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    
    if (gxt == TRUE & models$type[m] == 'normal'){
      ### good enough for now 
      gevnames[[m]] <- paste0(perm, '.gev')
      thresholds[[m]] <- summary(get(gevnames[[m]]))[sig.ix]
      num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]], lodcolumn = 3))
      
    } else if (!has_xchr_perms(get(perm))){ # no X-chr specific thresholds (e.g. backcross)
      
      # not using GEV 
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      ### Using GEV - copy from code on peanut GitHub 
      # gevnames[[m]] <- paste0(perm, '.gev') # GEV dist fit from perm object
      ### this method of defining length 20 thresholds won't work after Karl fixes bug 
      #thresholds[[m]] <- summary(get(gevnames[[m]]))[sig.ix]
      #num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      next 
      
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) == 1)) { # X-chr specific thresholds, 1 LOD 
      
      # not using GEV 
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      # using GEV but not handling X-chr specific thresholds (fine for virus research paper bc no X-chr QTL)
      gevnames[[m]] <- c(paste0(perm, 'A.gev'), paste0(perm, 'X.gev'))
      thresholds[[m]] <- summary(get(gevnames[[m]][1]))[sig.ix]
      num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      
      ### Handle X-chr specific thresholds 
      # this method of defining thresohlds won't work after Karl fixes bug 
      # thresholds[[m]] <- c(rep(summary(get(gevnames[[m]][1]))[sig.ix], 19), 
      #                    summary(HSpermX.gev)[sig.ix])
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) > 1)) { # X-chr specific thresholds >1 LOD (2part)
      
      # not using GEV 
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      gevnames[[m]] <- c(paste0(perm, c(rep('A.',3), rep('X.', 3)), colnames(get(perm)$A), '.gev'))
      ### Need to define thresholds[[m]] in order for later code in doc_peaks to work
      ### This threshold code doesn't work though 
      # thresholds[[m]] <- data.frame(lod.p.mu = c(rep(summary(get(gevnames[[m]][1]))[1], 19), summary(get(gevnames[[m]][4]))[1]),
      #                   lod.p = c(rep(summary(get(gevnames[[m]][2]))[1], 19), summary(get(gevnames[[m]][5]))[1]),
      #                   lod.mu = c(rep(summary(get(gevnames[[m]][3]))[1], 19), summary(get(gevnames[[m]][6]))[1]))
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
    } else {
      print("no conditions met")
    }
  }
  num.peaks <- sum(num.peaks.per.mod)
  
  
  #---------------------Get data for significant peaks-------------------------# 
  peak.data <- c('model', 'marker', 'chr', 'pos', 'lod', 'Bayes CI')
  peaks <- matrix(ncol = length(peak.data), nrow = num.peaks)
  colnames(peaks) <- peak.data
  
  row.ix = 0 # row in peaks table
  for (m in 1:nrow(models)){ # for each phenotype / model
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    #gevs <- gevnames[[m]]
    
    if (num.peaks.per.mod[m] == 0){
      print(paste("For model", models[m,'obj'], ", no LOD peaks above threshold"))
    } else {
      for (p in 1:num.peaks.per.mod[m]){ # for each peak 
        row.ix = row.ix + 1 
        
        # not using GEV 
        chrom = summary(get(mod), alpha = 0.10, perms = get(perm))[p,1]
        
        # using GEV
        # if (gxt == TRUE & models$type[m] == 'normal'){
        #   ### only looking at peaks significant for G + GxT at the moment 
        #   chrom = summary(get(mod), threshold = thresholds[[m]], lodcolumn = 3)[p,1]
        # } else if (ncol(get(perm)$A) == 1){ # not GxT normal, one threshold 
        #   chrom = summary(get(mod), threshold = thresholds[[m]])[p,1]
        # } else if (ncol(get(perm)$A) > 1){ # not GxT normal, multiple thresholds
        #   chrom = summary(get(mod), alpha = 0.10, perms = get(perm))[p,1]
        #   ### Fix and use below code to use GEV 
        #   #chrom = summary(get(mod), threshold = summary(get(gevs[[1]]))[sig.ix])[p,1]
        # }
        
        # bayesint() produces a table with 3 rows (lower interval, peak, upper interval), columns = chr, pos, lod 
        bayesCI <- bayesint(get(models$obj[m]), chr = chrom, prob = 0.95)
        # Populate peaks dataframe 
        peaks[row.ix,] <- c(mod, rownames(bayesCI)[2], chrom, bayesCI[2,2], bayesCI[2,3], 
                            paste(bayesCI[1,2],"-",bayesCI[3,2]))
      } 
    }
  }
  
  peaks <- as.data.frame(peaks)
  
  # Remove leading/trailing white space 
  for (i in 1:ncol(peaks)){
    peaks[,i] <- trimws(peaks[,i], which = c("both"))
  }
  
  return(peaks)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_unadj_pval: this function gets the unadjusted p-value from the linear regression
#           between a particular phenotype and the genotype at a particular marker.
# Input: 
#     cross: 
#     pheno.col: 
#     marker:
#     covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_unadj_pval <- function(cross, pheno.col, marker, covar = NULL){
  pheno <- pull.pheno(cross, pheno.col) ### COMMENT THIS OUT??
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  if (!is.null(covar)){
    covlist <- str_split(covar, ",")[[1]]
  }
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar), "", " + "),
                    paste(covlist, collapse = " + "))
  
  if (is.null(covar)){
    lmdata <- data.frame(cross$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(cross$pheno[,c(pheno.col,covlist)], marker.geno)
  }
  
  lmod <- lm(formula, lmdata)
  unadj_pval <- summary(lmod)$coefficients['marker.geno', 'Pr(>|t|)']
  return(unadj_pval)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_add_dom: given a cross and marker, gets additive and dominance terms at that marker 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_add_dom <- function(cross, marker){
  # calculate genotype probabilities if necessary 
  if (is.null(cross$geno[[1]]$prob)){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    message("Running calc.genoprob...")
  }
  # identify X chromosome 
  geno_class <- sapply(cross$geno, class)
  autosome_idx <- which(geno_class == "A")
  chr <- find.markerpos(cross, marker)[1,1]
  
  # calculate dosage if necessary
  if (is.null(cross$geno[[1]]$dos)){ 
    # calculate dosage of B allele (additive effect)
    calc.dosage <- function(cross){
      # autosomes only (construct explicit X design in marker loop)
      for (c in autosome_idx){ 
        prob <- cross$geno[[c]]$prob 
        cross$geno[[c]]$dos <- prob[,,2] + 2*prob[,,3] # 0*P(AA) + 1*P(AB) + 2*P(BB)
      }
      return(cross)
    }
    cross <- calc.dosage(cross) 
    message("Running calc.dosage...")
  }
  
  add <- cross$geno[[chr]]$dos[,marker]
  dom <- cross$geno[[chr]]$prob[,marker,2]
  
  return(list(add = add, dom = dom))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_Xchr_geno_terms: given a cross and marker, gets X chromosome genotype terms at that marker  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_Xchr_geno_terms <- function(cross, marker){
  # calculate genotype probabilities if necessary 
  if (is.null(cross$geno[[1]]$prob)){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    message("Running calc.genoprob...")
  }
  # identify X chromosome 
  geno_class <- sapply(cross$geno, class)
  autosome_idx <- which(geno_class == "A")
  chr <- find.markerpos(cross, marker)[1,1]
  
  # calculate dosage if necessary
  if (is.null(cross$geno[[1]]$dos)){ 
    # calculate dosage of B allele (additive effect)
    calc.dosage <- function(cross){
      # autosomes only (construct explicit X design in marker loop)
      for (c in autosome_idx){ 
        prob <- cross$geno[[c]]$prob 
        cross$geno[[c]]$dos <- prob[,,2] + 2*prob[,,3] # 0*P(AA) + 1*P(AB) + 2*P(BB)
      }
      return(cross)
    }
    cross <- calc.dosage(cross) 
    message("Running calc.dosage...")
  }
  
  p1 <- cross$geno[[chr]]$prob[,marker,1]
  p2 <- cross$geno[[chr]]$prob[,marker,2]
  sexM <- (cross$pheno$sex == "M")
  sexF <- (cross$pheno$sex == "F")
  fwdF <- sexF & (cross$pheno$pgm == 0) # female, forward direction (pgm = 0)
  revF <- sexF & (cross$pheno$pgm == 1) # female, reverse direction (pgm = 1)
  AA <- as.numeric(fwdF)*p1
  ABf <- as.numeric(fwdF)*p2
  ABr <- as.numeric(revF)*p2
  AB <- ABf + ABr
  BB <- as.numeric(revF)*p1
  AY <- as.numeric(sexM)*p1
  BY = as.numeric(sexM)*p2
  
  return(list(AA = AA, AB = AB, ABf = ABf, ABr = ABr, BB = BB, AY = AY, BY = BY))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_var_expl: this function gets the variance explained of a particular 
#               phenotype by the genotype at a particular marker.
# Input: 
#     cross: 
#     pheno.col: 
#     marker:
#     covar: 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_var_expl <- function(cross, pheno, marker, cov_list = c("sex"), gxt = FALSE){
  # calculate genotype probabilities if necessary 
  if (is.null(cross$geno[[1]]$prob)){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    message("Running calc.genoprob...")
  }
  # identify X chromosome 
  geno_class <- sapply(cross$geno, class)
  autosome_idx <- which(geno_class == "A")
  chr <- find.markerpos(cross, marker)[1,1]
  
  # calculate dosage if necessary
  if (is.null(cross$geno[[1]]$dos)){ 
    # calculate dosage of B allele (additive effect)
    calc.dosage <- function(cross){
      # autosomes only (construct explicit X design in marker loop)
      for (c in autosome_idx){ 
        prob <- cross$geno[[c]]$prob 
        cross$geno[[c]]$dos <- prob[,,2] + 2*prob[,,3] # 0*P(AA) + 1*P(AB) + 2*P(BB)
      }
      return(cross)
    }
    cross <- calc.dosage(cross) 
    message("Running calc.dosage...")
  }

  lmdata <- data.frame(pheno = pheno, cross$pheno[,cov_list, drop = FALSE])
  lmdata$geno <- cross$geno[[chr]]$data[,marker]
  if (chr %in% autosome_idx){
    lmdata$add <- cross$geno[[chr]]$dos[,marker]
    lmdata$dom <- cross$geno[[chr]]$prob[,marker,2]
    # define LM formulas 
    formula0 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "))
    if (gxt){ 
      formula1 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "), " + add + dom + add*infection + dom*infection")
    } else {
      formula1 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "), " + add + dom")    
    }
  } else {
    p1 <- cross$geno[[chr]]$prob[,marker,1]
    p2 <- cross$geno[[chr]]$prob[,marker,2]
    sexM <- (lmdata$sex == "M")
    sexF <- (lmdata$sex == "F")
    fwdF <- sexF & (lmdata$pgm == 0) # female, forward direction (pgm = 0)
    revF <- sexF & (lmdata$pgm == 1) # female, reverse direction (pgm = 1)
    lmdata$ABf <- as.numeric(fwdF)*p2
    lmdata$BB <- as.numeric(revF)*p1
    lmdata$BY = as.numeric(sexM)*p2
    # define LM formulas 
    formula0 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "))
    if (gxt){
      formula1 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "), " + ABf + BB + BY + ABf*infection + BB*infection + BY*infection")
    } else {
      formula1 <- paste0("pheno ~ ", paste(cov_list, collapse = " + "), " + ABf + BB + BY")
    }
  } 
  
  fit0 <- lm(formula0, lmdata)
  fit1 <- lm(formula1, lmdata)
  
  pheno_mean <- mean(pheno, na.rm = TRUE)
  SStotal <- sum((pheno - pheno_mean)^2, na.rm = TRUE) # equal to n-1 if data is centered/scaled 
  SSqtl <- deviance(fit0) - deviance(fit1)
  h2 <- SSqtl/SStotal
  
  # get effect size if on autosome and not GxT 
  if (chr != "X" & !gxt){ 
    add_effect <- coef(fit1)["add"] 
  } else {
    add_effect <- NA 
  }  
  
  return(list(h2 = h2,
              add_effect = add_effect))
}



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




#-----------------------------------Plotting-----------------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
#
# ### NEED TO UPDATE THIS CODE TO USE GEV AND PLOT SEPARATE LINES FOR X-CHR
# ### NOT IMPORTANT FOR PNUT PAPER BC GEV THRESHOLDS ARE NOT SIGNIFICANTLY 
# ### DIFFERENT EXCEPT IN THE CASE OF X-CHR-SPECIFIC PERMUTATIONS  
# 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans: this function plots genome scans for all models specified in the 
#             models table. Adds a solid line for 5 % significance threshold and 
#             a dotted line for 10 % significance threshold, if the permutation
#             object exists. 
# Input:
#     models: table containing information about phenotypes with single-QTL models 
#     ylim: y limits, if they need to be consistent across plots. If not set, will 
#           use the default. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_scans <- function(models, save = FALSE, save.dir = NULL, ...){
  
  if(!is.null(save.dir)) { ensure_directory(save.dir) }
  
  for (m in 1:nrow(models)){
    mod_name <- models$obj[m]
    mod <- get(mod_name)
    
    col = c('black', 'blue', 'red')
    lwd = c(2,2,2)
    if (ncol(mod)>3) { # if ncol>3, multi-LOD model
      lodcolumn=1:3
      if (names(mod)[3] == 'lod1'){ # GxT model
        # col <- c('#E66100', '#5D3A9B', 'black')
        col <- c("#D95F02","#1B9E77","black")
        #col <- c("#7570b3","#1B9E77","black")
        lwd <- c(3,3,1.5)
      }   
    } else { # single LOD model 
      lodcolumn=1
    }  
   
    if (save == TRUE){ # saving, so open connection to png and make axis labels bigger 
      fname <- paste0(save.dir, mod_name, '.png')
      png(fname, width = 750)
      plot_modX(mod, ylab = "", xlab = "", main = models$name[m], lodcolumn = lodcolumn, lwd = lwd,
           col = col, alternate.chrid = T, bandcol = "gray80", cex.main = 2, cex.axis = 2, ...) # cex.main was 1.8, cex.axis was 1.5
      title(ylab = "LOD", line = 2.5, cex.lab = 2)
      title(xlab = "Chromosome", cex.lab = 2.3, line = 3) # cex.lab was 2
    } else { # not saving, just printing
      plot(mod, ylab = "LOD", main = models$name[m], lodcolumn = lodcolumn, lwd = lwd, 
           col = col, alternate.chrid = T, bandcol = "gray90", ...)
    }
    
    # Add significance thresholds 
    perm <- models$perm.obj[m]

    if (exists(perm)){
      if (!has_xchr_perms(get(perm)) & ncol(mod)==3){ # no X-chr-specific perms, 1 lod column
        permgev <- paste0(perm, '.gev')
        perm.sum <- summary(get(permgev))
        abline(h = c(perm.sum[1], perm.sum[2]), lty=1:2)
      } else if (!has_xchr_perms(get(perm)) & names(mod)[3] == "lod.p.mu"){ # no X-chr-specific perms, 2part model
        perm.sum <- summary(get(perm))
        abline(h = perm.sum[,1], lty=1:2) # lod.p.mu
        abline(h = perm.sum[,2], lty=1:2, col = 'blue') # lod.p 
        abline(h = perm.sum[,3], lty=1:2, col = 'red') # lod.mu
      } else if (!has_xchr_perms(get(perm)) & names(mod)[3] == "lod1"){ # no X-chr-specific perms, GxT mod 
        ### update this when you have more specific thresholds 
        perm.sum <- summary(get(perm))
        abline(h = perm.sum, lty = 1:2)
      } else if (has_xchr_perms(get(perm)) & ncol(mod)==3){ # X-chr-specific perms, 1 LOD column
        ### Should handle X-chr specific thresholds, but not necessary for virus research paper bc no X-chr QTL
        # perm.sum = summary(get(perm))
        # abline(h = c(perm.sum$A[1], perm.sum$A[2]), lty = 1:2)
        permgevA <- paste0(perm, 'A.gev')
        perm.sumA <- summary(get(permgevA))
        abline(h = c(perm.sumA[1], perm.sumA[2]), lty=1:2)
        # permgevX <- paste0(perm, 'X.gev')
      } else if (has_xchr_perms(get(perm)) & ncol(mod)>3){ # X-chr-specific perms, >1 LOD column 
        perm.sum <- summary(get(perm))
        abline(h = perm.sum$A[,1], lty=1:2) # lod.p.mu
        abline(h = perm.sum$A[,2], lty=1:2, col = 'blue') # lod.p 
        abline(h = perm.sum$A[,3], lty=1:2, col = 'red') # lod.mu
      }
    }
    
    if (save == TRUE){ dev.off() } # close connection 
  }
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans: use different thresholds for X-chr 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans <- function(models, ylim = NULL){
#   for (i in 1:nrow(models)){
#     mod <- get(models[i, 'obj'])
#     
#     if (ncol(mod)>3) {lodcolumn=1:3} else {lodcolumn=1} # if ncol>3, multi-LOD model 
#     
#     if (is.null(ylim)){
#       plot(mod, ylab = "LOD", main = models[i,'name'], alternate.chrid = T, 
#            bandcol = "gray90", lodcolumn=lodcolumn)
#     } else {
#       plot(mod, ylab = "LOD", ylim = ylim, main = models[i,'name'], 
#            alternate.chrid = T, bandcol = "gray90", lodcolumn=lodcolumn)
#     }
#     
#     if (exists(models[i,'perm.obj'])){
#       perm.sum <- summary(get(models[i,'perm.obj']))
#       
#       if (length(perm.sum) == 10){ # separate thresholds for X chromosome 
#         ### CHANGE TO if (length(perm.sum == 2))
#         # Determine x-axis segments representing autosomes and X
#         temp <- mod
#         begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
#         rownames(begend) <- unique(mod[,1])
#         chr <- unique(as.character(mod[,1]))
#         begend <- begend[as.character(chr),,drop=FALSE]
#         len <- begend[,2]-begend[,1]
#         start <- c(0,cumsum(len))-c(begend[,1],0)
#         maxx <- sum(len)
#         
#         # end of chr1
#         x1=tail(mod[mod[,1] == chr[m],2], n=1)
#         
#         # beginning and end of chr2
#         
#         
#         # beginning and end of all chromosomes 
#         i <- c(1:20)
#         
#         segments(x0 = 0, x1 = 2600, y0 = perm.sum[[1]][1], lty=1)
#         segments(x0 = 0, x1 = 2600, y0 = perm.sum[[1]][2], lty=2)
#         segments(x0 = 2600, x1 = 3000, y0 = perm.sum[[2]][1], lty=1)
#         segments(x0 = 2600, x1 = 3000, y0 = perm.sum[[2]][2], lty=2)
#       } else { # only one set of thresholds 
#         ### REMOVE [[1]] from perm.sum calls 
#         abline(h = c(perm.sum[[1]][1], perm.sum[[1]][2]), lty=1:2)
#       }
#     }
#   }
# }


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots for significant peaks
# Input:
#     cross.obj: r/qtl cross object 
#     cross.type: bc or f2
#     raw.data: raw phenotypes (cross.obj$pheno before transformations)
#     geno.map: list mapping A and B allele to more informative labels
#     peaks: dataframe of significant peaks 
#     plot.type: type of plot to produce for each marker (can be "effect" or "pxg")
#     qtl_map: dataframe mapping chromosome to QTL name (for plot titles)
#     theme: ggplot theme to use for PxG
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_pxg <- function(cross.obj, cross.type, raw.data, peaks, plot.type, qtl_map, 
                     theme, type, save = FALSE, save.dir = NULL,...){
  for (k in 1:nrow(peaks)){
    
    if (save == TRUE){
      plotname <- paste0(peaks[k,'model'], '-chr', peaks[k,'chr'])
      fname <- paste0(save.dir, plotname, '.png')
      png(fname, width = 550)
    }
    
    if (plot.type == "effect"){
      effectplot(cross.obj, mname1 = peaks[k,'marker'], 
                 pheno.col = models[which(models$obj == peaks[k,'model']), 'colname'], 
                 ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
                 main = paste(peaks[k,'marker'], " (chr", peaks[k,'chr'], ")", sep=""))
    } else if (plot.type == "pxg"){
      pheno.col = as.character(models[which(models$obj == peaks[k,'model']), 'colname'])
      pheno = raw.data[,pheno.col]
      marker.name = peaks[k,'marker']
        
      p <- pxg(cross.obj, pheno = pheno, 
               marker = marker.name,  
               ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
               theme = theme, ...)
      print(p)
      
      # using built-in r/qtl function
      # plotPXG(cross.obj, marker = peaks[k,'marker'], 
      #         pheno.col = models[which(models$obj == peaks[k,'model']), 'colname'], 
      #         ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
      #         main = paste(peaks[k,'marker'], " (chr", peaks[k,'chr'], ")", sep=""))
    }
    
    if (save == TRUE){ dev.off() }
  }
}





#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots with ggplot. 
# Alternative to qtl::plotPXG, uses ggplot and more intuitive error bars  
# Input:
#     cross: r/qtl cross object 
#     pheno: vector of phenotype values for all mice 
#     marker: marker whose genotype to plot
#     geno.map: list mapping A and B allele to more informative labels 
#     qtl.map: list mapping chromosomes to QTL names for x-axis label
#     xlab: x axis label
#     ylab: y axis label 
#     title: title
#     theme: ggplot theme
#     type: scatter (default), boxplot or violin 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg <- function(cross, pheno, marker, geno.map, qtl.map = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                title = NULL, bestfit = TRUE, color = "black", theme = rmd_theme, type = 'scatter', ...){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if ((chr == 'X') & ('M' %in% levels(cross$pheno$sex))){ bestfit = FALSE } 
  
  if(is.null(xlab) & !is.null(qtl.map)){
    qtl.name <- qtl.map$qtl_name[qtl.map$chr == chr]
    xlab <- as.expression(bquote(italic(.(qtl.name))*" (chr"*.(chr)*") "*Genotype))
  } else if (is.null(xlab) & is.null(qtl.map)){
    xlab <- "Genotype"
  }
  
  cross.imp <- fill.geno(cross) 
  marker.genos <- cross.imp$geno[[chr]]$data[,marker]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      X.data <- data.frame(sex = sex, geno = df$geno, pgm = pgm)
      X.data$X.geno <- with(df, 
                            ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                            ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                            ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                            ifelse(sex==1 & geno==1, 'AY', 'BY'))))))
      
      df$geno <- mapvalues(X.data$X.geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  print(dfsum)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, colour = color)) + 
      geom_jitter(width = 0.1, height = 0) + 
      scale_y_continuous(limits = ylim) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean, ymax = mean), 
                    width = 0.2, size = 1, alpha = 0.5) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean-sd, ymax = mean+sd), 
                    width = 0.3, size = 1, alpha = 0.5) +
      labs(title = title, y = ylab, x = xlab) + 
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, colour = color)} + 
      theme
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
      geom_boxplot(fill = color, notch = TRUE) + 
      scale_y_continuous(limits = ylim) + 
      labs(title = title, y = ylab, x = xlab) + 
      theme
    # p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
    #   geom_boxplot(fill = "#1B9E77", notch = TRUE) + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
    #   labs(title = title, y = ylab, x = xlab) +
    #   theme
  } else if (type == 'violin'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_violin(fill = "#1B9E77") + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((type == 'scatter') & (bestfit == TRUE)){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# make_qtl_scan 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_qtl_scan <- function(flow_pheno, palette, incl_legend = TRUE, legend.x = 2300, 
                          legend.y = NULL, ylim = NULL, cex = 1){
  if (is.list(palette)) palette <- unlist(palette, use.names = FALSE)
  
  mod0 <- readRDS(paste0("results/qtl_mapping/modRDS/PBS/", flow_pheno, ".rds"))
  mod1 <- readRDS(paste0("results/qtl_mapping/modRDS/SARS/", flow_pheno, ".rds"))
  mod2 <- readRDS(paste0("results/qtl_mapping/modRDS/SARS2/", flow_pheno, ".rds"))
  mod0$p_SARS <- mod1$lod
  mod0$p_SARS2 <- mod2$lod
  
  perm0 <- readRDS(paste0("results/qtl_mapping/permRDS/PBS/", flow_pheno, ".rds"))
  perm1 <- readRDS(paste0("results/qtl_mapping/permRDS/SARS/", flow_pheno, ".rds"))
  perm2 <- readRDS(paste0("results/qtl_mapping/permRDS/SARS2/", flow_pheno, ".rds"))
  
  scan_plot <- ggplotify::as.grob(function(){
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(7,5,2,1) + 0.1, xpd = NA) # increase bottom margin 
    
    if (is.null(ylim)) {
      plot_modX(mod0, lodcolumn = 1:3, col = palette, lwd = 1.5,
                cex.axis = cex, cex.x = cex, cex.lab = cex, 
                alternate.chrid = TRUE, ylab = "", xlab = "")
    } else {
      plot_modX(mod0, lodcolumn = 1:3, col = palette, lwd = 1.5,
                cex.axis = cex, cex.x = cex, cex.lab = cex, ylim = ylim,
                alternate.chrid = TRUE, ylab = "", xlab = "")
    }
  
    mtext(expression(-log[10](P)), side = 2, line = 1.8, cex = cex)
    mtext("Chromosome", side = 1, line = 2.2, cex = cex)
    abline(h = summary(perm0)[1], lty = 2, col = palette[1])
    abline(h = summary(perm1)[1], lty = 2, col = palette[2])
    abline(h = summary(perm2)[1], lty = 2, col = palette[3])
    
    if (incl_legend){
      if (is.null(legend.y)) {
        legend.y <- if (!is.null(ylim)) ylim[2] - 2 else max(unlist(mod0[, 3:5])) - 2
      }
      usr <- par("usr")
      legend(x = legend.x, y = legend.y, legend = c("Control", "SARS-CoV", "SARS-CoV-2"), 
             col = palette, lwd = 2, lty = 1, bty = "n", xjust = 0, yjust = 0, cex = cex, y.intersp = 1.8)
    }
  })
  
  return(scan_plot)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# make_gxt_qtl_scan 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_gxt_qtl_scan <- function(flow_pheno, palette, legend.x = 2300, legend.y.diff = NULL, ylim = NULL, cex = 0.8){
  if (is.list(palette)) palette <- unlist(palette, use.names = FALSE)
  if (is.null(legend.y.diff)) legend.y.diff <- 1.5
  
  mod <- readRDS(paste0("results/qtl_mapping/modRDS/GxT/", flow_pheno, ".rds"))
  perm <- readRDS(paste0("results/qtl_mapping/permRDS/GxT/", flow_pheno, ".rds"))
  
  scanplot <- ggplotify::as.grob(function(){
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mar = c(7,5,2,1) + 0.1, xpd = NA) # increase bottom margin
    
    plot_modX(mod, lodcolumn = 1:3, col = palette, lwd = 1.5,
              cex.axis = cex, cex.x = cex, cex.lab = cex, alternate.chrid = TRUE,
              ylab = "", xlab = "", ylim = ylim) 
    
    mtext(expression(-log[10](P)), side = 2, line = 1.8, cex = cex)
    mtext("Chromosome", side = 1, line = 2.2, cex = cex)
    
    abline(h = summary.permgev(perm)[1], lty = 2)
    
    usr <- par("usr")
    legend.y <- if (!is.null(ylim)){ 
      ylim[2] - legend.y.diff
    } else {
      max(unlist(mod0[,3:5])) - legend.y.diff
    }  # par("usr")[4] - 0.5
    legend(x = legend.x, y = legend.y, legend = c("Overall", "Genetic", "GxT"), 
           col = palette, lwd = 2, lty = 1, bty = "n", xjust = 0, yjust = 0, cex = cex, y.intersp = 1.8)
  })
  
  return(scanplot)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# make_qtl_pxg_plot 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_qtl_pxg_plot <- function(flow_pheno, marker, chr, pos, ylab, labels, cex = cex){
  if (is.null(marker)){
    marker <- find.marker(cross, chr = chr, pos = pos)
  }
  pxg <- pxg_group(cross, pheno = cross_data[[flow_pheno]], marker = marker,
                   geno.map = list(A = "CC006", B = "CC044"), xlab = "Genotype", ylab = ylab, 
                   groupcolumn = "infection", grouplabels = grouplabels,
                   groupcolors = inf_pal, bestfit = TRUE) + theme_minimal()
  
  mod0 <- readRDS(paste0("results/qtl_mapping/modRDS/PBS/", flow_pheno, ".rds"))
  mod1 <- readRDS(paste0("results/qtl_mapping/modRDS/SARS/", flow_pheno, ".rds"))
  mod2 <- readRDS(paste0("results/qtl_mapping/modRDS/SARS2/", flow_pheno, ".rds"))
  mod0$p_SARS <- mod1$lod
  mod0$p_SARS2 <- mod2$lod
  
  perm0 <- readRDS(paste0("results/qtl_mapping/permRDS/PBS/", flow_pheno, ".rds"))
  perm1 <- readRDS(paste0("results/qtl_mapping/permRDS/SARS/", flow_pheno, ".rds"))
  perm2 <- readRDS(paste0("results/qtl_mapping/permRDS/SARS2/", flow_pheno, ".rds"))
  
  scan_plot <- ggplotify::as.grob(function(){
    plot_modX(mod0, lodcolumn = 1:3, col = brewer.pal(n=3,"Set2"), 
              cex.axis = cex, cex.x = cex, cex.lab = cex,
              alternate.chrid = TRUE, ylab = "", xlab = "")
    title(ylab = expression(-log[10](P)), line = 1.8, cex.lab = cex)
    title(xlab = "Chromosome", line = 2.2, cex.lab = cex)
    abline(h = summary(perm0)[1], lty = 2, col = brewer.pal(n=3,"Set2")[1])
    abline(h = summary(perm1)[1], lty = 2, col = brewer.pal(n=3,"Set2")[2])
    abline(h = summary(perm2)[1], lty = 2, col = brewer.pal(n=3,"Set2")[3])
  })
  
  pxg_wspacer <- plot_grid(NULL, pxg, nrow = 2, rel_heights = c(0.1,1)) 
  pxg_scan <- plot_grid(scan_plot, pxg_wspacer, rel_widths = c(1.2,1), labels = labels)
  return(pxg_scan)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# make_gxt_qtl_pxg_plot 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_gxt_qtl_pxg_plot <- function(flow_pheno, chr, pos, ylab, labels){
  marker <- find.marker(cross, chr = chr, pos = pos)
  pxg <- pxg_group(cross, pheno = cross_data[[flow_pheno]], marker = marker,
                   geno.map = list(A = "CC006", B = "CC044"), xlab = "Genotype", ylab = ylab, 
                   groupcolumn = "infection", grouplabels = grouplabels,
                   groupcolors = inf_pal, bestfit = TRUE) + theme_minimal()
  
  mod <- readRDS(paste0("results/qtl_mapping/modRDS/GxT/", flow_pheno, ".rds"))
  perm <- readRDS(paste0("results/qtl_mapping/permRDS/GxT/", flow_pheno, ".rds"))
  
  scanplot <- ggplotify::as.grob(function(){
    plot_modX(mod, lodcolumn = 1:3, col = c("black", rev(brewer.pal(n = 3, "Set1")[1:2])), 
              cex.axis = 0.8, cex.x = 0.8, cex.lab = 0.8, alternate.chrid = TRUE,
              ylab = "", xlab = "") 
    title(ylab = expression(-log[10](P)), line = 1.8, cex.lab = 0.8)
    title(xlab = "Chromosome", line = 2.2, cex.lab = 0.8)
    abline(h = summary.permgev(perm)[1], lty = 2)
    legend_x_pos <- 2000
    legend_y_pos <- par("usr")[4] - 0.5
    legend(legend_x_pos, legend_y_pos, legend = c("Overall", "Genetic", "GxT"), lwd = 2, cex = 1,
           col = c("black", rev(brewer.pal(n = 3, "Set1")[1:2])), bty = "y", bg = "white")
  })
  pxg_wspacer <- plot_grid(NULL, pxg, nrow = 2, rel_heights = c(0.1,1)) 
  pxg_scan <- plot_grid(scanplot, pxg_wspacer, rel_widths = c(1.2,1), labels = labels)
  return(pxg_scan)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for subsetting 
# Removed code that colored imputed genotypes differently 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_subset <- function(cross, pheno, marker, cross.type, geno.map, xlab = 'Genotype', ylab, 
                title = NULL, bestfit = TRUE, groupcol, groupselect, theme = rmd_theme){
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if (chr == 'X'){bestfit = FALSE} 
  
  # plot title 
  if (is.null(title)){
    title = paste('Genotype vs. phenotype at chr', chr, 'marker', marker)
  }
  
  cross <- fill.geno(cross) # calculate imputed genotypes, just for plotting
  marker.genos <- cross$geno[[chr]]$data[,marker]
  marker.genos <- marker.genos[pull.pheno(cross, groupcol)==groupselect]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos, pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
      }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$imp_geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$imp_geno), FUN = sd, na.rm = TRUE)$x)
  
  p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
      geom_jitter(width = 0.1, height = 0) + 
      geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) + 
      geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) +
      labs(title = title, y = ylab, x = xlab) + 
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
      theme
  
  if (bestfit == TRUE){
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for multiple groups (e.g. infection, control) to pxg()
# Trying to make it so it works with and without groups, so there's only one function
# Maybe only have one pxg function, and that function calls pxg_group if there's a group 
# and pxg_ind if not? 
### 1. Works for markers with no missing genotypes, need to validate with markers that 
###    have missing genotypes. 
### 2. Need to validate with no grouping (groupcolumn = NULL); I know that df doesn't 
###    group column when groupcolumn=NULL
### 3. Need to validate with bestfit = FALSE (validate on X chr)
###    groupcolumn = 'infection'/NULL with bestfit= TRUE/FALSE 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_group <- function(cross, pheno, marker, geno.map, 
                      xlab = 'Genotype', ylab, plot.title = NULL, 
                      groupcolumn = NULL, groupcolors = NULL, grouplabels = NULL,
                      bestfit = TRUE, theme = theme_minimal(), 
                      type = c("scatter", "boxplot", "violin"),
                      # box/violin overlay 
                      show.points = FALSE, point.alpha = 0.25, point.size = 0.8,
                      jitter.width = 0.12, show.n = TRUE, n.size = 3, n.vjust = -0.4, 
                      violin.scale = c("area", "count", "width"),
                      # X chr handling 
                      facet.x.sex = TRUE, sex.pheno = "sex", X_combine_AB = FALSE){
  
  type <- match.arg(type)
  violin.scale <- match.arg(violin.scale)
  
  cross.type <- class(cross)[1] # cross type 
  chr <- which_chr(cross, marker) # id chromosome marker is on 
  if (chr == 'X'){bestfit = FALSE} # best fit line only appropriate for autosomes 
  cross.imp <- fill.geno(cross) # impute genotypes 
  
  # handle grouping 
  if (is.null(groupcolumn)){
    group <- as.factor(rep("All", nind(cross)))
  } else {
    group <- factor(cross$pheno[,groupcolumn]) # category for grouping 
  }
  
  sex <- NULL
  if (chr == "X" && facet.x.sex) {
    sex <- factor(pull.pheno(cross, sex.pheno))
  }
  
  df <- data.frame(geno = cross$geno[[chr]]$data[,marker], 
                   imp_geno = cross.imp$geno[[chr]]$data[,marker],
                   pheno = pheno, 
                   group = group)
  if (!is.null(sex)) { df$sex <- sex }
  
  # ---------------------genotype labels from A/B mapping--------------------- # 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotype codes to labels 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
    df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2), to = geno.labels) 
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sexnum <- as.numeric(pull.pheno(cross, sex.pheno) == "M")
      geno_imp <- df$imp_geno
      X.geno <- ifelse(sexnum == 0 & geno_imp == 1 & pgm == 0, "AA",
                       ifelse(sexnum == 0 & geno_imp == 1 & pgm == 1, "BB", 
                              ifelse(sexnum == 0 & geno_imp == 2 & pgm == 0, "ABf",
                                     ifelse(sexnum == 0 & geno_imp == 2 & pgm == 1, "ABr", 
                                            ifelse(sexnum == 1 & geno_imp == 1, "AY", "BY")))))
      
      # X.geno with NA for missing genotypes 
      temp <- X.geno
      temp[which(is.na(df$geno))] <- NA 
      df$geno <- temp
      # full X.geno
      df$imp_geno <- X.geno
      
      # remap to genotype labels 
      df$geno <- mapvalues(df$geno, from = c("AA", "ABf", "ABr", "BB", "AY", "BY"), to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c("AA", "ABf", "ABr", "BB", "AY", "BY"), to = geno.labels) 
    }
  }
  df$geno <- factor(df$geno, levels = geno.labels)
  df$imp_geno <- factor(df$imp_geno, levels = geno.labels)
  df$x_num <- as.numeric(df$imp_geno)
  
  # if X chr, optionally combine ABf and ABr into a single AB level
  if (cross.type == "f2" && chr == "X" && X_combine_AB) {
    ABf_lab <- geno.labels[2]
    ABr_lab <- geno.labels[3]
    AB_lab  <- paste(geno.map$A, geno.map$B, sep = "/")
    
    df$imp_geno <- as.character(df$imp_geno)
    df$imp_geno[df$imp_geno %in% c(ABf_lab, ABr_lab)] <- AB_lab
    df$imp_geno <- factor(df$imp_geno,
                          levels = c(geno.labels[1], AB_lab, geno.labels[4], geno.labels[5], geno.labels[6]))
  }

  # palette, order, labels 
  breaks <- NULL
  labels <- NULL
  if (!is.null(groupcolors)){
    groupcolors <- unlist(groupcolors, use.names = TRUE)
    breaks <- names(groupcolors)
    df$group <- factor(df$group, levels = breaks)
    labels <- if (!is.null(grouplabels)) unname(grouplabels[breaks]) else breaks
  }
  
  # counts for n-labels (respect sex facet if present)
  if (show.n) {
    if (!is.null(sex)) {
      ndf <- aggregate(df$pheno, by = list(df$imp_geno, df$group, df$sex), FUN = function(x) sum(!is.na(x)))
      colnames(ndf) <- c("imp_geno", "group", "sex", "n")
    } else {
      ndf <- aggregate(df$pheno, by = list(df$imp_geno, df$group), FUN = function(x) sum(!is.na(x)))
      colnames(ndf) <- c("imp_geno", "group", "n")
    }
  }
  
  # --------------------------------build plot-------------------------------- #
  if (type == "scatter"){ # change dodge.width = 0.3 and position_dodge(width=0.3) when 2 groups
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno, color = group)) + 
      geom_point(shape = 20, position = position_jitterdodge(dodge.width = 0.7)) +
      labs(title = plot.title, y = ylab, x = xlab, color = tools::toTitleCase(groupcolumn)) + 
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", position = position_dodge(width = 0.7), width = 0.5) +
      stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.7)) +
      {if (!is.null(groupcolors)) scale_color_manual(values = groupcolors, breaks = breaks, labels = labels, name = groupcolumn)} + 
      {if (bestfit == TRUE) geom_smooth(data = df, aes(x = x_num, y = pheno, color = group),
                                        method = "lm", se = FALSE, inherit.aes = FALSE, linewidth = 0.5)} + 
      scale_x_discrete(limits = levels(df$imp_geno)) + 
      theme
  } else if (type %in% c("boxplot", "violin")){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno, fill = group)) +
      labs(title = plot.title, y = ylab, x = xlab, fill = tools::toTitleCase(groupcolumn)) +
      {if (!is.null(groupcolors)) scale_fill_manual(values = groupcolors, breaks = breaks, labels = labels, name = tools::toTitleCase(groupcolumn))} +
      theme
    
    if (type == "violin"){
      p <- p + 
        geom_violin(scale = violin.scale, trim = TRUE)
    } else {
      p <- p + 
        geom_boxplot(outliers = FALSE)
    }
    
    # optional faint points 
    if (show.points) {
      p <- p + 
        geom_point(aes(color = group), position = position_jitterdodge(jitter.width = jitter.width, dodge.width = 0.7),
                   alpha = point.alpha, size = point.size, inherit.aes = TRUE, show.legend = FALSE) +
        {if (!is.null(groupcolors)) scale_color_manual(values = groupcolors, breaks = breaks, labels = labels, name = tools::toTitleCase(groupcolumn))}
    }
    
    # optional n labels
    if (show.n) {
      p <- p + geom_text(data = ndf, aes(x = imp_geno, y = Inf, label = paste0("n=", n),group = group),
                         position = position_dodge(width = 0.7), vjust = n.vjust, size = n.size, inherit.aes = FALSE)
    }
  }
  
  # optional facet by sex 
  if (chr == "X" && facet.x.sex && !is.null(sex)) {
    facet_labs <- 
    p <- p + 
      ggh4x::facet_grid2(~sex, scales = "free_x", space = "free_x", labeller = as_labeller(c(F = "Female", M = "Male"))) + 
      theme(strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
            strip.text = element_text(face = "bold"))
  }
  
  return(p)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for coloring dots by 2nd marker genotype
# Remove color-by-imputed-or-not because we're coloring by the genotype at a 2nd marker 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_colbym <- function(cross, pheno, marker, geno.map, xlab = 'Genotype', ylab, 
                title = NULL, bestfit = FALSE, theme = rmd_theme, marker2,
                type = 'scatter'){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # what chromosome is marker2 on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker2)
  chr2 <- names(cross$geno)[o]
  
  cross <- fill.geno(cross)
  marker.genos <- cross$geno[[chr]]$data[,marker]
  marker2.genos <- cross$geno[[chr2]]$data[,marker2]
  
  df <- data.frame(geno = marker.genos,
                   geno2 = marker2.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
    df$geno2 <- mapvalues(df$geno2, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    
    ### NEED TO VALIDATE THIS
    
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$geno2 <- factor(df$geno2, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, col = geno2)) + 
      geom_jitter(width = 0.1, height = 0) + 
      labs(title = title, y = ylab, x = xlab) + 
      scale_color_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) +
      theme
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, fill = geno2)) + 
      geom_boxplot() +
      scale_fill_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) + 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((bestfit == TRUE) & (type == 'scatter')){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}






#-------------------------------Statistical model------------------------------#
# Calculates number of markers in a cross object 
calc.num.markers <- function(cross){
  p = 0
  for (c in 1:length(cross$geno)){
    p = p + ncol(cross$geno[[c]]$data)
  }
  return(p)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
### This function does random forest imputation for a phenotype?
#+++++++++++++++++++++++++++++++++++++++++++++++++++
imputeRF <- function(cross_data, pheno_name, pct_train, n_flow_start, plot_dir, 
                     surr_col_name){ # removed mod_perf, mod_perf_i, 
  #--------------------------------prep data-----------------------------------#
  null_data <- cross_data[is.na(cross_data[[pheno_name]]),]
  non_null_data <- cross_data[!is.na(cross_data[[pheno_name]]),]
  
  # select training data, with varying n per infection group 
  train_n <- round(table(non_null_data$infection) * pct_train)
  train_df <- non_null_data %>% 
    group_by(infection) %>% 
    arrange(infection) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(n = train_n) %>% 
    mutate(samp = map2(data, n, sample_n)) %>%
    dplyr::select(-c("data", "n")) %>% 
    unnest(samp)
  
  test_df <- non_null_data[!(non_null_data$mouse_ID %in% train_df$mouse_ID),]
  
  #-------------------------------train model----------------------------------#
  # define predictor variables / formula 
  cols_to_incl <- colnames(train_df)[n_flow_start:ncol(cross_data)] # all flow phenotypes and genotypes 
  cols_to_incl <- cols_to_incl[which(cols_to_incl != pheno_name)] # remove target phenotype 
  form <- paste(pheno_name, "~ infection + sex + batch + pgm + F1_dam_geno + d0 + d1 + d2 + d3 + d4 + HS + Titer +", 
                paste0(cols_to_incl, collapse = " + "))
  
  ranger_mod <- ranger(form, data = train_df) # train random forest model 
  pred_test <- predict(ranger_mod, test_df) # predict test data
  
  #---------------------------evaluate performance-----------------------------#
  # plot observed vs predicted 
  plot_fname <- paste0(plot_dir, pheno_name, ".png")
  png(plot_fname)
  plot(test_df[[pheno_name]], pred_test$predictions)
  abline(a=0, b=1)
  dev.off()
  
  # calculate RMSE
  rmse <- sqrt(sum((pred_test$predictions - test_df[[pheno_name]])^2) / nrow(test_df))
  # calculate R-squared 
  rss <- sum((pred_test$predictions - test_df[[pheno_name]])^2)
  tss <- sum((test_df[[pheno_name]] - mean(test_df[[pheno_name]]))^2)
  rsq <- 1 - (rss/tss)
  #mod_perf[mod_perf_i,] <<- c(pheno_name, rmse, Rsq)
  
  # combine test_df with null_data and predict all values based on random forest  
  pred_data <- rbind(null_data, test_df)
  pred <- predict(ranger_mod, pred_data)
  pred_data[[surr_col_name]] <- pred$predictions
  data_to_join <- pred_data[,c("mouse_ID", surr_col_name)]
  
  # data for qtl mapping 
  joined_data <- left_join(cross_data, data_to_join, by = "mouse_ID")
  attr(joined_data, "rmse") <- rmse 
  attr(joined_data, "rsq") <- rsq
  return(joined_data)
}


### What is the difference between these two fxns?
### COMMENTED OUT ALL, USING QTL_FUNCTIONS_SS.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++
### This function does univariate modeling on non-missing target, for comparison
### with SynSurr results 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# univar_scan <- function(cross, target, covar = NULL, model = "normal", gxt = FALSE, weights = NULL){
#   # Calculate genotype probabilities if necessary 
#   if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
#     cross <- calc.genoprob(cross)
#     print("Running calc.genoprob...")
#   }
#   cross <- calc.dosage(cross) # calculate dosage (of B allele)
#   p <- calc.num.markers(cross) # number of markers
#   n <- nrow(cross$pheno) # number of mice / observations
#   
#   # Pre-allocate vectors for results 
#   markers <- vector(length = p)
#   chrs <- vector(length = p)
#   positions <- vector(length = p)
#   lods <- vector(length = p)
#   
#   covar$geno <- rep(NA, nrow(covar)) # covariates data frame for FitBNR
#   marker.pos = 0
#   for (c in 1:length(cross$geno)){
#     for (k in 1:ncol(cross$geno[[c]]$data)){
#       marker.pos = marker.pos + 1
#       covar$geno <- cross$geno[[c]]$dos[,k] # add genotype dosage to covariate dataframe 
#       
#       if (model == "normal"){
#         fitlm <- lm(target ~ sexM + geno, data = covar)
#         pval <- summary(fitlm)$coefficients['geno', 'Pr(>|t|)']
#       } else if (model == "np"){ # doesn't support covariates 
#         pval <- kruskal.test(target ~ round(covar$geno))$p.value
#       } else if (model == "wls"){
#         covar_lm <- covar[!is.na(target),]
#         pheno_lm <- target[!is.na(target)]
#         weights_lm <- weights[!is.na(target)] 
#         
#         if (gxt == TRUE){
#           wls_lmGxT <- lm(pheno_lm ~ covar_lm$sexM + covar_lm$infectionSARSCoV + 
#                             covar_lm$infectionSARS2CoV + round(covar_lm$geno) + 
#                             covar_lm$infectionSARSCoV*round(covar_lm$geno) + 
#                             covar_lm$infectionSARS2CoV*round(covar_lm$geno),
#                           weights = weights_lm)
#           wls_lm0 <- lm(pheno_lm ~ covar_lm$sexM + covar_lm$infectionSARSCoV + 
#                           covar_lm$infectionSARS2CoV,
#                         weights = weights_lm)
#           aov <- anova(wls_lm0, wls_lmGxT)
#           pval <- aov$`Pr(>F)`[2]
#         } else { # no GxT 
#           wls_lm <- lm(pheno_lm ~ covar_lm$sexM + covar_lm$infectionSARSCoV + 
#                          covar_lm$infectionSARS2CoV + round(covar_lm$geno),
#                        weights = weights_lm)
#           pval <- summary(wls_lm)$coefficients[5,4]
#         } 
#       } 
#       
#       marker.name <- colnames(cross$geno[[c]]$data)[k]
#       markers[marker.pos] <- marker.name
#       chrs[marker.pos] <- c
#       positions[marker.pos] = cross$geno[[c]]$map[marker.name]
#       lods[marker.pos] = -log10(pval)
#     }
#   }
#   
#   # Create object 
#   model.df <- data.frame(chrs, positions, lods)
#   names(model.df) <- c("chr", "pos", "lod") # same column names as scanone object
#   rownames(model.df) <- markers
#   class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
#   
#   return(model.df)
# }


#+++++++++++++++++++++++++++++++++++++++++++++++++++
### This function does univariate modeling on non-missing target, for comparison
### with SynSurr results 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# univar_scan_wmv <- function(cross, target_mat, wls = TRUE){
#   # Calculate genotype probabilities if necessary 
#   if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
#     cross <- calc.genoprob(cross)
#     print("Running calc.genoprob...")
#   }
#   cross <- calc.dosage(cross) # calculate dosage (of B allele)
#   p <- calc.num.markers(cross) # number of markers
#   n <- nrow(cross$pheno) # number of mice / observations
#   
#   # Pre-allocate vectors for results 
#   markers <- vector(length = p)
#   chrs <- vector(length = p)
#   positions <- vector(length = p)
#   lods <- vector(length = p)
#   
#   ind_w_pheno <- !is.na(target_mat[,1]) # fully observed for first column of target mat
#   target_mat <- as.matrix(target_mat[ind_w_pheno,])
#   
#   df <- cross$pheno[ind_w_pheno, c("L_NKCells_pctCD25", "infection", "sex")]
#   df$group <- interaction(df$sex, df$infection, drop = TRUE) # define 6 groups (sex*infection)
#   ols <- lm(L_NKCells_pctCD25 ~ sex + infection, data = df) # get OLS estimates 
#   res <- resid(ols) # get residuals
#   vhat <- tapply(res^2, df$group[!is.na(df$L_NKCells_pctCD25)], mean) # estimate group variances 
#   weights <- 1 / vhat[df$group]
#   
#   diagW <- diag(sqrt(weights)) # get diagonal matrix containing the square roots of the weights
#   Ywt <- diagW %*% target_mat # pre-multiple outcome matrix by weights matrix 
#   
#   marker.pos = 0
#   for (c in 1:length(cross$geno)){
#     for (k in 1:ncol(cross$geno[[c]]$data)){
#       marker.pos = marker.pos + 1
#       df$geno <- cross$geno[[c]]$dos[ind_w_pheno,k] # add genotype dosage to covariate dataframe 
#       
#       fullX <- model.matrix(~ sex + infection + geno + sex*geno + infection*geno, data = df)
#       Xwt <- diagW %*% fullX # pre-multiply model matrix by weights matrix 
#       mlm1 <- lm(target_mat ~ .-1, data = as.data.frame(Xwt))
#       mlm0 <- lm(target_mat ~ `(Intercept)` + sexM + infectionSARSCoV + infectionSARS2CoV, data = as.data.frame(Xwt))
#       pval <- anova(mlm1, mlm0)$`Pr(>F)`[2]
#       
#       marker.name <- colnames(cross$geno[[c]]$data)[k]
#       markers[marker.pos] <- marker.name
#       chrs[marker.pos] <- c
#       positions[marker.pos] = cross$geno[[c]]$map[marker.name]
#       lods[marker.pos] = -log10(pval)
#     }
#   }
#   
#   # Create object 
#   model.df <- data.frame(chrs, positions, lods)
#   names(model.df) <- c("chr", "pos", "lod") # same column names as scanone object
#   rownames(model.df) <- markers
#   class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
#   
#   return(model.df)
# }


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Univariate permutation tests 
# Note that column names in covar dataframe are hard-coded
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# univar_perm <- function(cross, target, covar = NULL, num_perms = 1000, gxt = FALSE, weights = rep(1, length(target))){
#   # calculate dosage for hk regression 
#   if (!with(cross$geno[[1]], exists('dos'))){  # calc.dosage hasn't been run
#     print("Running calc.dosage...")
#     cross <- calc.dosage(cross)
#   }
#   
#   target_og <- target # data to permute 
#   ind_w_pheno <- !is.na(target_og)
#   target <- target_og[ind_w_pheno] # remove NAs; use ind_w_pheno to filter genotypes later 
#   
#   # create 1000 permutations of phenotype 
#   target_perms <- matrix(nrow = length(target), ncol = num_perms)
#   for (c in 1:num_perms){
#     target_perms[,c] <- target[sample(length(target))]
#   }
#   
#   lmdata <- as.data.frame(covar[ind_w_pheno,]) # data frame for lm()
#   lmdata$target <- target
#   lmdata$geno <- rep(NA, nrow(lmdata)) # to hold genotypes
#   
#   if (gxt == TRUE){
#     ### only need dummy gxt variable for SynSurr::FitBNR
#     #lmdata$gxt <- rep(NA, nrow(lmdata)) # to hold gxt interaction term 
#     #formula0 <- as.formula("target ~ sex + infection")
#     formula0 <- as.formula(target ~ sexM + infectionSARSCoV + infectionSARS2CoV)
#     #formula1 <- as.formula("target ~ sex + infection + geno + gxt")
#     formula1 <- as.formula(target ~ sexM + infectionSARSCoV + infectionSARS2CoV + 
#                              geno + infectionSARSCoV*geno + infectionSARS2CoV*geno)
#   } else {
#     formula0 <- as.formula("target ~ sex")
#     formula1 <- as.formula("target ~ sex + geno")
#   }
#   
#   # create empty matrix for storing results 
#   perm_logPs <- matrix(nrow = num_perms, ncol = calc.num.markers(cross))
#   marker.pos = 0
#   
#   for (c in 1:length(cross$geno)){ # for each chr
#     for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
#       lmdata$geno <- NA # clear geno col
#       marker.pos = marker.pos + 1
#       geno <- cross$geno[[c]]$dos[,k]
#       geno <- geno[ind_w_pheno] # remove genotypes for mice that don't have the phenotype 
#       lmdata$geno <- geno # update geno col with marker genotype data
#       ### only need dummy gxt col for SynSurr::FitBNR
#       #if (gxt == TRUE){ lmdata$gxt <- lmdata$infection * lmdata$geno } # update gxt col 
#       
#       #--------------------------Calculate -log(P)s----------------------------#
#       res <- lm.multiresponse(formula = formula1, 
#                               response.matrix = target_perms, 
#                               data = lmdata, 
#                               null.formula = formula0, 
#                               logP = TRUE,
#                               weights = weights)
#       
#       #----------------------------Save -log(P)s-------------------------------#
#       perm_logPs[,marker.pos] <- res$logP
#     }
#   }
#   
#   # find highest LOD from each genome scan and create distribution 
#   max_logPs <- apply(perm_logPs, 1, max)
#   # return - vector of 1000 LOD scores with some attributes 
#   attributes(max_logPs) <- list(method = 'hk', model = 'normal',
#                                 type = 'f2', class = c('scanoneperm', 'matrix'))
#   return(max_logPs) 
# }



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Univariate permutation tests NOT using lm.multiresponse 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# univar_perm2 <- function(cross, target, covar = NULL, num_perms = 1000, gxt = FALSE, wls = FALSE){
#   # calculate dosage for hk regression 
#   if (!with(cross$geno[[1]], exists('dos'))){  # calc.dosage hasn't been run
#     print("Running calc.dosage...")
#     cross <- calc.dosage(cross)
#   }
# 
#   ind_w_pheno <- !is.na(target)
#   lmdata_og <- as.data.frame(covar)
#   lmdata_og$target <- target
#   lmdata_og <- lmdata_og[ind_w_pheno,]
#   if (wls == TRUE){
#     df <- cross$pheno[ind_w_pheno, c("L_NKCells_pctCD25", "infection", "sex")]
#     df$group <- interaction(df$sex, df$infection, drop = TRUE) # define 6 groups (sex*infection)
#     ols <- lm(L_NKCells_pctCD25 ~ sex + infection, data = df) # get OLS estimates 
#     res <- resid(ols) # get residuals
#     vhat <- tapply(res^2, df$group[!is.na(df$L_NKCells_pctCD25)], mean) # estimate group variances 
#     lmdata_og$weights <- 1 / vhat[df$group]
#   } else {
#     lmdata_og$weights <- rep(1, nrow(lmdata_og))
#   }
#   
#   if (gxt == TRUE){
#     ### only need dummy gxt variable for SynSurr::FitBNR
#     #lmdata$gxt <- rep(NA, nrow(lmdata)) # to hold gxt interaction term 
#     #formula0 <- as.formula("target ~ sex + infection")
#     formula0 <- as.formula(target ~ sexM + infectionSARSCoV + infectionSARS2CoV)
#     #formula1 <- as.formula("target ~ sex + infection + geno + gxt")
#     formula1 <- as.formula(target ~ sexM + infectionSARSCoV + infectionSARS2CoV + 
#                              geno + infectionSARSCoV*geno + infectionSARS2CoV*geno)
#   } else {
#     formula0 <- as.formula(target ~ sex)
#     formula1 <- as.formula(target ~ sex + geno)
#   }
#   
#   # create empty matrix for storing results 
#   perm_logPs <- matrix(nrow = num_perms, ncol = calc.num.markers(cross))
#   marker.pos = 0
#   
#   for (i in 1:num_perms){
#     resamp <- sample(nrow(lmdata_og), replace = FALSE)
#     lmdata <- lmdata_og[resamp,]
#     wt <- lmdata$weights
#     lmdata$geno <- rep(NA, nrow(lmdata)) # to hold genotypes
#     marker.pos <- 0
#     
#     for (c in 1:length(cross$geno)){ # for each chr
#       for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
#         lmdata$geno <- NA # clear geno col
#         marker.pos = marker.pos + 1
#         geno <- cross$geno[[c]]$dos[ind_w_pheno,k] # remove genotypes for mice that don't have the phenotype
#         lmdata$geno <- geno # update geno col with marker genotype data
#         ### only need dummy gxt col for SynSurr::FitBNR
#         #if (gxt == TRUE){ lmdata$gxt <- lmdata$infection * lmdata$geno } # update gxt col 
#       
#         #--------------------------Calculate -log(P)s----------------------------#
#         fit1 <- lm(formula1, data = lmdata, weights = wt)
#         fit0 <- lm(formula0, data = lmdata, weights = wt)
#         aov <- anova(fit1, fit0)
#         pval <- aov$`Pr(>F)`[2]
#       
#         #----------------------------Save -log(P)s-------------------------------#
#         perm_logPs[i, marker.pos] <- -log10(pval)
#       }
#     }
#   }
#   
#   # find highest LOD from each genome scan and create distribution 
#   max_logPs <- apply(perm_logPs, 1, max)
#   # return - vector of 1000 LOD scores with some attributes 
#   attributes(max_logPs) <- list(method = 'hk', model = 'normal',
#                                 type = 'f2', class = c('scanoneperm', 'matrix'))
#   return(max_logPs) 
# }



#+++++++++++++++++++++++++++++++++++++++++++++++++++
### This function does mixed modeling (using Haley-Knott regression)
### To-do:
###   1. Figure out how to dynamically create random effect terms based on length(randoms).
###   For now, assume there are two (cage and batch) - line 44
###   2. Handle covariates in addition to random effects (also dynamically)
###   3. X-chr approach validated for BC but needs to be validated for F2 (i.e. compare 
###   with results from R/qtl) - line 24
###   4. Permutations 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
hk <- function(cross, pheno.col, covar = NULL){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  
  cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n = nrow(cross$pheno) # number of mice / observations
  
  # Pre-allocate vectors for results 
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods <- vector(length = p)
  
  df <- covar # data frame for lm()
  df$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  
  # Create null model: just covariates 
  # formula0 <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  # fit0 <- lm(formula0, df)
  
  # pheno <- cross$pheno[,pheno.col]
  # cage <- cross$pheno[,"Cage"]
  # batch <- cross$pheno[,"Batch"]
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      
      #-----------------------------Create models------------------------------#
      ### For dynamically creating random effect terms
      #fit <- lmer(pheno ~ cross$geno[[c]]$dos[,k] +
      #              (1|cross$pheno[,randoms[1]]) + (1|cross$pheno[,randoms[2]]))
      
      fit <- lmer(pheno ~ geno + (1|batch/cage))
      
      # fit0 <- lm(pheno ~ sex + batch)
      # fit1 <- lm(pheno ~ sex + batch + geno)
      
      #----------------------------Calculate LODs------------------------------#
      # aov <- anova(fit1, fit0)
      # Fval <- aov$F[2] 
      # df <- abs(aov$Df[2])
      # lod <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      t <- summary(fit)$coefficients['geno','t value']
      Fval <- t^2
      df <- 1 # (in backcross, 2 in f2)
      lod <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      #-------------------------------Save LODs--------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      lods[marker.pos] = lod
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods) 
  names(model.df) <- c("chr", "pos", "lod") # same column names as scanone object
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
mv_qtl <- function(cross, pheno.cols, covar = NULL){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  if (!with(cross$geno[[1]], exists('dos'))){
    cross <- calc.dosage(cross) # calculate dosage (of B allele)
    print("Running calc.dosage...")
  }
  
  # make phenotype matrix 
  phenos <- cross$pheno[,pheno.cols]
  mice_w_obs <- which(rowSums(is.na(phenos)) != ncol(phenos))
  phenos <- as.matrix(phenos[mice_w_obs,])
  sex <- covar[mice_w_obs, "sex"]
  
  n <- nrow(phenos) # number of mice / observations
  p <- calc.num.markers(cross) # number of markers
  # Pre-allocate vectors for results 
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods <- vector(length = p)
  pvals <- vector(length = p)
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[mice_w_obs,k]
      
      #-----------------------------Create models------------------------------#
      mlm <- lm(phenos ~ sex + geno)
      aov <- Anova(mlm)
      Fval <- get_summary_for_print(aov)['geno', 'approx F']
      df <- get_summary_for_print(aov)['geno', 'Df']
      
      #----------------------------Calculate LODs------------------------------#
      lod <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval <- get_summary_for_print(aov)['geno','Pr(>F)']
      pval <- -log10(pval)
      
      #-------------------------------Save LODs--------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      lods[marker.pos] = lod
      pvals[marker.pos] = pval
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods, pvals) 
  names(model.df) <- c("chr", "pos", "lod", "Pval") # same column names as scanone object
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
# Output: scanone object with three LOD columns  
### TO-DO: 
###   1. Handle X chromosome markers - shouldn't use dosage? 
###   2. Handle variance heterogeneity
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt <- function(cross, pheno.col, addcovar, trt){
  # Calculate genotype probabilities if necessary
  if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }

  cross <- calc.dosage(cross) # calculate dosage (of B allele)

  p <- calc.num.markers(cross) # number of markers
  n <- nrow(cross$pheno) # number of mice / observations

  # Pre-allocate vectors for results
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p)
  lods2 <- vector(length = p)
  lods3 <- vector(length = p)

  lmdata <- as.data.frame(addcovar) # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(addcovar)) # to hold genotypes

  # Create null model: just covariates
  formula0 <- paste0("pheno ~ ", paste(colnames(addcovar), collapse=" + "))
  fit0 <- lm(formula0, lmdata)
  
  formula1 <- paste0(formula0, " + geno")
  formula2 <- paste0(formula1, " + ", paste0("geno*", trt))

  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      lmdata$geno <- geno # update geno col with marker genotype data

      #-----------------------------Create models------------------------------#
      fit1 <- lm(formula1, lmdata)
      fit2 <- lm(formula2, lmdata)

      #----------------------------Calculate LODs------------------------------#
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in
      lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2])
      lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component.
      aov <- anova(fit0, fit2)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2])
      lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]

      lods1[marker.pos] = lod1
      lods2[marker.pos] = lod2
      lods3[marker.pos] = lod3
    }
  }

  # Create object
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3)

  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3")

  rownames(model.df) <- markers
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)

  return(model.df)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling handling variance heterogeneity with a dglm 
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
#   ### Note - handles variance heterogeneity between groups 
# Output: scanone object with three LOD columns  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt.dglm <- function(cross, pheno.col, covar, trt){
  # Calculate genotype probabilities if necessary
  # if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
  #   cross <- calc.genoprob(cross)
  #   print("Running calc.genoprob...")
  # }
  #cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n <- nrow(cross$pheno) # number of mice / observations
  
  # Pre-allocate vectors for results
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p) # marginal QTL
  lods2 <- vector(length = p) # GxT QTL
  lods3 <- vector(length = p) # marginal + GxT QTL 
  lods4 <- vector(length = p) # vQTL
  
  lmdata <- covar # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(covar)) # to hold genotypes
  
  # Create null model: just covariates
  formula0str <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  formula0 <- as.formula(formula0str)
  dformulastr <- paste0("~", trt)
  dformula <- as.formula(dformulastr)
  fit0 <- dglm(formula0, as.formula(dformula), data = lmdata, method = 'ml')
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      
      #geno <- cross$geno[[c]]$dos[,k]
      geno <- as.factor(cross$geno[[c]]$data[,k])
      lmdata$geno <- geno # update geno col with marker genotype data
      
      #-----------------------------Create models------------------------------#
      formula1str <- paste0(formula0str, " + geno")
      formula1 <- as.formula(formula1str)
      fit1 <- dglm(formula1, dformula, data = lmdata, method = 'ml')
      
      formula2str <- paste0(formula1str, " + ", paste0("geno*", trt))
      formula2 <- as.formula(formula2str)
      fit2 <- dglm(formula2, dformula, data = lmdata, method = 'ml')
      
      dformulavstr <- paste0(dformulastr, " + geno")
      dformulav <- as.formula(dformulastr)
      fitv <- dglm(formula1, dformulav, data = lmdata, method = 'ml')
      
      #----------------------------Calculate -logP-----------------------------#
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in
      # lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval1 <- -log10(aov$Seq.P[1])
      
      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2])
      # lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval2 <- -log10(aov$Seq.P[1])
      
      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component.
      aov <- anova(fit0, fit2)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2])
      # lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval3 <- -log10(aov$Seq.P[1])
      
      # vQTL 
      aov <- anova(fit1, fitv)
      pval4 <- -log10(aov$Seq.P[1])
      
      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      
      lods1[marker.pos] = pval1
      lods2[marker.pos] = pval2
      lods3[marker.pos] = pval3
      lods4[marker.pos] = pval4
    }
  }
  
  # Create object
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3, lods4)
  
  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3", "lod4")
  
  rownames(model.df) <- markers
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling with random effect modeling 
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
#   randoms: dataframe of covariates to model as random effects 
# Output: list of three scanone objects 
### TO-DO: 
###   1. Handle variance heterogeneity 
###   2. Handle random effect covariates? 
### NOT SURE HOW TO COMPARE FIT1 TO FIT0 WHEN BOTH HAVE RANDOM EFFECTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt_rand <- function(cross, pheno.col, covar, trt, randoms = NULL){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  
  cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n = nrow(cross$pheno) # Number of mice / observations
  
  # Pre-allocate vectors for results 
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p)
  lods2 <- vector(length = p)
  lods3 <- vector(length = p)
  
  lmdata <- cbind(covar,randoms) # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(covar)) # to hold genotypes 

  # Create null model: just covariates 
  fixedcovar <- paste(colnames(covar), collapse=" + ")
  randcovar <- paste0("(1|", colnames(randoms), ")", collapse=" + ")
  formula0 <- paste0("pheno ~ ", fixedcovar, " + ", randcovar)
  fit0 <- lmer(formula0, lmdata)
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr 
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker 
      lmdata$geno <- NA # clear geno col 
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      lmdata$geno <- geno # update geno col with marker genotype data 
      
      #-----------------------------Create models------------------------------#
      formula1 <- paste0(formula0, " + geno")
      fit1 <- lmer(formula1, lmdata)
      
      formula2 <- paste0(formula1, " + ", paste0("geno*", trt))
      fit2 <- lmer(formula2, lmdata)
      
      #----------------------------Calculate LODs------------------------------#
      ### Need to use MLE instead of LMER with ANOVA on LMMs 
      ### ANOVA gives you a Chisq value (not an F val) and p-val
      
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in 
      lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2])
      lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at 
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component. 
      aov <- anova(fit0, fit2)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2])
      lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      
      lods1[marker.pos] = lod1
      lods2[marker.pos] = lod2
      lods3[marker.pos] = lod3
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3) 
  
  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3") 
  
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}


#-------------------------------Permutation tests------------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to perform permutation tests more efficiently than R/qtl. Instead 
# of performing n.perm tests from start to finish and logging the highest LOD, 
# perform n.perm tests at each marker. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
scanone_perm <- function(cross, pheno.col, addcovar, n.perm){
  # create 1000 permutations of the phenotype, add to matrix 
  
  # lm multi response, add to matrix of results 
  
  # find highest LOD from each genome scan and create distribution 
  
  # return - vector of 1000 LOD scores with some attributes 
  # method = "em"
  # model = "normal" 
  # type = "bc"
  # class = scanoneperm, matrix 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Permutation tests for GxT 
# covar must contain trt 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt_perm <- function(cross, pheno.col, perm.strata=NULL, covar, n.perm, trt){
  
  pheno.og <- pull.pheno(combined, pheno.col)
  # remove NAs bc don't want to permute those 
  # use !is.na(pheno.og) to filter genotypes later on 
  pheno <- pheno.og[!is.na(pheno.og)]
  
  # create 1000 permutations of phenotype 
  pheno.perms <- matrix(nrow = length(pheno), ncol = n.perm)
  for (c in 1:n.perm){
    pheno.perms[,c] <- pheno[sample(length(pheno))]
  }
  
  # calculate dosage for hk regression 
  if (!with(cross$geno[[1]], exists('dos'))){  # calc.dosage hasn't been run
    print("Running calc.dosage...")
    cross <- calc.dosage(cross)
  }
  
  # at each marker, run lm.multiresponse(), add to matrix of results 
  # create empty matrix for storing results 
  permLODs <- matrix(nrow = n.perm, ncol = calc.num.markers(cross))
  
  lmdata <- as.data.frame(covar[!is.na(pheno.og),]) # data frame for lm()
  lmdata$geno <- rep(NA, nrow(lmdata)) # to hold genotypes
  
  # Create null model: just covariates
  formula0 <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  formula1 <- paste0(formula0, " + geno")
  formula2 <- paste0(formula1, " + ", paste0("geno*", trt))
  
  #fit0 <- lm(formula0, lmdata)
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      geno <- geno[!is.na(pheno.og)] # remove genotypes for mice that don't have the phenotype 
      
      lmdata$geno <- geno # update geno col with marker genotype data
      
      #-----------------------------Create models------------------------------#
      # fit1 <- lm(formula1, lmdata)
      # fit2 <- lm(formula2, lmdata)
      
      #----------------------------Calculate LODs------------------------------#
      res <- lm.multiresponse(formula = formula2, 
                       response.matrix = pheno.perms, 
                       data = lmdata, 
                       null.formula = formula0, 
                       LOD = TRUE)
      
      #------------------------------Save LODs---------------------------------#
      resLODs <- res$LOD
      permLODs[,marker.pos] <- resLODs
    }
  }
  
  # find highest LOD from each genome scan and create distribution 
  maxLODs <- apply(permLODs, 1, max)
  
  # return - vector of 1000 LOD scores with some attributes 
  attributes(maxLODs) <- list(method = 'hk', model = 'normal',
                              type = 'f2', class = c('scanoneperm', 'matrix'))
  
  return(maxLODs) 
}










#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots with ggplot
### THIS VERSION OF PXG() PLOTS IMPUTED GENOTYPES AS GRAY AND DIRECTLY GENOTYPED
### GENOTYPES AS BLACK. NO LONGER USED 
# Alternative to qtl::plotPXG, uses ggplot and more intuitive error bars  
# Input:
#     cross: r/qtl cross object 
#     pheno: vector of phenotype values for all mice 
#     marker: marker whose genotype to plot
#     geno.map: list mapping A and B allele to more informative labels 
#     xlab: x axis label
#     ylab: y axis label 
#     title: title
#     theme: ggplot theme
#     type: scatter (default), boxplot or violin 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_old <- function(cross, pheno, marker, geno.map, xlab = NULL, ylab, 
                title = NULL, bestfit = TRUE, theme = rmd_theme, type = 'scatter', ...){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if ((chr == 'X')&('M' %in% levels(cross$pheno$sex))){bestfit = FALSE} 
  
  # plot title 
  # if (is.null(title)){
  #   title = paste('Genotype vs. phenotype at chr', chr, 'marker', marker)
  # }
  
  if(is.null(xlab)){
    qtl_name <- qtl_map$qtl_name[qtl_map$chr == chr]
    xlab <- as.expression(bquote(italic(.(qtl_name))*" (chr"*.(chr)*") "*Genotype))
  }
  
  marker.genos <- cross$geno[[chr]]$data[,marker]
  which.missing <- which(is.na(marker.genos))
  # calculate imputed genotypes, for plotting and regression. Need to keep separate
  # from observed genotypes for plotting (imputed genotypes gray, observed genotypes black)
  cross.imp <- fill.geno(cross) 
  imp.marker.genos <- cross.imp$geno[[chr]]$data[,marker]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ imp.marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos,
                   imp_geno = imp.marker.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
    df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2), to = geno.labels) 
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$imp_geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$imp_geno <- X.data$X.geno
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                               to = geno.labels) 
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$imp_geno <- factor(df$imp_geno, levels = geno.labels)
  
  dfsum <- data.frame(imp_geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$imp_geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$imp_geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    if (length(which.missing) > 0){
      p <- ggplot(data = df[-which.missing,], mapping = aes(x = geno, y = pheno)) +
        geom_jitter(width = 0.1, height = 0) + scale_x_discrete(drop = FALSE) + 
        geom_jitter(data = df[which.missing,], mapping = aes(x = imp_geno, y = pheno), color = "gray", width = 0.1) +
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) +
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) + 
        labs(title = title, y = ylab, x = xlab) + 
        {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
        theme
    } else { 
      p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
        geom_jitter(width = 0.1, height = 0) + 
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) + 
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) +
        labs(title = title, y = ylab, x = xlab) + 
        {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
        theme
    }
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_boxplot(fill = "#1B9E77", notch = TRUE) + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  } else if (type == 'violin'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_violin(fill = "#1B9E77") + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((type == 'scatter') & (bestfit == TRUE)){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}



