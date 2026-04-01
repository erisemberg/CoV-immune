library(qtlpvl)
library(qtl)
library(tidyverse)
library(readxl)
library(extRemes)
library(ggplotify)
library(cowplot)
library(grid)
source("code-dependencies/cov_qtl_functions.R")

# ---------------------------------load data---------------------------------- #
cross_data <- read_csv("derived_data/cross_data.csv")

# load cross for QTL mapping 
cross <- read.cross(format = "csv", 
                    file = "derived_data/RqtlCC006xCC044_ctrlAndSARS.csv",
                    na.strings=c("-","NA", "na", "<NA>", "not tested"), 
                    genotypes=c("AA","AB","BB"),
                    alleles=c("A","B")) # A = CC006, B = CC044
cross <- jittermap(cross)
cross <- calc.genoprob(cross)

pheno_names <- read_excel("source_data/pheno_names.xlsx")
all_phenos <- pheno_names$flow_col_name
q <- length(all_phenos)
flow_phenos <- all_phenos[-1]

mice_w_any_flow <- rowSums(is.na(cross_data[,flow_phenos])) < length(flow_phenos)
flow <- cross_data[mice_w_any_flow,flow_phenos]
geno <- cross_data[mice_w_any_flow, 
                   c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]


groups <- c("PBS", "SARSCoV", "SARS2CoV", "GxT")
strat_groups <- c("PBS", "SARSCoV", "SARS2CoV")
group_dirs <- list("PBS" = "PBS", "SARSCoV" = "SARS", "SARS2CoV" = "SARS2", "GxT" = "GxT")
group_lookup <- c("C" = "PBS", "1" = "SARSCoV", "2" = "SARS2CoV", "A" = "GxT")

geno_map <- list(A = "CC006", B = "CC044")

qtl_map <- read_csv("results/qtl_map.csv")

# ------------------------------adjust p-values------------------------------- #
for (group in groups){
  group_dir <- group_dirs[[group]]
  minPs <- NULL
  for (k in 1:q){
    pheno <- all_phenos[k]
    if (group %in% c("PBS", "GxT") & pheno == "Titer"){ 
      minPs <- c(minPs, NA)
      next 
    }
    mod_fname <- paste0("results/qtl_mapping/modRDS/", group_dir, "/", pheno, ".rds")
    mod <- readRDS(mod_fname)
    perm_fname <- paste0("results/qtl_mapping/permRDS/", group_dir, "/", pheno, ".rds")
    perm <- readRDS(perm_fname)
    
    if (group %in% strat_groups & pheno != "Titer"){
      max_nlogP <- max(mod$lod)
      fevdsum <- summary.fevd(perm)
    } else if (group == "GxT"){ 
      max_nlogP <- max(unlist(mod[,3:5]))  
      fevdsum <- summary.fevd(perm)
    } else { # Titer 
      mod_df <- as.data.frame(unclass(mod))[,3:5]
      max_nlogP <- max(mod_df, na.rm = TRUE) # for simplicity, calling nlogP even though its an LOD for Titer 
      max_col <- which(mod_df == max_nlogP, arr.ind = TRUE)[2] 
      fevdsum <- summary(fevd(unclass(perm)[,max_col]))
    }
    minP <- pevd(q = max_nlogP, loc = fevdsum$par[1], scale = fevdsum$par[2],
                 shape = fevdsum$par[3], lower.tail = FALSE)
    minPs <- c(minPs, minP)
  }
  
  # run p.adjust on minimum p-values from this group, create df and save 
  df <- data.frame(phenotype = all_phenos, 
                   pval = minPs, 
                   qval = p.adjust(minPs, method = "BH"), 
                   rank = rank(minPs, ties.method = "last"))
  df <- df[order(df$rank),]
  write_csv(df, paste0("results/qtl_mapping/p_adjust/", group_dir, ".csv"))
}

# -------------------------------save all pvals------------------------------- #
all_pvals <- NULL
lod_col <- 1 # update when handling GxT, Titer 
for (group in groups){
  group_dir <- group_dirs[[group]]
  # function for calculating experiment-wide q-vals from genome-wide p-vals 
  fdr_df <- read_csv(paste0("results/qtl_mapping/p_adjust/", group_dir, ".csv"), show_col_types = FALSE)
  q_from_p <- make_interpolate_fn(fdr_df)
  
  for (k in seq_len(length(all_phenos))){
    pheno <- all_phenos[k]
    if (pheno == "Titer"){ next } # don't show titer in this plot, bc there are 2 p-vals
    
    # nominal (neg log10) p-vlaues 
    mod_fname <- paste0("results/qtl_mapping/modRDS/", group_dir, "/", pheno, ".rds")
    mod <- readRDS(mod_fname)
    df <- as.data.frame(mod)
    if (group == "GxT"){
      df$negLogP <- apply(df[,3:5], 1, max) ### SHOULD THIS USE OVERALL P-VAL INSTEAD OF MAX? 
      df <- df[,-c(3:5)]
    }
    colnames(df)[3] <- c("negLogP")
    
    # genome-wide adjusted p-values 
    perm_fname <- paste0("results/qtl_mapping/permRDS/", group_dir, "/", pheno, ".rds") 
    perm <- readRDS(perm_fname)
    fevdsum <- summary(fevd(mod[[2+lod_col]]))
    
    df$adjP <- pevd(q = df$negLogP, 
                    loc = fevdsum$par[1], scale = fevdsum$par[2], shape = fevdsum$par[3], 
                    lower.tail = FALSE)
    
    # experiment-wide adjusted q-values
    df$qval <- q_from_p(df$adjP)
    
    df <- df %>% 
      mutate(infection = group, phenotype = pheno, .before = "chr")
    
    all_pvals <- rbind(all_pvals, df)
  }
}
write_csv(all_pvals, "results/qtl_mapping/all_pvals.csv")

# ------------------------------QTL summary plot------------------------------ #
# create table with marker indices
marker_index <- all_pvals %>%
  distinct(chr, pos) %>%
  mutate(marker_idx = row_number())

# order by hclust based on phenotypic correlation 
# corr <- cor(flow, use = "pairwise")
# hc <- hclust(as.dist(1 - abs(corr)), method = "complete")
# pheno_order <- flow_cols[hc$order]

df_plot <- all_pvals %>%
  mutate(adjP_clip = pmin(adjP, 0.1)) %>%
  left_join(marker_index, by = c("chr", "pos")) 
#,phenotype = factor(phenotype, levels = pheno_order)

# table for labeling x-axis with chromosome boundaries 
marker_map <- df_plot %>%
  distinct(marker_idx, chr) %>%
  arrange(marker_idx)
chr_info <- marker_map %>%
  group_by(chr) %>%
  summarise(xmin = min(marker_idx) - 0.5,
            xmax = max(marker_idx) + 0.5,
            xmid = (xmin + xmax) / 2,
            .groups = "drop") %>%
  arrange(as.numeric(as.character(chr)) %>% ifelse(is.na(.), Inf, .))  # optional: numeric chr order
# boundaries between chromosomes (for vertical lines)
chr_boundaries <- chr_info %>%
  arrange(xmin) %>%
  slice(-1) %>%
  transmute(x = xmin)
chr_bands <- chr_info %>%
  mutate(band = as.integer(factor(chr)) %% 2)


# order by hclust based on similarity of QTL scan
mat <- df_plot %>%
  select(phenotype, infection, marker_idx, adjP_clip) %>%
  mutate(feature = paste(infection, marker_idx, sep = "_")) %>%
  select(-infection, -marker_idx) %>%
  pivot_wider(names_from = feature, values_from = adjP_clip)
rn <- mat$phenotype
dat <- mat %>%
  select(-phenotype) %>%
  as.matrix()
rownames(dat) <- rn
dat[is.na(dat)] <- 0.1
hc <- hclust(dist(dat), method = "complete")
pheno_order <- rownames(dat)[hc$order]
ord_df <- data.frame(pheno = rownames(dat),
                     order = hc$order)
write_csv(ord_df, "results/qtl_mapping/hclust_pheno_order.csv")

# order by hclust based on genetic correlations 
# library(sommer)
# library(AGHmatrix)
# grm <- Gmatrix(as.matrix(geno), method = "VanRaden")
# colnames(grm) <- rownames(grm) <- cross_data$mouse_ID[mice_w_any_flow]
# dat <- data.frame(mouse_ID = cross_data$mouse_ID[mice_w_any_flow])
# dat <- cbind(dat, as.data.frame(flow))
# # sanity checks 
# all(dat$mouse_ID == rownames(grm))
# all(setdiff(names(dat), "mouse_ID") == flow_phenos)
# fixed_form <- as.formula(paste0("cbind(", paste(flow_phenos, collapse = ","), ") ~ 1"))
# fit <- mmer(
#   fixed = fixed_form,
#   random = ~ vsr(mouse_ID, Gu = grm),
#   rcov = ~ units,
#   data = dat
# )
# G <- cov2cor(fit$sigma$`u:mouse_ID`)

df_plot2 <- df_plot %>%
  mutate(phenotype = factor(phenotype, levels = pheno_order),
         infection = factor(infection, levels = c("PBS", "SARSCoV", "SARS2CoV", "GxT")))
infection_labs <- c(PBS = "Control", SARSCoV = "SARS-CoV", SARS2CoV = "SARS-CoV-2", GxT = "Combined")

rug_data <- df_plot2 %>%
  filter(adjP_clip < 0.05) %>%  
  group_by(marker_idx) %>%
  summarise(qtl_count = n(), .groups = 'drop')

chr_info <- chr_info %>%
  mutate(label_y = ifelse(row_number() %% 2 == 1, 0, -1))

qtlsum <- ggplot(df_plot2, aes(x = marker_idx, y = phenotype, fill = adjP_clip)) +
  geom_tile() +
  facet_grid(infection ~ ., scales = "fixed", space = "fixed", 
             labeller = labeller(infection = infection_labs)) +
  geom_vline(data = chr_boundaries, aes(xintercept = x), inherit.aes = FALSE, linewidth = 0.2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(limits = c(0,0.1), direction = -1, name = "p-value") +
  labs(x = NULL, y = "Phenotype") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(vjust = 7),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(l = -5, t = 15),
        legend.justification = "top",
        legend.text = element_text(size = 8),
        strip.text = element_text(face = "bold"))
# histogram 
rug_plot <- ggplot(rug_data, aes(x = marker_idx, y = qtl_count)) +
  geom_col(fill = "gray30", width = 1) +
  geom_vline(data = chr_boundaries, aes(xintercept = x), inherit.aes = FALSE, linewidth = 0.2) +
  scale_x_continuous(limits = range(df_plot2$marker_idx), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::breaks_pretty(n = 3)) +
  labs(y = "QTL count", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
# chromosome label plot with alternating positions
chr_plot <- ggplot() +
  geom_blank() +
  scale_x_continuous(limits = range(df_plot2$marker_idx), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3, 0.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 10, 5)) +
  coord_cartesian(clip = "off")
# add alternating chromosome labels
for(i in 1:nrow(chr_info)) {
  chr_plot <- chr_plot + 
    annotation_custom(grob = textGrob(chr_info$chr[i], gp = gpar(fontsize = 8)),
                      xmin = chr_info$xmid[i], xmax = chr_info$xmid[i],
                      ymin = chr_info$label_y[i], ymax = chr_info$label_y[i])
}
# add "Chromosome" label
chr_plot <- chr_plot + 
  annotation_custom(grob = textGrob("Chromosome", gp = gpar(fontsize = 11)),
                    xmin = mean(range(df_plot2$marker_idx)), 
                    xmax = mean(range(df_plot2$marker_idx)),
                    ymin = -2.7, ymax = -2.7)

qtl_hm <- plot_grid(qtlsum, rug_plot, chr_plot, nrow = 3, rel_heights = c(10, 1, 0.5), align = "v", axis = "lr")
qtl_hm
# use this if background is not blue to distinguish chromosomes 
#   geom_rect(data = chr_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
#                                 alpha = factor(band)),
#           inherit.aes = FALSE, fill = "grey50") +
#   scale_alpha_manual(values = c(0.06, 0.00), guide = "none") +

# -------------------------------summarize QTL-------------------------------- #
#all_qtlv1 <- read_csv("results/qtl_mapping/all_qtl.csv")
all_qtl <- data.frame(infection = character(),
                      phenotype = character(),
                      chr = numeric(),
                      pos = numeric(),
                      CI_lwr = numeric(),
                      CI_upr = numeric(),
                      negLogP = numeric(),
                      adjP = numeric(),
                      qval = numeric(),
                      var_explained = numeric(),
                      beta = numeric())

ix <- 0; lod_col <- 1
for (group in groups){
  #method <- group_methods[[group]]
  group_dir <- group_dirs[[group]]
  
  if (group %in% strat_groups){ 
    cross_g <- subset(cross, ind = (cross$pheno$infection == group)) 
  } else { # GxT
    cross_g <- cross 
  }
  
  fdr_df <- read_csv(paste0("results/qtl_mapping/p_adjust/", group_dir, ".csv"), show_col_types = FALSE)
  q_from_p <- make_interpolate_fn(fdr_df)
  
  for (k in seq_len(length(all_phenos))){
    pheno <- all_phenos[k]
    if (pheno == "Titer" & group %in% c("PBS", "GxT")){ next } 
    
    if (group  %in% strat_groups){
      pheno_data <- cross_data[cross_data$infection == group, pheno][[1]]
    } else { # GxT
      pheno_data <- cross_data[[pheno]] 
    }
    
    mod_fname <- paste0("results/qtl_mapping/modRDS/", group_dir, "/", pheno, ".rds")
    mod <- readRDS(mod_fname)
    perm_fname <- paste0("results/qtl_mapping/permRDS/", group_dir, "/", pheno, ".rds") 
    perm <- readRDS(perm_fname)
    
    if (pheno == "Titer"){ # 2-part model; permutation object has 3 columns vs one for GxT QTL mapping 
      sum <- summary(mod, perms = perm, alpha = 0.05) # using Rqtl functionality instead of GEV, but good enough
    } else if (group == "GxT"){ # check if any of the 3 p-values are greater than the threshold 
      sum <- summary(mod, threshold = summary(perm)[1], format = "allpheno")
      if (sum(duplicated(sum$chr)) > 0){ # deduplicate bc each p-value counts 
        unique_chrs <- unique(sum$chr)
        dedupd_sum <- NULL
        for (j in 1:length(unique_chrs)){
          chr_sum <- sum[sum$chr == unique_chrs[j],]
          # choose a row based on highest p-value 
          max_p <- max(unlist(chr_sum[,3:5]))
          chosen_row <- chr_sum[(max_p == chr_sum$p_overall)|(max_p == chr_sum$p_G)|(max_p == chr_sum$p_GxT),]
          dedupd_sum <- rbind(dedupd_sum, chosen_row)
        }
        sum <- dedupd_sum
      }
    } else { # only need to check the one column (and there won't be duplicates)
      sum <- summary(mod, threshold = summary(perm)[1])
    }
    
    if (nrow(sum) == 0){ next } 
    if (pheno != "Titer"){ fevdsum <- summary.fevd(perm) }
    
    for (j in 1:nrow(sum)){
      qtl_driver <- NA
      ix <- ix + 1
      # unadjusted p-value 
      if (pheno == "Titer"){
        unadj_lod <- max(sum$lod.p.mu[j], sum$lod.p[j], sum$lod.mu[j])
        lod_col <- which.max(c(sum$lod.p[j], sum$lod.mu[j])) # lod.p.mu is the sum, so want to know which is the main driver of p and mu 
        qtl_driver <- ifelse(lod_col == 1, "binary", "continuous")
      } else if (group == "GxT"){
        unadj_negLogP <- max(sum$p_overall[j], sum$p_G[j], sum$p_GxT[j]) 
        lod_col <- which.max(c(sum$p_overall[j], sum$p_G[j], sum$p_GxT[j]))
        qtl_driver <- ifelse(lod_col == 2, "genetic", 
                             ifelse(lod_col == 3, "GxT", "overall"))
      } else { 
        unadj_negLogP <- sum$lod[j]
        lod_col <- 1
      }
      
      # chr:position(interval)
      chr <- sum$chr[j]
      bayesint_j <- bayesint(mod, chr = chr, lodcolumn = lod_col)
      lwr <- bayesint_j[1,"pos"]
      upr <- bayesint_j[3,"pos"]
      chr_pos <- paste0("Chr ", chr, ": ", 
                        format(round(sum$pos[j], 2), nsmall = 2), 
                        " (", format(round(lwr, 2), nsmall = 2), "-", 
                        format(round(upr,2), nsmall = 2), ")")
      
      # genome-wide adjusted p-value 
      if (pheno != "Titer"){
        adjP <- pevd(q = unadj_negLogP, loc = fevdsum$par[1], scale = fevdsum$par[2], shape = fevdsum$par[3], lower.tail = FALSE)
      } else {
        fevdsum <- summary(fevd(mod[[2+lod_col]]))
        adjP <- pevd(q = unadj_lod, loc = fevdsum$par[1], scale = fevdsum$par[2], shape = fevdsum$par[3], lower.tail = FALSE)
      }
      # experiment-wide (FDR) adjusted q-value 
      qval <- q_from_p(adjP)
      
      # variance explained by QTL 
      cov_list <- c("sex")
      if (chr == "X"){ cov_list <- c(cov_list, "pgm") }
      if (group == "GxT"){ cov_list <- c(cov_list, "infection") }
      var_expl <- get_var_expl(cross = cross_g, pheno = pheno_data, marker = rownames(sum)[j], cov_list = cov_list, gxt = (group == "GxT"))
      
      all_qtl[ix,] <- c(group, pheno, chr, sum$pos[j], lwr, upr, unadj_negLogP, adjP, qval, var_expl$h2, var_expl$add_effect)
    }
  }
}

write_csv(all_qtl, "results/qtl_mapping/all_qtl.csv")

# ----------------------------remove duplicates------------------------------- #
condensed_qtl <- all_qtl %>%
  mutate(CI_lwr = as.numeric(CI_lwr),
         CI_upr = as.numeric(CI_upr),
         qval = as.numeric(qval),
         infection = case_when(
           infection == "PBS" ~ "C",
           infection == "SARSCoV" ~ "1",
           infection == "SARS2CoV" ~ "2",
           infection == "GxT" ~ "A"),
         sig = case_when(
           qval <= 0.05 ~ "***",
           qval <= 0.1 ~ "**",
           adjP <= 0.05 ~ "*",
           TRUE ~ ""),
         inf_sig = paste0(infection, sig)) %>%
  group_by(chr, phenotype) %>%
  summarise(
    inf_sig1 = inf_sig[which.min(qval)],
    inf_sig_rest = paste(inf_sig[!(inf_sig %in% inf_sig1)], collapse = ", "),
    infection_sig = paste0(inf_sig1, ifelse(inf_sig_rest != "", ", ", ""), inf_sig_rest),
    infection_detect = {
      infs <- unique(infection[qval < 0.05])
      paste(infs, collapse = ",")
    },
    pos = as.numeric(pos[which.min(qval)]),               
    CI_lwr_conservative = min(CI_lwr, na.rm = TRUE),       
    CI_upr_conservative = max(CI_upr, na.rm = TRUE),       
    CI_lwr_liberal = max(CI_lwr, na.rm = TRUE),
    CI_upr_liberal = min(CI_upr, na.rm = TRUE),
    var_expl = as.numeric(var_explained[which.min(qval)]),
    qval = min(qval, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-inf_sig1, -inf_sig_rest) %>%
  filter(qval <= 0.10) %>% # only select QTL with at least one association where q < 0.1
  rowwise() %>%
  mutate(nominal_infs = test_nominal_sig(chr, pos, phenotype, infection_sig, cross_data, sig_thresh = 0.1),
         infection_sig = if_else(nominal_infs != "", 
                                 paste(infection_sig, nominal_infs, sep = ", "),
                                 infection_sig)) %>%
  ungroup() %>% 
  select(-nominal_infs) %>%
  mutate(sig_level = ifelse(qval <= 0.05, "significant", "suggestive")) %>%
  arrange(sig_level, as.numeric(chr))

write_csv(condensed_qtl, "results/qtl_mapping/condensed_qtl.csv")
#condensed_qtl <- read_csv("results/qtl_mapping/condensed_qtl.csv")

# how many control QTL?
condensed_qtl %>%
  filter(str_detect(infection_sig, "C\\*\\*")) %>%
  nrow()

# how many SARS-CoV QTL?
condensed_qtl %>%
  filter(str_detect(infection_sig, "1\\*\\*")) %>%
  nrow()

# how many SARS-CoV-2 QTL?
condensed_qtl %>%
  filter(str_detect(infection_sig, "2\\*\\*")) %>%
  nrow()

# do QTL from SARS2 mice explain greater variance on avg?
qtl2 <- condensed_qtl %>%
  filter(str_detect(infection_sig, "2\\*\\*"))
qtl1C <- condensed_qtl %>%
  filter(str_detect(infection_sig, "[1C]\\*\\*"))
var_expl_2 <- qtl2$var_expl
var_expl_1c <- qtl1C$var_expl
t.test(var_expl_2, var_expl_1c, var.equal = FALSE)
# yes, mean of SARS2 QTL = 35.5%, mean of control/SARS QTL = 12.4%





# -----------------------------create pub tables------------------------------ #
cond_qtl_disp <- condensed_qtl %>% 
  left_join(pheno_names, by = c("phenotype" = "flow_col_name")) %>%
  left_join(qtl_map, by = c("chr", "phenotype")) %>% 
  mutate(chr_pos = paste0("Chr ", chr, ": ", round(pos, 1), " (", round(CI_lwr_conservative, 1), "-", round(CI_upr_conservative, 1), ")"),
         qval = ifelse(qval == 0, "< 2.2e-308", signif(qval, 2)),
         var_expl = paste0(round(var_expl*100, 1), "%")) %>%
  select(sig_level, qtl_id, flow_display_name, infection_sig, chr_pos, qval, var_expl)

writexl::write_xlsx(cond_qtl_disp, "results/qtl_mapping/condensed_qtl.xlsx")


# -----------------------------QTL by infection------------------------------- #
qtl_by_infection <- cond_qtl_disp %>%
  distinct(qtl_id, infection_sig) %>%   # avoid phenotype duplication
  mutate(has_C = str_detect(infection_sig, "C"),
         has_1 = str_detect(infection_sig, "1"),
         has_2 = str_detect(infection_sig, "2")) %>%
  group_by(qtl_id) %>%
  summarise(Control = any(has_C),
            SARS = any(has_1),
            SARS2 = any(has_2),
            .groups = "drop") %>% 
  mutate(n_groups = Control + SARS + SARS2,
         class = case_when(n_groups == 1 ~ "infection_specific",
                           n_groups > 1 ~ "shared"))

shared_qtls <- qtl_by_infection %>%
  filter(class == "shared") %>%
  pull(qtl_id)

specific_qtls <- qtl_by_infection %>%
  filter(class == "infection_specific")

specific_by_group <- specific_qtls %>%
  mutate(infection = case_when(Control ~ "Control",
                               SARS ~ "SARS-CoV",
                               SARS2 ~ "SARS-CoV-2"))

# -------compare stratified to QTL ONLY identified by combined approach------- #
gxt_only_qtl <- condensed_qtl %>%
  mutate(
    has_A_sig = str_detect(infection_sig, "A\\*\\*+"),  # A has ** or ***
    has_strat_sig = str_detect(infection_sig, "[C12]\\*\\*+")  # C, 1, or 2 has ** or ***
  ) %>%
  filter(has_A_sig & !has_strat_sig) %>%
  select(-has_A_sig, -has_strat_sig)
dim(gxt_only_qtl)
sum(grepl("A\\*\\*\\*", gxt_only_qtl$infection_sig))

strat_qtl <- condensed_qtl %>%
  mutate(
    has_strat_sig = str_detect(infection_sig, "[C12]\\*\\*+"), # C, 1, or 2 has ** or ***
    has_C_sig = str_detect(infection_sig, "C\\*\\*"),
    has_1_sig = str_detect(infection_sig, "1\\*\\*"),
    has_2_sig = str_detect(infection_sig, "2\\*\\*"),
    num_subthreshold = str_count(infection_sig, "[C12]") - str_count(infection_sig, "[C12]\\*\\*") # sub-threshold associations, not incl. A 
  ) %>%
  filter(has_strat_sig) %>%
  select(-has_strat_sig)

strat_only_qtl <- condensed_qtl %>%
  mutate(
    has_A_sig = str_detect(infection_sig, "A\\*\\*+"),  # A has ** or ***
    has_strat_sig = str_detect(infection_sig, "[C12]\\*\\*+")  # C, 1, or 2 has ** or ***
  ) %>%
  filter(!has_A_sig & has_strat_sig) %>%
  select(-has_A_sig, -has_strat_sig)
dim(strat_only_qtl)

# variance explained lower for GxT QTL? Yes
var_expl_gxt <- gxt_only_qtl$var_expl
var_expl_strat <- strat_qtl$var_expl
t.test(var_expl_gxt, var_expl_strat, var.equal = FALSE)

# ----------------------------summaries for Fig 4----------------------------- #
# Create binary indicators for each group's presence (regardless of significance level)
qtl_presence <- condensed_qtl %>%
  mutate(
    has_C = str_detect(infection_sig, "C"),
    has_1 = str_detect(infection_sig, "1"),
    has_2 = str_detect(infection_sig, "2"),
    has_A = str_detect(infection_sig, "A"),
    qtl_id = paste(chr, phenotype, sep = ":")
  )

# create list format for upset plot - each element is a vector of qtl_ids present in that group
qtl_sets <- list(
  Control = qtl_presence$qtl_id[qtl_presence$has_C],
  `SARS-CoV` = qtl_presence$qtl_id[qtl_presence$has_1],
  `SARS-CoV-2` = qtl_presence$qtl_id[qtl_presence$has_2]
) # ,Combined = qtl_presence$qtl_id[qtl_presence$has_A]

# set order
set_order <- c("SARS-CoV-2", "SARS-CoV", "Control")
qtl_sets_ordered <- qtl_sets[set_order]

# Create upset plot
upset_plot <- UpSetR::upset(
  UpSetR::fromList(qtl_sets_ordered), 
  sets = set_order, 
  nsets = length(set_order), 
  keep.order = TRUE,
  sets.x.label = "Number of associations")
upset_plot

# create upset plot grob for Fig 3
grid.edit("arrange", name = "upset")
us_grob <- grid.grab()
upset_plot <- us_grob
grid.arrange(us_grob, var_hm)

# fit <- euler(qtl_sets)
# plot(fit, quantities = TRUE)  

# -------------------------are MHCII+ IM QTL distinct?------------------------ #
pheno <- cross_data$M_IntMacs_pctMHCII
pheno_rows <- which(condensed_qtl$phenotype == "M_IntMacs_pctMHCII")
df <- data.frame(pheno = pheno)
for (i in pheno_rows){
  chr <- condensed_qtl$chr[i]
  marker <- find.marker(cross, chr = chr, pos = condensed_qtl$pos[i])
  geno <- pull.geno(cross, chr = chr)[,marker]
  geno_colname <- paste0("geno", chr)
  df[[geno_colname]] <- geno
}
summary(lm(pheno ~ geno7 + geno10 + geno16, data = df))


# --------------------------------LM per QTL---------------------------------- #
# fit same linear model for all QTL, save beta, p-val, and variance explained by 
# additive, dominant, and interaction effects 
lm_df <- data.frame(chr = numeric(), 
                    phenotype = character(), 
                    infection = character(), 
                    addBeta = numeric(),
                    addP = numeric(),
                    domBeta = numeric(),
                    domP = numeric())
dom_qtl <- data.frame(chr = numeric(), 
                      pos = numeric(),
                      phenotype = character(),
                      group = character(), 
                      IDd_in_grp = logical(),
                      addP = numeric(), 
                      domP = numeric())
qtl_var_comp <- data.frame(chr = numeric(), 
                           pos = numeric(),
                           phenotype = character(), 
                           partR2_add = numeric(), 
                           partR2_dom = numeric(), 
                           partR2_gxt = numeric())

ix <- lm_ix <- dom_ix <- 0
for (i in 1:nrow(condensed_qtl)){
  # genotype data
  chr_i <- condensed_qtl$chr[i]
  if (chr_i == "X") { next } # skip X chromosome for now
  ix <- ix + 1
  pos_i <- as.numeric(condensed_qtl$pos[i])
  marker_i <- find.marker(cross, chr = chr_i, pos = pos_i)
  geno_terms <- get_add_dom(cross, marker_i) # for autosomes 
  # phenotype data
  pheno_name <- condensed_qtl$phenotype[i]
  pheno <- cross_data[[pheno_name]]
  # associated infections (at least nominal)
  detected_in_infections <- str_trim(strsplit(condensed_qtl$infection_detect[i], split = ",")[[1]])
  
  dat <- data.frame(pheno = pheno, 
                    sex = as.numeric(cross_data$sex), 
                    infection = cross_data$infection, 
                    add = geno_terms$add, 
                    dom = geno_terms$dom)
  
  for (group in strat_groups){  
    lm_ix = lm_ix + 1
    if (pheno_name == "Titer" & group == "PBS"){
      lm_df[lm_ix,] <- c(chr_i, pheno_name, group, rep(NA,4))
      next
    }
    
    dat_g <- dat[dat$infection == group,]
    fit_full <- lm(pheno ~ sex + add + dom, data = dat_g)
    summ_lm <- summary(fit_full)
    # get betas from full model
    add_effect <- signif(summ_lm$coefficients["add", "Estimate"], 2)
    dom_effect <- signif(summ_lm$coefficients["dom", "Estimate"], 2)
    # get p-values from type I (sequential) ANOVA
    summ_aov <- summary(aov(fit_full))  
    addP <- signif(summ_aov[[1]]["add", "Pr(>F)"], 2)
    domP <- signif(summ_aov[[1]]["dom", "Pr(>F)"], 2)
    lm_df[lm_ix,] <- c(chr_i, pheno_name, group, add_effect, addP, dom_effect, domP)
    
    # is this is a (probably rare) QTL where we could only have identified it with the dominant term? if so, additive effect not useful 
    if (addP > 0.05 & domP < 0.05){
      dom_ix <- dom_ix + 1
      IDd_in_grp <- if (group %in% detected_in_infections) { TRUE } else { FALSE }
      dom_qtl[dom_ix,] <- c(chr_i, pos_i, pheno_name, group, IDd_in_grp, signif(addP, 2), signif(domP, 2))
      # plot PxG for dominant QTL
      cross_g <- subset(cross, ind = (cross$pheno$infection == group))
      pxg_fname <- paste0("figures/qtl_pxg_dom/", group, "_chr", chr_i, "_", pheno_name, ".png")
      png(pxg_fname)
      print(pxg(cross_g, pheno = cross_data[cross_data$infection == group, pheno_name][[1]], marker = marker_i, geno.map = geno_map))
      dev.off()
    }
  }
  
  # variance decomposition in combined analysis - partial R-squared
  fit0 <- lm(pheno ~ sex + infection, data = dat)
  SSE0 <- summary(fit0)$sigma^2 * summary(fit0)$df[2]
  fit_add <- lm(pheno ~ sex + infection + add, data = dat)
  SSE_add <- summary(fit_add)$sigma^2 * summary(fit_add)$df[2]
  fit_dom <- lm(pheno ~ sex + infection + add + dom, data = dat)
  SSE_dom <- summary(fit_dom)$sigma^2 * summary(fit_dom)$df[2]
  fit_gxt <- lm(pheno ~ sex + infection + add + dom + add:infection + dom:infection, data = dat)
  SSE_gxt <- summary(fit_gxt)$sigma^2 * summary(fit_gxt)$df[2]
  
  # calculate partial R-squared for sequential ANOVA 
  partR2_add <- (SSE0 - SSE_add)/SSE0
  partR2_dom <- (SSE_add - SSE_dom)/SSE_add
  partR2_gxt <- (SSE_dom - SSE_gxt)/SSE_dom
  
  qtl_var_comp[ix,] <- c(chr_i, pos_i, pheno_name, partR2_add, partR2_dom, partR2_gxt)
}
#tmp <- lm_df

# variance decomposition stats (supp note)
mean(as.numeric(qtl_var_comp$partR2_add), na.rm = TRUE) # 5% 
range(as.numeric(qtl_var_comp$partR2_add), na.rm = TRUE) # 0-44%

# ---------------------------effect size heatmap------------------------------ #
# Create color variable: use addBeta if significant, otherwise NA (for gray)
sig_thresh <- 0.1
#lm_df <- tmp

lm_df <- lm_df %>% 
  left_join(condensed_qtl %>% select(chr, phenotype, discovery_infection = infection_detect), by = c("chr", "phenotype")) %>%
  left_join(qtl_map %>% select(chr, phenotype, qtl_name = qtl_id), by = c("chr", "phenotype")) %>%
  mutate(addBeta = as.numeric(addBeta), addP = as.numeric(addP),
         domBeta = as.numeric(domBeta), domP = as.numeric(domP)) %>%
  group_by(chr, phenotype) %>%
  mutate(any_add_sig = any(addP < sig_thresh),
         any_dom_sig = any(domP < sig_thresh),
         is_gxt_only = !grepl("C|1|2", discovery_infection[1]) & grepl("A", discovery_infection[1]),
         is_dom_qtl = !any_add_sig & any_dom_sig) %>%
  ungroup() %>%
  mutate(plot_value = case_when(
    (is_dom_qtl & domP < sig_thresh) ~ domBeta, 
    (!is_dom_qtl & addP < sig_thresh) ~ addBeta,
    TRUE ~ NA), 
    plot_value_scaled = sign(plot_value)*sqrt(abs(plot_value)),
    qtl_id = paste0(chr, "_", phenotype),
    qtl_type = case_when(is_dom_qtl ~ "dom",
                         is_gxt_only ~ "gxt",
                         TRUE ~ "add"),
    chr = as.numeric(chr),
    infection = factor(infection, levels = rev(strat_groups)),
    qtl_num = as.numeric(sub("^[^0-9]+", "", qtl_name))) %>%
  arrange(qtl_num) %>%
  mutate(qtl_id = factor(qtl_id, levels = unique(qtl_id))) %>%
  select(-qtl_num)


# Create effect heatmap WITHOUT labels
effect_hm <- ggplot(lm_df, aes(x = qtl_id, y = infection, fill = plot_value_scaled)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, 
                       oob = scales::squish,
                       na.value = "gray80",
                       name = "Genetic\neffect") +
  scale_x_discrete(drop = FALSE) + 
  scale_y_discrete(labels = c("PBS" = "Control", "SARSCoV" = "SARS-CoV", "SARS2CoV" = "SARS-CoV-2")) + 
  labs(x = NULL, y = "Infection", title = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.margin = margin(l = -5),
        panel.spacing = unit(0,"pt"))

# # Get dominance QTL info with positions
# dom_qtl_info <- lm_df %>%
#   filter(is_dom_qtl) %>%
#   mutate(position = match(qtl_id, levels(lm_df$qtl_id))) %>%
#   select(qtl_id, position) %>%
#   distinct()
# 
# # Add dashed boxes to effect_hm
# for (i in 1:nrow(dom_qtl_info)) {
#   effect_hm <- effect_hm +
#     annotate("rect",
#              xmin = dom_qtl_info$position[i] - 0.5, 
#              xmax = dom_qtl_info$position[i] + 0.5,
#              ymin = 0.5, 
#              ymax = 3.5,  # Adjust if you have different number of infection groups
#              fill = NA,
#              color = "black",
#              linewidth = 0.7)
# }

ggsave("figures/qtl_mapping/effect_heatmap.png", width = 8.5, height = 4)

# -----------------------variance decomposition per QTL------------------------ #
# Reshape data to long format
qtl_var_long <- qtl_var_comp %>%
  left_join(qtl_map, by = c("chr", "phenotype")) %>%
  left_join(pheno_names, by = c("phenotype" = "flow_col_name")) %>%
  select(-flow_display_name) %>%
  mutate(partR2_add = as.numeric(partR2_add),
         partR2_dom = as.numeric(partR2_dom),
         partR2_gxt = as.numeric(partR2_gxt),
         total_partR2 = as.numeric(partR2_add + partR2_dom + partR2_gxt),
         qtl_label = paste0(chr, "_", phenotype)) %>%
  pivot_longer(cols = starts_with("partR2_"),
               names_to = "component",
               values_to = "partial_R2",
               names_prefix = "partR2_") %>%
  mutate(component = factor(component, levels = c("gxt", "dom", "add"), labels = c("GxT", "Dominance", "Additive")),
         qtl_label = factor(qtl_label, levels = levels(lm_df$qtl_id)))
#   filter(!(qtl_id == "Irq16" & phenotype == "G_AlveolarMacs_pctMHCII") & !(qtl_id == "Irq6")) %>%
# qtl_label = paste0(qtl_id, " [", flow_display_name2, "]")


# Create variance heatmap WITHOUT labels
var_hm <- ggplot(qtl_var_long, aes(x = qtl_label, 
                                   y = component,
                                   fill = partial_R2*100)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    limits = c(0, 15),
    oob = scales::squish,
    name = "Partial\nR² (%)"
  ) +
  labs(x = NULL, 
       y = "Genetic Component",
       title = NULL) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.margin = margin(l = -5),
        panel.spacing = unit(0, "pt"))


# get QTL positions for labeling QTL
qtl_positions_multi <- lm_df %>%
  select(qtl_id, qtl_name) %>%
  distinct() %>%
  mutate(position = match(qtl_id, unique(lm_df$qtl_id))) %>%
  group_by(qtl_name) %>%
  mutate(positions = list(position),
         n_positions = n()) %>%
  filter(n_positions > 1) %>%
  ungroup() %>%
  select(qtl_name, positions, n_positions) %>%
  distinct()
# for QTL with one association 
qtl_positions_single <- lm_df %>%
  select(qtl_id, qtl_name) %>%
  distinct() %>%
  mutate(position = match(qtl_id, unique(lm_df$qtl_id))) %>%
  group_by(qtl_name) %>%
  filter(n() == 1) %>%
  ungroup()
# Create SEPARATE label plot
label_plot <- ggplot() +
  scale_x_continuous(limits = c(0.5, nlevels(lm_df$qtl_id) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
        panel.spacing = unit(0, "pt"),
        panel.background = element_blank(),
        plot.background = element_blank())
# Add single QTL labels (in the middle)
y_label_single <- 0.5
for (i in 1:nrow(qtl_positions_single)) {
  label_plot <- label_plot +
    annotate("text", 
             x = qtl_positions_single$position[i], 
             y = y_label_single, 
             label = qtl_positions_single$qtl_name[i], fontface = "italic",
             size = 2.5, angle = 90, hjust = 0.5, vjust = 0.5)
}
# Add DOWNWARD brackets with labels in between
y_bracket_down <- 0.35
y_label_down <- 0.5  # Same as single labels
for (i in 1:nrow(qtl_positions_multi)) {
  qtl_info <- qtl_positions_multi[i,]
  positions <- sort(unlist(qtl_info$positions))
  #shrink_label <- ((length(positions) < 3) & (nchar(qtl_info$qtl_name) > 4))
  shrink_label <- qtl_info$qtl_name %in% c("Irq13", "Irq20", "Irq21", "Irq22") # these are crowded
  label_size <- ifelse(shrink_label, 2, 2.5)
  chunk_min <- min(positions)
  chunk_max <- max(positions)
  label_plot <- label_plot +
    # Horizontal line
    annotate("segment", x = chunk_min - 0.3, xend = chunk_max + 0.3, 
             y = y_bracket_down, yend = y_bracket_down, linewidth = 0.5) +
    # Left vertical (pointing down)
    annotate("segment", x = chunk_min - 0.3, xend = chunk_min - 0.3, 
             y = y_bracket_down, yend = y_bracket_down - 0.15, linewidth = 0.5) +
    # Right vertical (pointing down)
    annotate("segment", x = chunk_max + 0.3, xend = chunk_max + 0.3, 
             y = y_bracket_down, yend = y_bracket_down - 0.15, linewidth = 0.5) +
    # Label in the middle
    annotate("text", x = (chunk_min + chunk_max) / 2, 
             y = y_label_down, 
             label = qtl_info$qtl_name, fontface = "italic",
             size = label_size, hjust = 0.5, vjust = 0.5)
}
# Add UPWARD brackets with labels in between
y_bracket_up <- 0.65
y_label_up <- 0.5  # Same as single labels
for (i in 1:nrow(qtl_positions_multi)) {
  qtl_info <- qtl_positions_multi[i,]
  positions <- sort(unlist(qtl_info$positions))
  chunk_min <- min(positions)
  chunk_max <- max(positions)
  label_plot <- label_plot +
    # Horizontal line
    annotate("segment", x = chunk_min - 0.3, xend = chunk_max + 0.3, 
             y = y_bracket_up, yend = y_bracket_up, linewidth = 0.5) +
    # Left vertical (pointing up)
    annotate("segment", x = chunk_min - 0.3, xend = chunk_min - 0.3, 
             y = y_bracket_up, yend = y_bracket_up + 0.15, linewidth = 0.5) +
    # Right vertical (pointing up)
    annotate("segment", x = chunk_max + 0.3, xend = chunk_max + 0.3, 
             y = y_bracket_up, yend = y_bracket_up + 0.15, linewidth = 0.5)
  # No label here - it's already added by the downward bracket loop
}

fig3 <- plot_grid(qtl_hm, NULL,
                  plot_grid(us_grob,
                            plot_grid(var_hm, NULL, label_plot, NULL, effect_hm, nrow = 5, 
                                      rel_heights = c(1,-0.04,0.2,-0.04,1), 
                                      labels = c("C", "", "D"), align = "v", axis = "lr"),
                            nrow = 2, rel_heights = c(1,2), labels = c("B", ""), hjust = -5),
                  ncol = 3, rel_widths = c(0.9, -0.08, 1), labels = c("A", "", ""))
fig3
ggsave("figures/Figure3.png", width = 12.5, height = 8, bg = "white")


# Create stacked bar plot
qtl_vc_barplot <- ggplot(qtl_var_long, aes(x = reorder(qtl_label, total_partR2), y = partial_R2 * 100, fill = component)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Additive" = "#2E86AB", 
                               "Dominance" = "#A23B72", 
                               "GxT" = "#F18F01")) +
  coord_flip() + 
  labs(x = "QTL", 
       y = "Partial R² (%)",
       fill = "Component",
       title = "Variance Explained by Genetic Components at Each QTL") +
  theme_minimal() +
  theme(legend.position = c(0.6,0.3),
        legend.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(reverse = TRUE))
qtl_vc_barplot
ggsave("figures/qtl_mapping/partialR2_by_component.png", height = 8, width = 6)

# histogram
qtl_vc_hist <- ggplot(qtl_var_long, aes(x = partial_R2*100, fill = component)) +
  geom_histogram() + # position = "identity", alpha = 0.6
  scale_fill_manual(values = c("Additive" = "#2E86AB", 
                               "Dominance" = "#A23B72", 
                               "GxT" = "#F18F01")) +
  labs(x = "Partial R² (%)", fill = "Component",
       title = "Variance explained by genetic components") +
  theme_minimal() + 
  guides(fill = guide_legend(reverse = TRUE))
qtl_vc_hist 


qtl_heat <- all_qtl %>%
  filter(infection %in% strat_groups) %>%
  mutate(qtl_id = paste0("chr", chr, ":", phenotype)) %>%
  select(qtl_id, infection, beta) %>%
  complete(qtl_id, infection = strat_groups) %>%
  mutate(beta = as.numeric(beta))
# fill in betas for nominally significant effects 
ggplot(qtl_heat, aes(x = infection, y = qtl_id, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "red",
    mid = "grey95",
    high = "blue",
    midpoint = 0,
    na.value = "grey80",
    name = "beta"
  ) +
  labs(x = NULL, y = NULL, title = "QTL effects across control and infected groups") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9)
  )




