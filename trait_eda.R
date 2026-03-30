library(qtl)
library(car)
library(readr)
library(readxl)
library(AGHmatrix)
library(lme4)
library(lme4qtl)
library(ordinal) # proportional odds model 
library(caret)
library(glmnet)
library(randomForest)
library(factoextra)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(r2glmm)
library(parallel)
library(doParallel)
library(foreach)
source("code-dependencies/cov_qtl_functions.R")

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

# load phenotypes, genotypes 
cross_data <- read_csv("derived_data/cross_data.csv")
cross_data$infection <- as.factor(cross_data$infection)
geno <- cross_data[,c(which(colnames(cross_data) == "gUNC145909"):which(colnames(cross_data) == "mJAX00187167"))]
cross_dataI <- read_csv("derived_data/cross_data_flow_imp.csv", show_col_types = FALSE) # processed and imputed data 

# figures 
ensure_directory("figures")
ensure_directory("figures/supplemental")
inf_pal <- list("PBS" = "#66C2A5", "SARSCoV" = "#FC8D62", "SARS2CoV" = "#8DA0CB")
grouplabels = list("PBS" = "Control", "SARSCoV" = "SARS-CoV", "SARS2CoV" = "SARS-CoV-2")

# ---------------------------------------------------------------------------- #
# ------------------------------------EDA------------------------------------- #
# ---------------------------------------------------------------------------- #
# sample sizes 
message(paste(sum(!is.na(cross_data$weight_aac)), "total mice"))
table(cross_data$infection[!is.na(cross_data$weight_aac)])

# does titer vary by infection group? 
message(paste("Titer range:", paste(round(range(cross_data$Titer, na.rm = TRUE), 2), collapse = " to ")))
t.test(Titer ~ infection, cross_data, var.equal = FALSE)

# does titer predict weight loss? only considers infected mice (don't have titer data for ctrl mice)
summary(lm(weight_aac ~ infection + sex + Titer + infection:sex, cross_data))
summary(lm(HS ~ infection + sex + Titer + infection:sex, cross_data))

# ----------------------------variance heterogeneity-------------------------- #
VHdf <- data.frame(phenotype = character(),
                   infectionVH = numeric(),
                   sexVH = numeric(),
                   infSexVH = numeric())
inf_var_g <- 0
for (i in 1:q){
  inf_vh_p <- leveneTest(flow[,i], as.factor(cross_dataI$infection))$`Pr(>F)`[1]
  sex_vh_p <- leveneTest(flow[,i], as.factor(cross_dataI$sex))$`Pr(>F)`[1]
  int_vh_p <- leveneTest(flow[,i], interaction(cross_dataI$infection, cross_dataI$sex))$`Pr(>F)`[1]
  VHdf[i,] <- c(flow_cols[i], inf_vh_p, sex_vh_p, int_vh_p)
  
  var0 <- var(flow[cross_dataI$infection == "PBS", i], na.rm = TRUE)
  var12 <- var(flow[cross_dataI$infection %in% c("SARSCoV", "SARS2CoV"), i], na.rm = TRUE)
  if (inf_vh_p <= 0.05 & var12 > var0){
    inf_var_g <- inf_var_g + 1
  }
}
message(paste("variance heterogeneity is present for", sum(VHdf$infectionVH < 0.05) , "phenotypes."))
message(paste("variance is greater in infected mice for", inf_var_g, "phenotypes."))


# -------------------------classification with glmnet------------------------- #
# split into train/test
flow_wInf <- cross_dataI[,c(paste0(flow_cols, "_imp"), "infection")]

train_ix <- createDataPartition(flow_wInf$infection, p = 0.8, list = FALSE)
train_data <- flow_wInf[train_ix,]
test_data <- flow_wInf[-train_ix,]

X_train <- model.matrix(infection ~ . - 1, data = train_data)
y_train <- train_data$infection
cv_fit <- cv.glmnet(X_train, y_train, family = "multinomial", type.measure = "class")

X_test <- model.matrix(infection ~ . -1, data = test_data)
glm_pred <- predict(cv_fit, newx = X_test, s = "lambda.min", type = "class")
# confusionMatrix(as.factor(glm_pred), as.factor(test_data$infection))

# extract all class-specific coefficient matrices
coefs <- coef(cv_fit, s = "lambda.min")

set.seed(123)
ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE, 
                     summaryFunction = multiClassSummary, allowParallel = FALSE)
tune_grid <- expand.grid(alpha = seq(0, 1, by = 0.1), lambda = seq(0.0001, 1, length = 10))
glmnet_cv <- train(infection ~ ., data = flow_wInf, method = "glmnet", trControl = ctrl, tuneGrid = tune_grid)
best <- glmnet_cv$bestTune
glmnet_cv$results[glmnet_cv$results$alpha == best$alpha & glmnet_cv$results$lambda == best$lambda,
                  c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]
var_imp <- varImp(glmnet_cv, scale = FALSE)
imp_df <- var_imp$importance
imp_df$overall <- rowSums(imp_df)
imp_df$phenotype <- rownames(imp_df)
imp_df <- imp_df[order(imp_df$overall, decreasing = TRUE),]

ensure_directory("results/var_selection")
ensure_directory("results/var_selection/frequentist")
write_csv(imp_df, "results/var_selection/frequentist/glmnet_var_importance.csv")


# ---------------------classification with random forest---------------------- #
rf_model <- randomForest(as.factor(infection) ~ ., data = train_data, ntree = 500, importance = TRUE)
rf_pred <- predict(rf_model, newdata = test_data)
#confusionMatrix(rf_pred, as.factor(test_data$infection))
#varImpPlot(rf_model)

### PREDICTION W RANDOM FOREST 
flow_wAAC <- cross_dataI[,c(paste0(flow_cols, "_imp"), "weight_aac")]
train_data2 <- flow_wAAC[train_ix,]
test_data2 <- flow_wAAC[-train_ix,]
rf_model2 <- randomForest(weight_aac ~ ., data = train_data2, ntree = 500, importance = TRUE)
rf_pred2 <- predict(rf_model2, newdata = test_data2)
#varImpPlot(rf_model2)

# -----------------------------T cell composition----------------------------- #
df <- data.frame(mouseID = cross$pheno$mouse_ID,
                 infection = cross$pheno$infection, 
                 L_CytoTCells_pctNaive = cross$pheno$L_CytoTCells_pctNaive,
                 L_CytoTCells_pctEffector = cross$pheno$L_CytoTCells_pctEffector,
                 L_CytoTCells_pctCentMem = cross$pheno$L_CytoTCells_pctCentMem, 
                 L_HelpTCells_pctNaive = cross$pheno$L_HelpTCells_pctNaive,
                 L_HelpTCells_pctEffector = cross$pheno$L_HelpTCells_pctEffector,
                 L_HelpTCells_pctCentMem = cross$pheno$L_HelpTCells_pctCentMem, 
                 L_DNTCells_pctNaive = cross$pheno$L_DNTCells_pctNaive,
                 L_DNTCells_pctEffector = cross$pheno$L_DNTCells_pctEffector,
                 L_DNTCells_pctCentMem = cross$pheno$L_DNTCells_pctCentMem)
df_long <- df %>%
  filter(!if_all(starts_with("L_"), ~ is.na(.x) | .x == 0)) %>%
  pivot_longer(cols = starts_with("L"), names_to = "phenotype", values_to = "proportion") %>%
  mutate(subset = case_when(str_detect(phenotype, "CytoTCells") ~ "Cytotoxic T cells",
                            str_detect(phenotype, "HelpTCells") ~ "Helper T cells",
                            str_detect(phenotype, "DNTCells")   ~ "DNT cells"),
         state = case_when(str_detect(phenotype, "Naive")     ~ "Naive",
                           str_detect(phenotype, "Effector")  ~ "Effector",
                           str_detect(phenotype, "CentMem")   ~ "Central memory")) %>%
  group_by(mouseID, infection, subset) %>%
  mutate(prop_norm = proportion / sum(proportion, na.rm = TRUE)) %>%
  ungroup()

# confirm they add up to ~100 per mouse 
df_long %>%
  group_by(mouseID, infection, subset) %>%
  summarise(total = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  filter(total < 90 | total > 110)
# nrow = 3036 = 1012*3 (before removing bad mice)

df_plot <- df_long %>% 
  group_by(infection, subset, state) %>%
  summarise(mean_prop = mean(prop_norm, na.rm = TRUE), .groups = "drop")

ggplot(df_plot, aes(x = infection, y = mean_prop, fill = state)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ subset, ncol = 1) +
  scale_x_discrete(limits = rev(names(inf_pal)), labels = grouplabels[names(inf_pal)]) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Normalized proportion", fill = "Differentiation state") +
  theme_minimal()
ggsave("figures/supplemental/Tcell_composition.png")

# ---------------------------------------------------------------------------- #
# ---------------------------------Supp Fig 1--------------------------------- #
# ---------------------------------------------------------------------------- #

# ---------------------------------weight loss-------------------------------- #
# cross_data[cross_data$mouse_ID == "CR_RB05_F_0094", c("d0", "d1", "d2", "d3", "d4")]
# remove control outliers 
cross_data$pd1[which(cross_data$mouse_ID=="CR_RB05_F_0815")] <- NA
cross_data$pd2[which(cross_data$mouse_ID=="CR_RB05_F_0094")] <- NA
# remove SARS2 outliers 
cross_data$pd1[which(cross_data$mouse_ID=="CR_RB05_F_1136")] <- NA
cross_data$pd2[which(cross_data$mouse_ID=="CR_RB05_M_1105")] <- NA
cross_data$pd2[which(cross_data$mouse_ID=="CR_RB05_M_0889")] <- NA
cross_data$pd4[which(cross_data$mouse_ID=="CR_RB05_M_0936")] <- NA
cross_data$pd4[which(cross_data$mouse_ID=="CR_RB05_M_0946")] <- NA
cross_data[which(cross_data$mouse_ID=="CR_RB05_F_0977"), c("pd3", "pd4")] <- NA
cross_data$pd4[which(cross_data$mouse_ID=="CR_RB05_M_1229")] <- NA

p1 <- cov_trajectory_plot(dat = cross_data[cross_data$infection == "PBS",], 
                          phenos = c("pd0", "pd1", "pd2", "pd3", "pd4"), 
                          incl.parents = FALSE, title = "Control", 
                          ylab = "% of starting weight", xlab = "", ylim = c(70,115)) +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- cov_trajectory_plot(dat = cross_data[cross_data$infection == "SARSCoV",], 
                          phenos = c("pd0", "pd1", "pd2", "pd3", "pd4"), 
                          incl.parents = TRUE, title = "SARS-CoV", parent.lty = 1,
                          ylab = "% of starting weight", ylim = c(70,115)) +
  theme_minimal() +
  theme(legend.position = "none", axis.title.y = element_blank()) 

p3 <- cov_trajectory_plot(dat = cross_data[cross_data$infection == "SARS2CoV",], 
                          phenos = c("pd0", "pd1", "pd2", "pd3", "pd4"), 
                          incl.parents = FALSE, title = "SARS-CoV-2", 
                          ylab = "% of starting weight", xlab = "", ylim = c(70,115)) +
  theme_minimal() + 
  theme(axis.title.y = element_blank()) 

legend_data <- data.frame(Strain = c("F2", "CC006", "CC044"),
                          color = c("gray63", "#E41A1C", "#377EB8"),
                          size = c(0.2, 1, 1))
color_values <- setNames(legend_data$color, legend_data$Strain)
size_values <- setNames(legend_data$size, legend_data$Strain)
legend_plot <- ggplot(legend_data, aes(x = 1, y = 1, color = Strain, size = Strain)) +
  geom_line() +
  scale_color_manual(name = "Strain", values = color_values) +
  scale_size_manual(name = "Strain",  values = size_values) +
  theme_void()

# --------------------------------viral titer--------------------------------- #
inf_pal <- brewer.pal(n = 3, name = "Set2")[2:3]

titer_df <- data.frame(Infection = cross_data$infection,
                       Titer = cross_data$Titer)
p4 <- ggplot(titer_df, aes(x = Titer)) + 
  geom_histogram(aes(fill = Infection), position = "dodge") + 
  scale_fill_manual(values = inf_pal, labels = c("SARS-CoV", "SARS-CoV-2")) +
  labs(x = "Viral titer", y = "Number of F2 mice") + 
  theme_minimal() + 
  theme(legend.position = "none")

# -------------------------------------HS------------------------------------- #
hs_df <- data.frame(Infection = cross_data$infection,
                    HS = cross_data$HS)
p5 <- ggplot(hs_df, aes(x = HS)) + 
  geom_histogram(aes(fill = Infection), position = "dodge", binwidth = 0.5) + 
  scale_fill_manual(values = inf_pal, labels = c("SARS-CoV", "SARS-CoV-2")) +
  labs(x = "Congestion score", y = "Number of F2 mice") + 
  theme_minimal()

# ---------------------------------Figure 1----------------------------------- #
f1_top <- plot_grid(p1, p2, p3, get_legend(legend_plot),
                    ncol = 4, rel_widths = c(0.3, 0.3, 0.3, 0.1), labels = c("A", "B", "C", ""))
f1_bottom <- plot_grid(p4, NULL, p5, ncol = 3, labels = c("D","", "E"), rel_widths = c(1,0.1,1.25))
fig1 <- plot_grid(f1_top, NULL, f1_bottom, nrow = 3, rel_heights = c(0.48, 0.04, 0.48))
ggsave("figures/Figure1.png", bg = "white", width = 11, height = 7.33)



# ---------------------------------------------------------------------------- #
# ---------------------------------Supp Fig 1--------------------------------- #
# ---------------------------------------------------------------------------- #

# ----------------------------illustrate AAC stat----------------------------- #
ex_mouse_id <- "CR_RB05_F_1004"
wt_phenos = c('pd0', 'pd1', 'pd2', 'pd3','pd4')  
dat <- cross_data %>%
  select(mouse_ID, all_of(wt_phenos)) %>%
  setNames(c("mouse_ID", as.character(0:4))) %>%
  pivot_longer(cols = -mouse_ID, names_to = "dpi", values_to = "weight") %>%
  mutate(dpi = as.numeric(dpi),
         weight = as.numeric(weight))

aac_plot <- dat %>%
  filter(mouse_ID == ex_mouse_id) %>%
  arrange(dpi) %>%
  mutate(baseline = first(weight),
         ymin = pmin(weight, baseline, na.rm = TRUE),
         ymax = pmax(weight, baseline, na.rm = TRUE)) %>%
  ggplot(aes(x = dpi, y = weight)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#6BD7AF", alpha = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_hline(aes(yintercept = baseline, linetype = "Baseline weight"), linewidth = 0.6, color = "black") +
  scale_linetype_manual(values = c("Baseline weight" = "22"), guide = "none") +
  scale_x_continuous(breaks = 0:4) + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(x = "Days post-infection", y = "% of baseline body weight") +
  theme_minimal()

# -----------------------------correlation plots------------------------------ #
df <- cross_data[,c("sex", "infection", "HS", "Titer", "weight_aac")]
df <- df[!is.na(df$Titer),]

# weight-titer correlation  
summary(lm(weight_aac ~ sex + Titer, data = df[df$infection == "SARSCoV",])) # p = 4.8e-08
summary(lm(weight_aac ~ sex + Titer, data = df[df$infection == "SARS2CoV",])) # p = 5.9e-07
sig_df <- sig_df <- data.frame(group = strat_groups[2:3],
                               label = c("p = 4.8e-08", "p = 5.9e-07"),
                               vjust = c(1.2, 2.6))
wt_titer <- ggplot(data = df, aes(x = Titer, y = weight_aac, color = infection)) +
  geom_point() +
  scale_color_manual(values = inf_pal) + 
  geom_smooth(method = "lm") +
  labs(x = "Viral titer", y = "Weight loss (AAC)") + 
  geom_text(data = sig_df, mapping = aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust),
            hjust = 1.1, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  theme_minimal()

# HS-titer correlation 
df$HS <- factor(df$HS, ordered = TRUE)
fit1 <- clm(HS ~ sex + Titer, data = df[df$infection == "SARSCoV",], link = "logit")
summary(fit1) # p = 0.05
fit2 <- clm(HS ~ sex + Titer, data = df[df$infection == "SARS2CoV",], link = "logit")
summary(fit2) # p = 3.8e-08
sig_df <- sig_df <- data.frame(group = strat_groups[2:3],
                               label = c("p = 0.05", "p = 3.8e-08"),
                               vjust = c(1.2, 2.6))
df_plot <- df[!is.na(df$HS),]
hs_titer <- ggplot(df_plot, aes(x = Titer, y = factor(HS, levels = seq(0, 4, by = 0.5), ordered = TRUE))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = infection), width = 0, height = 0.15, size = 1) +
  scale_color_manual(values = inf_pal, name = "Infection", breaks = names(inf_pal), labels = grouplabels[names(inf_pal)]) +
  geom_text(data = sig_df, 
            aes(x = Inf, y = Inf, label = label, color = group, vjust = vjust), 
            hjust = 1, size = 4, inherit.aes = FALSE, show.legend = FALSE) +
  labs(x = "Viral titer", y = "Lung congestion score") +
  theme_minimal()  

# -------------------------categorize immune traits--------------------------- #
pheno_annotated <- pheno_names %>%
  filter(flow_col_name != "Titer") %>%
  mutate(panel = case_when(
    str_detect(flow_col_name, "L_") ~ "Lymphocytes",
    str_detect(flow_col_name, "M_") ~ "Myeloid cells",
    str_detect(flow_col_name, "G_") ~ "Granulocytes"
  ), lineage = case_when(
    str_detect(flow_col_name, "BCells") ~ "B cells",
    str_detect(flow_col_name, "CytoTCells") ~ "Cytotoxic T cells",
    str_detect(flow_col_name, "HelpTCells") ~ "Helper T cells",
    str_detect(flow_col_name, "gdTCR") ~ "Gamma-delta T cells",
    str_detect(flow_col_name, "DNTCells") ~ "T cells (other)",
    str_detect(flow_col_name, "NKTCells") ~ "Natural killer T cells",
    str_detect(flow_col_name, "NKCells") ~ "Natural killer cells",
    str_detect(flow_col_name, "DC") ~ "Dendritic cells",
    str_detect(flow_col_name, "Mac") ~ "Monocytes/macrophages",
    str_detect(flow_col_name, "Mast") ~ "Mast cells",
    str_detect(flow_col_name, "Basophils") ~ "Basophils",
    str_detect(flow_col_name, "Eosinophils") ~ "Eosinophils",
    str_detect(flow_col_name, "Neutrophils") ~ "Neutrophils"
  ), subset = case_when(
    str_detect(flow_col_name, "Alv") ~ "Alveolar macrophages",
    str_detect(flow_col_name, "Int") ~ "Interstitial macrophages",
    str_detect(flow_col_name, "pDC") ~ "Plasmacytoid DCs",
    str_detect(flow_col_name, "Ly6Cp") ~ "Inflammatory (Ly6C+) macrophages",
    str_detect(flow_col_name, "Ly6Cn") ~ "Resident (Ly6C-) macrophages",
    str_detect(flow_col_name, "CD11Bp") ~ "CD11b+ DCs",
    str_detect(flow_col_name, "CD103p") ~ "CD103+ DCs"
  ), phenotype_class = case_when(
    str_detect(flow_col_name, "ofLive|ofCD45|ofCD3|ofDNT") ~ "Cell type frequency",
    str_detect(flow_col_name, "CCR7|CCR2") ~ "Migration",
    str_detect(flow_col_name, "CD25|CD69|PD1|CD40|CD80|CD86|MHCII") ~ "Activation",
    str_detect(flow_col_name, "Effector|Naive|CentMem") ~ "Differentiation",
    str_detect(flow_col_name, "FceR1|NK11") ~ "Other"
  )) 

pheno_desc <- ggplot(pheno_annotated, aes(y = fct_rev(fct_infreq(phenotype_class)))) +
  geom_bar(width = 0.75) +
  labs(x = "Number of phenotypes",
       y = "Phenotype category") + 
  theme_minimal() 

# -------------------------------heritability--------------------------------- #
varcomp_data <- read_csv("results/var_comp_res.csv", show_col_types = FALSE)
vc_long <- varcomp_data %>%
  select(pheno, sex_vp, inf_vp, h2, batch_vp) %>%
  mutate(resid = 1 - sex_vp - inf_vp - h2 - batch_vp) %>%
  pivot_longer(cols = -pheno, names_to = "component", values_to = "prop_var") %>%
  mutate(component = factor(component, levels = c("batch_vp", "inf_vp", "h2", "sex_vp", "resid")))
vc_plot <- ggplot(vc_long, aes(x = component, y = prop_var)) +
  geom_boxplot(fill = "gray") + 
  scale_x_discrete(labels = c("batch_vp" = "Batch", 
                              "inf_vp" = "Infection", 
                              "h2" = expression(h[SNP]^2), 
                              "sex_vp" = "Sex", 
                              "resid" = "Residual")) +
  labs(x = "Component", y = "Proportion of variation explained") + 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))

# -------------------------------Supp. Fig. 1--------------------------------- #
leg <- get_legend(hs_titer)
hs_titer <- hs_titer + theme(legend.position = "none")
wt_titer <- wt_titer + theme(legend.position = "none")

plot_grid(
  plot_grid(aac_plot, wt_titer, leg, hs_titer, ncol = 4, rel_widths = c(1.05,1,0.3,0.95), labels = c("A", "B", "C")),
  plot_grid(pheno_desc, NULL, vc_plot, NULL, ncol = 4, rel_widths = c(0.6, 0.05, 0.85, 0.45), labels = c("D", "E", "", "")), 
  nrow = 2
)
ggsave("figures/supplemental/SuppFig1.png", width = 11, height = 6, bg = "white")
  
  
