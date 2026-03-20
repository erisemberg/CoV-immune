library(readr)
library(ggplot2)
source("/Users/ellenrisemberg/Documents/ValdarFerris/scripts/qtl_functions.R")

# ---------------------------------Load data---------------------------------- #
cross_data <- read_csv("derived_data/cross_data.csv", show_col_types = FALSE)

# ----------------------------Figure 1: weight loss--------------------------- #
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

# ----------------------------Figure 1: viral titer--------------------------- #
inf_pal <- brewer.pal(n = 3, name = "Set2")[2:3]

titer_df <- data.frame(Infection = cross_data$infection,
                       Titer = cross_data$Titer)
p4 <- ggplot(titer_df, aes(x = Titer)) + 
  geom_histogram(aes(fill = Infection), position = "dodge") + 
  scale_fill_manual(values = inf_pal, labels = c("SARS-CoV", "SARS-CoV-2")) +
  labs(x = "Viral titer", y = "Number of F2 mice") + 
  theme_minimal() + 
  theme(legend.position = "none")

# --------------------------------Figure 1: HS-------------------------------- #
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
