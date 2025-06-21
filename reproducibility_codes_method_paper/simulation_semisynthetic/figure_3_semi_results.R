#######################################################
#' Generate Figure 3 in the manuscript using summaries of
#' semi-synthetic data simulation.
#######################################################

library(dplyr)
library(data.table)
library(ggtree)
library(ggplot2)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

## read in data
dat <- fread("simulation_semisynthetic/data/1_simu_semi_summary.csv")

# The color-blinded friendly palette with grey:
cbPalette <- c("#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
for (x in c("rmse", "tree_frobenius", "ari")) {
  dat[[x]] <- as.numeric(dat[[x]])
}
# keep only single linkage method results for BayesLCM
dat <- dat[(method_hclust == "single" & label == "BayesLCM (homogeneous var) + HC") |
             (method_hclust == "single" & label == "BayesLCM (heterogeneous var) + HC") | 
             (label %in% c("DDT-LCM", "DDT-LCM (true tree)", 
                           "DDT-LCM (misspecified tree)", "DDT-LCM (homogeneous var)")),] 

# a temporary data table to relabel the methods
method_labels <- data.table(
  label = c(
    "DDT-LCM", 
    "DDT-LCM (true tree)",
    "DDT-LCM (misspecified tree)",
    "DDT-LCM (homogeneous var)",
    "BayesLCM (heterogeneous var) + HC",
    "BayesLCM (homogeneous var) + HC"),
  labels_full = c(
    "(i) DDT-LCM", 
    "(ii) DDT-LCM \n(true tree)", 
    "(iii) DDT-LCM \n(misspecified tree)", 
    "(iv) DDT-LCM \n(homogeneous var)",
    "(v) BayesLCM \n(heterogeneous var) + HC",
    "(vi) BayesLCM \n(homogeneous var) + HC"
  ),
  labels_short = c("(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)")
)
dat <- merge(dat, method_labels, by = "label")

dat$labels_full <- factor(dat$labels_full, levels = c(
  "(i) DDT-LCM", 
  "(ii) DDT-LCM \n(true tree)", 
  "(iii) DDT-LCM \n(misspecified tree)", 
  "(iv) DDT-LCM \n(homogeneous var)",
  "(v) BayesLCM \n(heterogeneous var) + HC",
  "(vi) BayesLCM \n(homogeneous var) + HC"
))
dat$N <- paste0("N = ", dat$N)
dat$N <- factor(dat$N, levels = paste0("N = ", c(400, 800)))
dat$labels_short <- factor(dat$labels_short,
                           levels = c("(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)"))

### put all plots together -------------------------
dat$log_rmse <- log(dat$rmse)
dat_long <- melt(dat[, .(N, labels_full, labels_short, log_rmse, ari)], 
                 id.vars = c("N", "labels_full", "labels_short"))
dat2_long <- melt(dat[, .(N, labels_full, labels_short, tree_frobenius)],
                  id.vars = c("N", "labels_full", "labels_short"))
dat_all <- rbind(dat_long, dat2_long)
dat_all$variable <- factor(dat_all$variable,
                                levels = c("tree_frobenius", "log_rmse", "ari"), 
                                labels = c("Frobenius Norm of \nTree Covariance Matrix",
                                           "Log RMSE of Item \nResponse Probabilities",
                                           "Class Membership ARI")) 
p <- ggplot(data = dat_all, aes(x = labels_short, y = value, fill = labels_full)) +
  geom_boxplot() +
  facet_grid(variable~N, scales = "free_y", switch="y") +
  scale_fill_manual(values=cbPalette, drop = FALSE) +
  labs(y = "", x = "",
       fill = "Method") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.position = "right",
    legend.spacing.y = unit(0.3, 'cm')) +
  guides(fill = guide_legend(ncol=1, byrow=TRUE))
print(p)
if(!dir.exists("figures")){
  dir.create("figures")
}
ggsave("figures/combined_semi_boxplot.pdf", height = 5, width = 8)


