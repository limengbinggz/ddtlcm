#######################################################
#' Generate Figures S4.3 in the supplement using
#' summaries of sigma estimation.
#######################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(see)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)
options(warn = -1)

## specify parameters -----------------------------
ntree <- 4
ntheta <- 20
nY <- 5
K <- 6
J_g <- c(rep(10, 5), 15, 15)
J <- sum(J_g)
G <- length(J_g)
Ns <- c(100, 200, 400) 
# variance of logit response probabilities of items in each group
Sigma_by_group <- c(rep(0.6**2, 5), rep(2**2, 2))

dat <- fread("simulation_synthetic/data/1_simu_syn_sigma.csv")
dat[, seed := (seed_theta-1)*nY + seed_Y]
dat_summary <- dat[ , list(mean=mean(est), sd = sd(est), coverage = mean(cover)), 
                                 by = .(N, tree_idx, var_name)]
dat_summary[, true_value := rep(Sigma_by_group, length(Ns)*ntree)]
dat_summary[, `Coverage` := paste0(format(round(coverage, digits=2), nsmall = 2))]
dat_summary$N_label <- factor(dat_summary$N, levels = Ns, 
                        labels = paste0("N = ", Ns))
dat_summary$var_name <- gsub("Sigma", "sigma", dat_summary$var_name)
dat_summary$var_label <- paste0("$\\", dat_summary$var_name, "$")

dat_plot <- dat_summary[, .(tree_idx, N_label, var_name, Coverage)]
dat_plot$Coverage <- as.numeric(dat_plot$Coverage)

## Figure S4.3(a): coverage probability
ggplot(dat_plot, aes(y = `Coverage`, x = var_name)) +
  geom_point(aes(shape = N_label), size = 3) +
  geom_line(aes(group = N_label, linetype = N_label)) +
  scale_shape_manual(values=c(1, 16, 4)) +
  facet_grid(. ~ tree_idx) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Variable", y = "Coverage Probability", shape = "Sample Size", linetype = "Sample Size") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13)) 
if(!dir.exists("figures")){
  dir.create("figures")
}
ggsave("figures/sim_syn_sigma_coverage.pdf", height = 3.5, width = 10)



## Figure S4.3(b): distribution of posterior samples
dat_Sigma_by_group <- data.table(var_name = paste0("Sigma_", 1:G), true_value = Sigma_by_group)
dat2 <- merge(x = dat, y = dat_Sigma_by_group, by = "var_name", all.x = TRUE)
dat2$N_label <- paste0("N = ", dat2$N)
dat2$Type <- "True Value"
ggplot(dat2, aes(y = est, x = var_name)) +
  geom_violinhalf(aes(fill = factor(dat2$N_label))) +
  geom_point(aes(y = true_value, x = var_name, shape = Type), size = 3) +
  facet_wrap( ~ tree_idx, ncol = 2) +
  scale_shape_manual(values=c(4)) +
  scale_y_continuous(limits = c(0.3, 4.5)) +
  labs(x = "Variable", y = "Posterior Means", fill = "Sample Size", shape = "") +
  theme_bw() +
  scale_fill_material_d() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13)) 
if(!dir.exists("figures")){
  dir.create("figures")
}
ggsave("figures/sim_syn_sigma_distribution.pdf", height = 6, width = 8)



