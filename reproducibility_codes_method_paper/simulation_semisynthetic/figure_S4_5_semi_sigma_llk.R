#######################################################
#' Generate Figures S4.5 in the supplement using
#' summaries of sigma estimation.
#######################################################

library(data.table)
library(ggplot2)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)


dat_ddtlcm <- fread(file = "simulation_semisynthetic/data/2_simu_semi_selectK_summary.csv")
dat_ddtlcm <- dat_ddtlcm[K_candidate == 6, ]
colnames(dat_ddtlcm)[4] <- "K"
dat_ddtlcm[, mean_llk := mean(predictive_llk), by = .(N, seed_theta, seed_Y)]
dat_ddtlcm <- unique(dat_ddtlcm[, .(N, seed_theta, seed_Y, mean_llk)])

dat_ddtlcm_homovar <- fread(file = "simulation_semisynthetic/data/3_compare_homovar_llk_summary.csv")
dat_ddtlcm_homovar[, mean_llk_homovar := mean(predictive_llk), by = .(N, seed_theta, seed_Y)]
dat_ddtlcm_homovar <- unique(dat_ddtlcm_homovar[, .(N, seed_theta, seed_Y, mean_llk_homovar)])

dat <- merge(x = dat_ddtlcm, y = dat_ddtlcm_homovar,
             by = c("N", "seed_theta", "seed_Y"))
dat[, diff_llk := mean_llk - mean_llk_homovar]

## plot
p <- ggplot(dat, aes(x = factor(N), y = diff_llk)) +
  geom_boxplot() +
  labs(y = "Difference in Predictive Log Likelihood between \nDDT-LCM and DDT-LCM (homogeneous var)", 
       x = "Sample Size N") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12))
print(p)
if(!dir.exists("figures")){
  dir.create("figures")
}
ggsave("figures/sim_semisyn_sigma_llk_diff.pdf", height = 4, width = 4)










