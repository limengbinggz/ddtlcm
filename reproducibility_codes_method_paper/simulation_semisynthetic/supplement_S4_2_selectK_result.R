require(data.table)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

nseed_theta <- 20
nseed_Y <- 5
K_candidates <- c(4,5,6,7,8)
Ns <- c(400, 800)
nfold = 5

dat <- fread(file = "simulation_semisynthetic/data/2_simu_semi_selectK_summary.csv")
# compute average predictive loglikehood across folds
dat[, ave_predictive_llk := mean(predictive_llk), by=c("N", "seed_theta", "seed_Y", "K_candidate")]

# choose the K with the highest average log llk
dat_select <- dat[, .SD[which.max(ave_predictive_llk)], by=c("N", "seed_theta", "seed_Y")]

## uncomment if you want to save the output to a txt file
# if (!dir.exists("simulation_semisynthetic/output")){
#   dir.create("simulation_semisynthetic/output")
# }
# sink("simulation_semisynthetic/output/2_selectK_freq_table.txt")
cat("Frequency table of selecting K using predictive loglikelihood\n")
for (N in Ns) {
  cat("N =", N, ":\n")
  freq_table <- table(dat_select$K_candidate[dat_select$N == N])
  names(freq_table) <- paste("K =", names(freq_table))
  print(freq_table)
}
# sink()
