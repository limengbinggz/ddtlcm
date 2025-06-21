library(parallel)
library(data.table)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

nseed_theta <- 20
nseed_Y <- 5
Ns <- c(400, 800)
K_candidates <- c(4,5,6,7,8)
nfold = 5
J <- 78 # total number of items

### summarize cross validation results ---------------------------
colnames_dat <- c("N", "seed_theta", "seed_Y", "K_candidate", "fold", 
                  "predictive_llk", "predictive_llk_dist_mean")
run_one <- function(seed_theta){
  cat("seed_theta =", seed_theta, "\n")
  dat <- c()
  for (N in Ns) {
    cat("N =", Ns, "\n")
    for (seed_Y in 1:nseed_Y) {
      cat("seed_Y =", seed_Y, "\n")
      for (K_candidate in K_candidates) {
        cat("K_candidate =", K_candidate, "\n")
        for (fold_num in 1:nfold){
          cat("fold =", fold_num, "\n")
          save_dir <- "simulation_semisynthetic/data/2_semi_selectK"
          save_dir <- paste0(save_dir, "/N", N, "_J", J)
          save_dir <- paste0(save_dir, "/seedtheta", seed_theta, "_seedY", seed_Y)
          file_name <- paste0(save_dir, "/K", K_candidate, "_fold", fold_num, ".RData")
          tryCatch({
            load(file_name)
            res <- c(N, seed_theta, seed_Y, K_candidate, fold_num, 
                     predictive_llk[1], predictive_llk_dist_mean)
            dat <- rbind(dat, res)
          }, error=function(cond) {message(paste("File does not seem to exist:", file_name))})
        }
      }
    }
  }
  return(dat)
}
out <- mclapply(1:nseed_theta, run_one, mc.cores = 5L)
out <- do.call(rbind, out)
dat <- data.table(out)
colnames(dat) <- colnames_dat
write.csv(dat, file = "simulation_semisynthetic/data/2_simu_semi_selectK_summary.csv", 
          row.names = FALSE)
