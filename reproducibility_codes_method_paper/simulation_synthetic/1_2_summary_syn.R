library(data.table)
library(parallel)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)


K <- 3
nsim <- 20

## specify parameters -----------------------------
J_g <- num_items_per_group <- c(rep(10, 5), 15, 15)
G <- length(J_g)
class_probability <- c(0.4, 0.3, 0.3)
# variance of logit response probabilities of items in each group
Sigma_by_group <- c(rep(0.6**2, 5), rep(2**2, 2)) 
theta <- 0
alpha <- 0

# add locations and group selection -----------------------------------------
item_group_membership <- rep(1:G, J_g)
J <- sum(J_g) 

tree_idx_list <- 1:4
n_seed_theta <- 20
n_seed_Y <- 5
initialization_num_input <- 1
Ns <- c(100, 200, 400)
seed_theta_start <- 1

for (N in Ns) {

  ### classical LCM results with homogeneous variances -------------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (tree_idx in tree_idx_list) {
      for (seed_Y in 1:n_seed_Y) {
        file_name <- paste0("simulation_synthetic/data/1_syn_result", 
                            "/N", N, "_J", J, "_K", K, 
                            "/tree", tree_idx, "/seedtheta", seed_theta, 
                            "/bayeslcm_homovar/simdat_seed", seed_Y, ".RData")
        print(file_name)
        tryCatch({
          load(file_name)
          dat_result <- rbind(dat_result, summarized_result_bayeslcm_homovar)
        }, error=function(cond) {message(paste("File does not seem to exist:", file_name))}
        )
      }
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_bayeslcm_homovar <- data.table(out)
  

  ### classical LCM results with group variances -------------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (tree_idx in tree_idx_list) {
      for (seed_Y in 1:n_seed_Y) {
        file_name <- paste0("simulation_synthetic/data/1_syn_result", 
                            "/N", N, "_J", J, "_K", K, 
                            "/tree", tree_idx, "/seedtheta", seed_theta, 
                            "/bayeslcm_heterovar/simdat_seed", seed_Y, ".RData")
        print(file_name)
        tryCatch({
          load(file_name)
          dat_result <- rbind(dat_result, summarized_result_bayeslcm_heterovar)
        }, error=function(cond) {message(paste("File does not seem to exist:", file_name))}
        )
      }
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_bayeslcm_heterovar <- data.table(out)


  ### DDT-LCM with homogeneous variances results --------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (tree_idx in tree_idx_list) {
      for (seed_Y in 1:n_seed_Y) {
        file_name <- paste0("simulation_synthetic/data/1_syn_result", 
                            "/N", N, "_J", J, "_K", K, 
                            "/tree", tree_idx, "/seedtheta", seed_theta, 
                            "/ddtlcm_homovar/simdat_seed", seed_Y, ".RData")
        print(file_name)
        tryCatch({
          load(file_name)
          dat_result <- rbind(dat_result, summarized_result_ddtlcm_homovar)
        }, error=function(cond) {message(paste("File does not seem to exist:", file_name))})
      }
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_ddtlcm_homovar <- data.table(out)
  

  ### DDT-LCM results -------------------------------------
  fix_tree_list <- c("FALSE", "truth", "misspecified")
  run_one <- function(seed_theta){
    dat_result <- c()
    for (tree_idx in tree_idx_list) {
      for (seed_Y in 1:n_seed_Y) {
        for (fix_tree_at in fix_tree_list) {
          file_name <- paste0("simulation_synthetic/data/1_syn_result", 
                              "/N", N, "_J", J, "_K", K, 
                              "/tree", tree_idx, "/seedtheta", seed_theta, 
                              "/ddtlcm_fixtree_", fix_tree_at, 
                              "/simdat_seed", seed_Y, ".RData")
          print(file_name)
          tryCatch({
            load(file_name)
            dat_result <- rbind(dat_result, summarized_result_ddtlcm)
          }, error=function(cond) {
            message(paste("File does not seem to exist:", file_name))
          })
        }
      }
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_ddtlcm <- data.table(out)
  

  ## combine all subdatasets into an integrated one for directly plotting ---------
  dat_tmp <- rbind(dat_result_bayeslcm_homovar, dat_result_bayeslcm_groupvar, 
                   dat_result_ddtlcm_homovar, dat_result_ddtlcm)
  dat_all <- rbind(dat_all, dat_tmp)
}

write.csv(dat_all, 
          "simulation_synthetic/data/1_syn_result/1_simu_syn_summary.csv",
          quote = FALSE, row.names = FALSE)



### Part 2: summarize diffusion variance parameter estimates ------------------------------
fix_tree_at <- "FALSE"
dat_all <- c()
for (N in Ns) {
  run_one <- function(seed_theta){
    dat_result <- c()
    for (tree_idx in tree_idx_list) {
      for (seed_Y in 1:n_seed_Y) {
        for (fix_tree_at in fix_tree_list) {
          file_name <- paste0("simulation_synthetic/data/1_syn_result", 
                              "/N", N, "_J", J, "_K", K, 
                              "/tree", tree_idx, "/seedtheta", seed_theta, 
                              "/ddtlcm_fixtree_", fix_tree_at, 
                              "/simdat_seed", seed_Y, ".RData")
          print(file_name)
          tryCatch({
            load(file_name)
            dat_result <- rbind(dat_result, Sigma_dat)
          }, error=function(cond) {
            message(paste("File does not seem to exist:", file_name))
          })
        }
      }
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_tmp <- data.table(out)
  dat_all <- rbind(dat_all, dat_tmp)
}
write.csv(dat_all,
          file = "simulation_synthetic/data/1_syn_result/1_simu_syn_sigma.csv",
          quote = FALSE, row.names = FALSE)

