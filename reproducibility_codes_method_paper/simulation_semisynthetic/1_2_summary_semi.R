library(data.table)
library(parallel)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)
K <- 6
load("hchs_simulated_data/tree_real_nodata.RData")
tree <- tree_map

# add locations and group selection -----------------------------------------
num_items_per_group <- J_g
G <- length(J_g)
item_group_membership <- rep(1:G, J_g)
# total number of items
J <- sum(J_g) #78
class_probability_true <- class_probability
# group-specific diffusion variances
Sigma_by_group <- c(2.3, 1, 1, 2.3, 2.3, 1, 2.3)^2

n_seed_theta <- 20
n_seed_Y <- 5
initialization_num_input <- 1
Ns <- c(400, 800)
seed_theta_start <- 1

### Part 1: summarize tree and LCM model parameter estimates ------------------------------
dat_all <- c()
for (N in Ns) { 
  
  ### classical LCM results with homogeneous variances -------------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (seed_Y in 1:n_seed_Y) {
      file_name <- paste0("simulation_semisynthetic/data/1_semi_result", 
                          "/N", N, "_J", J, "_K", K, "/seedtheta", seed_theta, 
                          "/bayeslcm_homovar", "/simdat_seed", seed_Y, ".RData")
      print(file_name)
      tryCatch({
        load(file_name)
        dat_result <- rbind(dat_result, summarized_result_bayeslcm_homovar)
      }, error=function(cond) {message(paste("File does not seem to exist:", file_name))})
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_bayeslcm_homovar <- data.table(out)
  
  
  ### classical LCM results with group variances -------------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (seed_Y in 1:n_seed_Y) {
      file_name <- paste0("simulation_semisynthetic/data/1_semi_result", 
                          "/N", N, "_J", J, "_K", K, "/seedtheta", seed_theta,
                          "/bayeslcm_groupvar", "/simdat_seed", seed_Y, ".RData")
      print(file_name)
      tryCatch({
        load(file_name)
        dat_result <- rbind(dat_result, summarized_result_bayeslcm_groupvar)
      }, error=function(cond) {message(paste("File does not seem to exist:", file_name))})
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_bayeslcm_groupvar <- data.table(out)
  
  
  ### DDT-LCM results with homogeneous variances -------------------------------------
  run_one <- function(seed_theta){
    dat_result <- c()
    for (seed_Y in 1:n_seed_Y) {
      file_name <- paste0("simulation_semisynthetic/data/1_semi_result", 
                          "/N", N, "_J", J, "_K", K, "/seedtheta", seed_theta, 
                          "/ddtlcm_homovar", "/simdat_seed", seed_Y, ".RData")
      print(file_name)
      tryCatch({
        load(file_name)
        dat_result <- rbind(dat_result, summarized_result_ddtlcm_homovar)
      }, error=function(cond) {message(paste("File does not seem to exist:", file_name))})
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_ddtlcm_homovar <- data.table(out)
  
  
  ### DDT-LCM results -------------------------------------
  fix_tree_list <- c("FALSE", "misspecified", "truth") 
  run_one <- function(seed_theta){
    dat_result <- c()
    for (seed_Y in 1:n_seed_Y) {
      for (fix_tree_at in fix_tree_list) {
        file_name <- paste0("simulation_semisynthetic/data/1_semi_result", 
                            "/N", N, "_J", J, "_K", K, "/seedtheta", seed_theta, 
                            "/ddtlcm_fixtree_", fix_tree_at, 
                            "/simdat_seed", seed_Y, ".RData")
        print(file_name)
        tryCatch({
          load(file_name)
          dat_result <- rbind(dat_result, summarized_result_ddtlcm)
          rm(summarized_result_ddtlcm)
        }, error=function(cond) {
          message(paste("File does not seem to exist:", file_name))
        })
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
          "simulation_semisynthetic/data/1_semi_result/1_simu_semi_summary.csv",
          quote = FALSE, row.names = FALSE)




### Part 2: summarize diffusion variance parameter estimates ------------------------------
fix_tree_at <- "FALSE"
for (N in Ns) {
  run_one <- function(seed_theta){
    dat_result <- c()
    for (seed_Y in 1:n_seed_Y) {
      file_name <- paste0("simulation_semisynthetic/data/1_semi_result", 
                          "/N", N, "_J", J, "_K", K, "/seedtheta", seed_theta,
                          "/ddtlcm_fixtree_", fix_tree_at, 
                          "/simdat_seed", seed_Y, ".RData")
        print(file_name)
        tryCatch({
          load(file_name)
          dat_result <- rbind(dat_result, Sigma_dat)
          rm(Sigma_dat)
        }, error=function(cond) {
          message(paste("File does not seem to exist:", file_name))
        })
    }
    return(dat_result)
  }
  out <- mclapply(seed_theta_start:n_seed_theta, run_one, mc.cores = 7L)
  out <- do.call(rbind, out)
  dat_result_ddtlcm <- data.table(out)
  write.csv(dat_result_ddtlcm,
            file = "simulation_semisynthetic/data/1_semi_result/1_simu_semi_sigma.csv",
            quote = FALSE, row.names = FALSE)
}












