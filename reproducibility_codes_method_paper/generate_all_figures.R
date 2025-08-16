#######################################################
#' Master script to generate all figures.
#' Figures will be saved to the directory figures/
#######################################################

curr_dir <- "/Users/mengbing/Dropbox (University of Michigan)/from_box/research/tree_multivariate_binary/reproducibility_codes"
setwd(curr_dir)

### Part 1: fully synthetic data simulation ------------------

## Generate figure 2: compare DDT-LCM to baseline methods
source("simulation_synthetic/figure_2_syn_results.R")

## Generate figure S4.3: summarize estimation of diffusion variance parameters 
source("simulation_synthetic/figure_S4_3_syn_sigma_est.R")


### Part 2: semi-synthetic data simulation ------------------

## Generate figure 3: compare DDT-LCM to baseline methods
source("simulation_semisynthetic/figure_3_semi_results.R")

## Results reported in Supplementary Section S4.2
## Summarize selection probabilities of K \in {4,5,6,7,8}, by cross-validation
## results are printed on screen
source("simulation_semisynthetic/supplement_S4_2_selectK_result.R")

## Generate figure S4.5: compare 5-fold predictive loglikelihoods 
## of DDT-LCM and DDT-LCM (homogeneous var)
source("simulation_semisynthetic/figure_S4_5_semi_sigma_llk.R")





