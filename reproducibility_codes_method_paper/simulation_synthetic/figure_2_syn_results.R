#######################################################
#' Generate Figure 2 in the manuscript using summaries of
#' synthetic data simulation.
#######################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(grid)

# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

## read in data
dat <- fread("simulation_synthetic/data/1_simu_syn_summary.csv")
ntree <- 4
Ns <- c(100, 200, 400)

# The palette with grey:
cbPalette <- c("#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
for (x in c("rmse", "tree_frobenius", "ari")) {
  dat[[x]] <- as.numeric(dat[[x]])
}
# keep only single linkage method results for BayesLCM
dat <- dat[(method_hclust == "single" & label == "BayesLCM (homogeneous var) + HC") |
             (method_hclust == "single" & label == "BayesLCM (heterogeneous var) + HC") | 
             (label %in% c("DDT-LCM", "DDT-LCM (true tree)", 
                           "DDT-LCM (misspecified tree)", "DDT-LCM (homogeneous var)")),] 
dat$label <- factor(dat$label, levels = c(
  "DDT-LCM", 
  "DDT-LCM (true tree)",
  "DDT-LCM (misspecified tree)",
  "DDT-LCM (homogeneous var)",
  "BayesLCM (heterogeneous var) + HC",
  "BayesLCM (homogeneous var) + HC"
), labels = c(
  "(i) DDT-LCM", 
  "(ii) DDT-LCM \n(true tree)", 
  "(iii) DDT-LCM \n(misspecified tree)", 
  "(iv) DDT-LCM \n(homogeneous var)",
  "(v) BayesLCM \n(heterogeneous var) + HC",
  "(vi) BayesLCM \n(homogeneous var) + HC"
))
dat$N <- paste0("N = ", dat$N)
dat$N <- factor(dat$N, levels = paste0("N = ", Ns))
dat$tree <- paste0("Tree ", dat$tree)
dat$tree <- factor(dat$tree, levels = paste0("Tree ", 1:ntree))

### put all plots together -------------------------
dat$log_rmse <- log(dat$rmse)
dat_long <- melt(dat[, .(N, tree, label, log_rmse, ari)], 
                 id.vars = c("N", "tree", "label"))
dat2_long <- melt(dat[, .(N, tree, label, tree_frobenius)], 
                  id.vars = c("N", "tree", "label"))
dat_all <- rbind(dat_long, dat2_long)
dat_all$variable <- factor(dat_all$variable,
                           levels = c("tree_frobenius", "log_rmse", "ari"), 
                           labels = c("Frobenius Norm of \nTree Covariance Matrix",
                                      "Log RMSE of Item \nResponse Probabilities",
                                      "Class Membership ARI")) 


## add gaps between boxplots for unshown methods  -------------
# add pseudo data to create gaps in methods not shown in the figure
set.seed(1)
dat_pseudo <- c()
for (tree_idx in 1:ntree) {
  for (N in Ns) {
    for (method in c("(ii) DDT-LCM \n(true tree)", "(iii) DDT-LCM \n(misspecified tree)")) {
      dat_pseudo <- rbind(dat_pseudo,
                          c(N = paste0("N = ", N), 
                            tree = paste0("Tree ", tree_idx), 
                            label = method, 
                            variable = "Frobenius Norm of \nTree Covariance Matrix", 
                            value = runif(1, -2, -1)))
    }
  }
}
dat_pseudo <- data.table(dat_pseudo)
dat_pseudo$value <- as.numeric(dat_pseudo$value)
dat_all2 <- rbind(dat_all, dat_pseudo)

# We will truncate y axis to make the plot more appealing
axis_limits_discounted <- data.table(N = rep(c("N = 100", "N = 200", "N = 400"), 3),
                                     variable = rep(c("Frobenius Norm of \nTree Covariance Matrix", 
                                                      "Log RMSE of Item \nResponse Probabilities",
                                                      "Class Membership ARI"), each = 3),
                                     ymin = c(-1, 0.1, 0.1, 
                                              -3.3, -3.3, -3.3,
                                              0.5, 0.5, 0.5),
                                     ymax = c(1, 1, 1,       
                                              -2., -2., -2.,
                                              1.02, 1.02, 1.02))
lfun <- function(limits) {
  grp <- dat_all2$variable[which(abs(dat_all2$value - limits[1]) < 1e-10)[1]]
  grp <- as.character(grp)
  grp_N <- dat_all2$N[which(abs(dat_all2$value - limits[1]) < 1e-10)[1]]
  grp_N <- as.character(grp_N)
  lim_max <- axis_limits_discounted[variable == grp & N == grp_N,]$ymax
  lim_min <- axis_limits_discounted[variable == grp & N == grp_N,]$ymin
  return(c(lim_min, lim_max))
}


# solution found via https://stackoverflow.com/questions/12207419/r-how-do-i-use-coord-cartesian-on-facet-grid-with-free-ranging-axis
p1 <- ggplot(data = dat_all2, 
       aes(x = tree, y = value, fill = label, drop = FALSE)) +
  geom_boxplot(aes(middle = mean(value))) +
  facet_grid(variable~N, scales = "free_y", switch="y") +
  scale_fill_manual(values=cbPalette, drop = FALSE) +
  labs(y = "", x = "",
       fill = "Method") +
  theme_bw() +
  theme(
    axis.text.x = element_text(vjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "top") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) + 
  scale_y_continuous(limits = lfun, expand=expansion(0,0))

(p2 <- p1 + coord_cartesian(ylim = c(0, 1)))

#' Now we modify y-axis limits of top facet row. 
#' We generate two grobs g1 and g2 and replace the panels and 
#' y-axis in g1 with the corresponding elements from g2.
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

# Replace the upper panels and upper axis of p1 with that of p2
# i.e. replace the corresponding grobs in g1 with the versions of g2
# Panel numbering goes row-wise from top left to bottom right
panels_to_replace_with_g2 <- c(1, 2, 3)

# To get names of grobs run: lapply(g1[["grobs"]],function(x) x$name)
# Use names of grobs to find out indices of g1[["grobs"]] of the panels we want to replace
# as well as for the axis ticks.
pattern_for_specific_panels <- 
  paste0("^panel-((", paste0(panels_to_replace_with_g2, collapse = ")|("),"))")
pattern_for_axes_ticks <- "^GRID.absoluteGrob"
idx_panels_to_replace_from_g2 <- which(unlist(
  lapply(g1[["grobs"]], function(x) grepl(pattern_for_specific_panels, x$name))))
# > idx_panels_to_replace_from_g2
# [1] 2 5 8

idx_axesTicks <- which(unlist(
  lapply(g1[["grobs"]], function(x) grepl(pattern_for_axes_ticks, x$name))))

# Find out manually which of the defined axesTicks it is:
g_test <- g1
for (itr in idx_axesTicks) {
  g_test[["grobs"]][[itr]] <- g2[["grobs"]][[itr]]
  grid.newpage();grid.draw(g_test); grid.draw(textGrob(itr, just = "top"))
  Sys.sleep(0.1)
}
# We found out it is itr=19
idx_axesTicks_to_replace <- itr

#' Having now found out indices of the panels to be replaced 
#' idx_panels_to_replace_from_g2 as well as the index of the y-axis 
#' element idx_axesTicks_to_replace. We can replace them in the following.
# Replace panels
grid.newpage();grid.draw(g1)
for (iter in idx_panels_to_replace_from_g2) {
  g1[["grobs"]][[iter]] <- g2[["grobs"]][[iter]]
  grid.newpage();grid.draw(g1)
  Sys.sleep(0.1)
}
# Replace y-axis
g1[["grobs"]][[idx_axesTicks_to_replace]] <- g2[["grobs"]][[idx_axesTicks_to_replace]]

# Render plot
if(!dir.exists("figures")){
  dir.create("figures")
}
pdf("figures/combined_syn_boxplot.pdf", height = 7, width = 12)
grid.newpage()
grid.draw(g1)
dev.off()
