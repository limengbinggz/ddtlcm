#' Create trace plots of DDT-LCM parameters
#' @param x a "ddt_lcm" object
#' @param parameter_names a character vector to specify the parameters to be plotted. Each element can take the
#'  be 1) of format "parameter_index" to plot specific parameters with specific indices, 2) of format "parameter"
#'  to plot the parameters across all indices, or 3) equal to "all" to plot all parameters in the model. For 1),
#'  the item response probabilities should be named "responseprob_class,group,item"; the class probabilities should be
#'  named "classprob_class"; the divergence function parameter is "c"; the group-specific diffusion variances
#'  should be named "diffusionvar_group". For 2), "responseprob" to plot all item response probabilities;
#'  "classprob" to plot all class probabilities; "diffusionvar" to plot all diffusion variances.
#' @param burnin the number of posterior samples as burn-in, which will not be plotted.
#' @param \dots	Further arguments passed to each method
#' @method plot ddt_lcm
#' @importFrom graphics par
#' @return NULLs
#' @export
#' @examples
#' data(result_diet_1000iters)
#' # Plot "c" for the divergence function parameter; "diffusionvar_1" for diffusion variance of group 1
#' plot(x = result_diet_1000iters, parameter_names = c("c", "diffusionvar_1"), burnin = 500)
#' # Plot "responseprob_1,1,1" for the class 1 response probability of item 3 in major group 2
#' plot(x = result_diet_1000iters, parameter_names = "responseprob_1,1,1", burnin = 500)
#' # Plot "classprob_1" for the probability of being assigned to class 1
#' plot(x = result_diet_1000iters, parameter_names = "classprob_1", burnin = 500)
#' # plot all class probabilities
#' plot(x = result_diet_1000iters, parameter_names = "classprob", burnin = 500)
#' # plot all diffusion variances
#' plot(x = result_diet_1000iters, "diffusionvar", burnin = 500)
plot.ddt_lcm <- function(x, parameter_names = c("responseprob_1,1,1", "classprob_1", "c", "diffusionvar_1"), burnin = 50, ...){
  # number of classes
  K <- x$setting$K
  # number of major item groups
  G <- x$setting$G
  # extract list of item indices according to groups
  item_membership_list <- x$setting$item_membership_list
  # # get the number of items each group
  # num_items_per_group <- unlist(lapply(item_membership_list, length))
  # # cumulative index
  # cum_num_items_per_group <- c(0, cumsum(num_items_per_group))
  # number of iterations in total
  total_iters <- length(x$loglikelihood)

  ## input checking
  good_names <- grepl("^responseprob|^classprob|^c$|^diffusionvar", parameter_names)
  bad_names <- which(!good_names)
  if (length(bad_names) >= 1){
    stop(paste0("parameter_names[", bad_names, "] invalid. Please check the appropriate format in the function description."))
  }
  if (as.integer(burnin) != burnin | burnin < 0 | burnin >= total_iters){
    stop("The number of burn-in posterior samples, burnin, should be a positive integer smaller than the total number of samples.")
  }

  if (length(parameter_names) == 1 & parameter_names[1] == "all"){
    parameter_names <- c("responseprob", "classprob", "c", "diffusionvar")
  }
  parameter_names_all <- c()
  for (param in parameter_names) {
    param_split <- strsplit(param, "_")[[1]]
    # get which parameter of the model we are looking at
    param_type <- param_split[1]
    # if no indices are specified, then plot all parameters of that type
    if (length(param_split) == 1 & param_type != "c"){
      if (param_split == "responseprob") {
        index_group_item <- unlist(lapply(1:length(item_membership_list), function(x) paste0(x, ",", 1:length(item_membership_list[[x]]))))
        index_class_group_item <- as.vector(outer(1:K, index_group_item, paste, sep=","))
        parameter_names_all <- c(parameter_names_all, paste0("responseprob_", index_class_group_item))
      } else if (param_split == "classprob") {
        parameter_names_all <- c(parameter_names_all, paste0("classprob_", 1:K))
      } else if (param_split == "diffusionvar") {
        parameter_names_all <- c(parameter_names_all, paste0("diffusionvar_", 1:G))
      }
    } else{ # otherwise just keep the original parameter name
      parameter_names_all <- c(parameter_names_all, param)
    }
  }

  # now start plotting
  # if (length(parameter_names_all) > 1) par(ask = TRUE)
  for (param in parameter_names_all) {
    param_split <- strsplit(param, "_")[[1]]
    # get which parameter of the model we are looking at
    param_type <- param_split[1]
    if (param != "c"){
      # get index of the parameter as a string
      param_index <- param_split[2]
      # get the index of the parameter as a numeric vector
      param_index <- as.integer(strsplit(param_index, ",")[[1]])
    }
    if (param_type == "responseprob"){ # item response probability
      # find out the item index in all columns
      j <- item_membership_list[[param_index[2]]][param_index[3]]
      samples <- x$response_probs_samples[(burnin + 1):total_iters, param_index[1], j]
      plot(samples, type = 'l', xlab = 'Iteration',
           ylab = paste0('Positive response probability in class ', K, ' of item ', param_index[3], ' in group ', param_index[2]))
    } else if (param_type == "classprob"){ # item response probability
      samples <- x$class_probs_samples[param_index, (burnin + 1):total_iters]
      plot(samples, type = 'l', xlab = 'Iteration', ylab = paste0('Probability of class ', param_index))
    } else if (param_type == "c"){ # divergence function parameter
      samples <- x$c_samples[(burnin + 1):total_iters]
      plot(samples, type = 'l', xlab = 'Iteration', ylab = paste0('Divergence function parameter c'))
    } else if (param_type == "diffusionvar"){ # divergence function parameter
      samples <- x$Sigma_by_group_samples[param_index, (burnin + 1):total_iters]
      plot(samples, type = 'l', xlab = 'Iteration', ylab = paste0('Diffusion variance of group ', param_index))
    }
  }

}


#' Plot the MAP tree and class profiles of summarized DDT-LCM results
#' @param x a "summary.ddt_lcm" object
#' @param plot_option option to select which part of the plot to return. If "all", return
#'  the plot of MAP tree on the left and the plot of class profiles on the right. If "profile",
#'  only return the plot of class profiles. If "tree", only return the plot of MAP tree.
#' @param item_name_list a named list of G elements, where the g-th element contains a vector
#'  of item names for items in `item_membership_list[[g]]`. The name of the g-th element is
#'  the name of the major item group.
#' @param color_palette a vector of color names. Default is a color-blinded friendly palette.
#' @param log Default argument passed to plot(). Not used.
#' @param \dots	Further arguments passed to each method
#' @method plot summary.ddt_lcm
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggpubr ggarrange
#' @return a ggplot2 object. If plot_option is "all", then a plot with maximum a posterior
#'  tree structure on the left and a bar plot of item response probabilities (with 95%
#'  credible intervals and class probabilities) on the right. If plot_option is
#'  "profile", then only a bar plot of item response probabilities. If plot_option is
#'  "tree", then only a plot of the tree structure.
#' @export
#' @examples
#' data(result_diet_1000iters)
#' burnin <- 500
#' summarized_result <- summary(result_diet_1000iters, burnin, relabel = TRUE, be_quiet = TRUE)
#' plot(x = summarized_result, item_name_list = NULL, plot_option = "all")
plot.summary.ddt_lcm <- function(x, log=TRUE,
                                 plot_option = c("all", "profile", "tree"),
                                 item_name_list = NULL, color_palette = c("#E69F00", "#56B4E9", "#009E73", "#000000",
                                                                          "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999"),
                                 ...){
  plot_option = match.arg(plot_option, c("all", "profile", "tree"))
  tree_with_parameter <- x$tree_map
  K <- nrow(x$tree_Sigma)
  response_prob <- matrix(x$response_probs_summary[,"Mean"], nrow = K)
  item_membership_list <- x$setting$item_membership_list
  class_probability <- x$class_probs_summary[,"Mean"]
  class_probability_lower <- x$class_probs_summary[,"2.5%"]
  class_probability_higher <- x$class_probs_summary[,"97.5%"]
  plots <- plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list,
                         item_name_list, class_probability, class_probability_lower, class_probability_higher,
                         color_palette, return_separate_plots = TRUE)
  probs_lower <- x$response_probs_summary[, "2.5%"]
  probs_higher <- x$response_probs_summary[, "97.5%"]
  response_prob_dat <- plots[["response_prob_dat"]]
  response_prob_dat[, probs_lower := probs_lower]
  response_prob_dat[, probs_higher := probs_higher]
  if (plot_option == "all"){
    plot_lcm <- plots[['lcm_plot']] +
      geom_errorbar(data = response_prob_dat,
                    aes(ymin=probs_lower, ymax=probs_higher), linewidth=0.2)
    return(ggarrange(plots[['tree_plot']], plot_lcm, widths = c(1, 1.5), ncol = 2))
  } else if (plot_option == "profile"){
    plot_lcm <- plots[['lcm_plot']] +
      geom_errorbar(data = response_prob_dat,
                    aes(ymin=probs_lower, ymax=probs_higher), linewidth=0.2)
    return(plot_lcm)
  } else if (plot_option == "tree"){
    return(plots[['tree_plot']])
  } else {
    stop("Please select an option from 'all', 'profile', or 'tree'.")
  }
}


#' Plot the MAP tree and class profiles (bar plot) of summarized DDT-LCM results
#' @param tree_with_parameter a "phylo4d" tree with node parameters
#' @param response_prob a K by J matrix, where the k,j-th element is the response
#'  probability of item j for individuals in class k
#' @param item_membership_list a list of G elements, where the g-th element contains the
#'  column indices of the observed data matrix corresponding to items in major group g
#' @param item_name_list a named list of G elements, where the g-th element contains a vector
#'  of item names for items in `item_membership_list[[g]]`. The name of the g-th element is
#'  the name of the major item group.
#' @param class_probability a length K vector, where the k-th element is the
#'  probability of assigning an individual to class k. It does not have to sum up to 1
#' @param class_probability_lower a length K vector, 2.5% quantile of posterior
#'  the distribution.
#' @param class_probability_higher a length K vector, 97.5% quantile of posterior
#'  the distribution.
#' @param color_palette a vector of color names. Default is a color-blinded friendly palette.
#' @param return_separate_plots If FALSE (default), print the combined plot of MAP tree and
#'  class profiles. If TRUE, return the tree plot, class profile plot, and data.table
#'  used to create the plots in a list, without printing the combined plot.
#' @import ggplot2
#' @importFrom ggtext element_markdown
#' @importFrom ggpubr ggarrange
#' @return a ggplot2 object. A bar plot of item response probabilities.
#' @export
#' @examples
#' # load the MAP tree structure obtained from the real HCHS/SOL data
#' data(data_synthetic)
#' # extract elements into the global environment
#' list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
#' plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list)
plot_tree_with_barplot <- function(tree_with_parameter, response_prob, item_membership_list,
                                   item_name_list = NULL, class_probability = NULL,
                                   class_probability_lower = NULL, class_probability_higher = NULL,
                                   color_palette = c("#E69F00", "#56B4E9", "#009E73", "#000000",
                                                     "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999"),
                                   return_separate_plots = FALSE
                                           ){
  K <- nTips(tree_with_parameter)
  J <- length(unlist(item_membership_list))
  G <- length(item_membership_list)
  if (ncol(tree_with_parameter@data) != J){
    stop("The input list 'item_membership_list' should contain indices for all columns of
         node parameters on the tree.")
  }


  ## plot tree with branch length
  branch = branch.length = NULL # Setting the variables to NULL first to avoid NSE notes in R CMD check
  t1 <- ggtree(tree_with_parameter) +
    geom_tiplab() + #size = 8
    labs(title = "Tree") +
    coord_cartesian(clip = 'off') +
    theme(plot.title = element_text(hjust = 0.5)) + #, size = 16
    theme(plot.margin=unit(c(0,0,0,0), "cm")) +
    geom_text(aes(x=branch, label=signif(branch.length, 2)), #size=6,
              vjust=-.3, color = "firebrick")
  # need to change the tip label of the tree to match the class index from 1 to K
  class_order <- data.table(t1$data[t1$data$label %in% paste0("v", 1:K), c("label", "y")])
  class_order <- class_order[order(-y),]

  # due to NSE notes in R CMD check
  node_label_new = Class_index = Item = item_names = item_group = item_names_factor = NULL
  y = value = color = NULL

  class_order[, node_label_new := paste0("v", 1:K)]
  # class_order <- data.table(t1$data[t1$data$label %in% paste0("v", 1:K), c("label", "y")] %>% arrange(-y))
  # node_order_new <- order(as.numeric(gsub("v", "", class_order$label)))
  original_label <- t1$data$label[t1$data$label %in% paste0("v", 1:K)]
  node_order_new <- c()
  for (k in 1:K) {
    node_order_new <- c(node_order_new, class_order$node_label_new[class_order$label == original_label[k]])
  }
  t1$data$label[t1$data$label %in% paste0("v", 1:K)] <- node_order_new


  ### plot response profiles
  # response_prob_dat <- data.table(reshape2::melt(response_prob))
  response_prob_dat <- data.table(value = c(response_prob))
  response_prob_dat[, Class_index := rep(1:K, J)]
  response_prob_dat[, Item := rep(1:J, each = K)]
  # colnames(response_prob_dat) <- c("Class_index", "Item", "mean")
  if (!is.null(class_probability)){
    response_prob_dat[, class_probability := rep(round(class_probability, 3), J)]
  }

  item_order <- order(unlist(item_membership_list))
  if (!is.null(item_name_list)){
    item_names <- unlist(item_name_list)[unlist(item_membership_list)]
    item_names_ordered <- unlist(item_name_list)[item_order]
    item_groups <- names(item_name_list)
  } else{
    item_groups <- paste0("Group ", 1:G)
    item_names_ordered <- paste0("Item ", 1:J)
  }
  item_groups_vector <- c()
  for (g in 1:G) {
    item_groups_g <- rep(item_groups[g], length(item_membership_list[[g]]))
    item_groups_vector <- c(item_groups_vector, item_groups_g)
  }
  item_groups_vector <- item_groups_vector[item_order]
  response_prob_dat[, item_names := rep(item_names_ordered, each = K)]
  response_prob_dat[, item_group := rep(item_groups_vector, each = K)]
  response_prob_dat[, item_group := factor(item_group, levels = item_groups)]

  ## create color for major groups
  color_palette <- rep(color_palette, ceiling(length(color_palette) / G))[1:G]
  map_colors_to_nodes <- data.table(item_group = item_groups, color = color_palette)
  response_prob_dat <- merge(x = response_prob_dat, y = map_colors_to_nodes, by = "item_group", all.x = TRUE)
  response_prob_dat[, item_names_factor := factor(item_names, levels = unique(response_prob_dat$item_names))]
  # response_prob_dat$node_label <- paste0("v", response_prob_dat$Class_index)

  # add class prevalences
  if (!is.null(class_probability)) {
    class_label <- paste0("Class ", 1:K, ": ", round(class_probability * 100, 0), "%")
  } else{
    class_label <- paste0("Class ", 1:K)
  }
  if (!(is.null(class_probability_lower) && is.null(class_probability_higher))) {
    class_label <- paste0(
      class_label,
      " (", round(class_probability_lower * 100, 0), "%, ",
      round(class_probability_higher * 100, 0), "%)"
    )
  }
  response_prob_dat$Class <- factor(response_prob_dat$Class_index, levels = 1:K, labels = class_label)

  plot_lcm <- ggplot(response_prob_dat, aes(x = item_names_factor, y = value, fill = item_group)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    facet_wrap(~ Class, ncol = 1) +
    # geom_errorbar( aes(ymin=probs_lower, ymax=probs_higher), width=0.2, size=0.5) + #, color = root
    scale_fill_manual(values = unique(response_prob_dat[, .(item_group, color)])$color) +
    scale_color_manual(values = unique(response_prob_dat[, .(item_group, color)])$color) +
    scale_y_continuous(limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1, #size = 10,
                                         colour = unique(response_prob_dat[, .(item_names_factor, color)])$color),
          # axis.title = element_text(size = 16),
          # legend.text = element_text(size = 15),
          # legend.title = element_text(size = 15),
          # strip.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5), #size = 16,
          plot.margin=unit(c(0, 0, 0, 0.3), "cm")) +
    labs(x = "Items", y = "Item Response Probability", fill = "Group",
         title = "Class Profiles")

  if (!return_separate_plots){
    # grid.arrange(t1, plot_lcm, ncol=2, widths=c(1, 1.5))
    ggarrange(t1, plot_lcm, widths = c(1, 1.5), ncol = 2)
  } else{
    return(list(tree_plot = t1, lcm_plot = plot_lcm, response_prob_dat = response_prob_dat))
  }
}




#' Plot the MAP tree and class profiles (heatmap) of summarized DDT-LCM results
#' @param tree_with_parameter a "phylo4d" tree with node parameters
#' @param response_prob a K by J matrix, where the k,j-th element is the response
#'  probability of item j for individuals in class k
#' @param item_membership_list a list of G elements, where the g-th element contains the
#'  column indices of the observed data matrix corresponding to items in major group g
#' @import ggplot2
#' @importFrom ggtext element_markdown
#' @importFrom ggpubr ggarrange
#' @importFrom ggtree ggtree geom_tiplab geom_nodelab
#' @importFrom ggpubr ggarrange
#' @return a ggplot2 object. A plot with the tree structure on the left and a heatmap
#'  of item response probabilities on the right, with indication of item group
#'  memberships beneath the heatmap.
#' @export
#' @examples
#' # load the MAP tree structure obtained from the real HCHS/SOL data
#' data(data_synthetic)
#' # extract elements into the global environment
#' list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
#' plot_tree_with_heatmap(tree_with_parameter, response_prob, item_membership_list)
plot_tree_with_heatmap <- function(tree_with_parameter, response_prob, item_membership_list){
  K <- nTips(tree_with_parameter)
  J <- length(unlist(item_membership_list))
  if (ncol(tree_with_parameter@data) != J){
    stop("The input list 'item_membership_list' should contain indices for all columns of
         node parameters on the tree.")
  }
  branch = branch.length = NULL # Setting the variables to NULL first to avoid NSE notes in R CMD check

  x = NULL
  t1 <- ggtree(tree_with_parameter) +
    geom_tiplab(size = 4) +
    geom_nodelab(size = 4) +
    labs(title = paste0("Tree"), x = "Item", y = "Class") +
    theme(plot.title = element_text(hjust = 0.5, margin=margin(b=0), vjust = -1, size = 10)) +
    coord_cartesian(clip="off") +
    theme(plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
          axis.title = element_text(size = 10, color = "white")) +
    geom_text(aes(x=branch, label=signif(branch.length, 2)), #size=6,
              vjust=-.3, color = "firebrick") #+

  dat_response_prob <- data.table(value = c(response_prob))
  item = y = value = group_label = NULL # due to NSE notes in R CMD check
  dat_response_prob[, class := rep(1:K, J)]
  dat_response_prob[, item := rep(1:J, each = K)]
  # reorder rows of the heatmap to match the tree leaves
  class_order <- data.table(t1$data[t1$data$label %in% paste0("v", 1:K), c("label", "y")])
  class_order <- class_order[order(-y),]
  class_order <- as.integer(gsub("v", "", class_order$label))
  dat_response_prob$class <- factor(dat_response_prob$class, levels = class_order)

  plot_prob <- ggplot(data = dat_response_prob, aes(x = item, y = class, fill = value)) +
    geom_tile() +
    labs(x = "Item", y = "Class", fill = "Item Response Probabilities", title = "blank") +
    scale_fill_gradientn(limits = c(0,1), colours=c("gold", "springgreen4", "blue2") ) +
    theme(plot.title = element_text(hjust = 0.5, margin=margin(b=0), vjust = -1, size = 10, color = "white")) +
    coord_cartesian(clip = "off") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position="bottom",
          plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"))
  G <- length(item_membership_list)
  for (g in 1:G) {
    x1 = item_membership_list[[g]][1]; x2 = item_membership_list[[g]][length(item_membership_list[[g]])]
    y1 = 0.2; y2 = 0.3
    dat_bracket <- data.frame(
      x = c(x1, x1, x2, x2),
      y = c(y2, y1, y1, y2)
    )
    plot_prob <- plot_prob +
      annotate("text", x = mean(c(x1, x2)), y = y1 - 0.2, label = paste0(g), size = 3.5) +
      geom_line(data = dat_bracket, aes(x = x, y = y), inherit.aes = FALSE)
  }
  ggarrange(t1, plot_prob, widths = c(0.2, 0.3), ncol = 2, common.legend = TRUE)
}


