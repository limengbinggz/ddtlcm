
#' Plot the MAP tree and class profiles of summarized DDT-LCM results
#' @param x a "summary.ddt_lcm" object
#' @param plot_option option to select which part of the plot to return. If "all", return
#'  the plot of MAP tree on the left and the plot of class profiles on the right. If "profile",
#'  only return the plot of class profiles. If "tree", only return the plot of MAP tree.
#' @param item_name_list a named list of G elements, where the g-th element contains a vector 
#'  of item names for items in `item_membership_list[[g]]`. The name of the g-th element is 
#'  the name of the major item group.
#' @param color_palette a vector of color names. Default is a color-blinded friendly palette.
#' @method plot summary.ddt_lcm
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggpubr ggarrange
#' @export
plot.summary.ddt_lcm <- function(x, plot_option = c("all", "profile", "tree"),
                                 item_name_list = NULL, color_palette = c("#E69F00", "#56B4E9", "#009E73", "#000000", 
                                                                          "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999")){
  plot_option = match.arg(plot_option, c("all", "profile", "tree"))
  tree_with_parameter <- x$tree_map
  response_prob <- matrix(x$response_probs_summary[,"Mean"], nrow = K)
  item_membership_list <- x$setting$item_membership_list
  class_probability <- x$class_probs_summary[,"Mean"]
  plots <- plot_tree_with_barplot(tree_with_parameter, response_prob, item_membership_list, 
                         item_name_list, class_probability, color_palette,
                         return_separate_plots = T)
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
#' @param color_palette a vector of color names. Default is a color-blinded friendly palette.
#' @param return_separate_plots If FALSE (default), print the combined plot of MAP tree and 
#'  class profiles. If TRUE, return the tree plot, class profile plot, and data.table
#'  used to create the plots in a list, without printing the combined plot. 
#' @import ggplot2
#' @importFrom ggtext element_markdown
#' @importFrom ggpubr ggarrange
#' @export
plot_tree_with_barplot <- function(tree_with_parameter, response_prob, item_membership_list, 
                                   item_name_list = NULL, class_probability = NULL, 
                                   color_palette = c("#E69F00", "#56B4E9", "#009E73", "#000000", 
                                                     "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#999999"),
                                   return_separate_plots = F
                                           ){
  K <- nTips(tree_with_parameter)
  J <- length(unlist(item_membership_list))
  G <- length(item_membership_list)
  if (ncol(tree_with_parameter@data) != J){
    stop("The input list 'item_membership_list' should contain indices for all columns of
         node parameters on the tree.")
  }
  
  ## plot tree with branch length
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
  response_prob_dat <- merge(x = response_prob_dat, y = map_colors_to_nodes, by = "item_group", all.x = T)
  response_prob_dat[, item_names_factor := factor(item_names, levels = unique(response_prob_dat$item_names))]
  # response_prob_dat$node_label <- paste0("v", response_prob_dat$Class_index)
 
  # add class prevalences
  if (!is.null(class_probability)) {
    class_label <- paste0("Class ", 1:K, ": ", round(class_probability, 2))
  } else{
    class_label <- paste0("Class ", 1:K)
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
#' @importFrom ggbrace geom_brace
#' @importFrom ggtree ggtree geom_tiplab geom_nodelab
#' @importFrom ggpubr ggarrange
#' @export
plot_tree_with_heatmap <- function(tree_with_parameter, response_prob, item_membership_list){
  K <- nTips(tree_with_parameter)
  J <- length(unlist(item_membership_list))
  if (ncol(tree_with_parameter@data) != J){
    stop("The input list 'item_membership_list' should contain indices for all columns of
         node parameters on the tree.")
  }
  t1 <- ggtree(tree_with_parameter) +
    geom_tiplab(size = 4) +
    geom_nodelab(size = 4) +
    labs(title = paste0("Tree")) +
    theme(plot.title = element_text(hjust = 0.5, margin=margin(b=0), vjust = -1, size = 10)) +
    coord_cartesian(clip="off") +
    geom_brace(aes_(x=c(0, 1), y=c(0.4, 0), label = paste("Group", 1)),
               inherit.data=F, rotate=180, labelsize=3.5, color = "white") +
    theme(plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm"))
  # leaf_data <- expit(as.matrix(tree_with_parameter@data[as.character(1:K), 1:J]))
  dat_response_prob <- data.table(value = c(response_prob))
  dat_response_prob[, class := rep(1:K, J)]
  dat_response_prob[, item := rep(1:J, each = K)]
  
  # reorder rows of the heatmap to match the tree leaves
  class_order <- data.table(t1$data[t1$data$label %in% paste0("v", 1:K), c("label", "y")])
  class_order <- class_order[order(-y),]

  class_order <- as.integer(gsub("v", "", class_order$label))

  dat_response_prob$class <- factor(dat_response_prob$class, levels = class_order)
  
  plot_prob <- ggplot(data = dat_response_prob, aes(x = item, y = class, fill = value)) +
    # geom_tile(hjust = 0.5, vjust = 0.5, interpolate = FALSE, na.rm = TRUE) +
    geom_tile() +
    labs(x = "Item", y = "Class", fill = "Item Response Probabilities") +
    scale_fill_gradientn(limits = c(0,1), colours=c("gold", "springgreen4", "blue2") ) + 
    theme(plot.title = element_text(hjust = 0.5, vjust = 10, size = 15),
          panel.grid.major = element_blank(),
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
          plot.margin=unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    coord_cartesian(clip = "off")
  
  G <- length(item_membership_list)
  dat_brace <- c()
  for (g in 1:G) {
    dat_brace <- c(dat_brace, item_membership_list[[g]][1], item_membership_list[[g]][length(item_membership_list[[g]])])
  }
  dat_brace <- data.table(x = dat_brace)
  dat_brace[, y := rep(c(0.4, 0), G)]
  dat_brace[, group_label := rep(1:G, each = 2)]
  plot_prob + geom_brace(aes(dat_brace$x, dat_brace$y, label = group_label), inherit.data=F)
  for (g in 1:G) {
    suppressWarnings({
      plot_prob <- plot_prob +
        geom_brace(aes_(x=c(item_membership_list[[g]][1], item_membership_list[[g]][length(item_membership_list[[g]])]),
                      y=c(0.4, 0), label = g), 
                   inherit.data=F, rotate=180, labelsize=3.5)
    })
  }
  ggarrange(t1, plot_prob, widths = c(0.2, 0.3), ncol = 2, common.legend = T)
}


