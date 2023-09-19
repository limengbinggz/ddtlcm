## Notes

### Second submission

>   We still see:
    Please write TRUE and FALSE instead of T and F.
    'T' and 'F' instead of TRUE and FALSE:
      man/quiet.Rd:
        quiet(x, be_quiet = T)

    Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
    Missing Rd-tags:
         plot.summary.ddt_lcm.Rd: \value
         plot_tree_with_barplot.Rd: \value
         plot_tree_with_heatmap.Rd: \value
         predict.ddt_lcm.Rd: \value
         predict.summary.ddt_lcm.Rd: \value
         print.ddt_lcm.Rd: \value
         print.summary.ddt_lcm.Rd: \value

Sorry that I did not realize \value was missing because @return were not written in documentation. Everything has been corrected. Thanks!


### First submission

>   Please write TRUE and FALSE instead of T and F.
    'T' and 'F' instead of TRUE and FALSE:
      man/plot_tree_with_barplot.Rd:
        plot_tree_with_barplot(
          tree_with_parameter,
          response_prob,
          item_membership_list,
          item_name_list = NULL,
          class_probability = NULL,
          class_probability_lower = NULL,
          class_probability_higher = NULL,
          color_palette = c("#E69F00", "#56B4E9", "#009E73", "#000000", "#0072B2", "#D55E00",
            "#CC79A7", "#F0E442", "#999999"),
          return_separate_plots = F
        )
      man/quiet.Rd:
        quiet(x, be_quiet = T)
      man/summary.ddt_lcm.Rd:
        {summary}{ddt_lcm}(object, burnin = 3000, relabel = T, be_quiet = F, ...)

    Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
    Missing Rd-tags in up to 13 .Rd files, e.g.:
         expit.Rd: \value
         logit.Rd: \value
         plot.Rd: \value
         plot.summary.ddt_lcm.Rd: \value
         plot_tree_with_barplot.Rd: \value
         plot_tree_with_heatmap.Rd: \value
         ...

    Please add a few more small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.

Everything has been corrected. Thanks!



## Test environments
* local MacOS 13.3.1 (a), R 4.2.3
* winbuilder
* rhub check_for_cran

## R CMD check results

── R CMD check results  ddtlcm 0.1.2 ────

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


## Reverse dependencies

None