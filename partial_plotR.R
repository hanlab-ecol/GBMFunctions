#####Function to make partial dependence plots of all variables wanted
##
##uses the packages patchwork, tidyverse, gtable, magrittr
##arguments required: 
#data - data.frame of full output from xgboost
#hist.data - data.frame of frequencies for plotting histograms (column names currently are variable.name, varaible.value for the x axis and value for the y axis frequency)
#vars - a vector of the variables of interest. It may be simplest to simply do unique(data$variable.name) to grab all variables
#type - whether you want the plots to show the mean and 95% confidence interval (as derived from the t-distribution) or mean and all bootstrap results
#histogram - a logical argument with TRUE meaning that histogram plots are made and overlayed with the partial dependency plots. 
#cleaned names is an optional argument to give a vector of panel names (in the same order as the variables) that will replace the variable names. 
#This is more useful for publication purposes

##example: partial_plot(pd_out, out, vars = levels(out$variable.name), type = "mean", histogram = T)
partial_plot <- function(data, hist.data, vars, type = c("mean", "all"), histogram = T, cleaned_names = vars,...) {
  if(is.null(cleaned_names)) {
    cleaned_names <- vars
    }
    PLT <- lapply(vars, function(vars) {
    BOOT <- data %>%
      filter(variable.name == vars) %$%
      bootstrap_run %>%
      unique
    if(length(BOOT) > 1){
      ROWN <- sapply(BOOT, function(j) data %>%
                       filter(variable.name == vars & bootstrap_run == j) %$%
                       x %>%
                       length)
      ROWN <- do.call(c, lapply(ROWN, function(i) 1:i))
      UNIQ <- data %>%
        filter(variable.name == vars) %$%
        x %>%
        unique %>%
        as.character %>%
        unique %>%
        as.numeric
      TAB <- data %>%
        filter(variable.name == vars) %$%
        x %>%
        table
      NUMS <- sort(UNIQ[which(UNIQ %in% names(sort(TAB, decreasing = T)[1:max(ROWN)]))])
      if(length(NUMS) != max(ROWN)) {
        NUMA <- data %>%
          filter(variable.name == vars) %$%
          x %>%
          table/max(BOOT)
        if(length(NUMA) == 1) NUMS <- rep(names(NUMA), NUMA) else {
        NUMS <- do.call(c, sapply(1:length(NUMA), function(x) rep(names(NUMA)[x], NUMA[x])))
        }
      }
      if(type == "all") {
        data %>%
          filter(variable.name == vars) %>%
          mutate(rown = ROWN) %>%
          group_by(rown) %>%
          summarize(mean = mean(yhat),
                    conf25 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[1]),
                    conf95 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[2])) %>%
          mutate(var = NUMS) %>%
          ggplot() +
          geom_line(data = data %>%
                      filter(variable.name == vars),
                    aes(x = x,
                        y = yhat,
                        group = bootstrap_run),
                    color = "grey50") +
          #geom_ribbon(aes(x = var, ymin = conf25, ymax = conf95), alpha = 0.5) +
          geom_line(aes(x = var,
                        y = mean),
                    color = "black",
                    size = 1.5) +
          labs(x = NULL,
               y = NULL) +
          theme(panel.background = element_blank(),
                panel.border = element_rect(fill = NA,
                                            color = "black",
                                            size = 0.8),
                panel.grid = element_line(color = "transparent"),
                axis.title = element_text(face = "bold"))} 
      else {
            data %>%
          filter(variable.name == vars) %>%
          mutate(rown = ROWN) %>%
          group_by(rown) %>%
          summarise(mean = mean(yhat),
                    conf25 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[1]),
                    conf95 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[2])) %>%
          mutate(var = NUMS) %>%
          ggplot() +
              geom_ribbon(aes(x = var,
                              ymin = conf25,
                              ymax = conf95),
                          alpha = 0.5) +
              geom_line(aes(x = var,
                            y = mean),
                        color = "black",
                        size = 1.5) +
              labs(x = NULL,
                   y = NULL) +
              theme(panel.background = element_blank(),
                    panel.border = element_rect(fill = NA,
                                                color = "black",
                                                size = 0.8),
                    panel.grid = element_line(color = "transparent"),
                    axis.text = element_text(size = 7))
          }} else {
            data %>%
              filter(variable.name == vars) %>%
              ggplot() +
              geom_line(data = data %>%
                          filter(variable.name == vars),
                        aes(x = x,
                            y = yhat),
                        color = "black",
                        size = 1.5) +
              labs(x = NULL,
                   y = NULL) +
              theme(panel.background = element_blank(),
                    panel.border = element_rect(fill = NA,
                                                color = "black",
                                                size = 0.8),
                    panel.grid = element_line(color = "transparent"))  
          }})
  if(isTRUE(histogram)) {
    PLT2 <- lapply(vars, function(vars) {
      hist.data %>%
        filter(variable.name == vars) %>%
        ggplot(aes(x = variable.value,
                   y = value)) +
        geom_bar(stat = "identity",
                 color = "black",
                 fill = "grey85") +
        scale_y_continuous(position = "right") +
        labs(x = cleaned_names[i],
             y = NULL) +
        theme(panel.background = element_blank(),
              panel.border = element_rect(fill = NA,
                                          color = "black",
                                          size = 0.8),
              panel.grid.major = element_line(color = "grey90"),
              axis.title = element_text(size = 9),
              axis.text = element_text(size = 7))
    })
    PLTS <- lapply(1:length(PLT2), function(j) {
      base_g <- ggplot_gtable(ggplot_build(PLT2[[j]]))
      overlay_g <- ggplot_gtable(ggplot_build(PLT[[j]]))
      plt_panel = c(subset(base_g$layout, name == "panel", se = t:r))
      pnl_ind = which(overlay_g$layout$name == "panel")
      leg_ind = which(overlay_g$layout$name == "axis-l")
      lab_ind <- which(overlay_g$layout$name == "ylab-l")
      final_grob = gtable_add_grob(base_g,
                                   overlay_g$grobs[[pnl_ind]],
                                   plt_panel$t,
                                   plt_panel$l,
                                   plt_panel$b,
                                   plt_panel$r, name = "a")
      final_grob <- gtable_add_cols(final_grob, widths = unit(20, "points"), pos = 0)
      final_grob <- gtable_add_cols(final_grob, widths = unit(20, "points"), pos = 0)
      final_grob = gtable_add_grob(final_grob,
                                   overlay_g$grobs[[leg_ind]],
                                   7, #plt_axis$t,
                                   3, #plt_axis$l,
                                   7, #plt_axis$b,
                                   3, #plt_axis$r,
                                   name = "b")
      final_grob <- gtable_add_grob(final_grob,
                                    overlay_g$grob[[lab_ind]],
                                    7,
                                    1,
                                    7,
                                    1,
                                    name = "c")
      final_grob
    })
    LAB1 <- cowplot::ggdraw() +
    cowplot::draw_label("Marginal Effect on Prediction",
                        angle = 90,
                        fontface = "bold",
                        size = 18)
    LAB2 <- cowplot::ggdraw() +
    cowplot::draw_label("Frequency",
                        angle = 270,
                        fontface = "bold",
                        size = 18)
    cowplot::plot_grid(LAB1, wrap_plots(PLTS), LAB2,
                       nrow = 1,
                       rel_widths = c(0.1, 1, 0.1))} else wrap_plots(PLT)}
