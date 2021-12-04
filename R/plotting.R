append_avg = function(res, col) {
  avg_res = res %>% 
    # group_by(iter, method) %>% summarise("{colname}" := mean(.data[[col]])) %>%
    group_by(iter, method) %>% summarise("{{col}}" := mean({{col}})) %>% 
    add_column(rep = "avg", alpha = 1, size = 1)
  
  res %>% select(iter, method, rep, {{col}}) %>% 
    add_column(alpha = .5, size = 0.5) %>% 
    rbind(avg_res) %>% 
    mutate(rep = factor(rep)) # %>% mutate(alpha = factor(alpha))
}

plot_paths_and_avg = function(res, col, alpha_range = c(0.25,1)) {
  res %>% 
    append_avg({{col}}) %>% 
    ggplot(aes(x = iter, y = {{col}}, color = method, alpha = alpha)) +
    # scale_alpha_discrete(range = c(0.2, .9), guide="none") + 
    # geom_line(aes(size = method), alpha = 0.5) +
    geom_line(aes(group = interaction(rep, method), size = size)) + 
    scale_alpha( range = alpha_range, guide = 'none') +
    scale_size( range = c(0.5,1.2), guide = 'none') + 
    theme_minimal() +
    # scale_colour_manual(values = c(1,1.5)) +
    ggplot2::theme(
      plot.title = element_text(size=10),
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.8, 0.2),
      # legend.text = ggplot2::element_text(size=18),
    ) + 
    ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
    xlab("Iteration") 
}
