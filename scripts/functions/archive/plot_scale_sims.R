plot_scale_sims <- function(results) {
  
  #bias by scale sample size
  p1 <- ggplot(results$summary, aes(x = scales_n, y = mean_bias, color = factor(age))) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = mean_bias - sd_bias, 
                    ymax = mean_bias + sd_bias, 
                    fill = factor(age)), alpha = 0.2) +
    labs(title = "Bias in natural proportion estimates",
         x = "Number of scale samples", 
         y = "Mean bias",
         color = "Age", fill = "Age") +
    theme_minimal()
  
  #RMSE by scale sample size
  p2 <- ggplot(results$summary, aes(x = scales_n, y = rmse, color = factor(age))) +
    geom_line() +
    geom_point() +
    labs(title = "RMSE by sample size",
         x = "Number of scale samples", 
         y = "RMSE",
         color = "Age") +
    theme_minimal()
  
  #CI width (precision)
  p3 <- ggplot(results$summary, aes(x = scales_n, y = mean_ci_width, color = factor(age))) +
    geom_line() +
    geom_point() +
    labs(title = "95% CI width (precision)",
         x = "Number of scale samples", 
         y = "Confidence interval width",
         color = "Age") +
    theme_minimal()
  
  #coverage probability
  p4 <- ggplot(results$summary, aes(x = scales_n, y = coverage_rate, color = factor(age))) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    labs(title = "Coverage probability of 95% CIs",
         x = "Number of scale samples", 
         y = "Coverage rate",
         color = "Age") +
    ylim(0, 1) +
    theme_minimal()
  
  p_grid<-gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  ggsave(paste("outputs/sim_plots",Sys.Date(),".png"),p_grid,scale=2)
}
