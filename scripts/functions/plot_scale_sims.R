plot_scale_sims <- function(results,N,n_replicates) {
  ages_results<-results$age_summary_stats
  totals_results<-results$totals_summary_stats
  
  #bias by scale sample size
  p1a <- ggplot(ages_results, aes(x = scales_n, y = mean_bias, color = factor(age))) +
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
    theme_bw()
  
  #MOE by scale sample size
  p2a <- ggplot(ages_results, aes(x = scales_n, y = mean_moe_prop, color = factor(age))) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_moe_prop - sd_moe_prop, 
                    ymax = mean_moe_prop + sd_moe_prop, 
                    fill = factor(age)), alpha = 0.2) +
    labs(title = "mean MOE by sample size",
         x = "Number of scale samples", 
         y = "mean MOE",
         color = "Age", fill = "Age") +
    theme_bw()
  
  #CI width (precision)
  p3a <- ggplot(ages_results, aes(x = scales_n, y = mean_ci_width_prop, color = factor(age))) +
    geom_line() +
    geom_point() +
    labs(title = "95% CI width (precision)",
         x = "Number of scale samples", 
         y = "Confidence interval width",
         color = "Age") +
    theme_bw()
  
  #number of iterations within target moe
  p4a <- ggplot(ages_results, aes(x = scales_n, y = pct_meeting_target_moe_prop, color = factor(age))) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(from=0,to=100,by=10))+
    scale_x_continuous(limits=c(min(ages_results$scales_n),max(ages_results$scales_n)),
                       breaks=seq(from=0,to=max(ages_results$scales_n),by=200))+
    labs(title = paste("% iterations meeting target MOE (",target_moe*100,"%)",sep=""),
         x = "Number of scale samples", 
         y = "% iterations meeting target MOE",
         color = "Age") +
    theme_bw()
  
  #p4 solo
  #number of iterations within target moe
  p4a_solo <- ggplot(ages_results, aes(x = scales_n, y = pct_meeting_target_moe_prop, color = factor(age))) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits=c(0,100),
                       breaks=seq(from=0,to=100,by=10))+
    scale_x_continuous(limits=c(min(ages_results$scales_n),max(ages_results$scales_n)),
                       breaks=seq(from=0,to=max(ages_results$scales_n),by=100))+
    labs(title = paste("% iterations meeting target MOE (",target_moe*100,"%)",sep=""),
         x = "Number of scale samples", 
         y = "% iterations meeting target MOE",
         color = "Age") +
    theme_bw()
  
  main_title <- grid::textGrob(paste("N=",N,"boot iterations",n_replicates), 
                         gp = grid::gpar(fontsize = 16, fontface = "bold"))
  
  p_ages<-gridExtra::grid.arrange(p1a, p2a, p3a, p4a, ncol = 2,
                                  top = gridExtra::arrangeGrob(main_title,
                                                               nrow = 1
                                                               , heights = unit(c(1),"cm")))
  
  
  #bias by scale sample size
  p1 <- ggplot(totals_results, aes(x = scales_n, y = mean_bias)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = mean_bias - sd_bias, 
                    ymax = mean_bias + sd_bias), alpha = 0.2) +
    labs(title = "Bias in natural proportion estimates",
         x = "Number of scale samples", 
         y = "Mean bias") +
    theme_bw()
  
  #MOE by scale sample size
  p2 <- ggplot(totals_results, aes(x = scales_n, y = mean_moe_prop)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = mean_moe_prop - sd_moe_prop, 
                    ymax = mean_moe_prop + sd_moe_prop), alpha = 0.2) +
    labs(title = "mean MOE by sample size",
         x = "Number of scale samples", 
         y = "mean MOE") +
    theme_bw()
  
  #CI width (precision)
  p3 <- ggplot(totals_results, aes(x = scales_n, y = mean_ci_width_prop)) +
    geom_line() +
    geom_point() +
    labs(title = "95% CI width (precision)",
         x = "Number of scale samples", 
         y = "Confidence interval width") +
    theme_bw()
  
  #% number of iterations within target moe
  p4 <- ggplot(totals_results, aes(x = scales_n, y = pct_meeting_target_moe_prop)) +
    geom_line() +
    geom_point() +
    labs(title = paste("% iterations meeting target MOE (",target_moe*100,"%)",sep=""),
         x = "Number of scale samples", 
         y = "% iterations meeting target MOE") +
    theme_bw()
  
  p_ages<-gridExtra::grid.arrange(p1a, p2a, p3a, p4a, ncol = 2)
  p_totals<-gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  return(list("ages_plot"=p_ages,
              "totals_plot"=p_totals))
  
}
