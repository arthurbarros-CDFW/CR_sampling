#archive scripts

p <- ggplot(tag_estimates_df, aes(x = factor(fish_id), y = mean_k)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    title = "Estimated unrecovered tags for each CWT recovery",
    subtitle = paste("Sampling fraction theta =", round(theta, 3)),
    x = "Fish ID",
    y = "Estimated unrecovered tags (k)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)