#==============================================================================
# Sampling for Detection of Rare Natural-Origin Fish (<2%)
#==============================================================================

# Function to calculate detection probability
detection_probability <- function(p, n) {
  # p: true proportion of natural-origin fish (0-1)
  # n: sample size
  return(1 - (1 - p)^n)
}

# Function to find sample size for target detection probability
find_sample_size <- function(p, target_prob = 0.95) {
  # Find n such that detection probability ≥ target_prob
  n <- 1
  while (detection_probability(p, n) < target_prob) {
    n <- n + 1
  }
  return(n)
}

#------------------------------------------------------------------------------
# 1. Calculate for different proportions
#------------------------------------------------------------------------------

proportions <- c(0.02, 0.015, 0.01, 0.005, 0.002)  # 2%, 1.5%, 1%, 0.5%, 0.2%

results <- data.frame(
  True_Proportion = proportions,
  Proportion_Percent = proportions * 100,
  Sample_Size_95 = sapply(proportions, find_sample_size, target_prob = 0.95),
  Sample_Size_80 = sapply(proportions, find_sample_size, target_prob = 0.80),
  Detection_Prob_150 = sapply(proportions, function(p) detection_probability(p, 150)),
  Detection_Prob_300 = sapply(proportions, function(p) detection_probability(p, 300))
)

print("Sample Size Requirements for Detection:")
print(results)

#------------------------------------------------------------------------------
# 2. Create detection probability curves
#------------------------------------------------------------------------------

library(ggplot2)

# Generate data for plotting
plot_data <- expand.grid(
  p = c(0.005, 0.01, 0.02),  # 0.5%, 1%, 2%
  n = seq(1, 500, by = 5)
)

plot_data$detection_prob <- mapply(detection_probability, 
                                   p = plot_data$p, 
                                   n = plot_data$n)

# Plot
ggplot(plot_data, aes(x = n, y = detection_prob, color = as.factor(p*100))) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(title = "Probability of Detecting Natural-Origin Fish",
       subtitle = "For tributaries with <2% natural-origin fish",
       x = "Sample Size (number of fish)",
       y = "Probability of Detecting ≥1 Natural-Origin Fish",
       color = "True Proportion (%)") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(values = c("0.5" = "red", "1" = "blue", "2" = "darkgreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save plot
ggsave("rare_fish_detection_curve.png", width = 8, height = 6, dpi = 300)

#------------------------------------------------------------------------------
# 3. Rule of Three Analysis
#------------------------------------------------------------------------------

rule_of_three <- function(n) {
  # If you sample n fish and find 0 natural-origin,
  # 95% CI upper bound is approximately 3/n
  upper_bound <- 3/n
  return(upper_bound)
}

# Create table
rule_table <- data.frame(
  Sample_Size = c(50, 100, 150, 200, 300),
  Upper_Bound_95 = sapply(c(50, 100, 150, 200, 300), rule_of_three),
  Interpretation = c(
    "≤6% if 0 detected",
    "≤3% if 0 detected", 
    "≤2% if 0 detected",
    "≤1.5% if 0 detected",
    "≤1% if 0 detected"
  )
)

print("\nRule of Three Analysis (if NO natural-origin fish detected):")
print(rule_table)

#------------------------------------------------------------------------------
# 4. Practical sampling recommendations
#------------------------------------------------------------------------------

cat("\n=== PRACTICAL SAMPLING RECOMMENDATIONS ===\n")
cat("For tributaries with <2% natural-origin fish:\n\n")

cat("SCENARIO 1: Presence/Absence Monitoring\n")
cat("  • Sample 150 fish annually\n")
cat("  • If 0 natural-origin detected → ≤2% likely (95% confidence)\n")
cat("  • Cost-effective for routine monitoring\n\n")

cat("SCENARIO 2: Positive Detection Required\n")
cat("  • Sample 300 fish annually\n")
cat("  • 95% chance to detect if true proportion = 1%\n")
cat("  • 78% chance to detect if true proportion = 0.5%\n\n")

cat("SCENARIO 3: High-Certainty Detection\n")
cat("  • Sample 600 fish (or 300 for 2 years)\n")
cat("  • 95% chance to detect if true proportion = 0.5%\n")
cat("  • Consider multi-year sampling to reduce annual effort\n")

#------------------------------------------------------------------------------
# 5. Multi-year sampling strategy
#------------------------------------------------------------------------------

multi_year_detection <- function(p, n_per_year, years) {
  # Probability of detection over multiple years
  annual_miss_prob <- (1 - p)^n_per_year
  multi_year_miss_prob <- annual_miss_prob^years
  return(1 - multi_year_miss_prob)
}

multi_year_results <- data.frame(
  Annual_Sample = c(100, 100, 150, 150),
  Years = c(2, 3, 2, 3),
  P_1percent = c(
    multi_year_detection(0.01, 100, 2),
    multi_year_detection(0.01, 100, 3),
    multi_year_detection(0.01, 150, 2),
    multi_year_detection(0.01, 150, 3)
  ),
  P_0.5percent = c(
    multi_year_detection(0.005, 100, 2),
    multi_year_detection(0.005, 100, 3),
    multi_year_detection(0.005, 150, 2),
    multi_year_detection(0.005, 150, 3)
  )
)

print("\nMulti-Year Sampling Strategies:")
print(multi_year_results)

#==============================================================================
# Summary Table for Management
#==============================================================================

summary_table <- data.frame(
  Objective = c(
    "Confirm ≤2% if absent",
    "Detect 2% population",
    "Detect 1% population", 
    "Detect 0.5% population"
  ),
  Confidence_Level = c("95%", "95%", "95%", "95%"),
  Annual_Sample_Size = c(150, 150, 300, 600),
  Detection_Probability = c("NA (absence)", "95%", "95%", "95%"),
  Alternative_2yr_Plan = c("75/year", "75/year", "150/year", "300/year"),
  Notes = c(
    "Rule of Three upper bound",
    "Adequate for most monitoring",
    "For sensitive populations",
    "Extremely rare - consider genetics"
  )
)

print("\n=== MANAGEMENT DECISION TABLE ===")
print(summary_table)

# Save results
write.csv(results, "detection_sample_sizes.csv", row.names = FALSE)
write.csv(summary_table, "management_recommendations.csv", row.names = FALSE)