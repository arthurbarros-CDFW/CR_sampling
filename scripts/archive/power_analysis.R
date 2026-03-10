#==============================================================================
# Power Analysis for Estimating Proportion of Natural-Origin Fish
#==============================================================================

# Load required packages
library(pwr)
library(tidyverse)

#------------------------------------------------------------------------------
# 1. Binomial Sample Size Formula (Simple Case)
#------------------------------------------------------------------------------

calculate_sample_size_binomial <- function(p, margin_error, confidence = 0.95) {
  # p: expected proportion of natural-origin fish (0-1)
  # margin_error: desired margin of error (e.g., 0.03 for ±3%)
  # confidence: confidence level (e.g., 0.95 for 95%)
  
  z <- qnorm(1 - (1 - confidence) / 2)  # Z-score for confidence level
  n <- (z^2 * p * (1 - p)) / (margin_error^2)
  return(ceiling(n))
}

# Example: For p (proportion) = 0.05 (5%) and ±3% margin of error
sample_size <- calculate_sample_size_binomial(p = 0.5, margin_error = 0.05, confidence = 0.95)
cat("Required sample size (binomial formula):", sample_size, "\n")

#------------------------------------------------------------------------------
# 2. Power Analysis for Detecting a Difference from a Threshold
#------------------------------------------------------------------------------

# If you want to test whether natural-origin proportion is > a certain threshold
# (e.g., >2% when management requires monitoring of low proportions)

# p1: proportion under H0 (null hypothesis, e.g., 0.02)
# p2: proportion under H1 (alternative, e.g., 0.05)
# power: desired power (e.g., 0.8)
# sig.level: significance level (e.g., 0.05)

# Calculate effect size h manually
h <- 2 * asin(sqrt(0.05)) - 2 * asin(sqrt(0.02))
cat("Effect size h =", h, "\n")

#------------------------------------------------------------------------------
# 3. Simulation-Based Approach (More Realistic)
#------------------------------------------------------------------------------

# This accounts for:
# - Finite population size
# - Sampling variability
# - Potential bias in aging/marking

simulate_sample_size <- function(true_prop = 0.05, 
                                 population_size = 1000,
                                 margin_error = 0.03,
                                 conf_level = 0.95,
                                 n_simulations = 1000) {
  
  # Vector of possible sample sizes to test
  sample_sizes <- seq(50, 500, by = 50)
  results <- data.frame()
  
  for (n in sample_sizes) {
    coverages <- numeric(n_simulations)
    
    for (sim in 1:n_simulations) {
      # Simulate population
      pop <- c(rep(1, round(true_prop * population_size)),  # 1 = natural-origin
               rep(0, population_size - round(true_prop * population_size)))  # 0 = hatchery
      
      # Take sample
      sample_indices <- sample(1:population_size, n, replace = FALSE)
      sample_fish <- pop[sample_indices]
      
      # Calculate sample proportion and confidence interval
      p_hat <- mean(sample_fish)
      se <- sqrt(p_hat * (1 - p_hat) / n)
      z <- qnorm(1 - (1 - conf_level) / 2)
      ci_lower <- p_hat - z * se
      ci_upper <- p_hat + z * se
      
      # Check if true proportion is within CI
      coverages[sim] <- (ci_lower <= true_prop) & (true_prop <= ci_upper)
    }
    
    # Calculate coverage probability
    coverage_prob <- mean(coverages)
    results <- rbind(results, data.frame(
      sample_size = n,
      coverage = coverage_prob,
      margin_achieved = z * sqrt((true_prop * (1 - true_prop)) / n)
    ))
  }
  
  return(results)
}

# Run simulation
#set.seed(123)
sim_results <- simulate_sample_size(true_prop = 0.05, 
                                    population_size = 2000,
                                    margin_error = 0.03,
                                    conf_level = 0.95,
                                    n_simulations = 10000)

# Find minimum sample size achieving at least 95% coverage
adequate_samples <- sim_results %>%
  filter(coverage >= 0.95)

min_sample <- ifelse(nrow(adequate_samples) > 0, 
                     min(adequate_samples$sample_size),
                     NA)

cat("\nMinimum sample size from simulation:", min_sample, "\n")

#------------------------------------------------------------------------------
# 4. Create Sample Size Table for Different Scenarios
#------------------------------------------------------------------------------

scenarios <- expand.grid(
  prop = c(0.01, 0.05, 0.10, 0.20),  # Natural-origin proportions
  error = c(0.02, 0.03, 0.05),       # Margin of error (±%)
  conf = c(0.90, 0.95)               # Confidence level
)

scenarios$sample_size <- mapply(calculate_sample_size_binomial,
                                p = scenarios$prop,
                                margin_error = scenarios$error,
                                confidence = scenarios$conf)

# Format table
sample_size_table <- scenarios %>%
  pivot_wider(names_from = conf, 
              values_from = sample_size,
              names_prefix = "Conf_") %>%
  arrange(prop, error)

print("Sample Size Recommendations:")
print(sample_size_table)

#------------------------------------------------------------------------------
# 5. Visualization
#------------------------------------------------------------------------------

# Plot sample size vs. proportion for different margins of error
plot_data <- scenarios %>%
  filter(conf == 0.95)

ggplot(plot_data, aes(x = prop, y = sample_size, color = as.factor(error))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "Sample Size Requirements for Estimating Natural-Origin Proportion",
       x = "True Proportion of Natural-Origin Fish",
       y = "Required Sample Size",
       color = "Margin of Error (±)") +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save plot
ggsave("sample_size_power_analysis.png", width = 8, height = 6, dpi = 300)

#------------------------------------------------------------------------------
# 6. Adjust for Finite Population (if tributary population is small)
#------------------------------------------------------------------------------

adjust_for_finite_population <- function(n, N) {
  # n: sample size from infinite population formula
  # N: total number of spawners in the tributary
  n_adj <- n / (1 + (n - 1) / N)
  return(ceiling(n_adj))
}

# Example: If formula says n=400 but tributary has only 2000 spawners
n_adj <- adjust_for_finite_population(n = 400, N = 2000)
cat("Adjusted sample size for finite population (N=2000):", n_adj, "\n")

#==============================================================================
# Output Summary
#==============================================================================

cat("\n=== SAMPLE SIZE RECOMMENDATIONS ===\n")
cat("For tributaries with 5% natural-origin fish:\n")
cat("- To estimate within ±3% with 95% confidence:", 
    calculate_sample_size_binomial(0.05, 0.03, 0.95), "fish\n")
cat("- To detect increase from 2% to 5% with 80% power:", 
    ceiling(power_analysis$n), "fish\n")
cat("\nConsiderations:\n")
cat("1. Increase sample size if aging/marking error is high\n")
cat("2. Use stratified sampling if fish distribution is uneven\n")
cat("3. Monitor annually and adjust sampling based on prior years' estimates\n")

# Save results
write.csv(sample_size_table, "sample_size_recommendations.csv", row.names = FALSE)