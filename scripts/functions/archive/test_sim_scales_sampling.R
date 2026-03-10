#' @title sim_wrapper
#'
#' @description fill
#' 
#' @param N
#' 
#' @param theta 
#' 
#' @param k_iterations 
#' 
#' @param tag_rate
#' 
#' @param survey "hatchery" or "trib"
#' 
#' @param scales_n
#' 
#' @return fill

test_sim_scales_sampling <- function(pop,
                                spawning_results,
                                theta,
                                k_iterations,
                                scale_iterations,
                                tag_rate,
                                survey,
                                scales_n,
                                CI) {
  
  # Convert pop to data.table once at the beginning if it's not already
  if(!data.table::is.data.table(pop)) {
    pop_dt <- as.data.table(pop)
  } else {
    pop_dt <- pop
  }
  
  # Get untagged recovered fish - use data.table for speed
  recovered_fish_dt <- as.data.table(spawning_results$recovered_fish)
  untagged_recovered <- recovered_fish_dt[has_cwt == 0]
  
  # total number of untagged recoveries
  untagged_total <- nrow(untagged_recovered)
  
  time_start <- Sys.time()
  
  # Handle scale sampling
  if(untagged_total > scales_n) {
    # Use fast random sampling with data.table
    scale_samples <- untagged_recovered[sample(.N, scales_n, replace = FALSE)]
  } else {
    scale_samples <- untagged_recovered
    scales_n <- untagged_total
  }
  
  # Pre-compute population total once
  pop_total <- nrow(pop_dt)
  
  # Process hatchery ages - use data.table and split efficiently
  hatchery_ages_dt <- as.data.table(spawning_results$hatchery_ages)
  
  # If hatchery_ages_dt is empty, create template once
  if(nrow(hatchery_ages_dt) == 0) {
    empty_tagged_template <- data.table(
      age = 2:4,
      total_tags = 0,
      total_hatchery = 0,
      N = 0,
      k_iteration = 1:scale_iterations
    )
    # Replicate for all iterations efficiently
    hatchery_ages_split <- replicate(scale_iterations, 
                                     empty_tagged_template[age %in% 2:4, 
                                                           .(age, total_tags, total_hatchery, N)], 
                                     simplify = FALSE)
  } else {
    # Split by k_iteration - ensure all iterations present
    all_iterations <- 1:scale_iterations
    hatchery_ages_split <- vector("list", scale_iterations)
    
    # Use fast subsetting with keys
    setkey(hatchery_ages_dt, k_iteration)
    for(j in all_iterations) {
      # Get data for this iteration, or create empty if missing
      iter_data <- hatchery_ages_dt[.(j)]
      if(nrow(iter_data) == 0) {
        hatchery_ages_split[[j]] <- data.table(
          age = 2:4,
          total_tags = 0,
          total_hatchery = 0,
          N = 0
        )
      } else {
        hatchery_ages_split[[j]] <- iter_data[, .(age, total_tags, total_hatchery, N)]
      }
    }
  }
  
  # Pre-allocate results list
  boot_results <- vector("list", scale_iterations)
  
  # Pre-compute all ages template
  all_ages <- data.table(age = 2:4)
  
  # Main bootstrap loop - optimized
  for(j in 1:scale_iterations) {
    # Bootstrap sampling - much faster with integer indices
    sample_idx <- sample.int(scales_n, scales_n, replace = TRUE)
    scales_boot <- scale_samples[sample_idx]
    
    # Get tagged ages for this iteration (already pre-split)
    tagged_ages_boot <- hatchery_ages_split[[j]]
    
    # Calculate totals once
    total_tagged_boot <- sum(tagged_ages_boot$total_tags, na.rm = TRUE)
    untagged_total_boot <- pop_total - total_tagged_boot
    
    # Process untagged scales - optimized aggregation
    if(nrow(scales_boot) > 0) {
      # Fast counting by age
      untagged_counts <- scales_boot[, .(sample_untagged_count = .N), by = age]
      total_samples <- sum(untagged_counts$sample_untagged_count)
      
      # Calculate proportions and expand to all ages in one step
      untagged_scale_ages_complete <- all_ages[untagged_counts, on = "age"]
      untagged_scale_ages_complete[
        , `:=`(
          sample_untagged_count = nafill(sample_untagged_count, fill = 0),
          untagged_age_prop = nafill(sample_untagged_count / total_samples, fill = 0),
          untagged_age_total = nafill(sample_untagged_count / total_samples * untagged_total_boot, fill = 0)
        )
      ]
    } else {
      # No scales sampled
      untagged_scale_ages_complete <- copy(all_ages)
      untagged_scale_ages_complete[
        , `:=`(
          sample_untagged_count = 0,
          untagged_age_prop = 0,
          untagged_age_total = 0
        )
      ]
    }
    
    # Fast merge and calculations
    age_counts <- untagged_scale_ages_complete[tagged_ages_boot, on = "age"]
    
    # Fill NAs efficiently
    for(col in c("total_tags", "total_hatchery", "N")) {
      if(col %in% names(age_counts)) {
        set(age_counts, which(is.na(age_counts[[col]])), col, 0)
      }
    }
    
    # Vectorized calculations
    age_counts[, `:=`(
      untagged_hatchery = total_hatchery - total_tags,
      total_natural = pmax(untagged_age_total - (total_hatchery - total_tags), 0)
    )]
    
    # Calculate proportions
    total_natural_sum <- sum(age_counts$total_natural)
    if(total_natural_sum > 0) {
      age_counts[, natural_proportions := total_natural / total_natural_sum]
    } else {
      age_counts[, natural_proportions := 0]
    }
    
    # Store results - keep as data.table
    boot_results[[j]] <- age_counts[, .(age, total_hatchery, total_natural, natural_proportions)]
  }
  
  time_end <- Sys.time()
  boot_time <- time_end - time_start
  
  # Combine all bootstrap results efficiently
  all_boot_dt <- rbindlist(boot_results, idcol = "iteration")
  
  # Calculate age-specific summary statistics - optimized
  age_summary_stats <- all_boot_dt[
    !is.na(natural_proportions),
    .(
      est_count_hatchery = mean(total_hatchery),
      est_count_natural = mean(total_natural),
      mean_proportion_natural = mean(natural_proportions),
      median_proportion_natural = median(natural_proportions),
      sd_proportion_natural = sd(natural_proportions),
      se_proportion_natural = sd(natural_proportions) / sqrt(.N),
      lower_CI_prop = quantile(natural_proportions, probs = (1 - CI) / 2),
      upper_CI_prop = quantile(natural_proportions, probs = 1 - (1 - CI) / 2),
      lower_CI_count = quantile(total_natural, probs = (1 - CI) / 2),
      upper_CI_count = quantile(total_natural, probs = 1 - (1 - CI) / 2),
      n_iterations = .N
    ),
    by = age
  ][order(age)]
  
  # Calculate overall statistics - optimized
  iteration_overall <- all_boot_dt[
    !is.na(natural_proportions),
    .(
      total_hatchery_iter = sum(total_hatchery),
      total_natural_iter = sum(total_natural)
    ),
    by = iteration
  ][
    , overall_prop_natural := fifelse(
      total_hatchery_iter + total_natural_iter > 0,
      total_natural_iter / (total_hatchery_iter + total_natural_iter),
      0
    )
  ]
  
  overall_summary_stats <- iteration_overall[
    , .(
      est_total_hatchery = mean(total_hatchery_iter),
      est_total_natural = mean(total_natural_iter),
      mean_proportion_natural = mean(overall_prop_natural),
      median_proportion_natural = median(overall_prop_natural),
      sd_proportion_natural = sd(overall_prop_natural),
      se_proportion_natural = sd(overall_prop_natural) / sqrt(.N),
      lower_CI_prop = quantile(overall_prop_natural, probs = (1 - CI) / 2),
      upper_CI_prop = quantile(overall_prop_natural, probs = 1 - (1 - CI) / 2),
      lower_CI_total = quantile(total_natural_iter, probs = (1 - CI) / 2),
      upper_CI_total = quantile(total_natural_iter, probs = 1 - (1 - CI) / 2),
      cv_proportion = sd(overall_prop_natural) / mean(overall_prop_natural),
      n_iterations = .N
    )
  ]
  
  # Calculate true population statistics - optimized with data.table
  true_summary_stats <- pop_dt[, .(true_count = .N), by = .(age, origin)]
  
  # Complete missing combinations efficiently
  all_combinations <- CJ(age = 2:4, origin = c("hatchery", "natural"))
  true_summary_stats <- true_summary_stats[all_combinations, on = .(age, origin)]
  true_summary_stats[is.na(true_count), true_count := 0]
  
  # Pivot wider efficiently
  true_summary_stats <- dcast(true_summary_stats, age ~ origin, value.var = "true_count")
  setnames(true_summary_stats, c("hatchery", "natural"), c("true_count_hatchery", "true_count_natural"))
  
  # Calculate proportions
  total_hatchery_val <- sum(true_summary_stats$true_count_hatchery)
  total_natural_val <- sum(true_summary_stats$true_count_natural)
  
  true_summary_stats[, `:=`(
    true_proportion_hatchery = if (total_hatchery_val > 0) true_count_hatchery / total_hatchery_val else 0,
    true_proportion_natural = if (total_natural_val > 0) true_count_natural / total_natural_val else 0
  )]
  
  results <- list(
    "age_summary_stats" = age_summary_stats,
    "true_summary_stats" = true_summary_stats,
    "overall_summary_stats" = overall_summary_stats,
    "scales_collected" = scales_n
  )
  
  return(results)
}
