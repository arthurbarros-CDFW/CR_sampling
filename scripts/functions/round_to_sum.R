round_to_sum <- function(proportions, target_sum) {
  # Calculate unrounded values
  unrounded <- target_sum * proportions
  
  # Get floor values and remainders
  floor_vals <- floor(unrounded)
  remainders <- unrounded - floor_vals
  
  # Initial sum
  current_sum <- sum(floor_vals)
  remaining <- target_sum - current_sum
  
  # Allocate remaining units to largest remainders
  if(remaining > 0) {
    # Get indices to increment (largest remainders)
    idx_to_increment <- order(remainders, decreasing = TRUE)[1:remaining]
    floor_vals[idx_to_increment] <- floor_vals[idx_to_increment] + 1
  }
  
  return(floor_vals)
}
