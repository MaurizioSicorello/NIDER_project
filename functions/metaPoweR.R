metaPower <- function(effectSize, sampleSizes, tau = 0, alpha = 0.05, numTruePositives = 1) {
  
  # Fisher transformation
  effectSize_z <- atanh(effectSize)
  
  # Calculate standard error
  # Variances
  variances <- 1 / (sampleSizes - 3)
  
  # Fixed effect weights
  weights <- 1 / variances
  
  # Random effect weights
  weights_star <- 1 / (variances + tau^2)
  
  # Standard error
  SE_random <- sqrt(1 / sum(weights_star))
  SE_fixed <- sqrt(1 / sum(weights))
  
  # Non-centrality parameter
  z_random <- effectSize_z / SE_random
  z_fixed <- effectSize_z / SE_fixed
  
  # Power
  # Critical z-value
  z_crit <- qnorm(1 - alpha / 2)
  
  randomPower <- 1 + pnorm(-z_crit + z_random) - pnorm(z_crit + z_random)
  fixedPower <- 1 + pnorm(-z_crit + z_fixed) - pnorm(z_crit + z_fixed)
  
  # Probability of at least 1 true positive test
  randomPower <- 1 - (1 - randomPower)^numTruePositives
  fixedPower <- 1 - (1 - fixedPower)^numTruePositives
  
  # Return results
  results <- list(random = randomPower, fixed = fixedPower)
  
  return(results)
}
