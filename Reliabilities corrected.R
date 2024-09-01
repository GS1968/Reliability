# Install necessary packages if they are not already installed
required_packages <- c("psych", "lavaan", "GPArotation", "MBESS", "boot", "ggplot2", "tidyr", "ggpattern")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load necessary packages
library(psych)
library(lavaan)
library(GPArotation)
library(MBESS)
library(boot)
library(ggplot2)
library(tidyr)
library(ggpattern) # For adding patterns to the bars

# Avoid scientific notation
options(scipen = 999)

# Function to calculate reliability metrics for a given dataset
calculate_reliability_metrics <- function(data) {
  # Calculate Cronbach's alpha
  alpha <- tryCatch({
    psych::alpha(data)$total$raw_alpha
  }, warning = function(w) {
    NA
  }, error = function(e) {
    NA
  })
  
  # Calculate Omega total
  omega_total <- tryCatch({
    fa_result <- psych::fa(data, nfactors = 1, fm = "ml")
    lambda <- fa_result$loadings
    psi <- diag(fa_result$residual)
    omega_tot <- sum(lambda)^2 / (sum(lambda)^2 + sum(psi))
    omega_tot
  }, warning = function(w) {
    NA
  }, error = function(e) {
    NA
  })
  
  # Calculate Greatest Lower Bound (GLB)
  glb <- tryCatch({
    psych::glb.fa(data)$glb
  }, error = function(e) {
    NA
  })
  
  # Calculate maximal reliability using TLI from factor analysis
  maximal_reliability <- tryCatch({
    fa_result <- psych::fa(data, nfactors = 1, fm = "ml")
    maximal_reliability <- fa_result$TLI
    maximal_reliability
  }, warning = function(w) {
    NA
  }, error = function(e) {
    NA
  })
  
  return(list(alpha = alpha, omega_total = omega_total, glb = glb, maximal_reliability = maximal_reliability))
}

# Bootstrapping function for reliability metrics
bootstrap_reliability <- function(data, R = 1000) {
  boot_results <- boot(data, statistic = function(data, indices) {
    sample_data <- data[indices, ]
    metrics <- calculate_reliability_metrics(sample_data)
    # Return NA for any component if it fails during calculation
    sapply(metrics, function(x) ifelse(is.na(x), NA, x))
  }, R = R)
  
  # Calculate means of bootstrapped results, ignoring NA values
  boot_means <- apply(boot_results$t, 2, function(x) mean(x, na.rm = TRUE))
  
  return(list(
    alpha = boot_means[1],
    omega_total = boot_means[2],
    glb = boot_means[3],
    maximal_reliability = boot_means[4]
  ))
}

# Main function
calculate_reliability_corrected <- function(file_path, R = 1000) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Convert data to matrix
  data <- as.matrix(data)
  
  # Check that data is numeric
  if(!is.numeric(data)) {
    stop("Data must be numeric.")
  }
  
  # Calculate reliability metrics for the original sample
  sample_metrics <- calculate_reliability_metrics(data)
  
  # Bootstrap reliability metrics
  boot_metrics <- bootstrap_reliability(data, R)
  
  # Calculate corrected reliability coefficients for the original sample
  corrected_alpha_sample <- if(!is.na(sample_metrics$alpha) && sample_metrics$alpha > 0) sample_metrics$alpha / sqrt(sample_metrics$alpha) else NA
  corrected_omega_sample <- if(!is.na(sample_metrics$omega_total) && sample_metrics$omega_total > 0) sample_metrics$omega_total / sqrt(sample_metrics$omega_total) else NA
  corrected_glb_sample <- if(!is.na(sample_metrics$glb) && sample_metrics$glb > 0) sample_metrics$glb / sqrt(sample_metrics$glb) else NA
  corrected_maximal_sample <- if(!is.na(sample_metrics$maximal_reliability) && sample_metrics$maximal_reliability > 0) sample_metrics$maximal_reliability / sqrt(sample_metrics$maximal_reliability) else NA
  
  # Calculate corrected reliability coefficients for the bootstrapped means
  corrected_alpha_boot <- if(!is.na(boot_metrics$alpha) && boot_metrics$alpha > 0) boot_metrics$alpha / sqrt(boot_metrics$alpha) else NA
  corrected_omega_boot <- if(!is.na(boot_metrics$omega_total) && boot_metrics$omega_total > 0) boot_metrics$omega_total / sqrt(boot_metrics$omega_total) else NA
  corrected_glb_boot <- if(!is.na(boot_metrics$glb) && boot_metrics$glb > 0) boot_metrics$glb / sqrt(boot_metrics$glb) else NA
  corrected_maximal_boot <- if(!is.na(boot_metrics$maximal_reliability) && boot_metrics$maximal_reliability > 0) boot_metrics$maximal_reliability / sqrt(boot_metrics$maximal_reliability) else NA
  
  # Create a data frame to present the results in a table
  results_table <- data.frame(
    Reliability_Coefficient = c("Cronbach's Alpha", "Omega Total", "GLB", "Maximal Reliability"),
    Sample_Based = c(sample_metrics$alpha, sample_metrics$omega_total, sample_metrics$glb, sample_metrics$maximal_reliability),
    Bootstrap_Mean_Based = c(boot_metrics$alpha, boot_metrics$omega_total, boot_metrics$glb, boot_metrics$maximal_reliability),
    Corrected_for_Attenuation = c(corrected_alpha_sample, corrected_omega_sample, corrected_glb_sample, corrected_maximal_sample),
    Corrected_Bootstrap_Attenuation = c(corrected_alpha_boot, corrected_omega_boot, corrected_glb_boot, corrected_maximal_boot)
  )
  
  # Print the results table
  print(results_table)
  
  # Create a plot to visualize the reliability estimates
  results_long <- tidyr::pivot_longer(results_table, cols = -Reliability_Coefficient, names_to = "Method", values_to = "Estimate")
  
  # Reorder the Method factor levels to desired order in the plot
  results_long$Method <- factor(results_long$Method, levels = c("Sample_Based", "Bootstrap_Mean_Based", "Corrected_for_Attenuation", "Corrected_Bootstrap_Attenuation"))
  
  reliability_plot <- ggplot(results_long, aes(x = Reliability_Coefficient, y = Estimate, fill = Method, pattern = Method)) +
    geom_bar_pattern(stat = "identity", position = position_dodge(), color = "black", 
                     pattern_density = 0.1, pattern_fill = "black", pattern_spacing = 0.05) +
    geom_text(aes(label = round(Estimate, 2)), vjust = -0.3, position = position_dodge(0.9), size = 4.5) +
    scale_fill_manual(values = c("Sample_Based" = "black", 
                                 "Bootstrap_Mean_Based" = "orange", 
                                 "Corrected_for_Attenuation" = "lightgreen", 
                                 "Corrected_Bootstrap_Attenuation" = "grey")) +
    scale_pattern_manual(values = c("Sample_Based" = "none", 
                                    "Bootstrap_Mean_Based" = "stripe", 
                                    "Corrected_for_Attenuation" = "crosshatch", 
                                    "Corrected_Bootstrap_Attenuation" = "circle")) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "red", size = 1) +
    labs(title = "Reliability Estimates Comparison",
         x = "Reliability Coefficient",
         y = "Estimate") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Enlarge title font
          axis.title = element_text(size = 14),  # Enlarge axis titles
          axis.text = element_text(size = 14),  # Enlarge axis text
          legend.text = element_text(size = 14))  # Enlarge legend text
  
  print(reliability_plot)
  
  return(results_table)
}

# Example usage with seed for replication:
#set.seed(123)  # Set seed for reproducibility
#file_path <- "data.csv"
#reliability_results <- calculate_reliability_corrected(file_path, R = 1000)
