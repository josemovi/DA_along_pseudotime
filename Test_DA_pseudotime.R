
# Subset data fram by lineage
dfl1<-df_sample[df_sample$lineages == "Lineage2",]
# Define the range of pseudotime values
min_pt <- min(dfl1$pseudotime)
max_pt <- max(dfl1$pseudotime)

# Create a sequence of pseudotime values at intervals of 0.5 units
pseudotime_seq <- seq(min_pt, max_pt, by = 0.5)

# Initialize empty dataframe to store density values for control and PE samples
density_df <- data.frame(Pseudotime = pseudotime_seq,
                         Density_Control = rep(0, length(pseudotime_seq)),
                         Density_PE = rep(0, length(pseudotime_seq)),
                         Wilcox_Pvalue = rep(0, length(pseudotime_seq)))

# Loop through pseudotime sequence
for (i in seq_along(pseudotime_seq)) {
  # Subset data for each pseudotime interval
  subset_data <- dfl1[dfl1$pseudotime >= pseudotime_seq[i] & dfl1$pseudotime < pseudotime_seq[i] + 0.5, ]
  
  # Count the number of control and PE samples in the interval
  control_count <- sum(subset_data$conditions == "Control")
  PE_count <- sum(subset_data$conditions == "PE")
  
  # Update the density dataframe
  density_df$Density_Control[i] <- control_count
  density_df$Density_PE[i] <- PE_count
  
  # Check if both conditions have enough observations for the test
  if (control_count >= 3 & PE_count >= 3) {
    # Contingency table for Wilcoxon rank-sum test
    contingency_table <- table(subset_data$conditions, subset_data$Sample_id)
    
    # Extract counts for Control and PE conditions
    control_v <- as.vector(contingency_table["Control", ])
    pe_v <- as.vector(contingency_table["PE", ])
    
    # Remove zeros from the vectors
    control_v <- control_v[control_v != 0]
    pe_v <- pe_v[pe_v != 0]
    
    # Perform Wilcoxon rank-sum test
    test_result <- wilcox.test(pe_v, control_v)
    
    # Update Wilcox_Pvalue column in density dataframe
    density_df$Wilcox_Pvalue[i] <- test_result$p.value
  } else {
    # If there are not enough observations for the test, set p-value to 1
    density_df$Wilcox_Pvalue[i] <- 1
  }
}



# Determine which condition is more abundant in each interval
density_df$More_Abundant <- ifelse(density_df$Density_Control > density_df$Density_PE, "Control", "PE")

density_df$Difference <- density_df$Density_Control - density_df$Density_PE

# Display significant DA intervals
density_df[density_df$Wilcox_Pvalue <= 0.05,]
