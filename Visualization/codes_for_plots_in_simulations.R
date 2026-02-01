# ==============================================================================
# Visualization: General Plotting Functions
# ==============================================================================
#
# This script contains code for generating visualizations including:
#   - Kaplan-Meier survival curves by cluster
#   - MDS (Multi-Dimensional Scaling) plots for cluster visualization
#   - Covariate distribution plots
#
# Usage:
#   - Load simulation or real-world results
#   - Generate KM curves and MDS plots
#   - Customize plots for publication
#
# ==============================================================================

# Load required libraries
library(survival)      
library(flexsurv)    
library(MASS)              
library(dplyr)         
library(survminer)    

# Clustering packages
library(cluster)        
library(clustMixType)  
library(DisimForMixed) 
library(mclust)         
library(factoextra)     

# ==============================================================================
# Load Simulation Results
# ==============================================================================

# Load simulation results 
all_results_s1_n_500 <- readRDS("all_results_s1_n_500.rds")

full_data <- all_results_s1_n_500[[1]]$full_data

# ==============================================================================
# Kaplan-Meier Survival Curves
# ==============================================================================

# Create survival object (time to event, event indicator)
surv_obj <- Surv(time = full_data$time, event = full_data$status)

# Fit Kaplan-Meier curves stratified by true cluster
fit <- survfit(surv_obj ~ cluster, data = full_data)

# Plot using ggsurvplot (publication-ready plots)
ggsurvplot(
  fit,
  data = full_data,
  risk.table = TRUE,        # Show number at risk table
  pval = TRUE,              # Display log-rank test p-value
  conf.int = TRUE,          # Show confidence intervals
  legend.title = "True Cluster",
  # Adjust legend labels based on number of clusters
  # For 2 clusters:
  legend.labs = c("Cluster 1", "Cluster 2"),
  # For 3 clusters, use:
  # legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  palette = "Set1",         # Color palette
  title = "Scenario 1",     # Plot title
  xlab = "Time",            # X-axis label
  ylab = "Survival Probability"  # Y-axis label
)

# ==============================================================================
# Multi-Dimensional Scaling (MDS) Plot
# ==============================================================================
# MDS reduces high-dimensional data to 2D for visualization
# Useful for visualizing cluster separation in covariate space

# Prepare mixed-type covariates for clustering
mixed_data <- full_data %>%
  dplyr::select(age, hemoglobin, tumor_burden, gender, race, kps, disease_sites) %>%
  mutate(across(c(gender, race, kps, disease_sites), as.factor))

# Compute Gower distance matrix for mixed-type data
# Gower distance handles both continuous and categorical variables
gower_dist <- daisy(
  mixed_data,
  metric = "gower"
)

# Perform classical MDS (2 dimensions)
mds_coords <- cmdscale(gower_dist, k = 2)

# Prepare data frame for plotting
mds_df <- as.data.frame(mds_coords)
colnames(mds_df) <- c("Dim1", "Dim2")

# Add cluster labels (assuming covs_all exists from simulation)
# Note: This requires covs_all to be available in the environment
# If not available, use: mds_df$cluster <- factor(full_data$cluster)
mds_df$cluster <- factor(full_data$cluster)

# Visualize clusters in 2D MDS space
library(ggplot2)
ggplot(mds_df, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 2) +
  labs(
    title = "Scenario 1 MDS Dimension Reduction",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2",
    color = "Cluster"
  ) +
  theme_minimal()
