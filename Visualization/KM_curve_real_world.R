# ==============================================================================
# Visualization: Kaplan-Meier Curves
# ==============================================================================
#
# This script generates Kaplan-Meier survival curves for real-world data
# with cluster assignments.
#
# Input:
#   - real_world_full_data.rds: Full clinical trial dataset
#   - real_world_mixed_data.rds: Mixed-type covariates
#   - [n]_[k]clus_cluster_labels.csv: Cluster assignments
#
# Output:
#   - KM survival plots
#   - Summary tables
#
# ==============================================================================

# Load required libraries
library(dplyr)     
library(survival)   
library(survminer)  
library(gtsummary)  
library(RColorBrewer)  

# ==============================================================================
# Load Data
# ==============================================================================

# Load real-world clinical trial data
real_world_full_data <- readRDS("real_world_full_data.rds")
real_world_mixed_data <- readRDS("real_world_mixed_data.rds")

# Load cluster assignments (adjust filename as needed)
# Example: "214_3clus_cluster_labels.csv" for 214 patients, 3 clusters
X214_3clus_cluster_labels <- read.csv("214_3clus_cluster_labels.csv")

# Extract PAM clustering results
pam <- X214_3clus_cluster_labels[
  X214_3clus_cluster_labels$method == "PAM",
]

# Merge cluster assignments with full data
data <- merge(real_world_full_data, pam, by = "id")

# ==============================================================================
# Prepare Data for KM Curves
# ==============================================================================

# Convert cluster to character for easier manipulation
data <- data %>%
  mutate(cluster = as.character(cluster))

# Create duplicate dataset with "Overall" label for overall survival curve
overall_data <- data %>%
  mutate(cluster = "Overall")

# Combine cluster-specific and overall data
combined_data <- bind_rows(data, overall_data)

# Set factor levels with desired order (Overall first, then clusters)
combined_data$cluster <- factor(
  combined_data$cluster,
  levels = c("Overall", "1", "2", "3"),
  labels = c("Overall", "Cluster 1", "Cluster 2", "Cluster 3")
)

# ==============================================================================
# Fit Kaplan-Meier Curves
# ==============================================================================

# Fit KM model stratified by cluster (including Overall)
fit_combined <- survfit(Surv(time, status) ~ cluster, data = combined_data)

# ==============================================================================
# Plot Kaplan-Meier Curves
# ==============================================================================

ggsurvplot(
  fit_combined,
  data = combined_data,
  risk.table = TRUE,        # Display number at risk table
  pval = TRUE,              # Show log-rank test p-value
  conf.int = TRUE,          # Display confidence intervals
  legend.title = "Group",
  # Color palette: black for Overall, Set1 for clusters
  palette = c("black", RColorBrewer::brewer.pal(3, "Set1")),
  title = "KM Curve in 214 Study",
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 10        # Time axis breaks every 10 units
)

# ==============================================================================
# Summary Table of Baseline Characteristics
# ==============================================================================

# Merge mixed data with cluster assignments for summary table
mixed_data <- merge(real_world_mixed_data, pam, by = "id")

# Remove "NOT REPORTED" level from IMDC if present
mixed_data$IMDC <- droplevels(
  mixed_data$IMDC[mixed_data$IMDC != "NOT REPORTED"]
)

# Create summary table of baseline characteristics by cluster
# Shows median (min - max) for continuous variables
data %>%
  dplyr::select(age, sex, region, IMDC, ECOG, PDL1, cluster) %>%
  tbl_summary(
    by = cluster,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})"
    )
  ) %>%
  add_overall()  # Add overall column
