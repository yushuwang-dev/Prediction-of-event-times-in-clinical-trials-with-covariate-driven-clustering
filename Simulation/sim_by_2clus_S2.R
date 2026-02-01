# ==============================================================================
# Simulation Study: Scenario 2 (2 Clusters, Exponential Distribution)
# ==============================================================================
#
# This script performs a simulation study to evaluate milestone event time
# prediction methods using covariate-driven clustering.
# Output:
#   - all_results_s2_v2[n]_n_[n_total].rds: Complete simulation results
#   - predict_timepoints_s2_v2[n]_n_[n_total].csv: Results by landmark time
#   - predict_IF_s2_v2[n]_n_[n_total].csv: Results by information fraction
#
# ==============================================================================

library(simsurv)
library(flexsurv)
library(dplyr)
library(cluster)
library(clustMixType)
library(kamila)
library(parallel)
library(future.apply)
library(MASS)

# Source functions
source("simulation_functions.R")

# ==============================================================================
# Simulation Parameters
# ==============================================================================

scenario <- "less"
n_total <- 500
n1 <- n2 <- n_total / 2
n_simulations <- 1000
FU_time <- 24
target_event <- 0.7 * n_total
n_cores <- detectCores()

clustering_methods <- c("None", "PAM", "Hierarchical", "K-prototypes", 
                        "Kamila", "Truth")
model_types <- c("Exponential", "Weibull")

# Hazard parameters: covariate effects (log-hazard ratios) and baseline hazards
betas1 <- c(age = 0.03, gender = 0.2, hemoglobin = -0.3, tumor_burden = 0.03)
betas2 <- c(age = 0.02, gender = 0.1, hemoglobin = -0.1, tumor_burden = 0.01)
lambda1 <- 0.15
gamma1 <- 1
lambda2 <- 0.05
gamma2 <- 1

# Generate cluster-specific patient covariates (Scenario 2)
simulate_covariates <- function(n, cluster_id) {
  if (cluster_id == 1) {
    data.frame(
      id = seq_len(n),
      cluster = cluster_id,
      age = rnorm(n, mean = 60, sd = 5),
      hemoglobin = rnorm(n, mean = 13, sd = 1),
      tumor_burden = rnorm(n, mean = 55, sd = 10),
      gender = rbinom(n, 1, 0.55),
      race = sample(c("White", "Black", "Other"), n, replace = TRUE, 
                   prob = c(0.6, 0.2, 0.2)),
      kps = sample(1:4, n, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2)),
      disease_sites = sample(1:5, n, replace = TRUE, 
                            prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
    )
  } else {
    data.frame(
      id = seq_len(n) + n1,
      cluster = cluster_id,
      age = rnorm(n, mean = 58, sd = 5),
      hemoglobin = rnorm(n, mean = 13.5, sd = 1),
      tumor_burden = rnorm(n, mean = 50, sd = 10),
      gender = rbinom(n, 1, 0.45),
      race = sample(c("White", "Black", "Other"), n, replace = TRUE,
                   prob = c(0.5, 0.25, 0.25)),
      kps = sample(1:4, n, replace = TRUE, prob = c(0.2, 0.2, 0.3, 0.3)),
      disease_sites = sample(1:5, n, replace = TRUE,
                            prob = c(0.25, 0.25, 0.2, 0.15, 0.15))
    )
  }
}


# Main simulation function: generate one complete iteration
sim_func <- function(i, clustering_methods, model_types, targetev, sims, n_clus) {
  if (n_clus == 2) {
    covs1 <- simulate_covariates(n1, 1)
    covs2 <- simulate_covariates(n2, 2)
    covs_all <- bind_rows(covs1, covs2)
    
    # Simulate survival times
    sim1 <- simsurv(dist = "exponential", lambdas = lambda1, 
                    betas = betas1, x = covs1)
    sim2 <- simsurv(dist = "exponential", lambdas = lambda2, 
                    betas = betas2, x = covs2)
    sim_all <- bind_rows(sim1, sim2)
    
  } else if (n_clus == 3) {
    covs1 <- simulate_covariates(n1, 1)
    covs2 <- simulate_covariates(n2, 2)
    covs3 <- simulate_covariates(n3, 3)
    covs_all <- bind_rows(covs1, covs2, covs3)
    
    sim1 <- simsurv(dist = "weibull", lambdas = lambda1, gammas = gamma1, betas = betas1,
                    x = covs1, maxt = 5)
    sim2 <- simsurv(dist = "weibull", lambdas = lambda2, gammas = gamma2, betas = betas2,
                    x = covs2, maxt = 5)
    sim3 <- simsurv(dist = "weibull", lambdas = lambda3, gammas = gamma3, betas = betas3,
                    x = covs3, maxt = 5)
    sim_all <- bind_rows(sim1, sim2, sim3)
  }
  
  
  # Apply dropout and administrative censoring
  drop_time <- rexp(n_total, 0.01)
  sim_all$dropout <- ifelse(drop_time < sim_all$eventtime, 1, 0)
  sim_all$status <- ifelse(sim_all$dropout == 1, 0, sim_all$status)
  sim_all$eventtime <- pmin(sim_all$eventtime, drop_time)
  sim_all$id <- 1:n_total
  colnames(sim_all)[2] <- "time"
  
  # Staggered entry and administrative censoring
  sim_all$entry_time <- runif(n_total, 0, 2)
  sim_all$ttf <- sim_all$entry_time + sim_all$time
  sim_all$status <- ifelse(sim_all$ttf > FU_time, 0, sim_all$status)
  sim_all$ttf <- ifelse(sim_all$ttf < FU_time, sim_all$ttf, FU_time)
  sim_all$time <- ifelse(sim_all$ttf < FU_time, sim_all$time, FU_time)
  
  
  # Prepare data for analysis
  full_data <- left_join(sim_all, covs_all, by = "id") %>%
    mutate(
      gender = factor(gender, labels = c("Female", "Male")),
      race = factor(race),
      kps = factor(kps, ordered = TRUE),
      disease_sites = factor(disease_sites, ordered = TRUE)
    )
  mixed_data <- full_data %>% 
    dplyr::select(age, hemoglobin, tumor_burden, gender, race, kps, disease_sites)
  sim_all <- merge(sim_all, covs_all[, 1:2], by = "id")
  full_data$sim_id <- i
  
  # Apply clustering methods and predict
  rolling_list <- list()
  rolling_predtime_list <- list()
  if_list <- list()
  if_predtimes_list <- list()
  
  for (cm in clustering_methods) {
    res <- tryCatch({
      rolling_result <- rolling_milestone_prediction(
        full_data = full_data,
        mixed_data = mixed_data,
        clustering_method = cm,
        modeltypes = model_types,
        criterion = "BIC",
        targetev = targetev,
        time_points = as.integer(seq.int(2, FU_time, length.out = 6)),
        sims = sims,
        n_clus = n_clus,
        sim_id = i)
      
      info_result <- info_frac_prediction(
        full_data = full_data,
        mixed_data = mixed_data,
        clustering_method = cm,
        modeltypes = model_types,
        criterion = "BIC",
        targetev = targetev,
        info_frac = c(0.4, 0.6, 0.75, 0.8),
        sims = sims,
        n_clus = n_clus,
        sim_id = i)
      
      rolling_list[[length(rolling_list) + 1]] <- rolling_result$pred_results
      rolling_predtime_list[[length(rolling_predtime_list) + 1]] <- rolling_result$pred_times
      if_list[[length(if_list) + 1]] <- info_result$pred_results
      if_predtimes_list[[length(if_list) + 1]] <- info_result$pred_times
    }, error = function(e) {
      message(paste("Error in", cm, ":", e$message))
      NULL
    })
  }
  
  list(
    rolling = do.call(rbind, rolling_list),
    rolling_predtimes = do.call(rbind, rolling_predtime_list),
    info = do.call(rbind, if_list),
    info_predtimes = do.call(rbind, if_predtimes_list),
    full_data = full_data
  )
}



# Run simulations in parallel
plan(multisession, workers = n_cores) 
all_results <- future_lapply(
  1:n_simulations, 
  sim_func,
                             clustering_methods = clustering_methods,
  model_types = model_types, 
  targetev = target_event, 
  sims = 1000, 
  n_clus = 2, 
  future.seed = TRUE
)

# Extract and combine results
rolling_results <- bind_rows(lapply(all_results, `[[`, "rolling"))
info_results <- bind_rows(lapply(all_results, `[[`, "info"))
saveRDS(all_results, paste0("all_results_s2_v2", scenario, "_n_", n_total, ".rds"))
rm(all_results)
gc()

# Process and save rolling landmark results
summary_rolling <- rolling_results %>%
  mutate(
    #bias      = as.numeric(bias),
    #abs_err   = as.numeric(abs_err),
    #mse       = as.numeric(mse),
    mean_rae  = as.numeric(mean_rae),
    mean_err  = as.numeric(mean_err),
    mean_ae   = as.numeric(mean_ae)
  )

# Remove simulations with infinite values if < 5% of total
inf_sim_df <- summary_rolling %>%
  group_by(Sim) %>%
  summarise(has_inf = any(is.infinite(c_across(c(
    mean_rae, mean_err, mean_ae
  )))), .groups = "drop") %>%
  filter(has_inf)

n_total_sim <- summary_rolling$Sim %>% unique() %>% length()
n_inf_sim <- nrow(inf_sim_df)

if (n_inf_sim / n_total_sim < 0.05) {
  summary_rolling <- summary_rolling %>%
    filter(!(Sim %in% inf_sim_df$Sim))
}

summary_by_time <- summary_rolling %>%
  group_by(Landmark_Time, Clustering_Method) %>%
  summarise(
    Mean_RAE    = mean(mean_rae, na.rm = TRUE),
    Mean_ERR    = mean(mean_err, na.rm = TRUE),
    MEAN_AE     = mean(mean_ae, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_by_time, 
          paste0("predict_timepoints_s2_v2", scenario, "_n_", n_total, ".csv"))

# Process and save information fraction results
summary_if <- info_results %>%
  mutate(
    #bias      = as.numeric(bias),
    #abs_err   = as.numeric(abs_err),
    #mse       = as.numeric(mse),
    mean_rae  = as.numeric(mean_rae),
    mean_err  = as.numeric(mean_err),
    mean_ae   = as.numeric(mean_ae)
  )

inf_sim_df <- summary_if %>%
  group_by(Sim) %>%
  summarise(has_inf = any(is.infinite(c_across(c(
    mean_rae, mean_err, mean_ae
  )))), .groups = "drop") %>%
  filter(has_inf)

n_total_sim <- summary_if$Sim %>% unique() %>% length()
n_inf_sim <- nrow(inf_sim_df)

if (n_inf_sim / n_total_sim < 0.05) {
  summary_if <- summary_if %>%
    filter(!(Sim %in% inf_sim_df$Sim))
}

summary_by_IF <- summary_if %>%
  group_by(IF, Clustering_Method) %>%
  summarise(
    Mean_RAE    = mean(mean_rae, na.rm = TRUE),
    Mean_ERR    = mean(mean_err, na.rm = TRUE),
    MEAN_AE     = mean(mean_ae, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_by_IF, paste0("predict_IF_s2_v2", scenario, "_n_", n_total, ".csv"))


