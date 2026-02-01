# ==============================================================================
# Real-World Data Analysis: Milestone Event Time Prediction
# ==============================================================================
# 
# This script analyzes real clinical trial data to predict milestone event times
# using covariate-driven clustering methods.
#
# Input: 
#   - real_world_full_data.rds: Full clinical trial dataset
#   - real_world_mixed_data.rds: Mixed-type covariates for clustering
#
# Output:
#   - CSV files with prediction metrics and cluster assignments
#
# ==============================================================================

# Load required libraries
library(flexsurv)     
library(dplyr)         
library(cluster)       
library(clustMixType)  
library(kamila)        
library(MASS)          
library(caret)         

# ==============================================================================
# Data Loading and Preparation
# ==============================================================================

# Load pre-processed data files
full_data <- readRDS("real_world_full_data.rds")
mixed_data <- readRDS("real_world_mixed_data.rds")

# Convert dropout indicator to numeric (0/1)
full_data$dropout <- as.numeric(full_data$dropout) - 1

# ==============================================================================
# Configuration Parameters
# ==============================================================================

# Clustering methods to evaluate
clustering_methods <- c("None", "PAM", "Hierarchical", "K-prototypes", 
                        "Kamila", "region", "imdc")

# Survival model types to fit per cluster
model_types <- c("Exponential", "Weibull", "Gompertz")

# Method for determining optimal number of clusters
n_clus_method <- "silhouette_pam"

# Target event number for milestone prediction
target_event <- 763

# Known true milestone values (for evaluation)
true_milestone <- 105.9406
true_date <- as.Date("2023-08-23")
entry_date <- as.Date("2014-10-26")

# ==============================================================================
# Optimal Number of Clusters Selection
# ==============================================================================

gower_dist <- daisy(mixed_data[, -1], metric = "gower")

# Method 1: Silhouette score with PAM
sil_width <- c()
for (k in 2:6) {
  pam_fit <- pam(gower_dist, diss = TRUE, k = k)
  sil_width[k] <- pam_fit$silinfo$avg.width
}
best_noc_pam <- which.max(sil_width)

# Method 2: Silhouette score with Hierarchical Clustering
hc <- hclust(gower_dist, method = "complete")
sil_scores <- c()
for (k in 2:6) {
  cluster_labels <- cutree(hc, k = k)
  sil <- silhouette(cluster_labels, dist = gower_dist)
  sil_scores[k] <- mean(sil[, 3])
}
best_noc_hc <- which.max(sil_scores)

# Method 3: Elbow method with K-prototypes
wss <- c()
for (k in 2:6) {
  kpres <- kproto(mixed_data[, -1], k = k, verbose = FALSE)
  wss[k] <- kpres$tot.withinss
}
plot(2:6, wss[2:6], type = "b", pch = 19,
     xlab = "Number of clusters (k)", 
     ylab = "Total within-cluster SS",
     main = "Elbow Method using K-prototypes")

# Select best number of clusters based on chosen method
if (n_clus_method == "silhouette_pam") {
  best_noc <- 3
} else if (n_clus_method == "silhouette_hc") {
  best_noc <- 2
} else if (n_clus_method == "wss_kprototype") {
  best_noc <- 4
} else {
  stop("Unsupported method for selecting number of clusters")
}

# ==============================================================================
# Main Prediction Function
# ==============================================================================

#' Predict milestone event times using cluster-specific survival models
#'
#' @param full_data Full dataset with survival times and covariates
#' @param mixed_data Mixed-type covariates for clustering (exclude ID)
#' @param clustering_method Clustering method to use
#' @param modeltypes Vector of survival model types to try
#' @param criterion Model selection criterion ("AIC" or "BIC")
#' @param targetev Target event number for milestone
#' @param sims Number of Monte Carlo simulation iterations
#' @param n_clus Number of clusters
#'
#' @return List containing:
#'   - metrics: Data frame with prediction metrics
#'   - clusters: Data frame with cluster assignments
predict_milestone <- function(full_data, mixed_data, clustering_method = "None",
                              modeltypes, criterion = c("AIC", "BIC"), 
                              targetev = target_event, sims = 1000, 
                              n_clus = best_noc) {
  
  # Apply clustering method
  if (clustering_method == "PAM") {
    gower_dist <- daisy(mixed_data[, -1], metric = "gower")
    pam_res <- pam(gower_dist, k = n_clus, diss = TRUE)
    mixed_clusters <- pam_res$clustering
    largest_cluster <- which.max(table(mixed_clusters))
    full_data$pred_cluster <- largest_cluster
    full_data$pred_cluster[full_data$id %in% mixed_data$id] <- mixed_clusters
    
  } else if (clustering_method == "Hierarchical") {
    gower_dist <- daisy(mixed_data[, -1], metric = "gower")
    hc_res <- hclust(gower_dist, method = "complete")
    mixed_clusters <- cutree(hc_res, k = n_clus)
    largest_cluster <- which.max(table(mixed_clusters))
    full_data$pred_cluster <- largest_cluster
    full_data$pred_cluster[full_data$id %in% mixed_data$id] <- mixed_clusters
    
  } else if (clustering_method == "K-prototypes") {
    kproto_res <- kproto(mixed_data[, -1], k = n_clus, verbose = FALSE)
    mixed_clusters <- kproto_res$cluster
    largest_cluster <- which.max(table(mixed_clusters))
    full_data$pred_cluster <- largest_cluster
    full_data$pred_cluster[full_data$id %in% mixed_data$id] <- mixed_clusters
    
  } else if (clustering_method == "Kamila") {
    cat_cols <- names(mixed_data)[sapply(mixed_data, function(x) is.factor(x))]
    con_cols <- names(mixed_data)[sapply(mixed_data, is.numeric)]
    kamila_res <- kamila(catFactor = mixed_data[, cat_cols], 
                         conVar = as.data.frame(mixed_data[, con_cols]), 
                         numClust = n_clus, numInit = 10)
    mixed_clusters <- kamila_res$finalMemb
    largest_cluster <- which.max(table(mixed_clusters))
    full_data$pred_cluster <- largest_cluster
    full_data$pred_cluster[full_data$id %in% mixed_data$id] <- mixed_clusters
    
  } else if (clustering_method == "None") {
    full_data$pred_cluster <- 1
    
  } else if (clustering_method == "region") {
    full_data$pred_cluster <- as.numeric(full_data$region)
    
  } else if (clustering_method == "imdc") {
    full_data$pred_cluster <- as.numeric(full_data$IMDC)
    full_data$pred_cluster[full_data$pred_cluster == "4"] <- 2
    
  } else {
    stop("Unsupported clustering method.")
  }
  
  # Fit dropout model
  data <- full_data
  dmod <- tryCatch(
    flexsurvreg(Surv(time, dropout) ~ 1, data = data, dist = "exp"),
    error = function(e) NULL
  )
  drop_mu  <- dmod$res.t[[1]]
  drop_var <- as.numeric(dmod$cov)
  
  # Fit survival models per cluster and simulate
  cluster_levels <- sort(unique(data$pred_cluster))
  sim_matrix <- matrix(NA, nrow = sims, ncol = nrow(data))
  model_per_cluster <- list()
  
  for (cl in cluster_levels) {
    cluster_data <- data[data$pred_cluster == cl, ]
    row_idx <- which(data$pred_cluster == cl)
    
    fits <- lapply(modeltypes, function(mtype) {
      tryCatch({
        flexsurvreg(Surv(time, status) ~ 1, data = cluster_data, 
                    dist = tolower(mtype))
      }, error = function(e) NULL)
    })
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0) next
    
    if (criterion == "AIC") {
      best_fit <- fits[[which.min(sapply(fits, AIC))]]
    } else {
      best_fit <- fits[[which.min(sapply(fits, BIC))]]
    }
    best_type <- best_fit$dlist[[1]]
    model_per_cluster[[as.character(cl)]] <- best_type
    t1 <- best_fit$res.t[[1]]
    cov <- best_fit$cov
    
    # Define event time simulation function
    sim_event_abs <- switch(best_type,
      "exp" = {
        function() {
          simles <- mvrnorm(1, mu = c(t1), Sigma = cov)
          rate <- exp(simles[1])
          abs_time <- -log(runif(nrow(cluster_data)) * 
                          exp(-rate * cluster_data$time)) / rate
          ifelse(cluster_data$status == 1, cluster_data$time, abs_time)
        }
      },
      "weibull.quiet" = {
        t2 <- best_fit$res.t[[2]]
        function() {
          simles <- mvrnorm(1, mu = c(t1, t2), Sigma = cov)
          shape <- exp(simles[1])
          scale <- exp(simles[2])
          S_t0 <- exp(-(cluster_data$time / scale)^shape)
          abs_time <- scale * (-log(runif(nrow(cluster_data)) * S_t0))^(1 / shape)
          ifelse(cluster_data$status == 1, cluster_data$time, abs_time)
        }
      },
      "gompertz" = {
        t2 <- best_fit$res.t[[2]]
        function() {
          simles <- mvrnorm(1, mu = c(t1, t2), Sigma = cov)
          shape <- simles[1]
          rate <- exp(simles[2])
          S_t0 <- exp(-rate / shape * (exp(shape * cluster_data$time) - 1))
          U <- runif(nrow(cluster_data)) * S_t0
          abs_time <- log(1 - (shape / rate) * log(U)) / shape
          ifelse(cluster_data$status == 1, cluster_data$time, abs_time)
        }
      },
      "llogis" = {
        t2 <- best_fit$res.t[[2]]
        function() {
          simles <- mvrnorm(1, mu = c(t1, t2), Sigma = cov)
          shape <- exp(simles[1])
          scale <- exp(simles[2])
          S_t0 <- 1 / (1 + (cluster_data$time / scale)^shape)
          U <- runif(nrow(cluster_data)) * S_t0
          abs_time <- scale * ((1 / U - 1)^(1 / shape))
          ifelse(cluster_data$status == 1, cluster_data$time, abs_time)
        }
      },
      stop("Unsupported model type")
    )
    
    # Simulate with dropout censoring
    simfun <- function() {
      lambda_drop <- exp(rnorm(1, mean = drop_mu, sd = sqrt(drop_var)))
      t_event_abs <- sim_event_abs()
      t_drop_abs <- cluster_data$time + rexp(nrow(cluster_data), rate = lambda_drop)
      t_pred_abs <- ifelse(
        cluster_data$status == 1,
        cluster_data$time,
        ifelse(t_drop_abs < t_event_abs, NA, t_event_abs)
      )
      t_pred_abs + cluster_data$entry_time
    }
    sim_cluster <- replicate(sims, simfun())
    sim_matrix[, row_idx] <- t(sim_cluster)
  }
  
  # Compute milestone predictions and metrics
  model_vec <- rep(NA, n_clus)
  model_fitted <- unlist(model_per_cluster)
  cluster_ids <- as.numeric(names(model_fitted))
  model_vec[cluster_ids] <- model_fitted
  model_cols <- setNames(as.list(model_vec), 
                        paste0("Model_Type_", seq_len(n_clus)))
  
  milestone_times <- apply(sim_matrix, 1, function(x) {
    y <- sort(x, na.last = NA)
    if (length(y) >= targetev) y[targetev] else NA_real_
  })
  pred_times <- colMeans(sim_matrix, na.rm = TRUE)
  
  rae <- abs(milestone_times - true_milestone) / true_milestone
  err <- milestone_times - true_milestone
  ae <- abs(milestone_times - true_milestone)
  
  cluster_assignments <- data.frame(
    id = full_data$id,
    cluster = full_data$pred_cluster,
    method = clustering_method
  )
  
  return(list(
    metrics = data.frame(
      Clustering_Method = clustering_method,
      model_cols,
      mean_rae = mean(rae), 
      mean_err = mean(err),
      mean_ae = mean(ae)
    ),
    clusters = cluster_assignments
  ))
}

# Run analysis for all clustering methods
result_list <- list()
cluster_label_list <- list()

for (cm in clustering_methods) {
  res <- tryCatch({
    pred_result <- predict_milestone(
      full_data = full_data,
      mixed_data = mixed_data,
      clustering_method = cm,
      modeltypes = model_types,
      criterion = "BIC",
      targetev = target_event,
      sims = 10000,
      n_clus = best_noc
    )
    result_list[[length(result_list) + 1]] <- pred_result$metrics
    cluster_label_list[[cm]] <- pred_result$clusters
  }, error = function(e) {
    message(paste("Error in", cm, ":", e$message))
    NULL
  })
}

# Save initial results
all_cols <- unique(unlist(lapply(result_list, names)))
standardized_list <- lapply(result_list, function(df) {
  missing_cols <- setdiff(all_cols, names(df))
  df[missing_cols] <- NA
  df <- df[all_cols]
  return(df)
})
result <- bind_rows(standardized_list)
write.csv(result, paste0("214_", best_noc, "clus_BIC.csv"))

all_clusters <- bind_rows(cluster_label_list)
write.csv(all_clusters, 
          paste0("214_", best_noc, "clus_cluster_labels.csv"), 
          row.names = FALSE)

# Retry logic for failed methods (optional)
max_attempts <- 500
sims_base    <- 1000
sims_step    <- 0

result_list <- setNames(vector("list", length(clustering_methods)), 
                       clustering_methods)
cluster_label_list <- setNames(vector("list", length(clustering_methods)), 
                               clustering_methods)

attempt <- 0
repeat {
  attempt <- attempt + 1
  sims_now <- sims_base + sims_step * (attempt - 1)
  set.seed(sample.int(.Machine$integer.max, 1))
  
  for (cm in clustering_methods) {
    if (!is.null(result_list[[cm]]) && 
        !is.na(result_list[[cm]]$mean_rae)) {
      next
    }
    
    pred_result <- tryCatch(
      predict_milestone(
        full_data = full_data,
        mixed_data = mixed_data,
        clustering_method = cm,
        modeltypes = model_types,
        criterion = "BIC",
        targetev = target_event,
        sims = sims_now,
        n_clus = best_noc
      ),
      error = function(e) {
        message(sprintf("Attempt %d for %s failed: %s", 
                       attempt, cm, e$message))
        NULL
      }
    )
    
    if (!is.null(pred_result)) {
      result_list[[cm]] <- transform(pred_result$metrics, 
                                     attempt = attempt, 
                                     sims = sims_now)
      cluster_label_list[[cm]] <- pred_result$clusters
    }
  }
  
  all_ok <- all(vapply(result_list, 
                     function(df) !is.null(df) && !is.na(df$mean_rae), 
                     logical(1)))
  if (all_ok || attempt >= max_attempts) break
}

failed <- names(result_list)[vapply(result_list, 
                                   function(df) is.null(df) || 
                                              is.na(df$mean_rae), 
                                   logical(1))]
if (length(failed) > 0) {
  warning(sprintf("mean_rae still NA after %d attempts for: %s",
                 attempt, paste(failed, collapse = ", ")))
}

# Save final results
metrics_df <- do.call(rbind, result_list)
row.names(metrics_df) <- NULL
write.csv(metrics_df, paste0("214_", best_noc, "clus_BIC_v2.csv"))
