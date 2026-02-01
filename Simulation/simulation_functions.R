# ==============================================================================
# Shared Simulation Functions
# ==============================================================================
# 
# This file contains common functions used across all simulation scenarios.
# Source this file at the beginning of each simulation script.
#
# Functions:
#   - apply_clustering(): Apply clustering method to data
#   - predict_milestone(): Predict milestone event times
#   - rolling_milestone_prediction(): Rolling landmark analysis
#   - info_frac_prediction(): Information fraction analysis
#
# ==============================================================================

# Apply clustering method to mixed-type data
apply_clustering <- function(full_data, mixed_data, clustering_method, n_clus) {
  if (clustering_method == "PAM") {
    gower_dist <- daisy(mixed_data, metric = "gower")
    pam_res <- pam(gower_dist, k = n_clus, diss = TRUE)
    full_data$pred_cluster <- pam_res$clustering
  } else if (clustering_method == "Hierarchical") {
    gower_dist <- daisy(mixed_data, metric = "gower")
    hc_res <- hclust(gower_dist, method = "complete")
    full_data$pred_cluster <- cutree(hc_res, k = n_clus)
  } else if (clustering_method == "K-prototypes") {
    kproto_res <- kproto(mixed_data, k = n_clus, verbose = FALSE)
    full_data$pred_cluster <- kproto_res$cluster
  } else if (clustering_method == "Kamila") {
    kamila_res <- kamila(catFactor = mixed_data[, 4:7], conVar = mixed_data[, 1:3],
                         numClust = n_clus, numInit = 10)
    full_data$pred_cluster <- kamila_res$finalMemb
  } else if (clustering_method == "None") {
    full_data$pred_cluster <- 1
  } else if (clustering_method == "Truth") {
    full_data$pred_cluster <- full_data$cluster
  } else {
    stop("Unsupported clustering method.")
  }
  return(full_data)
}

# Predict milestone event times using cluster-specific survival models
predict_milestone <- function(data, dropout_data, modeltypes, 
                              criterion = c("AIC", "BIC"), targetev = 270, 
                              sims = 1000, t_cut = NULL) {
  
  # Fit dropout model
  dmod <- tryCatch(
    flexsurvreg(Surv(time, dropout) ~ 1, data = dropout_data, dist = "exp"),
    error = function(e) NULL
  )
  drop_mu  <- dmod$res.t[[1]]
  drop_var <- as.numeric(dmod$cov)
  
  # Fit survival models per cluster and simulate
  cluster_levels <- sort(unique(data$pred_cluster))
  sim_matrix <- matrix(NA_real_, nrow = sims, ncol = nrow(data))
  model_per_cluster <- list()
  
  for (cl in cluster_levels) {
    cluster_data <- data[data$pred_cluster == cl, ]
    row_idx <- which(data$pred_cluster == cl)
    
    fits <- lapply(modeltypes, function(mtype) {
      tryCatch(
        flexsurvreg(Surv(time, status) ~ 1, data = cluster_data, 
                    dist = tolower(mtype)),
        error = function(e) NULL
      )
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
    t1  <- best_fit$res.t[[1]]
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
  
  milestone_times <- apply(sim_matrix, 1, function(x) {
    y <- sort(x, na.last = NA)
    if (length(y) >= targetev) y[targetev] else NA_real_
  })
  pred_times <- colMeans(sim_matrix, na.rm = TRUE)
  
  list(finalpred = milestone_times,
       sim_matrix = sim_matrix,
       modeltype  = model_per_cluster,
       pred_times = pred_times)
}

# Rolling landmark analysis: predictions at multiple time points
rolling_milestone_prediction <- function(full_data, mixed_data, 
                                         clustering_method = "None", 
                                         modeltypes = c("Exponential", "Weibull", 
                                                        "Gompertz", "llogis"), 
                                         criterion = c("AIC", "BIC"),
                                         targetev = 270, time_points = 1:5, 
                                         sims = 1000, n_clus = 2, sim_id = NA) {
  results <- data.frame()
  pred_df <- data.frame()
  
  # Apply clustering
  full_data <- apply_clustering(full_data, mixed_data, clustering_method, n_clus)
  
  # Rolling landmark analysis
  for (t_landmark in time_points) {
    data <- full_data
    
    # Censor observations beyond landmark time
    data$status <- ifelse(data$ttf <= t_landmark, data$status, 0)
    data$dropout <- ifelse(data$ttf <= t_landmark, data$dropout, 0)
    data$ttf <- pmin(data$ttf, t_landmark)
    data$time <- data$ttf - data$entry_time
    data_full <- data
    data <- subset(data_full, !(dropout == 1 & ttf < t_landmark))
    data <- subset(data, !(ttf < entry_time))
    
    pred_res <- predict_milestone(
      data = data,
      dropout_data = data_full,
      modeltypes = modeltypes, 
      criterion = criterion,
      targetev = targetev,
      sims = sims,
      t_cut = t_landmark
    )
    
    pred <- pred_res$finalpred
    
    model_vec <- rep(NA, n_clus)
    model_fitted <- unlist(pred_res$modeltype)
    cluster_ids <- as.numeric(names(model_fitted))
    model_vec[cluster_ids] <- model_fitted
    model_cols <- setNames(as.list(model_vec), paste0("Model_Type_", seq_len(n_clus)))
    
    true_milestone <- sort(full_data$ttf[full_data$status == 1])[targetev]
    
    if (sum(full_data$status == 1) >= targetev) {
      rae <- abs(pred - true_milestone) / true_milestone
      err <- pred - true_milestone
      ae <- abs(pred - true_milestone)
    } else {
      rae <- NA
      err <- NA
      ae <- NA
    }
    
    results <- rbind(results, data.frame(
      Sim = sim_id,
      Clustering_Method = clustering_method,
      model_cols,
      Landmark_Time = t_landmark,
      mean_rae = mean(rae, na.rm = T), 
      mean_err = mean(err, na.rm = T),
      mean_ae = mean(ae, na.rm = T)
    ))
    
    pred_df <- rbind(pred_df, data.frame(
      Sim = sim_id,
      id = subset(full_data, !(dropout == 1 & ttf < t_landmark))$id,
      Clustering_Method = clustering_method,
      model_cols,
      Landmark_Time = t_landmark,
      pred_mtimes = mean(pred),
      true_mtime = true_milestone, 
      pred_time = pred_res$pred_times,
      pred_cluster = subset(full_data, !(dropout == 1 & ttf < t_landmark))$pred_cluster
    ))
  }
  
  return(list(pred_results = results, pred_times = pred_df))
}

# Information fraction analysis: predictions at different event counts
info_frac_prediction <- function(full_data, mixed_data, clustering_method = "None", 
                                 modeltypes = c("Exponential", "Weibull", "Gompertz", "llogis"), 
                                 criterion = c("AIC", "BIC"),
                                 targetev = 270, info_frac = c(0.4, 0.6, 0.8), sims = 1000, 
                                 n_clus = 2, sim_id = NA) {
  results <- data.frame()
  pred_df <- data.frame()
  
  # Apply clustering
  full_data <- apply_clustering(full_data, mixed_data, clustering_method, n_clus)
  
  for (percent in info_frac) {
    data <- full_data
    event_count <- floor(targetev * percent)
    
    # Censor those beyond the event count
    target_time <- sort(data$time[data$status == 1])[event_count]
    data$status <- ifelse(data$time <= target_time, data$status, 0)
    data$time <- pmin(data$time, target_time)
    data$dropout <- ifelse(data$time <= target_time, data$dropout, 0)
    data_full <- data
    
    data <- subset(data, !(dropout == 1 & time < target_time))
    
    pred_res <- predict_milestone(
      data = data,
      dropout_data = data_full,
      modeltypes = modeltypes, 
      criterion = criterion,
      targetev = targetev,
      sims = sims,
      t_cut = target_time
    )
    
    pred <- pred_res$finalpred
    
    model_vec <- rep(NA, n_clus)
    model_fitted <- unlist(pred_res$modeltype)
    cluster_ids <- as.numeric(names(model_fitted))
    model_vec[cluster_ids] <- model_fitted
    model_cols <- setNames(as.list(model_vec), paste0("Model_Type_", seq_len(n_clus)))
    
    true_milestone <- sort(full_data$ttf[full_data$status == 1])[targetev]
    
    if (sum(full_data$status == 1) >= targetev) {
      rae <- abs(pred - true_milestone) / true_milestone
      err <- pred - true_milestone
      ae <- abs(pred - true_milestone)
    } else {
      rae <- NA
      err <- NA
      ae <- NA
    }
    
    results <- rbind(results, data.frame(
      Sim = sim_id,
      Clustering_Method = clustering_method,
      model_cols,
      IF = percent,
      mean_rae = mean(rae, na.rm = T), 
      mean_err = mean(err, na.rm = T),
      mean_ae = mean(ae, na.rm = T)
    ))
    pred_df <- rbind(pred_df, data.frame(
      Sim = sim_id,
      id = subset(full_data, !(dropout == 1 & time < target_time))$id,
      Clustering_Method = clustering_method,
      model_cols,
      IF = percent,
      pred_mtimes = mean(pred),
      true_mtime = true_milestone, 
      pred_time = pred_res$pred_times
    ))
  }
  
  return(list(pred_results = results, pred_times = pred_df))
}
