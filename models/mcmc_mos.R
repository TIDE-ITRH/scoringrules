
#' Gibbs sampler for Model Output Statistics regression
#'
#' @param measured Tibble of measured values
#' @param vars_meas A character vector of measured variable names
#' @param forecast Tibble of forecast values
#' @param vars_fore A character vector of forecast variable names
#' @param n_mcmc MCMC iterations
#'
#' @return A list of MCMC samples for each parameter
#' @export
#'
#' @examples
#' 
sample_mcmc_mos <- function(
  measured,
  vars_meas, 
  forecast,
  vars_fore,
  n_mcmc
){

  horizons <- forecast$horizon %>% unique()
  H <- length(horizons)

  K <- length(vars_meas)

  beta_sample <- array(NA, c(K + 1, K, H, n_mcmc))
  sigma_sample <- array(NA, c(K, K, H, n_mcmc))

  V0 <- diag(0.05, K)
  nu0 <- 2

  beta0 <- matrix(c(0, 1, 0, 0, 0, 1), ncol = K)
  A0 <- diag(1, K + 1)

  for (h in 1:H){
    beta_sample[ , , h, 1] <- beta0
    sigma_sample[ , , h, 1] <- V0
  }

  data_combined <- measured %>%
      left_join(forecast, by = c("time" = "time_predict"))

  for (i in 2:n_mcmc){

    cat(sprintf("\rFitting MCMC MOS: Iteration %i/%i", i, n_mcmc))

    for (h in horizons){

      h_idx <- which(h == horizons)

      Y_mat <- data_combined %>%
        filter(horizon == h) %>%
        dplyr::select(all_of({{ vars_meas }})) %>%
        as.matrix()

      X_tbl <- data_combined %>%
        filter(horizon == h) %>%
        dplyr::select(all_of({{ vars_fore }}))

      X_mat <- cbind(
        const = rep(1, nrow(X_tbl)), 
        X_tbl %>% dplyr::select(all_of({{ vars_fore }}))
      ) %>% as.matrix()

      Y_mat[which(is.na(Y_mat))] <- mean(Y_mat, na.rm = TRUE)
      X_mat[which(is.na(X_mat))] <- mean(X_mat, na.rm = TRUE)

      Bn <- solve(t(X_mat) %*% X_mat + A0) %*% 
        (t(X_mat) %*% (Y_mat) + A0 %*% beta0)
      An <- t(X_mat) %*% X_mat + A0

      Vn <- V0 + t(Y_mat - X_mat %*% Bn) %*% 
        (Y_mat - X_mat %*% Bn) + t(Bn - beta0) %*% A0 %*% (Bn - beta0)
      nun <- nu0 + nrow(X_mat)

      sigma_sample[ , , h_idx, i] <- rinvwishart(nun, Vn)

      beta_sample[ , , h_idx, i] <- rmatrixnorm(
        Bn, 
        round(solve(An), 9), 
        sigma_sample[ , , h_idx, i]
      )  

    }

  }

  mcmc_fit <- list(
    beta_samples = beta_sample,
    sigma_samples = sigma_sample,
    K = K,
    horizons = horizons,
    n_mcmc = n_mcmc
  )
  
}

#' Generate predictions for MOS forecasts from MCMC samples
#'
#' @param mcmc_list MCMC results from sample_mcmc_mos
#' @param measured Tibble of measured values
#' @param vars_meas A character vector of measured variable names
#' @param forecast Tibble of forecast values
#' @param vars_fore A character vector of forecast variable names
#'
#' @return A list with an array of predictions and their times
#' @export
#'
#' @examples
#' 
predict_mcmc_mos <- function(
  mcmc_list,
  measured,
  vars_meas, 
  forecast,
  vars_fore
){

  horizons <- mcmc_list$horizons
  K <- mcmc_list$K
  n_mcmc <- mcmc_list$n_mcmc

  H <- length(horizons)

  times_issued <- forecast$time_issued %>% unique()

  y_predict <- array(NA, dim = c(length(times_issued), K, H, n_mcmc))

  for (i in 1:n_mcmc){

    for (h in horizons){

      cat(sprintf("\rPredicting MCMC MOS: Horizon %i Iteration %i/%i     ", h, i, n_mcmc))

      h_idx <- which(h == horizons)

      data_tmp <- tibble(
        time_issued = times_issued,
        time_predict = times_issued + hours(h)
      ) %>%
        left_join(measured, by = c("time_predict" = "time")) %>%
        left_join(forecast, by = c("time_predict", "time_issued"))

      X_mat <- cbind(
        const = rep(1, nrow(data_tmp)), 
        data_tmp %>% dplyr::select(all_of({{ vars_fore }}))
      ) %>% as.matrix()

      beta <- mcmc_list$beta_sample[ , , h_idx, i]
      sigma <- mcmc_list$sigma_sample[ , , h_idx, i]

      Y_pred <- X_mat %*% beta + rmvnorm(nrow(data_tmp), sigma = sigma)

      y_predict[ , , h_idx, i] <- Y_pred

    }

  }

  mcmc_predictions <- list(
    y_predict = y_predict,
    beta_samples = mcmc_list$beta_sample,
    sigma_samples = mcmc_list$sigma_sample,
    K = K,
    horizons = horizons,
    n_mcmc = n_mcmc
  )
  
}

#' Calculate diagnostics of MOS model from MCMC samples
#'
#' @param mcmc_predictions MCMC results from predict_mcmc_mos
#' @param measured Tibble of measured values
#' @param vars_meas A character vector of measured variable names
#' @param forecast Tibble of forecast values
#' @param vars_fore A character vector of forecast variable names
#' @param burnin A numeric value of burn-in index
#'
#' @return A tibble with diagnostic results
#' @export
#'
#' @examples
#' 
calc_diagnostics_mos <- function(
  mcmc_predictions,
  measured,
  vars_meas, 
  forecast,
  vars_fore,
  burnin
){

  y_predict <- mcmc_predictions$y_predict

  beta_sample <- mcmc_predictions$beta_samples
  sigma_sample <- mcmc_predictions$sigma_samples

  n_mcmc <- mcmc_predictions$n_mcmc

  K <- mcmc_predictions$K

  horizons <- mcmc_predictions$horizons
  H <- length(horizons)

  times_issued <- forecast$time_issued %>% unique()

  rmse_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating RMS: Horizon %i     ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    mean <- apply(y_predict[ , , h_idx, burnin:n_mcmc], c(1,2), mean, na.rm = TRUE)

    resid <- mean - Y_mat
    na_idx <- c(which(is.nan(resid)), which(is.na(resid)))
    resid[na_idx] <- mean(resid, na.rm = TRUE)

    rms <- mean(sqrt(resid[,1]^2 + resid[,2]^2))

    tibble(
      horizon = h,
      score = rms,
      metric = "RMSE",
      model = "MOS"
    )

  }) %>% bind_rows()

  dss_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating DSS: Horizon %i    ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    mean <- apply(y_predict[ , , h_idx, burnin:n_mcmc], c(1,2), mean, na.rm = TRUE)

    resid <- mean - Y_mat
    na_idx <- c(which(is.nan(resid)), which(is.na(resid)))
    resid[na_idx] <- mean(resid, na.rm = TRUE)

    dss <- array(NA, dim = nrow(resid))

    for(i in 1:nrow(resid)){

      if(i %in% na_idx){
        dss[i] <- NA
      } else {
        var_calc <- cov(t(y_predict[i, , h_idx, burnin:n_mcmc]))
        dss[i] <- -dmvnorm(resid[i, ], sigma = var_calc, log = TRUE)
      }

    }
    
    tibble(
      horizon = h,
      score = mean(dss, na.rm = TRUE),
      metric = "DSS",
      model = "MOS"
    )

  }) %>% bind_rows()

  es_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating ES: Horizon %i   ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    es <- array(NA, dim = nrow(Y_mat))

    for(i in 1:nrow(Y_mat)){

      es[i] <- es_sample(Y_mat[i,], y_predict[i, , h_idx, burnin:n_mcmc])

    }
    
    tibble(
      horizon = h,
      score = mean(es, na.rm = TRUE),
      metric = "ES",
      model = "MOS"
    )

  }) %>% bind_rows()
 
  return(
    bind_rows(rmse_tbl, dss_tbl, es_tbl)
  )

}


#' Calculate diagnostics of MOS model from MCMC samples
#'
#' @param mcmc_predictions MCMC results from predict_mcmc_mos
#' @param measured Tibble of measured values
#' @param vars_meas A character vector of measured variable names
#' @param forecast Tibble of forecast values
#' @param vars_fore A character vector of forecast variable names
#' @param burnin A numeric value of burn-in index
#'
#' @return A tibble with diagnostic results
#' @export
#'
#' @examples
#' 
calc_long_diagnostics_mos <- function(
  mcmc_predictions,
  measured,
  vars_meas, 
  forecast,
  vars_fore,
  burnin
){

  y_predict <- mcmc_predictions$y_predict

  beta_sample <- mcmc_predictions$beta_samples
  sigma_sample <- mcmc_predictions$sigma_samples

  n_mcmc <- mcmc_predictions$n_mcmc

  K <- mcmc_predictions$K

  horizons <- mcmc_predictions$horizons[which(mcmc_predictions$horizons <= 48)]
  H <- length(horizons)

  times_issued <- forecast$time_issued %>% unique()

  rmse_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating RMS: Horizon %i     ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    mean <- apply(y_predict[ , , h_idx, burnin:n_mcmc], c(1,2), mean, na.rm = TRUE)

    resid <- mean - Y_mat
    na_idx <- c(which(is.nan(resid)), which(is.na(resid)))
    resid[na_idx] <- mean(resid, na.rm = TRUE)

    tibble(
      times_issued = times_issued,
      horizon = h,
      SE = resid[,1]^2 + resid[,2]^2,
      metric = "SE",
      model = "MOS"
    )

  }) %>% bind_rows()

  dss_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating DSS: Horizon %i    ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    mean <- apply(y_predict[ , , h_idx, burnin:n_mcmc], c(1,2), mean, na.rm = TRUE)

    resid <- mean - Y_mat
    na_idx <- c(which(is.nan(resid)), which(is.na(resid)))
    resid[na_idx] <- mean(resid, na.rm = TRUE)

    dss <- array(NA, dim = nrow(resid))

    for(i in 1:nrow(resid)){

      if(i %in% na_idx){
        dss[i] <- NA
      } else {
        var_calc <- cov(t(y_predict[i, , h_idx, burnin:n_mcmc]))
        dss[i] <- -dmvnorm(resid[i, ], sigma = var_calc, log = TRUE)
      }

    }

    tibble(
      times_issued = times_issued,
      horizon = h,
      DSS = dss,
      metric = "DSS",
      model = "MOS"
    )

  }) %>% bind_rows()

  es_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating ES: Horizon %i   ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time")) %>%
      left_join(forecast, by = c("time_predict", "time_issued"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    es <- array(NA, dim = nrow(Y_mat))

    for(i in 1:nrow(Y_mat)){

      es[i] <- es_sample(Y_mat[i,], y_predict[i, , h_idx, burnin:n_mcmc])

    }
    
    tibble(
      times_issued = times_issued,
      horizon = h,
      ES = es,
      metric = "ES",
      model = "MOS"
    )

  }) %>% bind_rows()
 
  return(
    list(rmse_tbl, dss_tbl, es_tbl)
  )

}
