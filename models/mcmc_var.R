
#' Gibbs sampler for auto-regression
#'
#' @param measured Tibble of measured values
#' @param vars_meas A character vector of measured variable names
#' @param n_mcmc MCMC iterations
#' @param n_mcmc VAR length
#'
#' @return A list of MCMC samples for each parameter
#' @export
#'
#' @examples
#' 
sample_mcmc_var <- function(
  measured,
  vars_meas,
  n_mcmc,
  p = 24
){

  K <- length(vars_meas)

  phi_sample <- array(NA, c(K * p, K, n_mcmc))
  sigma_sample <- array(NA, c(K, K, n_mcmc))

  V0 <- diag(0.05, K)
  nu0 <- 2

  phi0 <- matrix(c(c(0.9, 0, 0, 0.9), rep(0, K * K * (p-1))), ncol = 2, byrow = TRUE)
  A0_phi <- diag(1, K*p)

  R0 <- matrix(0, nrow = p, ncol = K)

  sigma_sample[ , , 1] <- V0
  phi_sample[ , , 1] <- phi0

  Y <- measured %>%
    dplyr::select(all_of({{vars_meas}})) %>%
    as.matrix()

  Y[which(is.na(Y))] <- mean(Y, na.rm = TRUE)

  for (i in 2:n_mcmc){

    cat(sprintf("\rFitting MCMC AR: Iteration %i/%i", i, n_mcmc))

    E_lag <- embed(rbind(R0, Y), p + 1)[ ,-(1:K)]

    Phin <- solve(t(E_lag) %*% E_lag + A0_phi) %*% (t(E_lag) %*% Y + A0_phi %*% phi0)
    An_phi <- t(E_lag) %*% E_lag + A0_phi

    phi_sample[ , , i] <- rmatrixnorm(
      Phin, 
      round(solve(An_phi), 9), 
      sigma_sample[ , , i-1]
    )

    A_bold <- rbind(
      t(phi_sample[ , , i]), 
      diag(1, K*(p-1)) %>% cbind(matrix(0, ncol = K, nrow = K*(p-1)))
    )

    if(eigen(A_bold)$values %>% abs() %>% max() >= 1) {
      phi_sample[ , , i] <- phi_sample[ , , i-1]
      warning(sprintf("Stability Conditions Violated: Iteration %i/%i", i, n))
    }

    var_comp <- E_lag %*% phi_sample[ , , i]

    Vn <- V0 + t(Y - var_comp) %*% (Y - var_comp)
    nun <- nu0 + nrow(Y)

    sigma_sample[ , , i] <- rinvwishart(nun, Vn)

  }

  mcmc_fit <- list(
    phi_sample = phi_sample,
    sigma_sample = sigma_sample,
    K = K,
    p = p,
    n_mcmc = n_mcmc
  )
  
}

#' Generate predictions for VAR forecasts from MCMC samples
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
predict_mcmc_var <- function(
  mcmc_list,
  measured,
  vars_meas, 
  forecast,
  max_horizon
){

  horizons <- 1:max_horizon
  H <- max_horizon
  K <- mcmc_list$K
  p <- mcmc_list$p
  n_mcmc <- mcmc_list$n_mcmc

  R0 <- matrix(0, nrow = p, ncol = K)

  times_issued <- forecast %>%
    filter(time_issued >= min(measured$time)+hours(p)) %>%
    pull(time_issued) %>%
    unique()

  times_issued <- times_issued[-which(!times_issued %in% measured$time)]

  y_predict <- array(NA, dim = c(length(times_issued), K, H, n_mcmc))

  for (i in 1:n_mcmc){

    for (t in times_issued){

      idx <- which(measured$time == t)

      cat(
        sprintf("\rIteration %i/%i Time %i/%i     ", i, n_mcmc,
        which(t == times_issued), length(times_issued))
      )

      pred <- array(NA, c(K, H))

      E <- as.matrix(measured[(idx-p+1):idx, c(2,3)])

      pred[ ,1] <- t(mcmc_list$phi_sample[ , , i]) %*% t(embed(E, p)) + 
        t(rmvnorm(1, sigma = mcmc_list$sigma_sample[ , , i]))

      for(h in 2:H){

        E <- rbind(E[-1, ], pred[, h-1])

        pred[ ,h] <- t(mcmc_list$phi_sample[ , , i]) %*% t(embed(E, p)) + 
        t(rmvnorm(1, sigma = mcmc_list$sigma_sample[ , , i]))  

      }

      y_predict[which(t == times_issued), , , i] <- pred

    }

  }

  mcmc_predictions <- list(
    y_predict = y_predict,
    times_issued = times_issued,
    phi_samples = mcmc_list$phi_sample,
    sigma_samples = mcmc_list$sigma_sample,
    K = K,
    p = p,
    horizons = horizons,
    n_mcmc = n_mcmc
  )

  return(mcmc_predictions)
  
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
calc_diagnostics_var <- function(
  mcmc_predictions,
  measured,
  vars_meas, 
  burnin
){

  y_predict <- mcmc_predictions$y_predict

  n_mcmc <- mcmc_predictions$n_mcmc

  K <- mcmc_predictions$K

  horizons <- mcmc_predictions$horizons
  H <- length(horizons)

  times_issued <- mcmc_predictions$times_issued

  rmse_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating RMS: Horizon %i     ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

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
      model = "VAR"
    )

  }) %>% bind_rows()

  dss_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating DSS: Horizon %i    ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

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
      model = "VAR"
    )

  }) %>% bind_rows()

  es_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating ES: Horizon %i   ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

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
      model = "VAR"
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
calc_long_diagnostics_var <- function(
  mcmc_predictions,
  measured,
  vars_meas, 
  burnin
){

  y_predict <- mcmc_predictions$y_predict

  # beta_sample <- mcmc_predictions$beta_samples
  # sigma_sample <- mcmc_predictions$sigma_samples

  n_mcmc <- mcmc_predictions$n_mcmc

  K <- mcmc_predictions$K

  horizons <- mcmc_predictions$horizons
  H <- length(horizons)

  times_issued <- mcmc_predictions$times_issued

  rmse_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating RMS: Horizon %i     ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

    Y_mat <- data_tmp %>%
      dplyr::select(all_of({{ vars_meas }})) %>%
      as.matrix()

    mean <- apply(y_predict[ , , h_idx, burnin:n_mcmc], c(1,2), mean, na.rm = TRUE)

    resid <- mean - Y_mat
    na_idx <- c(which(is.nan(resid)), which(is.na(resid)))
    resid[na_idx] <- mean(resid, na.rm = TRUE)

    # rms <- mean(sqrt(resid[,1]^2 + resid[,2]^2))

    tibble(
      times_issued = times_issued,
      horizon = h,
      SE = resid[,1]^2 + resid[,2]^2,
      metric = "SE",
      model = "VAR"
    )

  }) %>% bind_rows()

  dss_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating DSS: Horizon %i    ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

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
      model = "VAR"
    )

  }) %>% bind_rows()

  es_tbl <- lapply(horizons, function(h){

    cat(sprintf("\rCalculating ES: Horizon %i   ", h))

    h_idx <- which(h == horizons)

    data_tmp <- tibble(
      time_issued = times_issued,
      time_predict = times_issued + hours(h)
    ) %>%
      left_join(measured, by = c("time_predict" = "time"))

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
      model = "VAR"
    )

  }) %>% bind_rows()
 
  return(
    list(rmse_tbl, dss_tbl, es_tbl)
  )

}