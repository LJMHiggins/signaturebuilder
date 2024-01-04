#' Sample skewed distribution
#'
#' @param cp_params String vector of centred parameters for sn package.
#' @param sample_no Integer. Number of samples to draw from simulated distribution.
#'
#' @return List with 2 data.frames, simulated_dists contains sampled values.
#'  df_quants contains quantiles and median.
#' @export
#'
#' @examples
simulate_distribution <- function(cp_params, sample_no = 500){
  if (length(cp_params) == 3){
    params <- sn::cp2dp(cp_params, family = "SN")
    boot.dists <- sapply(1:1000, function(i) {
      sim <- sn::rsn(sample_no, dp = params)
    }
    )
  } else if (length(cp_params) == 4){
    params <- sn::cp2dp(cp_params, family = "ST")
    boot.dists <- sapply(1:1000, function(i) {
      sim <- sn::rst(sample_no, dp = params)
    }
    )
  }
  simulated_dists <- data.frame(sim.x = as.vector(boot.dists),
                                iteration = rep(1:ncol(boot.dens), each = sample_no))

  quants <- apply(boot.dists, 1, quantile, c(0.025, 0.5, 0.975))
  df_quants <- data.frame(
    lower = quants[1, ],
    upper = quants[3, ],
    middle = quants[2, ]
  )
  df_quants <- df_quants[rep(seq_len(nrow(df_quants)), each = sample_no), ]
  return(list(boot.dists = simulated_dists,
              quants = df_quants))
}

#' Simulate skewed density
#'
#' @param cp_params String vector of centred parameters for sn package.
#' @param sample_no Integer. Number of samples to draw from simulated distribution.
#'
#' @return data.frame with sampled values, respective probability density and iteration id.
#' @export
#'
#' @examples
simulate_density <- function(cp_params, sample_no = 500){

  .skew_norm <- function(params, sample_no){

    boot.dists <- sapply(1:1000, function(i) {
      sim <- sn::rsn(sample_no, dp = params)
    })
    boot.dens <- apply(boot.dists, 2, function(col){
      sim.dens <- sn::dsn(col, dp = params)
    })
    simulated_pdfs <- data.frame(x = as.vector(boot.dists),
                                 density = as.vector(boot.dens),
                                 iteration = rep(1:ncol(boot.dens), each = sample_no))
    return(simulated_pdfs)
  }
  .skew_t <- function(params, sample_no){

    boot.dists <- sapply(1:1000, function(i) {
      sim <- sn::rst(sample_no, dp = params)
    })
    boot.dens <- apply(boot.dists, 2, function(col){
      sim.dens <- sn::dst(col, dp = params)
    })
    simulated_pdfs <- data.frame(x = as.vector(boot.dists),
                                 density = as.vector(boot.dens),
                                 iteration = rep(1:ncol(boot.dens), each = sample_no))
    return(simulated_pdfs)
  }
  if (length(cp_params) == 3){
    params <- sn::cp2dp(cp_params, family = "SN")
    simulated_pdfs <- .skew_norm(params = params, sample_no = sample_no)
  } else if (length(cp_params) == 4){
    params <- sn::cp2dp(cp_params, family = "ST")
    simulated_pdfs <- .skew_t(params = params, sample_no = sample_no)
  }
  return(simulated_pdfs)
}

#' Fit skewed distributions and estimate parameters
#'
#' @param data_df data.frame of data to be modelled. Each column represents a
#'  separate distribution.
#' @param distributions String vector with distributions to fit.
#'
#' @return data.frame with estimated distribution parameters and AIC metric.
#' @export
#'
#' @examples
fit_skewed_dists <- function(data_df,
                             distributions = c("normal", "skew-normal", "skew-t")){

  grid_exec <- expand.grid(distributions,
                           colnames(data_df),
                           stringsAsFactors = FALSE)

  dist_fit_results <- foreach::foreach(i = 1:nrow(grid_exec), .combine = "rbind")%do%{
    dist_type <- grid_exec[i, 1]
    rnai_i <- grid_exec[i, 2]
    rnai_data <- d2_sample[,rnai_i] %>% na.omit()

    if (dist_type == "normal"){
      fit_output <- tryCatch({

        fit <- MASS::fitdistr(rnai_data, densfun = dist_type)
        AIC_test <- AIC(fit)
        mean_est <- fit$estimate[[1]]
        sd_est <- fit$estimate[[2]]
        return(data.frame(rnai = rnai_i,
                          distribution = dist_type,
                          est_mean = mean_est,
                          est_sd = sd_est,
                          est_skew = NA,
                          est_kurtosis = NA,
                          AIC_score = AIC_test))

      }, error = function(e){
        print(paste("Could not fit ", dist_type, " distribution to ", rnai_i))
        return(NULL)
      }
      )
    } else if (dist_type == "skew-normal"){
      fit_output <- tryCatch({
        fit <- sn::selm(rnai_data ~ 1, family = "SN")
        AIC_test <- AIC(fit)
        mean_est <- fit@param$cp[[1]]
        sd_est <- fit@param$cp[[2]]
        gamma1_est <- fit@param$cp[[3]]
        return(data.frame(rnai = rnai_i,
                          distribution = dist_type,
                          est_mean = mean_est,
                          est_sd = sd_est,
                          est_skew = gamma1_est,
                          est_kurtosis = NA,
                          AIC_score = AIC_test))

      }, error = function(e){
        print(paste("Could not fit ", dist_type, " distribution to ", rnai_i))
        return(NULL)
      }
      )
    } else if (dist_type == "skew-t"){
      fit_output <- tryCatch({
        fit <- sn::selm(rnai_data ~ 1, family = "ST")
        AIC_test <- AIC(fit)
        mean_est <- fit@param$cp[[1]]
        sd_est <- fit@param$cp[[2]]
        gamma1_est <- fit@param$cp[[3]]
        gamma2_est <- fit@param$cp[[4]]
        return(data.frame(rnai = rnai_i,
                          distribution = dist_type,
                          est_mean = mean_est,
                          est_sd = sd_est,
                          est_skew = gamma1_est,
                          est_kurtosis = gamma2_est,
                          AIC_score = AIC_test))

      }, error = function(e){
        print(paste("Could not fit ", dist_type, " distribution to ", rnai_i))
        return(NULL)
      }
      )
    }
    return(fit_output)
  }
  # Final output df with parameters
  return(dist_fit_results)
}
