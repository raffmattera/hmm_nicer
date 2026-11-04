# Original HMM mathematical functions, mechanically consolidated.
# Source: the previous NICER_AOAS_repro project; function bodies are unchanged.
# Per-file SHA-256 values are recorded in README.md.

# ===== VERBATIM SOURCE: R/hmm_original/utils.R =====
# Internal numerical utilities -------------------------------------------------

logSumExp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

safe_log <- function(x, eps = 1e-12) log(pmax(x, eps))

`%^%` <- function(A, m) {
  if (m == 0) return(diag(nrow(A)))
  if (m == 1) return(A)
  res <- diag(nrow(A))
  base <- A
  mm <- as.integer(m)
  while (mm > 0) {
    if (mm %% 2 == 1) res <- res %*% base
    mm <- mm %/% 2
    if (mm > 0) base <- base %*% base
  }
  res
}

make_A_power_cache <- function(A) {
  cache <- new.env(parent = emptyenv())
  cache[["0"]] <- diag(nrow(A))
  cache[["1"]] <- A
  function(m) {
    key <- as.character(m)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, A %^% m, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
}

is_missing_skellam <- function(y, expo) (expo <= 0) | is.na(y)
is_missing_pois2 <- function(H, S, expo) (expo <= 0) | is.na(H) | is.na(S)

find_missing_blocks <- function(miss) {
  T <- length(miss)
  blocks <- list()
  i <- 1L
  b <- 0L
  while (i <= T) {
    if (!miss[i]) {
      i <- i + 1L
    } else {
      s <- i
      while (i <= T && miss[i]) i <- i + 1L
      e <- i - 1L
      b <- b + 1L
      blocks[[b]] <- c(s, e)
    }
  }
  blocks
}

#' Bayesian Information Criterion for HMMs
#'
#' @param loglik Observed log-likelihood.
#' @param K Number of latent states.
#' @param n_eff Effective sample size, usually the number of bins with positive exposure.
#' @param emission Emission family.
#' @export
bic_hmm <- function(loglik, K, n_eff,
                    emission = c("skellam", "pois2", "bivpois")) {
  
  emission <- match.arg(emission)
  
  n_transition <- (K - 1) + K * (K - 1)
  
  n_emission <- switch(
    emission,
    skellam = 2 * K,
    pois2 = 2 * K,
    bivpois = 3 * K
  )
  
  p <- n_transition + n_emission
  
  -2 * loglik + p * log(n_eff)
}

#' Integrated Completed Likelihood criterion for HMMs
#'
#' @param loglik Observed log-likelihood.
#' @param gamma Posterior state-probability matrix.
#' @param K Number of latent states.
#' @param n_eff Effective sample size.
#' @export
icl_hmm <- function(loglik, gamma, K, n_eff,
                    emission = c("skellam", "pois2", "bivpois")) {
  
  emission <- match.arg(emission)
  
  n_transition <- (K - 1) + K * (K - 1)
  
  n_emission <- switch(
    emission,
    skellam = 2 * K,
    pois2 = 2 * K,
    bivpois = 3 * K
  )
  
  p <- n_transition + n_emission
  
  bic <- -2 * loglik + p * log(n_eff)
  
  g <- pmax(gamma, 1e-12)
  ent <- -sum(g * log(g), na.rm = TRUE)
  
  bic + 2 * ent
}
#' Sort Skellam states by increasing differential mean
#' @export
sort_skellam_states <- function(fit) {
  ord <- order(fit$mu1 - fit$mu2)
  fit$mu1 <- fit$mu1[ord]
  fit$mu2 <- fit$mu2[ord]
  fit$pi1 <- fit$pi1[ord]
  fit$A <- fit$A[ord, ord, drop = FALSE]
  if (!is.null(fit$gamma)) fit$gamma <- fit$gamma[, ord, drop = FALSE]
  fit
}

#' Sort bivariate Poisson states by increasing first-channel rate
#' @export
sort_pois2_states <- function(fit) {
  ord <- order(fit$lamH)
  fit$lamH <- fit$lamH[ord]
  fit$lamS <- fit$lamS[ord]
  fit$pi1 <- fit$pi1[ord]
  fit$A <- fit$A[ord, ord, drop = FALSE]
  if (!is.null(fit$gamma)) fit$gamma <- fit$gamma[, ord, drop = FALSE]
  fit
}

sort_pois2_states_by_contrast <- function(model) {
  ord <- order(model$lamH - model$lamS)
  
  model$lamH <- model$lamH[ord]
  model$lamS <- model$lamS[ord]
  model$pi1  <- model$pi1[ord]
  model$A    <- model$A[ord, ord]
  
  if (!is.null(model$gamma)) {
    model$gamma <- model$gamma[, ord, drop = FALSE]
  }
  
  model
}

sort_skellam_states_by_contrast <- function(model) {
  ord <- order(model$mu1 - model$mu2)
  
  model$mu1 <- model$mu1[ord]
  model$mu2 <- model$mu2[ord]
  model$pi1 <- model$pi1[ord]
  model$A   <- model$A[ord, ord]
  
  if (!is.null(model$gamma)) {
    model$gamma <- model$gamma[, ord, drop = FALSE]
  }
  
  model
}

sort_bivpois_states_by_contrast <- function(model) {
  ord <- order(model$lamH - model$lamS)
  
  model$lamH <- model$lamH[ord]
  model$lamS <- model$lamS[ord]
  model$lamC <- model$lamC[ord]
  model$pi1  <- model$pi1[ord]
  model$A    <- model$A[ord, ord]
  
  if (!is.null(model$gamma)) {
    model$gamma <- model$gamma[, ord, drop = FALSE]
  }
  
  model
}
# Utils simulations


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

safe_select_K <- function(criteria, K_grid) {
  if (all(is.na(criteria))) return(NA_integer_)
  K_grid[which.min(criteria)]
}

safe_ARI <- function(z_true, z_hat) {
  ok <- !is.na(z_true) & !is.na(z_hat)
  if (sum(ok) < 2) return(NA_real_)
  adjustedRandIndex(z_true[ok], z_hat[ok])
}

safe_sojourn_error <- function(A_hat, A_true) {
  if (is.null(A_hat) || is.null(A_true)) return(NA_real_)
  
  soj_hat <- 1 / pmax(1 - diag(A_hat), 1e-8)
  soj_true <- 1 / pmax(1 - diag(A_true), 1e-8)
  
  sqrt(mean((soj_hat - soj_true)^2))
}

get_iter <- function(fit) {
  if (is.null(fit)) return(NA_integer_)
  if (!is.null(fit$iter)) return(fit$iter)
  NA_integer_
}

get_converged <- function(fit) {
  if (is.null(fit)) return(NA)
  if (!is.null(fit$converged)) return(fit$converged)
  NA
}

# ------------------------------------------------------------
# Retrieve fitted model corresponding to selected K
# ------------------------------------------------------------

get_fit_by_selected_K <- function(fits, K_selected, K_fit_grid) {
  
  if (is.null(K_selected) || is.na(K_selected)) {
    return(NULL)
  }
  
  idx <- match(K_selected, K_fit_grid)
  
  if (is.na(idx)) {
    return(NULL)
  }
  
  fits[[idx]]
}

contrast_trajectory_rmse <- function(z_true, z_hat, true_contrast, est_contrast) {
  
  ok <- !is.na(z_true) & !is.na(z_hat)
  
  if (sum(ok) < 2) return(NA_real_)
  
  z_true <- as.integer(z_true[ok])
  z_hat <- as.integer(z_hat[ok])
  
  if (any(z_true < 1 | z_true > length(true_contrast))) return(NA_real_)
  if (any(z_hat < 1 | z_hat > length(est_contrast))) return(NA_real_)
  
  d_true <- true_contrast[z_true]
  d_hat <- est_contrast[z_hat]
  
  sqrt(mean((d_true - d_hat)^2))
}

nearest_contrast_state <- function(x, reference_contrast) {
  vapply(
    x,
    function(xx) which.min(abs(xx - reference_contrast)),
    integer(1)
  )
}

contrast_state_accuracy <- function(z_true, z_hat, true_contrast, est_contrast) {
  
  ok <- !is.na(z_true) & !is.na(z_hat)
  
  if (sum(ok) < 2) return(NA_real_)
  
  z_true <- as.integer(z_true[ok])
  z_hat <- as.integer(z_hat[ok])
  
  if (any(z_true < 1 | z_true > length(true_contrast))) return(NA_real_)
  if (any(z_hat < 1 | z_hat > length(est_contrast))) return(NA_real_)
  
  true_state <- z_true
  
  est_state <- nearest_contrast_state(
    x = est_contrast[z_hat],
    reference_contrast = true_contrast
  )
  
  mean(true_state == est_state)
}

contrast_sign_accuracy <- function(z_true, z_hat, true_contrast, est_contrast) {
  
  ok <- !is.na(z_true) & !is.na(z_hat)
  
  if (sum(ok) < 2) return(NA_real_)
  
  z_true <- as.integer(z_true[ok])
  z_hat <- as.integer(z_hat[ok])
  
  if (any(z_true < 1 | z_true > length(true_contrast))) return(NA_real_)
  if (any(z_hat < 1 | z_hat > length(est_contrast))) return(NA_real_)
  
  s_true <- sign(true_contrast[z_true])
  s_hat <- sign(est_contrast[z_hat])
  
  mean(s_true == s_hat)
}

# ------------------------------------------------------------
# Block-gap utilities
# ------------------------------------------------------------

generate_block_lengths <- function(T_obs,
                                   D,
                                   min_block_length = 20,
                                   jitter = 0.20,
                                   seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(D >= 1)
  stopifnot(T_obs >= D * min_block_length)
  
  base <- rep(floor(T_obs / D), D)
  remainder <- T_obs - sum(base)
  
  if (remainder > 0) {
    base[seq_len(remainder)] <- base[seq_len(remainder)] + 1
  }
  
  block_lengths <- base
  
  if (jitter > 0 && D > 1) {
    
    noise <- round(rnorm(D, mean = 0, sd = jitter * mean(base)))
    block_lengths <- base + noise
    block_lengths <- pmax(min_block_length, block_lengths)
    
    diff_total <- T_obs - sum(block_lengths)
    
    while (diff_total != 0) {
      
      if (diff_total > 0) {
        idx <- sample.int(D, 1)
        block_lengths[idx] <- block_lengths[idx] + 1
        diff_total <- diff_total - 1
      } else {
        eligible <- which(block_lengths > min_block_length)
        
        if (length(eligible) == 0) {
          stop("Cannot reduce block lengths further while preserving min_block_length.")
        }
        
        idx <- sample(eligible, 1)
        block_lengths[idx] <- block_lengths[idx] - 1
        diff_total <- diff_total + 1
      }
    }
  }
  
  as.integer(block_lengths)
}

generate_gap_lengths <- function(D,
                                 block_lengths,
                                 gap_scenario = c("small_gaps", "large_gaps"),
                                 small_gap_factor = 0.15,
                                 large_gap_factor = 1.00,
                                 gap_cv = 0.50,
                                 gap_min = 1,
                                 seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  gap_scenario <- match.arg(gap_scenario)
  
  if (D < 2) return(integer(0))
  
  mean_block_length <- mean(block_lengths)
  
  if (gap_scenario == "small_gaps") {
    gap_mean <- small_gap_factor * mean_block_length
  } else {
    gap_mean <- large_gap_factor * mean_block_length
  }
  
  shape <- 1 / (gap_cv^2)
  rate <- shape / gap_mean
  
  gap_lengths <- ceiling(rgamma(
    n = D - 1,
    shape = shape,
    rate = rate
  ))
  
  gap_lengths <- pmax(gap_min, gap_lengths)
  
  as.integer(gap_lengths)
}

summarise_gap_design <- function(block_lengths, gap_lengths) {
  
  observed_total <- sum(block_lengths)
  gap_total <- sum(gap_lengths)
  full_total <- observed_total + gap_total
  
  data.frame(
    D = length(block_lengths),
    mean_block_length = mean(block_lengths),
    sd_block_length = ifelse(length(block_lengths) > 1, sd(block_lengths), 0),
    mean_gap_length = ifelse(length(gap_lengths) > 0, mean(gap_lengths), 0),
    sd_gap_length = ifelse(length(gap_lengths) > 1, sd(gap_lengths), 0),
    gap_to_block_ratio = ifelse(
      length(gap_lengths) > 0,
      mean(gap_lengths) / mean(block_lengths),
      0
    ),
    observed_total = observed_total,
    gap_total = gap_total,
    full_total = full_total,
    missing_fraction_time = ifelse(full_total > 0, gap_total / full_total, 0)
  )
}

generate_gap_lengths <- function(D,
                                 block_lengths,
                                 gap_scenario = c("small_gaps", "large_gaps"),
                                 small_gap_factor = 0.15,
                                 large_gap_factor = 1.00,
                                 gap_cv = 0.50,
                                 gap_min = 1,
                                 seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  gap_scenario <- match.arg(gap_scenario)
  
  if (D < 2) {
    return(integer(0))
  }
  
  mean_block_length <- mean(block_lengths)
  
  if (gap_scenario == "small_gaps") {
    gap_mean <- small_gap_factor * mean_block_length
  } else {
    gap_mean <- large_gap_factor * mean_block_length
  }
  
  shape <- 1 / (gap_cv^2)
  rate <- shape / gap_mean
  
  gap_lengths <- ceiling(rgamma(
    n = D - 1,
    shape = shape,
    rate = rate
  ))
  
  gap_lengths <- pmax(gap_min, gap_lengths)
  
  as.integer(gap_lengths)
}

summarise_gap_design <- function(block_lengths, gap_lengths) {
  
  observed_total <- sum(block_lengths)
  gap_total <- sum(gap_lengths)
  full_total <- observed_total + gap_total
  
  data.frame(
    D = length(block_lengths),
    mean_block_length = mean(block_lengths),
    sd_block_length = ifelse(length(block_lengths) > 1, sd(block_lengths), 0),
    mean_gap_length = ifelse(length(gap_lengths) > 0, mean(gap_lengths), 0),
    sd_gap_length = ifelse(length(gap_lengths) > 1, sd(gap_lengths), 0),
    gap_to_block_ratio = ifelse(
      length(gap_lengths) > 0,
      mean(gap_lengths) / mean(block_lengths),
      0
    ),
    observed_total = observed_total,
    gap_total = gap_total,
    full_total = full_total,
    missing_fraction_time = ifelse(full_total > 0, gap_total / full_total, 0)
  )
}


# ------------------------------------------------------------
# Cross-gap negative log-score: Skellam emission
# ------------------------------------------------------------

crossgap_nll_skellam <- function(y_full,
                                 expo_full,
                                 observed,
                                 z_pre_full,
                                 A_hat,
                                 mu1,
                                 mu2,
                                 use_power = TRUE,
                                 eps = 1e-12) {
  
  obs_idx <- which(observed)
  
  if (length(obs_idx) < 2) {
    return(NA_real_)
  }
  
  pre_idx <- obs_idx[-length(obs_idx)]
  post_idx <- obs_idx[-1]
  
  d <- post_idx - pre_idx
  is_gap <- d > 1
  
  if (!any(is_gap)) {
    return(NA_real_)
  }
  
  pre_idx <- pre_idx[is_gap]
  post_idx <- post_idx[is_gap]
  d <- d[is_gap]
  
  K <- nrow(A_hat)
  nll <- rep(NA_real_, length(d))
  
  for (g in seq_along(d)) {
    
    z_pre <- z_pre_full[pre_idx[g]]
    
    if (is.na(z_pre)) {
      next
    }
    
    z_pre <- as.integer(z_pre)
    
    if (z_pre < 1 || z_pre > K) {
      next
    }
    
    A_prop <- if (use_power) {
      A_hat %^% d[g]
    } else {
      A_hat
    }
    
    pred_state_prob <- A_prop[z_pre, ]
    
    y_post <- y_full[post_idx[g]]
    expo_post <- expo_full[post_idx[g]]
    
    if (is.na(y_post) || is.na(expo_post) || expo_post <= 0) {
      next
    }
    
    # Use the same Skellam emission log-probabilities used by the EM code.
    log_em <- emission_loglik_vec_skellam(
      y_t = y_post,
      expo_t = expo_post,
      mu1 = mu1,
      mu2 = mu2
    )
    
    valid <- pred_state_prob > 0 & !is.na(log_em) & is.finite(log_em)
    
    if (!any(valid)) {
      next
    }
    
    log_terms <- log(pred_state_prob[valid]) + log_em[valid]
    m <- max(log_terms)
    log_pred_prob <- m + log(sum(exp(log_terms - m)))
    
    nll[g] <- -log_pred_prob
  }
  
  mean(nll, na.rm = TRUE)
}


# ------------------------------------------------------------
# Cross-gap negative log-score: bivariate Poisson emission
# ------------------------------------------------------------

crossgap_nll_pois2 <- function(H_full,
                               S_full,
                               expo_full,
                               observed,
                               z_pre_full,
                               A_hat,
                               lamH,
                               lamS,
                               use_power = TRUE,
                               eps = 1e-12) {
  
  obs_idx <- which(observed)
  
  if (length(obs_idx) < 2) {
    return(NA_real_)
  }
  
  pre_idx <- obs_idx[-length(obs_idx)]
  post_idx <- obs_idx[-1]
  
  d <- post_idx - pre_idx
  is_gap <- d > 1
  
  if (!any(is_gap)) {
    return(NA_real_)
  }
  
  pre_idx <- pre_idx[is_gap]
  post_idx <- post_idx[is_gap]
  d <- d[is_gap]
  
  K <- nrow(A_hat)
  nll <- rep(NA_real_, length(d))
  
  for (g in seq_along(d)) {
    
    z_pre <- z_pre_full[pre_idx[g]]
    
    if (is.na(z_pre)) {
      next
    }
    
    z_pre <- as.integer(z_pre)
    
    if (z_pre < 1 || z_pre > K) {
      next
    }
    
    A_prop <- if (use_power) {
      A_hat %^% d[g]
    } else {
      A_hat
    }
    
    pred_state_prob <- A_prop[z_pre, ]
    
    H_post <- H_full[post_idx[g]]
    S_post <- S_full[post_idx[g]]
    expo_post <- expo_full[post_idx[g]]
    
    if (
      is.na(H_post) || is.na(S_post) ||
      is.na(expo_post) || expo_post <= 0
    ) {
      next
    }
    
    log_em <- rep(NA_real_, K)
    
    for (k in seq_len(K)) {
      log_em[k] <-
        dpois(
          H_post,
          lambda = pmax(expo_post * lamH[k], eps),
          log = TRUE
        ) +
        dpois(
          S_post,
          lambda = pmax(expo_post * lamS[k], eps),
          log = TRUE
        )
    }
    
    # Stable mixture probability:
    # log sum_k pred_state_prob[k] * exp(log_em[k])
    valid <- pred_state_prob > 0 & !is.na(log_em)
    
    if (!any(valid)) {
      next
    }
    
    log_terms <- log(pred_state_prob[valid]) + log_em[valid]
    m <- max(log_terms)
    log_pred_prob <- m + log(sum(exp(log_terms - m)))
    
    nll[g] <- -log_pred_prob
  }
  
  mean(nll, na.rm = TRUE)
}


# ------------------------------------------------------------
# Cross-gap switch Brier score
# ------------------------------------------------------------

crossgap_switch_brier <- function(Z_true_full,
                                  observed,
                                  z_pre_full,
                                  A_hat,
                                  use_power = TRUE) {
  
  obs_idx <- which(observed)
  
  if (length(obs_idx) < 2) {
    return(NA_real_)
  }
  
  pre_idx <- obs_idx[-length(obs_idx)]
  post_idx <- obs_idx[-1]
  
  d <- post_idx - pre_idx
  is_gap <- d > 1
  
  if (!any(is_gap)) {
    return(NA_real_)
  }
  
  pre_idx <- pre_idx[is_gap]
  post_idx <- post_idx[is_gap]
  d <- d[is_gap]
  
  K <- nrow(A_hat)
  brier <- rep(NA_real_, length(d))
  
  for (g in seq_along(d)) {
    
    z_pre_hat <- z_pre_full[pre_idx[g]]
    
    if (is.na(z_pre_hat)) {
      next
    }
    
    z_pre_hat <- as.integer(z_pre_hat)
    
    if (z_pre_hat < 1 || z_pre_hat > K) {
      next
    }
    
    z_true_pre <- Z_true_full[pre_idx[g]]
    z_true_post <- Z_true_full[post_idx[g]]
    
    if (is.na(z_true_pre) || is.na(z_true_post)) {
      next
    }
    
    A_prop <- if (use_power) {
      A_hat %^% d[g]
    } else {
      A_hat
    }
    
    # Predicted probability that the state changes across the gap,
    # conditional on the estimated pre-gap state.
    p_switch <- 1 - A_prop[z_pre_hat, z_pre_hat]
    
    # Observed switch event in the simulated truth.
    switch_obs <- as.numeric(z_true_pre != z_true_post)
    
    brier[g] <- (p_switch - switch_obs)^2
  }
  
  mean(brier, na.rm = TRUE)
}

############ SIMULATIONS

run_one_S1 <- function(job,
                       K_fit_grid,
                       max_iter_sk,
                       max_iter_ps,
                       tol_sk,
                       tol_ps,
                       common_background_type,
                       background_cycles,
                       rate_scale,
                       D_blocks,
                       block_jitter,
                       small_gap_factor,
                       large_gap_factor,
                       gap_cv,
                       gap_min) {
  
  # ------------------------------------------------------------
  # Job settings
  # ------------------------------------------------------------
  
  T_obs <- job$T_length
  K_true <- job$K_true
  replication <- job$replication
  gap_scenario <- as.character(job$gap_scenario)
  background_amplitude <- job$background_amplitude
  job_id <- job$job_id
  
  # If design_id exists, use it to keep the same gap design
  # across different background amplitudes.
  if ("design_id" %in% names(job)) {
    design_id <- job$design_id
  } else {
    design_id <- job_id
  }
  
  # Comparable number of background cycles across T_length.
  background_period_job <- T_obs / background_cycles
  
  set.seed(100000 + design_id)
  
  true_par <- make_true_params_S1(
    K_true = K_true,
    scale = rate_scale
  )
  
  min_block_length <- max(
    10,
    floor(0.40 * T_obs / D_blocks)
  )
  
  block_lengths <- generate_block_lengths(
    T_obs = T_obs,
    D = D_blocks,
    min_block_length = min_block_length,
    jitter = block_jitter
  )
  
  gap_lengths <- generate_gap_lengths(
    D = D_blocks,
    block_lengths = block_lengths,
    gap_scenario = gap_scenario,
    small_gap_factor = small_gap_factor,
    large_gap_factor = large_gap_factor,
    gap_cv = gap_cv,
    gap_min = gap_min
  )
  
  gap_summary <- summarise_gap_design(
    block_lengths = block_lengths,
    gap_lengths = gap_lengths
  )
  
  dat <- simulate_xrb_hmm_blocks(
    D = D_blocks,
    block_lengths = block_lengths,
    gap_lengths = gap_lengths,
    K = K_true,
    muH = true_par$muH,
    muS = true_par$muS,
    A = true_par$A,
    pi1 = true_par$pi1,
    partial_exposure_prob = 0,
    common_background = common_background_type,
    background_amplitude = background_amplitude,
    background_period = background_period_job,
    seed = 200000 + job_id
  )
  
  # ------------------------------------------------------------
  # Full and observed data
  # ------------------------------------------------------------
  
  y <- dat$Omega_full
  H <- dat$H_full
  S <- dat$S_full
  expo <- dat$expo_full
  Z_true <- dat$Z_true_full
  
  obs <- which(dat$observed)
  
  y_obs <- dat$Omega
  H_obs <- dat$H
  S_obs <- dat$S
  expo_obs <- dat$expo
  Z_obs <- dat$Z_true_obs
  
  n_eff <- length(obs)
  missing_fraction <- gap_summary$missing_fraction_time
  mean_exposure <- mean(expo_obs, na.rm = TRUE)
  
  # ------------------------------------------------------------
  # Fit gap-aware models for K in K_fit_grid
  # ------------------------------------------------------------
  
  fits_sk_gap <- vector("list", length(K_fit_grid))
  fits_ps_gap <- vector("list", length(K_fit_grid))
  
  BIC_sk_gap <- rep(NA_real_, length(K_fit_grid))
  ICL_sk_gap <- rep(NA_real_, length(K_fit_grid))
  BIC_ps_gap <- rep(NA_real_, length(K_fit_grid))
  ICL_ps_gap <- rep(NA_real_, length(K_fit_grid))
  
  for (ii in seq_along(K_fit_grid)) {
    
    Kfit <- K_fit_grid[ii]
    
    fits_sk_gap[[ii]] <- tryCatch(
      hmm_skellam_EM_gapaware(
        y = y,
        expo = expo,
        K = Kfit,
        max_iter = max_iter_sk,
        tol = tol_sk,
        verbose = FALSE
      ),
      error = function(e) {
        message(
          "Skellam gap-aware failed | job_id=", job_id,
          " | Kfit=", Kfit,
          " | error=", conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(fits_sk_gap[[ii]])) {
      
      BIC_sk_gap[ii] <- bic_hmm(
        fits_sk_gap[[ii]]$loglik,
        K = Kfit,
        n_eff = n_eff,
        emission = "skellam"
      )
      
      ICL_sk_gap[ii] <- icl_hmm(
        fits_sk_gap[[ii]]$loglik,
        fits_sk_gap[[ii]]$gamma[obs, , drop = FALSE],
        K = Kfit,
        n_eff = n_eff
      )
    }
    
    fits_ps_gap[[ii]] <- tryCatch(
      hmm_pois2_EM_gapaware(
        H = H,
        S = S,
        expo = expo,
        K = Kfit,
        max_iter = max_iter_ps,
        tol = tol_ps,
        verbose = FALSE
      ),
      error = function(e) {
        message(
          "Poisson2 gap-aware failed | job_id=", job_id,
          " | Kfit=", Kfit,
          " | error=", conditionMessage(e)
        )
        NULL
      }
    )
    
    if (!is.null(fits_ps_gap[[ii]])) {
      
      BIC_ps_gap[ii] <- bic_hmm(
        fits_ps_gap[[ii]]$loglik,
        K = Kfit,
        n_eff = n_eff,
        emission = "pois2"
      )
      
      ICL_ps_gap[ii] <- icl_hmm(
        fits_ps_gap[[ii]]$loglik,
        fits_ps_gap[[ii]]$gamma[obs, , drop = FALSE],
        K = Kfit,
        n_eff = n_eff
      )
    }
  }
  
  # ------------------------------------------------------------
  # Select K
  # ------------------------------------------------------------
  
  K_BIC_sk_gap <- safe_select_K(BIC_sk_gap, K_fit_grid)
  K_ICL_sk_gap <- safe_select_K(ICL_sk_gap, K_fit_grid)
  K_BIC_ps_gap <- safe_select_K(BIC_ps_gap, K_fit_grid)
  K_ICL_ps_gap <- safe_select_K(ICL_ps_gap, K_fit_grid)
  
  # BIC-selected fits for main analysis
  fit_sk_BIC <- get_fit_by_selected_K(
    fits = fits_sk_gap,
    K_selected = K_BIC_sk_gap,
    K_fit_grid = K_fit_grid
  )
  
  fit_ps_BIC <- get_fit_by_selected_K(
    fits = fits_ps_gap,
    K_selected = K_BIC_ps_gap,
    K_fit_grid = K_fit_grid
  )
  
  # Oracle-K fits from the fitted gap-aware lists
  idx_true <- match(K_true, K_fit_grid)
  
  if (is.na(idx_true)) {
    fit_sk_gap <- NULL
    fit_ps_gap <- NULL
  } else {
    fit_sk_gap <- fits_sk_gap[[idx_true]]
    fit_ps_gap <- fits_ps_gap[[idx_true]]
  }
  
  # ------------------------------------------------------------
  # Fit compressed/no-gap models with K = K_true
  # ------------------------------------------------------------
  
  fit_sk_nogap <- tryCatch(
    hmm_skellam_EM_gapaware(
      y = y_obs,
      expo = expo_obs,
      K = K_true,
      max_iter = max_iter_sk,
      tol = tol_sk,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  fit_ps_nogap <- tryCatch(
    hmm_pois2_EM_gapaware(
      H = H_obs,
      S = S_obs,
      expo = expo_obs,
      K = K_true,
      max_iter = max_iter_ps,
      tol = tol_ps,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  # ------------------------------------------------------------
  # Initialise oracle-K metrics
  # ------------------------------------------------------------
  
  A_error_sk_gap <- A_error_ps_gap <- NA_real_
  A_error_sk_nogap <- A_error_ps_nogap <- NA_real_
  
  sojourn_error_sk_gap <- sojourn_error_ps_gap <- NA_real_
  sojourn_error_sk_nogap <- sojourn_error_ps_nogap <- NA_real_
  
  contrast_error_sk_gap <- contrast_error_ps_gap <- NA_real_
  contrast_error_sk_nogap <- contrast_error_ps_nogap <- NA_real_
  
  contrast_traj_rmse_sk_gap <- contrast_traj_rmse_ps_gap <- NA_real_
  contrast_traj_rmse_sk_nogap <- contrast_traj_rmse_ps_nogap <- NA_real_
  
  contrast_state_acc_sk_gap <- contrast_state_acc_ps_gap <- NA_real_
  contrast_state_acc_sk_nogap <- contrast_state_acc_ps_nogap <- NA_real_
  
  contrast_sign_acc_sk_gap <- contrast_sign_acc_ps_gap <- NA_real_
  contrast_sign_acc_sk_nogap <- contrast_sign_acc_ps_nogap <- NA_real_
  
  ARI_sk_gap <- ARI_ps_gap <- NA_real_
  ARI_sk_nogap <- ARI_ps_nogap <- NA_real_
  
  crossgap_nll_sk_gap <- crossgap_nll_ps_gap <- NA_real_
  crossgap_nll_sk_nogap <- crossgap_nll_ps_nogap <- NA_real_
  
  crossgap_brier_sk_gap <- crossgap_brier_ps_gap <- NA_real_
  crossgap_brier_sk_nogap <- crossgap_brier_ps_nogap <- NA_real_
  
  # ------------------------------------------------------------
  # Initialise BIC-selected metrics
  # ------------------------------------------------------------
  
  ARI_sk_BIC <- ARI_ps_BIC <- NA_real_
  
  contrast_traj_rmse_sk_BIC <- contrast_traj_rmse_ps_BIC <- NA_real_
  contrast_state_acc_sk_BIC <- contrast_state_acc_ps_BIC <- NA_real_
  contrast_sign_acc_sk_BIC <- contrast_sign_acc_ps_BIC <- NA_real_
  
  crossgap_nll_sk_BIC <- crossgap_nll_ps_BIC <- NA_real_
  crossgap_brier_sk_BIC <- crossgap_brier_ps_BIC <- NA_real_
  
  # ------------------------------------------------------------
  # Oracle-K Skellam gap-aware metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_sk_gap)) {
    
    fit_sk_gap <- sort_skellam_states_by_contrast(fit_sk_gap)
    est_contrast_sk_gap <- fit_sk_gap$mu1 - fit_sk_gap$mu2
    
    A_error_sk_gap <- sqrt(sum((fit_sk_gap$A - true_par$A)^2))
    sojourn_error_sk_gap <- safe_sojourn_error(fit_sk_gap$A, true_par$A)
    contrast_error_sk_gap <- sqrt(mean((est_contrast_sk_gap - true_par$contrast)^2))
    
    z_sk_gap <- viterbi_skellam_maxplus(
      y = y,
      expo = expo,
      pi1 = fit_sk_gap$pi1,
      A = fit_sk_gap$A,
      mu1 = fit_sk_gap$mu1,
      mu2 = fit_sk_gap$mu2,
      fill_gaps = "NA"
    )
    
    ARI_sk_gap <- safe_ARI(
      Z_true[obs],
      z_sk_gap[obs]
    )
    
    contrast_traj_rmse_sk_gap <- contrast_trajectory_rmse(
      Z_true[obs],
      z_sk_gap[obs],
      true_par$contrast,
      est_contrast_sk_gap
    )
    
    contrast_state_acc_sk_gap <- contrast_state_accuracy(
      Z_true[obs],
      z_sk_gap[obs],
      true_par$contrast,
      est_contrast_sk_gap
    )
    
    contrast_sign_acc_sk_gap <- contrast_sign_accuracy(
      Z_true[obs],
      z_sk_gap[obs],
      true_par$contrast,
      est_contrast_sk_gap
    )
    
    crossgap_nll_sk_gap <- crossgap_nll_skellam(
      y_full = y,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_sk_gap,
      A_hat = fit_sk_gap$A,
      mu1 = fit_sk_gap$mu1,
      mu2 = fit_sk_gap$mu2,
      use_power = TRUE
    )
    
    crossgap_brier_sk_gap <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_sk_gap,
      A_hat = fit_sk_gap$A,
      use_power = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # Oracle-K Poisson gap-aware metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_ps_gap)) {
    
    fit_ps_gap <- sort_pois2_states_by_contrast(fit_ps_gap)
    est_contrast_ps_gap <- fit_ps_gap$lamH - fit_ps_gap$lamS
    
    A_error_ps_gap <- sqrt(sum((fit_ps_gap$A - true_par$A)^2))
    sojourn_error_ps_gap <- safe_sojourn_error(fit_ps_gap$A, true_par$A)
    contrast_error_ps_gap <- sqrt(mean((est_contrast_ps_gap - true_par$contrast)^2))
    
    z_ps_gap <- viterbi_pois2_maxplus(
      H = H,
      S = S,
      expo = expo,
      pi1 = fit_ps_gap$pi1,
      A = fit_ps_gap$A,
      lamH = fit_ps_gap$lamH,
      lamS = fit_ps_gap$lamS,
      fill_gaps = "NA"
    )
    
    ARI_ps_gap <- safe_ARI(
      Z_true[obs],
      z_ps_gap[obs]
    )
    
    contrast_traj_rmse_ps_gap <- contrast_trajectory_rmse(
      Z_true[obs],
      z_ps_gap[obs],
      true_par$contrast,
      est_contrast_ps_gap
    )
    
    contrast_state_acc_ps_gap <- contrast_state_accuracy(
      Z_true[obs],
      z_ps_gap[obs],
      true_par$contrast,
      est_contrast_ps_gap
    )
    
    contrast_sign_acc_ps_gap <- contrast_sign_accuracy(
      Z_true[obs],
      z_ps_gap[obs],
      true_par$contrast,
      est_contrast_ps_gap
    )
    
    crossgap_nll_ps_gap <- crossgap_nll_pois2(
      H_full = H,
      S_full = S,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_ps_gap,
      A_hat = fit_ps_gap$A,
      lamH = fit_ps_gap$lamH,
      lamS = fit_ps_gap$lamS,
      use_power = TRUE
    )
    
    crossgap_brier_ps_gap <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_ps_gap,
      A_hat = fit_ps_gap$A,
      use_power = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # BIC-selected Skellam metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_sk_BIC)) {
    
    fit_sk_BIC <- sort_skellam_states_by_contrast(fit_sk_BIC)
    est_contrast_sk_BIC <- fit_sk_BIC$mu1 - fit_sk_BIC$mu2
    
    z_sk_BIC <- viterbi_skellam_maxplus(
      y = y,
      expo = expo,
      pi1 = fit_sk_BIC$pi1,
      A = fit_sk_BIC$A,
      mu1 = fit_sk_BIC$mu1,
      mu2 = fit_sk_BIC$mu2,
      fill_gaps = "NA"
    )
    
    ARI_sk_BIC <- safe_ARI(
      Z_true[obs],
      z_sk_BIC[obs]
    )
    
    contrast_traj_rmse_sk_BIC <- contrast_trajectory_rmse(
      Z_true[obs],
      z_sk_BIC[obs],
      true_par$contrast,
      est_contrast_sk_BIC
    )
    
    contrast_state_acc_sk_BIC <- contrast_state_accuracy(
      Z_true[obs],
      z_sk_BIC[obs],
      true_par$contrast,
      est_contrast_sk_BIC
    )
    
    contrast_sign_acc_sk_BIC <- contrast_sign_accuracy(
      Z_true[obs],
      z_sk_BIC[obs],
      true_par$contrast,
      est_contrast_sk_BIC
    )
    
    crossgap_nll_sk_BIC <- crossgap_nll_skellam(
      y_full = y,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_sk_BIC,
      A_hat = fit_sk_BIC$A,
      mu1 = fit_sk_BIC$mu1,
      mu2 = fit_sk_BIC$mu2,
      use_power = TRUE
    )
    
    crossgap_brier_sk_BIC <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_sk_BIC,
      A_hat = fit_sk_BIC$A,
      use_power = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # BIC-selected Poisson metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_ps_BIC)) {
    
    fit_ps_BIC <- sort_pois2_states_by_contrast(fit_ps_BIC)
    est_contrast_ps_BIC <- fit_ps_BIC$lamH - fit_ps_BIC$lamS
    
    z_ps_BIC <- viterbi_pois2_maxplus(
      H = H,
      S = S,
      expo = expo,
      pi1 = fit_ps_BIC$pi1,
      A = fit_ps_BIC$A,
      lamH = fit_ps_BIC$lamH,
      lamS = fit_ps_BIC$lamS,
      fill_gaps = "NA"
    )
    
    ARI_ps_BIC <- safe_ARI(
      Z_true[obs],
      z_ps_BIC[obs]
    )
    
    contrast_traj_rmse_ps_BIC <- contrast_trajectory_rmse(
      Z_true[obs],
      z_ps_BIC[obs],
      true_par$contrast,
      est_contrast_ps_BIC
    )
    
    contrast_state_acc_ps_BIC <- contrast_state_accuracy(
      Z_true[obs],
      z_ps_BIC[obs],
      true_par$contrast,
      est_contrast_ps_BIC
    )
    
    contrast_sign_acc_ps_BIC <- contrast_sign_accuracy(
      Z_true[obs],
      z_ps_BIC[obs],
      true_par$contrast,
      est_contrast_ps_BIC
    )
    
    crossgap_nll_ps_BIC <- crossgap_nll_pois2(
      H_full = H,
      S_full = S,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_ps_BIC,
      A_hat = fit_ps_BIC$A,
      lamH = fit_ps_BIC$lamH,
      lamS = fit_ps_BIC$lamS,
      use_power = TRUE
    )
    
    crossgap_brier_ps_BIC <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_ps_BIC,
      A_hat = fit_ps_BIC$A,
      use_power = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # Oracle-K compressed/no-gap Skellam metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_sk_nogap)) {
    
    fit_sk_nogap <- sort_skellam_states_by_contrast(fit_sk_nogap)
    est_contrast_sk_nogap <- fit_sk_nogap$mu1 - fit_sk_nogap$mu2
    
    A_error_sk_nogap <- sqrt(sum((fit_sk_nogap$A - true_par$A)^2))
    sojourn_error_sk_nogap <- safe_sojourn_error(fit_sk_nogap$A, true_par$A)
    contrast_error_sk_nogap <- sqrt(mean((est_contrast_sk_nogap - true_par$contrast)^2))
    
    z_sk_nogap <- viterbi_skellam_maxplus(
      y = y_obs,
      expo = expo_obs,
      pi1 = fit_sk_nogap$pi1,
      A = fit_sk_nogap$A,
      mu1 = fit_sk_nogap$mu1,
      mu2 = fit_sk_nogap$mu2,
      fill_gaps = "NA"
    )
    
    ARI_sk_nogap <- safe_ARI(
      Z_obs,
      z_sk_nogap
    )
    
    contrast_traj_rmse_sk_nogap <- contrast_trajectory_rmse(
      Z_obs,
      z_sk_nogap,
      true_par$contrast,
      est_contrast_sk_nogap
    )
    
    contrast_state_acc_sk_nogap <- contrast_state_accuracy(
      Z_obs,
      z_sk_nogap,
      true_par$contrast,
      est_contrast_sk_nogap
    )
    
    contrast_sign_acc_sk_nogap <- contrast_sign_accuracy(
      Z_obs,
      z_sk_nogap,
      true_par$contrast,
      est_contrast_sk_nogap
    )
    
    z_sk_nogap_full <- rep(NA_integer_, length(y))
    z_sk_nogap_full[obs] <- z_sk_nogap
    
    crossgap_nll_sk_nogap <- crossgap_nll_skellam(
      y_full = y,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_sk_nogap_full,
      A_hat = fit_sk_nogap$A,
      mu1 = fit_sk_nogap$mu1,
      mu2 = fit_sk_nogap$mu2,
      use_power = FALSE
    )
    
    crossgap_brier_sk_nogap <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_sk_nogap_full,
      A_hat = fit_sk_nogap$A,
      use_power = FALSE
    )
  }
  
  # ------------------------------------------------------------
  # Oracle-K compressed/no-gap Poisson metrics
  # ------------------------------------------------------------
  
  if (!is.null(fit_ps_nogap)) {
    
    fit_ps_nogap <- sort_pois2_states_by_contrast(fit_ps_nogap)
    est_contrast_ps_nogap <- fit_ps_nogap$lamH - fit_ps_nogap$lamS
    
    A_error_ps_nogap <- sqrt(sum((fit_ps_nogap$A - true_par$A)^2))
    sojourn_error_ps_nogap <- safe_sojourn_error(fit_ps_nogap$A, true_par$A)
    contrast_error_ps_nogap <- sqrt(mean((est_contrast_ps_nogap - true_par$contrast)^2))
    
    z_ps_nogap <- viterbi_pois2_maxplus(
      H = H_obs,
      S = S_obs,
      expo = expo_obs,
      pi1 = fit_ps_nogap$pi1,
      A = fit_ps_nogap$A,
      lamH = fit_ps_nogap$lamH,
      lamS = fit_ps_nogap$lamS,
      fill_gaps = "NA"
    )
    
    ARI_ps_nogap <- safe_ARI(
      Z_obs,
      z_ps_nogap
    )
    
    contrast_traj_rmse_ps_nogap <- contrast_trajectory_rmse(
      Z_obs,
      z_ps_nogap,
      true_par$contrast,
      est_contrast_ps_nogap
    )
    
    contrast_state_acc_ps_nogap <- contrast_state_accuracy(
      Z_obs,
      z_ps_nogap,
      true_par$contrast,
      est_contrast_ps_nogap
    )
    
    contrast_sign_acc_ps_nogap <- contrast_sign_accuracy(
      Z_obs,
      z_ps_nogap,
      true_par$contrast,
      est_contrast_ps_nogap
    )
    
    z_ps_nogap_full <- rep(NA_integer_, length(H))
    z_ps_nogap_full[obs] <- z_ps_nogap
    
    crossgap_nll_ps_nogap <- crossgap_nll_pois2(
      H_full = H,
      S_full = S,
      expo_full = expo,
      observed = dat$observed,
      z_pre_full = z_ps_nogap_full,
      A_hat = fit_ps_nogap$A,
      lamH = fit_ps_nogap$lamH,
      lamS = fit_ps_nogap$lamS,
      use_power = FALSE
    )
    
    crossgap_brier_ps_nogap <- crossgap_switch_brier(
      Z_true_full = Z_true,
      observed = dat$observed,
      z_pre_full = z_ps_nogap_full,
      A_hat = fit_ps_nogap$A,
      use_power = FALSE
    )
  }
  
  # ------------------------------------------------------------
  # Return one-row data.frame
  # ------------------------------------------------------------
  
  data.frame(
    scenario = "S1_block_gaps_common_mode_grid",
    gap_scenario = gap_scenario,
    background_type = common_background_type,
    background_amplitude = background_amplitude,
    background_period = background_period_job,
    background_cycles = background_cycles,
    
    T_length = T_obs,
    T_full = dat$T_full,
    K_true = K_true,
    replication = replication,
    job_id = job_id,
    design_id = design_id,
    
    rate_scale = rate_scale,
    D_blocks = D_blocks,
    
    mean_block_length = gap_summary$mean_block_length,
    sd_block_length = gap_summary$sd_block_length,
    mean_gap_length = gap_summary$mean_gap_length,
    sd_gap_length = gap_summary$sd_gap_length,
    gap_to_block_ratio = gap_summary$gap_to_block_ratio,
    observed_total = gap_summary$observed_total,
    gap_total = gap_summary$gap_total,
    
    n_eff = n_eff,
    missing_fraction = missing_fraction,
    mean_exposure = mean_exposure,
    
    K_BIC_sk_gap = K_BIC_sk_gap,
    K_ICL_sk_gap = K_ICL_sk_gap,
    K_BIC_ps_gap = K_BIC_ps_gap,
    K_ICL_ps_gap = K_ICL_ps_gap,
    
    BIC_sk_gap_over = K_BIC_sk_gap > K_true,
    ICL_sk_gap_over = K_ICL_sk_gap > K_true,
    BIC_ps_gap_over = K_BIC_ps_gap > K_true,
    ICL_ps_gap_over = K_ICL_ps_gap > K_true,
    
    BIC_sk_gap_under = K_BIC_sk_gap < K_true,
    ICL_sk_gap_under = K_ICL_sk_gap < K_true,
    BIC_ps_gap_under = K_BIC_ps_gap < K_true,
    ICL_ps_gap_under = K_ICL_ps_gap < K_true,
    
    # Oracle-K gap-aware/no-gap parameter recovery
    A_error_sk_gap = A_error_sk_gap,
    A_error_ps_gap = A_error_ps_gap,
    A_error_sk_nogap = A_error_sk_nogap,
    A_error_ps_nogap = A_error_ps_nogap,
    
    sojourn_error_sk_gap = sojourn_error_sk_gap,
    sojourn_error_ps_gap = sojourn_error_ps_gap,
    sojourn_error_sk_nogap = sojourn_error_sk_nogap,
    sojourn_error_ps_nogap = sojourn_error_ps_nogap,
    
    contrast_error_sk_gap = contrast_error_sk_gap,
    contrast_error_ps_gap = contrast_error_ps_gap,
    contrast_error_sk_nogap = contrast_error_sk_nogap,
    contrast_error_ps_nogap = contrast_error_ps_nogap,
    
    # Oracle-K recovery
    contrast_traj_rmse_sk_gap = contrast_traj_rmse_sk_gap,
    contrast_traj_rmse_ps_gap = contrast_traj_rmse_ps_gap,
    contrast_traj_rmse_sk_nogap = contrast_traj_rmse_sk_nogap,
    contrast_traj_rmse_ps_nogap = contrast_traj_rmse_ps_nogap,
    
    contrast_state_acc_sk_gap = contrast_state_acc_sk_gap,
    contrast_state_acc_ps_gap = contrast_state_acc_ps_gap,
    contrast_state_acc_sk_nogap = contrast_state_acc_sk_nogap,
    contrast_state_acc_ps_nogap = contrast_state_acc_ps_nogap,
    
    contrast_sign_acc_sk_gap = contrast_sign_acc_sk_gap,
    contrast_sign_acc_ps_gap = contrast_sign_acc_ps_gap,
    contrast_sign_acc_sk_nogap = contrast_sign_acc_sk_nogap,
    contrast_sign_acc_ps_nogap = contrast_sign_acc_ps_nogap,
    
    ARI_sk_gap = ARI_sk_gap,
    ARI_ps_gap = ARI_ps_gap,
    ARI_sk_nogap = ARI_sk_nogap,
    ARI_ps_nogap = ARI_ps_nogap,
    
    delta_contrast_state_acc_sk = contrast_state_acc_sk_gap - contrast_state_acc_sk_nogap,
    delta_contrast_state_acc_ps = contrast_state_acc_ps_gap - contrast_state_acc_ps_nogap,
    
    delta_contrast_traj_rmse_sk = contrast_traj_rmse_sk_nogap - contrast_traj_rmse_sk_gap,
    delta_contrast_traj_rmse_ps = contrast_traj_rmse_ps_nogap - contrast_traj_rmse_ps_gap,
    
    delta_ARI_sk = ARI_sk_gap - ARI_sk_nogap,
    delta_ARI_ps = ARI_ps_gap - ARI_ps_nogap,
    
    # Oracle-K cross-gap gap-aware/no-gap
    crossgap_nll_sk_gap = crossgap_nll_sk_gap,
    crossgap_nll_ps_gap = crossgap_nll_ps_gap,
    crossgap_nll_sk_nogap = crossgap_nll_sk_nogap,
    crossgap_nll_ps_nogap = crossgap_nll_ps_nogap,
    
    crossgap_brier_sk_gap = crossgap_brier_sk_gap,
    crossgap_brier_ps_gap = crossgap_brier_ps_gap,
    crossgap_brier_sk_nogap = crossgap_brier_sk_nogap,
    crossgap_brier_ps_nogap = crossgap_brier_ps_nogap,
    
    delta_crossgap_nll_sk = crossgap_nll_sk_nogap - crossgap_nll_sk_gap,
    delta_crossgap_nll_ps = crossgap_nll_ps_nogap - crossgap_nll_ps_gap,
    
    delta_crossgap_brier_sk = crossgap_brier_sk_nogap - crossgap_brier_sk_gap,
    delta_crossgap_brier_ps = crossgap_brier_ps_nogap - crossgap_brier_ps_gap,
    
    # BIC-selected recovery metrics for main workflow
    ARI_sk_BIC = ARI_sk_BIC,
    ARI_ps_BIC = ARI_ps_BIC,
    
    contrast_traj_rmse_sk_BIC = contrast_traj_rmse_sk_BIC,
    contrast_traj_rmse_ps_BIC = contrast_traj_rmse_ps_BIC,
    
    contrast_state_acc_sk_BIC = contrast_state_acc_sk_BIC,
    contrast_state_acc_ps_BIC = contrast_state_acc_ps_BIC,
    
    contrast_sign_acc_sk_BIC = contrast_sign_acc_sk_BIC,
    contrast_sign_acc_ps_BIC = contrast_sign_acc_ps_BIC,
    
    crossgap_nll_sk_BIC = crossgap_nll_sk_BIC,
    crossgap_nll_ps_BIC = crossgap_nll_ps_BIC,
    
    crossgap_brier_sk_BIC = crossgap_brier_sk_BIC,
    crossgap_brier_ps_BIC = crossgap_brier_ps_BIC,
    
    delta_ARI_BIC_sk_minus_ps = ARI_sk_BIC - ARI_ps_BIC,
    
    delta_contrast_acc_BIC_sk_minus_ps =
      contrast_state_acc_sk_BIC - contrast_state_acc_ps_BIC,
    
    delta_sign_acc_BIC_sk_minus_ps =
      contrast_sign_acc_sk_BIC - contrast_sign_acc_ps_BIC,
    
    delta_contrast_rmse_BIC_ps_minus_sk =
      contrast_traj_rmse_ps_BIC - contrast_traj_rmse_sk_BIC,
    
    delta_crossgap_nll_BIC_ps_minus_sk =
      crossgap_nll_ps_BIC - crossgap_nll_sk_BIC,
    
    delta_crossgap_brier_BIC_ps_minus_sk =
      crossgap_brier_ps_BIC - crossgap_brier_sk_BIC,
    
    iter_sk_gap_Ktrue = get_iter(fit_sk_gap),
    iter_ps_gap_Ktrue = get_iter(fit_ps_gap),
    iter_sk_nogap_Ktrue = get_iter(fit_sk_nogap),
    iter_ps_nogap_Ktrue = get_iter(fit_ps_nogap),
    
    conv_sk_gap_Ktrue = get_converged(fit_sk_gap),
    conv_ps_gap_Ktrue = get_converged(fit_ps_gap),
    conv_sk_nogap_Ktrue = get_converged(fit_sk_nogap),
    conv_ps_nogap_Ktrue = get_converged(fit_ps_nogap)
  )
}

# ===== VERBATIM SOURCE: R/hmm_original/emissions.R =====
# Emission models --------------------------------------------------------------

#' Skellam log-density with exposure offset
#'
#' @param y Integer difference count.
#' @param mu1 First intrinsic rate.
#' @param mu2 Second intrinsic rate.
#' @param expo Effective exposure in the bin.
#' @export
#' 
# log_skellam_offset <- function(y, mu1, mu2, expo) {
#   if (expo <= 0 || is.na(y)) return(0.0)
#   m1 <- mu1 * expo
#   m2 <- mu2 * expo
#   if (!is.finite(m1) || !is.finite(m2) || m1 <= 0 || m2 <= 0) return(-1e300)
#   ll <- extraDistr::dskellam(y, mu1 = m1, mu2 = m2, log = TRUE)
#   if (!is.finite(ll)) -1e300 else ll
# }
log_skellam_offset <- function(y, mu1, mu2, expo) {
  if (expo <= 0 || is.na(y)) return(0.0)
  
  m1 <- mu1 * expo
  m2 <- mu2 * expo
  
  if (!is.finite(m1) || !is.finite(m2) || m1 <= 0 || m2 <= 0) {
    return(-1e300)
  }
  
  ll <- extraDistr::dskellam(
    y,
    mu1 = m1,
    mu2 = m2,
    log = TRUE
  )
  
  if (!is.finite(ll) || ll <= -1e290) {
    ll_norm <- dnorm(
      y,
      mean = m1 - m2,
      sd = sqrt(m1 + m2),
      log = TRUE
    )
    
    if (is.finite(ll_norm)) {
      return(ll_norm)
    } else {
      return(-1e300)
    }
  }
  
  ll
}

emission_loglik_vec_skellam <- function(y_t, expo_t, mu1, mu2) {
  vapply(seq_along(mu1), function(k) log_skellam_offset(y_t, mu1[k], mu2[k], expo_t), numeric(1))
}

#' Bivariate independent Poisson log-density with exposure offset
#'
#' @param H First count series, usually hard-band counts.
#' @param S Second count series, usually soft-band counts.
#' @param lamH First intrinsic rate.
#' @param lamS Second intrinsic rate.
#' @param expo Effective exposure in the bin.
#' @export
log_pois2_offset <- function(H, S, lamH, lamS, expo) {
  if (expo <= 0 || is.na(H) || is.na(S)) return(0.0)
  mH <- lamH * expo
  mS <- lamS * expo
  if (!is.finite(mH) || !is.finite(mS) || mH <= 0 || mS <= 0) return(-1e300)
  ll <- stats::dpois(H, mH, log = TRUE) + stats::dpois(S, mS, log = TRUE)
  if (!is.finite(ll)) -1e300 else ll
}

emission_loglik_vec_pois2 <- function(H_t, S_t, expo_t, lamH, lamS) {
  vapply(seq_along(lamH), function(k) log_pois2_offset(H_t, S_t, lamH[k], lamS[k], expo_t), numeric(1))
}

# Correlated bivariate Poisson emission ----------------------------------------

log_bivpois_offset <- function(H, S, lamH, lamS, lamC, expo) {
  if (expo <= 0 || is.na(H) || is.na(S)) return(0.0)
  
  mH <- lamH * expo
  mS <- lamS * expo
  mC <- lamC * expo
  
  if (!is.finite(mH) || !is.finite(mS) || !is.finite(mC) ||
      mH <= 0 || mS <= 0 || mC <= 0) {
    return(-1e300)
  }
  
  h <- as.integer(H)
  s <- as.integer(S)
  
  if (h < 0 || s < 0) return(-1e300)
  
  x <- 0:min(h, s)
  
  log_terms <- 
    x * log(mC) - lgamma(x + 1) +
    (h - x) * log(mH) - lgamma(h - x + 1) +
    (s - x) * log(mS) - lgamma(s - x + 1)
  
  ll <- -(mH + mS + mC) + logSumExp(log_terms)
  
  if (!is.finite(ll)) -1e300 else ll
}

emission_loglik_vec_bivpois <- function(H_t, S_t, expo_t, lamH, lamS, lamC) {
  vapply(
    seq_along(lamH),
    function(k) log_bivpois_offset(
      H = H_t,
      S = S_t,
      lamH = lamH[k],
      lamS = lamS[k],
      lamC = lamC[k],
      expo = expo_t
    ),
    numeric(1)
  )
}

expected_common_bivpois <- function(H, S, lamH, lamS, lamC, expo) {
  if (expo <= 0 || is.na(H) || is.na(S)) return(NA_real_)
  
  mH <- lamH * expo
  mS <- lamS * expo
  mC <- lamC * expo
  
  if (!is.finite(mH) || !is.finite(mS) || !is.finite(mC) ||
      mH <= 0 || mS <= 0 || mC <= 0) {
    return(0.0)
  }
  
  h <- as.integer(H)
  s <- as.integer(S)
  
  if (h < 0 || s < 0) return(NA_real_)
  
  x <- 0:min(h, s)
  
  logw <- 
    x * log(mC) - lgamma(x + 1) +
    (h - x) * log(mH) - lgamma(h - x + 1) +
    (s - x) * log(mS) - lgamma(s - x + 1)
  
  lw <- logw - logSumExp(logw)
  p <- exp(lw)
  
  sum(x * p)
}

# ===== VERBATIM SOURCE: R/hmm_original/forward_backward.R =====
# Gap-aware forward-backward recursions ---------------------------------------

forward_hybrid <- function(logB, pi1, A, miss) {
  T <- nrow(logB); K <- ncol(logB)
  A_pow <- make_A_power_cache(A)
  logA_pow <- function(m) safe_log(A_pow(m))
  logalpha <- matrix(NA_real_, T, K)
  logalpha[1, ] <- safe_log(pi1) + logB[1, ]
  t <- 1L
  while (t < T) {
    if (miss[t + 1L]) {
      s <- t + 1L
      e <- s
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      L <- e - t
      logAt <- logA_pow(L)
      for (j in 1:K) logalpha[e, j] <- logSumExp(logalpha[t, ] + logAt[, j])
      t <- e
    } else {
      tt <- t + 1L
      logA <- safe_log(A)
      for (j in 1:K) logalpha[tt, j] <- logB[tt, j] + logSumExp(logalpha[t, ] + logA[, j])
      t <- tt
    }
  }
  if (any(is.na(logalpha[T, ]))) {
    last <- max(which(!is.na(logalpha[, 1])))
    if (last < T) {
      L <- T - last
      logAt <- logA_pow(L)
      for (j in 1:K) logalpha[T, j] <- logSumExp(logalpha[last, ] + logAt[, j])
    }
  }
  list(logalpha = logalpha, loglik = logSumExp(logalpha[T, ]))
}

backward_hybrid <- function(logB, A, miss) {
  T <- nrow(logB); K <- ncol(logB)
  A_pow <- make_A_power_cache(A)
  logA_pow <- function(m) safe_log(A_pow(m))
  logbeta <- matrix(NA_real_, T, K)
  logbeta[T, ] <- 0
  t <- T
  while (t > 1L) {
    if (miss[t]) {
      e <- t
      s <- e
      while (s >= 1L && miss[s]) s <- s - 1L
      s <- s + 1L
      prev <- s - 1L
      if (prev < 1L) break
      L <- e - prev
      logAt <- logA_pow(L)
      for (i in 1:K) logbeta[prev, i] <- logSumExp(logAt[i, ] + logbeta[e, ])
      t <- prev
    } else {
      prev <- t - 1L
      logA <- safe_log(A)
      for (i in 1:K) logbeta[prev, i] <- logSumExp(logA[i, ] + logB[t, ] + logbeta[t, ])
      t <- prev
    }
  }
  if (any(is.na(logbeta[1, ]))) logbeta[1, ] <- 0
  logbeta
}

expected_transitions_in_block <- function(alpha_prev, beta_e, A, L) {
  K <- length(alpha_prev)
  if (L <= 0) return(matrix(0, K, K))
  AL <- A %^% L
  denom <- as.numeric(alpha_prev %*% AL %*% beta_e)
  if (!is.finite(denom) || denom <= 0) return(matrix(0, K, K))
  b_list <- vector("list", L)
  b <- beta_e
  for (u in 1:L) {
    b_list[[u]] <- b
    b <- A %*% b
  }
  A_cnt <- matrix(0, K, K)
  a <- alpha_prev
  for (r in 0:(L - 1L)) {
    bR <- b_list[[L - r]]
    A_cnt <- A_cnt + (A * (a %o% as.vector(bR))) / denom
    a <- as.vector(a %*% A)
  }
  A_cnt
}

# ===== VERBATIM SOURCE: R/hmm_original/maxplus.R =====
# Max-Plus utilities -----------------------------------------------------------

max_plus_mult <- function(X, Y) {
  N <- nrow(X); M <- ncol(Y)
  res <- matrix(-Inf, N, M)
  for (i in 1:N) {
    for (j in 1:M) res[i, j] <- max(X[i, ] + Y[, j])
  }
  res
}

#' Fast matrix power in the Max-Plus semiring
#'
#' @param logA Logarithmic transition matrix.
#' @param m Positive integer power.
#' @export
max_plus_power <- function(logA, m) {
  if (m == 0) {
    K <- nrow(logA)
    res <- matrix(-Inf, K, K)
    diag(res) <- 0
    return(res)
  }
  if (m == 1) return(logA)
  K <- nrow(logA)
  res <- matrix(-Inf, K, K)
  diag(res) <- 0
  base <- logA
  mm <- as.integer(m)
  while (mm > 0) {
    if (mm %% 2 == 1) res <- max_plus_mult(res, base)
    mm <- mm %/% 2
    if (mm > 0) base <- max_plus_mult(base, base)
  }
  res
}

# ===== VERBATIM SOURCE: R/hmm_original/viterbi.R =====
# Viterbi decoding -------------------------------------------------------------

#' Max-Plus Viterbi decoding for the exposure-adjusted Skellam-HMM
#'
#' @param y Integer difference series, e.g. H - S.
#' @param expo Effective exposure per bin.
#' @param pi1 Initial distribution.
#' @param A Transition matrix.
#' @param mu1 First Skellam intrinsic rate.
#' @param mu2 Second Skellam intrinsic rate.
#' @param fill_gaps How to fill internal unobserved bins. Currently `constant` fills them with the left boundary state for plotting.
#' @export
viterbi_skellam_maxplus <- function(y, expo, pi1, A, mu1, mu2,
                                    fill_gaps = c("constant", "NA")) {
  fill_gaps <- match.arg(fill_gaps)
  stopifnot(length(y) == length(expo))
  T <- length(y); K <- length(pi1)
  miss <- is_missing_skellam(y, expo)
  logA <- safe_log(A)

  logB <- matrix(0, nrow = T, ncol = K)
  for (t in 1:T) if (!miss[t]) logB[t, ] <- emission_loglik_vec_skellam(y[t], expo[t], mu1, mu2)

  delta <- matrix(NA_real_, T, K)
  psi <- matrix(0L, T, K)
  delta[1, ] <- safe_log(pi1) + logB[1, ]

  t <- 1L
  while (t < T) {
    if (miss[t + 1L]) {
      s <- t + 1L
      e <- s
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      L <- e - t
      logAt <- max_plus_power(logA, L)
      for (j in 1:K) {
        tmp <- delta[t, ] + logAt[, j]
        psi[e, j] <- which.max(tmp)
        delta[e, j] <- max(tmp)
      }
      t <- e
    } else {
      tt <- t + 1L
      for (j in 1:K) {
        tmp <- delta[t, ] + logA[, j]
        psi[tt, j] <- which.max(tmp)
        delta[tt, j] <- max(tmp) + logB[tt, j]
      }
      t <- tt
    }
  }

  z <- integer(T)
  if (fill_gaps == "NA") z[] <- NA_integer_
  last <- max(which(!is.na(delta[, 1])))
  z[last] <- which.max(delta[last, ])
  tt <- last
  while (tt > 1L) {
    if (miss[tt]) {
      e <- tt
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      s <- tt
      while (s >= 1L && miss[s]) s <- s - 1L
      s <- s + 1L
      prev <- s - 1L
      if (prev < 1L) break
      z[prev] <- psi[e, z[e]]
      if (fill_gaps == "constant") z[s:e] <- z[prev]
      tt <- prev
    } else {
      z[tt - 1L] <- psi[tt, z[tt]]
      tt <- tt - 1L
    }
  }
  z
}

#' Max-Plus Viterbi decoding for the exposure-adjusted bivariate Poisson HMM
#' @export
viterbi_pois2_maxplus <- function(H, S, expo, pi1, A, lamH, lamS,
                                  fill_gaps = c("constant", "NA")) {
  fill_gaps <- match.arg(fill_gaps)
  stopifnot(length(H) == length(S), length(H) == length(expo))
  T <- length(H); K <- length(pi1)
  miss <- is_missing_pois2(H, S, expo)
  logA <- safe_log(A)

  logB <- matrix(0, nrow = T, ncol = K)
  for (t in 1:T) if (!miss[t]) logB[t, ] <- emission_loglik_vec_pois2(H[t], S[t], expo[t], lamH, lamS)

  delta <- matrix(NA_real_, T, K)
  psi <- matrix(0L, T, K)
  delta[1, ] <- safe_log(pi1) + logB[1, ]
  t <- 1L
  while (t < T) {
    if (miss[t + 1L]) {
      s <- t + 1L
      e <- s
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      L <- e - t
      logAt <- max_plus_power(logA, L)
      for (j in 1:K) {
        tmp <- delta[t, ] + logAt[, j]
        psi[e, j] <- which.max(tmp)
        delta[e, j] <- max(tmp)
      }
      t <- e
    } else {
      tt <- t + 1L
      for (j in 1:K) {
        tmp <- delta[t, ] + logA[, j]
        psi[tt, j] <- which.max(tmp)
        delta[tt, j] <- max(tmp) + logB[tt, j]
      }
      t <- tt
    }
  }

  z <- integer(T)
  if (fill_gaps == "NA") z[] <- NA_integer_
  last <- max(which(!is.na(delta[, 1])))
  z[last] <- which.max(delta[last, ])
  tt <- last
  while (tt > 1L) {
    if (miss[tt]) {
      e <- tt
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      s <- tt
      while (s >= 1L && miss[s]) s <- s - 1L
      s <- s + 1L
      prev <- s - 1L
      if (prev < 1L) break
      z[prev] <- psi[e, z[e]]
      if (fill_gaps == "constant") z[s:e] <- z[prev]
      tt <- prev
    } else {
      z[tt - 1L] <- psi[tt, z[tt]]
      tt <- tt - 1L
    }
  }
  z
}

viterbi_bivpois_maxplus <- function(H, S, expo, pi1, A, lamH, lamS, lamC,
                                    fill_gaps = c("constant", "NA")) {
  
  fill_gaps <- match.arg(fill_gaps)
  
  stopifnot(length(H) == length(S), length(H) == length(expo))
  
  T <- length(H)
  K <- length(pi1)
  
  miss <- is_missing_pois2(H, S, expo)
  logA <- safe_log(A)
  
  logB <- matrix(0, nrow = T, ncol = K)
  
  for (t in 1:T) {
    if (!miss[t]) {
      logB[t, ] <- emission_loglik_vec_bivpois(
        H_t = H[t],
        S_t = S[t],
        expo_t = expo[t],
        lamH = lamH,
        lamS = lamS,
        lamC = lamC
      )
    }
  }
  
  delta <- matrix(NA_real_, T, K)
  psi <- matrix(0L, T, K)
  
  delta[1, ] <- safe_log(pi1) + logB[1, ]
  
  t <- 1L
  
  while (t < T) {
    
    if (miss[t + 1L]) {
      
      s <- t + 1L
      e <- s
      
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      
      L <- e - t
      logAt <- max_plus_power(logA, L)
      
      for (j in 1:K) {
        tmp <- delta[t, ] + logAt[, j]
        psi[e, j] <- which.max(tmp)
        delta[e, j] <- max(tmp)
      }
      
      t <- e
      
    } else {
      
      tt <- t + 1L
      
      for (j in 1:K) {
        tmp <- delta[t, ] + logA[, j]
        psi[tt, j] <- which.max(tmp)
        delta[tt, j] <- max(tmp) + logB[tt, j]
      }
      
      t <- tt
    }
  }
  
  z <- integer(T)
  
  if (fill_gaps == "NA") {
    z[] <- NA_integer_
  }
  
  last <- max(which(!is.na(delta[, 1])))
  z[last] <- which.max(delta[last, ])
  
  tt <- last
  
  while (tt > 1L) {
    
    if (miss[tt]) {
      
      e <- tt
      
      while (e <= T && miss[e]) e <- e + 1L
      e <- e - 1L
      
      s <- tt
      
      while (s >= 1L && miss[s]) s <- s - 1L
      s <- s + 1L
      
      prev <- s - 1L
      
      if (prev < 1L) break
      
      z[prev] <- psi[e, z[e]]
      
      if (fill_gaps == "constant") {
        z[s:e] <- z[prev]
      }
      
      tt <- prev
      
    } else {
      
      z[tt - 1L] <- psi[tt, z[tt]]
      tt <- tt - 1L
    }
  }
  
  z
}

#' Observed-only naive Viterbi baseline for Skellam-HMM
#'
#' This baseline drops zero-exposure bins and treats consecutive observed bins as consecutive in time.
#' @export
viterbi_skellam_observed_only <- function(y, expo, pi1, A, mu1, mu2) {
  obs <- which(expo > 0 & !is.na(y))
  yo <- y[obs]
  eo <- expo[obs]
  K <- length(pi1); To <- length(yo)
  logA <- safe_log(A)
  logB <- matrix(0, To, K)
  for (t in 1:To) logB[t, ] <- emission_loglik_vec_skellam(yo[t], eo[t], mu1, mu2)
  delta <- matrix(-Inf, To, K)
  psi <- matrix(0L, To, K)
  delta[1, ] <- safe_log(pi1) + logB[1, ]
  for (t in 2:To) {
    for (j in 1:K) {
      tmp <- delta[t - 1, ] + logA[, j]
      psi[t, j] <- which.max(tmp)
      delta[t, j] <- max(tmp) + logB[t, j]
    }
  }
  zobs <- integer(To)
  zobs[To] <- which.max(delta[To, ])
  for (t in (To - 1):1) zobs[t] <- psi[t + 1, zobs[t + 1]]
  z <- rep(NA_integer_, length(y))
  z[obs] <- zobs
  z
}

# ===== VERBATIM SOURCE: R/hmm_original/em_skellam.R =====
# Skellam HMM estimation -------------------------------------------------------

#' Gap-aware EM estimation for an exposure-adjusted Skellam-HMM
#'
#' @param y Integer difference series, e.g. H - S.
#' @param expo Effective exposure per bin.
#' @param K Number of latent states.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Log-likelihood convergence tolerance.
#' @param mu1_init Optional initial first rates.
#' @param mu2_init Optional initial second rates.
#' @param A_init Optional initial transition matrix.
#' @param pi1_init Optional initial distribution.
#' @param verbose Print iteration log-likelihood.
#' @export
hmm_skellam_EM_gapaware <- function(y, expo, K,
                                    max_iter = 200, tol = 1e-8,
                                    mu1_init = NULL, mu2_init = NULL,
                                    A_init = NULL, pi1_init = NULL,
                                    verbose = TRUE) {
  stopifnot(length(y) == length(expo))
  T <- length(y)
  miss <- is_missing_skellam(y, expo)
  blocks <- find_missing_blocks(miss)
  pi1 <- if (is.null(pi1_init)) rep(1 / K, K) else pi1_init
  A <- if (is.null(A_init)) matrix(1 / K, K, K) else A_init

  if (is.null(mu1_init) || is.null(mu2_init)) {
    idx <- which(!miss)
    yy <- y[idx]
    if (length(yy) < K) stop("Not enough observed points to initialise.")
    km <- stats::kmeans(yy, centers = K, nstart = 10)
    mu1 <- mu2 <- numeric(K)
    for (k in 1:K) {
      yk <- yy[km$cluster == k]
      m <- mean(yk)
      v <- stats::var(yk)
      if (!is.finite(v)) v <- abs(m) + 1e-2
      v <- max(v, abs(m) + 1e-2)
      ebar <- mean(expo[idx][km$cluster == k])
      ebar <- max(ebar, 1e-6)
      m_rate <- m / ebar
      v_rate <- v / ebar
      v_rate <- max(v_rate, abs(m_rate) + 1e-3)
      mu1[k] <- (v_rate + m_rate) / 2
      mu2[k] <- (v_rate - m_rate) / 2
    }
    mu1 <- pmax(mu1, 1e-6)
    mu2 <- pmax(mu2, 1e-6)
  } else {
    mu1 <- mu1_init
    mu2 <- mu2_init
  }

  loglik_old <- -Inf
  gamma <- NULL
  loglik <- NA_real_
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    logB <- matrix(0, nrow = T, ncol = K)
    for (t in 1:T) if (!miss[t]) logB[t, ] <- emission_loglik_vec_skellam(y[t], expo[t], mu1, mu2)

    fw <- forward_hybrid(logB, pi1, A, miss)
    logalpha <- fw$logalpha
    logbeta <- backward_hybrid(logB, A, miss)
    loglik <- fw$loglik

    loggamma <- logalpha + logbeta
    gamma <- exp(loggamma - matrixStats::rowLogSumExps(loggamma))
    pi1 <- gamma[1, ] / sum(gamma[1, ])

    A_cnt <- matrix(0, K, K)
    logA <- safe_log(A)
    for (t in 1:(T - 1L)) {
      if (miss[t + 1L]) next
      if (any(is.na(logalpha[t, ])) || any(is.na(logbeta[t + 1L, ]))) next
      logxi <- outer(logalpha[t, ], logbeta[t + 1L, ] + logB[t + 1L, ], "+") + logA
      xi <- exp(logxi - logSumExp(logxi))
      A_cnt <- A_cnt + xi
    }
    for (blk in blocks) {
      s <- blk[1]; e <- blk[2]
      prev <- s - 1L
      if (prev < 1L) next
      L <- e - prev
      if (any(is.na(logalpha[prev, ])) || any(is.na(logbeta[e, ]))) next
      alpha_prev <- exp(logalpha[prev, ] - logSumExp(logalpha[prev, ]))
      beta_e <- matrix(exp(logbeta[e, ] - logSumExp(logbeta[e, ])), ncol = 1)
      A_cnt <- A_cnt + expected_transitions_in_block(alpha_prev, beta_e, A, L)
    }
    A <- pmax(A_cnt, 1e-12)
    A <- A / rowSums(A)

    obs <- which(!miss)
    
    for (k in 1:K) {
      w <- gamma[obs, k]
      
      if (sum(w, na.rm = TRUE) < 1e-10) {
        mu1[k] <- pmax(mu1[k], 1e-6)
        mu2[k] <- pmax(mu2[k], 1e-6)
        next
      }
      
      yk <- y[obs]
      ek <- expo[obs]
      
      fit <- stats::optim(
        par = log(pmax(c(mu1[k], mu2[k]), 1e-6)),
        fn = function(p) {
          m1 <- exp(p[1])
          m2 <- exp(p[2])
          
          ll <- mapply(
            log_skellam_offset,
            y = yk,
            expo = ek,
            MoreArgs = list(mu1 = m1, mu2 = m2)
          )
          
          if (any(!is.finite(ll))) return(1e200)
          
          val <- -sum(w * ll)
          if (!is.finite(val)) 1e200 else val
        },
        method = "L-BFGS-B",
        lower = log(1e-6),
        upper = log(1e6)
      )
      
      mu1[k] <- pmax(exp(fit$par[1]), 1e-6)
      mu2[k] <- pmax(exp(fit$par[2]), 1e-6)
    }

    if (verbose) cat("Iter:", iter, "LogLik:", loglik, "\n")
    if (iter > 1 && abs(loglik - loglik_old) < tol) {
      converged <- TRUE
      break
    }
    loglik_old <- loglik
  }
  list(pi1 = pi1, A = A, mu1 = mu1, mu2 = mu2, gamma = gamma, loglik = loglik, iter = iter, converged = converged)
}

# ===== VERBATIM SOURCE: R/hmm_original/em_pois2.R =====
# Bivariate Poisson HMM estimation --------------------------------------------

#' Gap-aware EM estimation for an exposure-adjusted bivariate Poisson HMM
#'
#' @param H First count series, usually hard-band counts.
#' @param S Second count series, usually soft-band counts.
#' @param expo Effective exposure per bin.
#' @param K Number of latent states.
#' @export
hmm_pois2_EM_gapaware <- function(H, S, expo, K,
                                  max_iter = 200, tol = 1e-8,
                                  lamH_init = NULL, lamS_init = NULL,
                                  A_init = NULL, pi1_init = NULL,
                                  verbose = TRUE) {
  stopifnot(length(H) == length(S), length(H) == length(expo))
  T <- length(H)
  miss <- is_missing_pois2(H, S, expo)
  blocks <- find_missing_blocks(miss)
  pi1 <- if (is.null(pi1_init)) rep(1 / K, K) else pi1_init
  A <- if (is.null(A_init)) matrix(1 / K, K, K) else A_init

  if (is.null(lamH_init) || is.null(lamS_init)) {
    idx <- which(!miss)
    if (length(idx) < K) stop("Not enough observed points to initialise.")
    proxy <- H[idx] - S[idx]
    km <- stats::kmeans(proxy, centers = K, nstart = 10)
    lamH <- lamS <- numeric(K)
    for (k in 1:K) {
      ii <- idx[km$cluster == k]
      denom <- max(sum(expo[ii]), 1e-12)
      lamH[k] <- sum(H[ii]) / denom
      lamS[k] <- sum(S[ii]) / denom
    }
    lamH <- pmax(lamH, 1e-12)
    lamS <- pmax(lamS, 1e-12)
  } else {
    lamH <- lamH_init
    lamS <- lamS_init
  }

  loglik_old <- -Inf
  gamma <- NULL
  loglik <- NA_real_
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    logB <- matrix(0, nrow = T, ncol = K)
    for (t in 1:T) if (!miss[t]) logB[t, ] <- emission_loglik_vec_pois2(H[t], S[t], expo[t], lamH, lamS)

    fw <- forward_hybrid(logB, pi1, A, miss)
    logalpha <- fw$logalpha
    logbeta <- backward_hybrid(logB, A, miss)
    loglik <- fw$loglik

    loggamma <- logalpha + logbeta
    gamma <- exp(loggamma - matrixStats::rowLogSumExps(loggamma))
    pi1 <- gamma[1, ] / sum(gamma[1, ])

    A_cnt <- matrix(0, K, K)
    logA <- safe_log(A)
    for (t in 1:(T - 1L)) {
      if (miss[t + 1L]) next
      if (any(is.na(logalpha[t, ])) || any(is.na(logbeta[t + 1L, ]))) next
      logxi <- outer(logalpha[t, ], logbeta[t + 1L, ] + logB[t + 1L, ], "+") + logA
      xi <- exp(logxi - logSumExp(logxi))
      A_cnt <- A_cnt + xi
    }
    for (blk in blocks) {
      s <- blk[1]; e <- blk[2]
      prev <- s - 1L
      if (prev < 1L) next
      L <- e - prev
      if (any(is.na(logalpha[prev, ])) || any(is.na(logbeta[e, ]))) next
      alpha_prev <- exp(logalpha[prev, ] - logSumExp(logalpha[prev, ]))
      beta_e <- matrix(exp(logbeta[e, ] - logSumExp(logbeta[e, ])), ncol = 1)
      A_cnt <- A_cnt + expected_transitions_in_block(alpha_prev, beta_e, A, L)
    }
    A <- pmax(A_cnt, 1e-12)
    A <- A / rowSums(A)

    obs <- which(!miss)
    for (k in 1:K) {
      w <- gamma[obs, k]
      denom <- max(sum(w * expo[obs]), 1e-12)
      lamH[k] <- max(sum(w * H[obs]) / denom, 1e-12)
      lamS[k] <- max(sum(w * S[obs]) / denom, 1e-12)
    }

    if (verbose) cat("Iter:", iter, "LogLik:", loglik, "\n")
    if (iter > 1 && abs(loglik - loglik_old) < tol) {
      converged <- TRUE
      break
    }
    loglik_old <- loglik
  }
  list(pi1 = pi1, A = A, lamH = lamH, lamS = lamS, gamma = gamma, loglik = loglik, iter = iter, converged = converged)
}

# ===== VERBATIM SOURCE: R/hmm_original/em_bivpois.R =====
# Correlated bivariate Poisson HMM estimation ----------------------------------

hmm_bivpois_EM_gapaware <- function(H, S, expo, K,
                                    max_iter = 200, tol = 1e-8,
                                    lamH_init = NULL,
                                    lamS_init = NULL,
                                    lamC_init = NULL,
                                    A_init = NULL,
                                    pi1_init = NULL,
                                    verbose = TRUE) {
  
  stopifnot(length(H) == length(S), length(H) == length(expo))
  
  T <- length(H)
  
  eps_rate <- 1e-8
  eps_A <- 1e-12
  
  miss <- is_missing_pois2(H, S, expo)
  blocks <- find_missing_blocks(miss)
  
  pi1 <- if (is.null(pi1_init)) rep(1 / K, K) else pi1_init
  A <- if (is.null(A_init)) matrix(1 / K, K, K) else A_init
  
  # ------------------------------------------------------------------
  # Initialisation
  # ------------------------------------------------------------------
  
  if (is.null(lamH_init) || is.null(lamS_init) || is.null(lamC_init)) {
    
    idx <- which(!miss)
    
    if (length(idx) < K) {
      stop("Not enough observed points to initialise.")
    }
    
    proxy <- H[idx] - S[idx]
    km <- stats::kmeans(proxy, centers = K, nstart = 10)
    
    lamH <- lamS <- lamC <- numeric(K)
    
    for (k in 1:K) {
      ii <- idx[km$cluster == k]
      denom <- max(sum(expo[ii]), eps_rate)
      
      margH <- sum(H[ii]) / denom
      margS <- sum(S[ii]) / denom
      
      # Conservative common-component initialisation.
      common0 <- 0.10 * min(margH, margS)
      common0 <- min(common0, 0.80 * min(margH, margS))
      common0 <- max(common0, eps_rate)
      
      lamC[k] <- common0
      lamH[k] <- max(margH - common0, eps_rate)
      lamS[k] <- max(margS - common0, eps_rate)
    }
    
  } else {
    
    lamH <- lamH_init
    lamS <- lamS_init
    lamC <- lamC_init
  }
  
  lamH <- pmax(lamH, eps_rate)
  lamS <- pmax(lamS, eps_rate)
  lamC <- pmax(lamC, eps_rate)
  
  loglik_old <- -Inf
  gamma <- NULL
  loglik <- NA_real_
  converged <- FALSE
  
  # ------------------------------------------------------------------
  # EM iterations
  # ------------------------------------------------------------------
  
  for (iter in 1:max_iter) {
    
    # E-step: emission matrix
    logB <- matrix(0, nrow = T, ncol = K)
    
    for (t in 1:T) {
      if (!miss[t]) {
        logB[t, ] <- emission_loglik_vec_bivpois(
          H_t = H[t],
          S_t = S[t],
          expo_t = expo[t],
          lamH = lamH,
          lamS = lamS,
          lamC = lamC
        )
      }
    }
    
    fw <- forward_hybrid(logB, pi1, A, miss)
    logalpha <- fw$logalpha
    logbeta <- backward_hybrid(logB, A, miss)
    loglik <- fw$loglik
    
    loggamma <- logalpha + logbeta
    gamma <- exp(loggamma - matrixStats::rowLogSumExps(loggamma))
    
    pi1 <- gamma[1, ] / sum(gamma[1, ])
    
    # ----------------------------------------------------------------
    # Transition update, identical to the other models
    # ----------------------------------------------------------------
    
    A_cnt <- matrix(0, K, K)
    logA <- safe_log(A)
    
    for (t in 1:(T - 1L)) {
      if (miss[t + 1L]) next
      if (any(is.na(logalpha[t, ])) || any(is.na(logbeta[t + 1L, ]))) next
      
      logxi <- outer(
        logalpha[t, ],
        logbeta[t + 1L, ] + logB[t + 1L, ],
        "+"
      ) + logA
      
      xi <- exp(logxi - logSumExp(logxi))
      A_cnt <- A_cnt + xi
    }
    
    for (blk in blocks) {
      s <- blk[1]
      e <- blk[2]
      prev <- s - 1L
      
      if (prev < 1L) next
      
      L <- e - prev
      
      if (any(is.na(logalpha[prev, ])) || any(is.na(logbeta[e, ]))) next
      
      alpha_prev <- exp(logalpha[prev, ] - logSumExp(logalpha[prev, ]))
      beta_e <- matrix(exp(logbeta[e, ] - logSumExp(logbeta[e, ])), ncol = 1)
      
      A_cnt <- A_cnt + expected_transitions_in_block(
        alpha_prev,
        beta_e,
        A,
        L
      )
    }
    
    A <- pmax(A_cnt, eps_A)
    A <- A / rowSums(A)
    
    # ----------------------------------------------------------------
    # Emission update using expected common counts
    # ----------------------------------------------------------------
    
    obs <- which(!miss)
    
    for (k in 1:K) {
      
      w <- gamma[obs, k]
      denom <- max(sum(w * expo[obs]), eps_rate)
      
      c_hat <- vapply(
        obs,
        function(t) expected_common_bivpois(
          H = H[t],
          S = S[t],
          lamH = lamH[k],
          lamS = lamS[k],
          lamC = lamC[k],
          expo = expo[t]
        ),
        numeric(1)
      )
      
      c_hat[!is.finite(c_hat)] <- 0
      
      lamC[k] <- sum(w * c_hat, na.rm = TRUE) / denom
      lamH[k] <- sum(w * (H[obs] - c_hat), na.rm = TRUE) / denom
      lamS[k] <- sum(w * (S[obs] - c_hat), na.rm = TRUE) / denom
      
      lamH[k] <- max(lamH[k], eps_rate)
      lamS[k] <- max(lamS[k], eps_rate)
      lamC[k] <- max(lamC[k], eps_rate)
    }
    
    if (verbose) {
      cat("Iter:", iter, "LogLik:", loglik, "\n")
    }
    
    if (iter > 1 && abs(loglik - loglik_old) < tol) {
      converged <- TRUE
      break
    }
    
    loglik_old <- loglik
  }
  
  list(
    pi1 = pi1,
    A = A,
    lamH = lamH,
    lamS = lamS,
    lamC = lamC,
    gamma = gamma,
    loglik = loglik,
    iter = iter,
    converged = converged
  )
}

# ===== VERBATIM SOURCE: R/legacy_original/simulate.R =====
# Simulation utilities ---------------------------------------------------------

#' Generate an exposure sequence with complete gaps and partial exposures
#' @export
generate_exposure <- function(T_length, num_gaps = 5, max_gap_len = 60, num_partials = 50) {
  expo <- rep(1.0, T_length)
  if (num_gaps > 0) {
    for (i in 1:num_gaps) {
      gap_len <- sample(1:max_gap_len, 1)
      if (gap_len >= T_length) gap_len <- T_length - 1L
      start_idx <- sample(1:(T_length - gap_len + 1), 1)
      expo[start_idx:(start_idx + gap_len - 1)] <- 0.0
    }
  }
  if (num_partials > 0) {
    valid_idx <- which(expo > 0)
    n <- min(length(valid_idx), num_partials)
    if (n > 0) {
      partial_idx <- sample(valid_idx, n)
      expo[partial_idx] <- stats::runif(n, min = 0.1, max = 0.9)
    }
  }
  expo
}

#' Simulate an exposure-adjusted XRB HMM with optional common-mode background
#' @export
simulate_xrb_hmm <- function(T_length, K, mu1, mu2, A, pi1, expo,
                             common_background = c("none", "sinusoidal"),
                             background_amplitude = 10,
                             background_period = 500) {
  common_background <- match.arg(common_background)
  Z <- integer(T_length)
  H <- rep(NA_integer_, T_length)
  S <- rep(NA_integer_, T_length)
  Z[1] <- sample.int(K, size = 1, prob = pi1)
  for (t in 2:T_length) Z[t] <- sample.int(K, size = 1, prob = A[Z[t - 1], ])
  for (t in 1:T_length) {
    if (expo[t] > 0) {
      bg <- 0
      if (common_background == "sinusoidal") {
        bg <- background_amplitude * sin(2 * pi * t / background_period) + background_amplitude
      }
      rate_H <- (mu1[Z[t]] + bg) * expo[t]
      rate_S <- (mu2[Z[t]] + bg) * expo[t]
      H[t] <- stats::rpois(1, lambda = rate_H)
      S[t] <- stats::rpois(1, lambda = rate_S)
    }
  }
  list(H = H, S = S, Omega = H - S, Z_true = Z, expo = expo)
}

simulate_xrb_hmm_blocks <- function(D,
                                    block_lengths,
                                    gap_lengths,
                                    K,
                                    muH,
                                    muS,
                                    A,
                                    pi1,
                                    partial_exposure_prob = 0,
                                    partial_range = c(0.3, 0.9),
                                    common_background = c("none", "sinusoidal"),
                                    background_amplitude = 0,
                                    background_period = 500,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  common_background <- match.arg(common_background)
  
  stopifnot(length(block_lengths) == D)
  stopifnot(length(gap_lengths) == D - 1)
  stopifnot(nrow(A) == K, ncol(A) == K)
  stopifnot(length(muH) == K, length(muS) == K)
  stopifnot(length(pi1) == K)
  
  # normalise, just in case R decides to tolerate nonsense silently
  A <- A / rowSums(A)
  pi1 <- pi1 / sum(pi1)
  
  T_full <- sum(block_lengths) + sum(gap_lengths)
  
  Z_full <- integer(T_full)
  H_full <- rep(NA_integer_, T_full)
  S_full <- rep(NA_integer_, T_full)
  Omega_full <- rep(NA_integer_, T_full)
  expo_full <- rep(0, T_full)
  observed <- rep(FALSE, T_full)
  
  # ------------------------------------------------------------
  # 1. Generate latent Markov chain over the full timeline
  #    The chain evolves also during unobserved gaps.
  # ------------------------------------------------------------
  
  Z_full[1] <- sample.int(K, size = 1, prob = pi1)
  
  for (t in 2:T_full) {
    previous_state <- Z_full[t - 1]
    Z_full[t] <- sample.int(K, size = 1, prob = A[previous_state, ])
  }
  
  # ------------------------------------------------------------
  # 2. Define observed blocks separated by complete gaps
  # ------------------------------------------------------------
  
  block_start <- integer(D)
  block_end <- integer(D)
  
  current_pos <- 1
  
  for (d in seq_len(D)) {
    
    idx <- current_pos:(current_pos + block_lengths[d] - 1)
    
    block_start[d] <- idx[1]
    block_end[d] <- idx[length(idx)]
    
    observed[idx] <- TRUE
    expo_full[idx] <- 1
    
    # optional partial exposure inside observed blocks
    if (partial_exposure_prob > 0) {
      partial_flag <- runif(length(idx)) < partial_exposure_prob
      
      if (any(partial_flag)) {
        expo_full[idx[partial_flag]] <- runif(
          sum(partial_flag),
          min = partial_range[1],
          max = partial_range[2]
        )
      }
    }
    
    current_pos <- idx[length(idx)] + 1
    
    # skip the unobserved gap after the block
    if (d < D) {
      current_pos <- current_pos + gap_lengths[d]
    }
  }
  
  # ------------------------------------------------------------
  # 3. Generate counts only in observed positions
  # ------------------------------------------------------------
  
  for (t in seq_len(T_full)) {
    
    if (observed[t] && expo_full[t] > 0) {
      
      bg <- 0
      
      if (common_background == "sinusoidal") {
        bg <- background_amplitude * sin(2 * pi * t / background_period) +
          background_amplitude
      }
      
      k <- Z_full[t]
      
      rate_H <- (muH[k] + bg) * expo_full[t]
      rate_S <- (muS[k] + bg) * expo_full[t]
      
      H_full[t] <- rpois(1, lambda = rate_H)
      S_full[t] <- rpois(1, lambda = rate_S)
      Omega_full[t] <- H_full[t] - S_full[t]
    }
  }
  
  # ------------------------------------------------------------
  # 4. Extract observed sequence
  # ------------------------------------------------------------
  
  obs_idx <- which(observed)
  
  H_obs <- H_full[obs_idx]
  S_obs <- S_full[obs_idx]
  Omega_obs <- Omega_full[obs_idx]
  expo_obs <- expo_full[obs_idx]
  Z_obs <- Z_full[obs_idx]
  
  # elapsed time between observed records
  # within a block: delta = 1
  # across a gap: delta = gap_length + 1
  delta_obs <- c(NA_integer_, diff(obs_idx))
  
  # mark starts of new observed blocks
  is_new_block <- c(TRUE, diff(obs_idx) > 1)
  
  list(
    H = H_obs,
    S = S_obs,
    Omega = Omega_obs,
    expo = expo_obs,
    Z_true_obs = Z_obs,
    time_obs = obs_idx,
    delta_obs = delta_obs,
    is_new_block = is_new_block,
    
    Z_true_full = Z_full,
    H_full = H_full,
    S_full = S_full,
    Omega_full = Omega_full,
    expo_full = expo_full,
    observed = observed,
    
    block_start = block_start,
    block_end = block_end,
    gap_lengths = gap_lengths,
    block_lengths = block_lengths,
    T_full = T_full,
    A_true = A,
    pi1_true = pi1,
    muH_true = muH,
    muS_true = muS
  )
}

# ===== VERBATIM SOURCE: R/transition_metrics.R =====
# Utility for simulation update requested by statistical referee.
# State label alignment is essential and must be done BEFORE calling this function.
# For nontrivial K, supply an independently justified label mapping.
transition_metric_rows <- function(A_true, A_est, epsilon=NULL) {
  stopifnot(is.matrix(A_true),is.matrix(A_est),
            identical(dim(A_true),dim(A_est)),nrow(A_true)==ncol(A_true))
  K <- nrow(A_true)
  if(any(!is.finite(A_true)) || any(!is.finite(A_est)) ||
     any(A_true<0) || any(A_est<0)) stop('Invalid transition probability matrix')
  if(any(abs(rowSums(A_true)-1)>1e-7) || any(abs(rowSums(A_est)-1)>1e-7))
     stop('Rows do not sum to one')
  if(!is.null(epsilon)) stop('No silent epsilon smoothing: KL infinity is meaningful')
  tv <- rowSums(abs(A_true-A_est))/2
  kl <- vapply(seq_len(K), function(i) {
    t <- A_true[i,]; e <- A_est[i,]
    if(any(t>0 & e==0)) return(Inf)
    pos <- t>0
    sum(t[pos]*log(t[pos]/e[pos]))
  }, numeric(1))
  data.frame(from_state=seq_len(K),TV=tv,KL=kl)
}
mean_transition_metrics <- function(A_true,A_est) {
  v <- transition_metric_rows(A_true,A_est)
  c(mean_TV=mean(v$TV),mean_KL=mean(v$KL))
}

aoas_check_A <- function(A) {
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A) ||
      nrow(A) < 1L || any(!is.finite(A)) || any(A < -1e-10) ||
      any(abs(rowSums(A) - 1) > 1e-6)) {
    stop("A must be a finite row-stochastic square transition matrix.")
  }
  invisible(TRUE)
}

aoas_A_power <- function(A, L) {
  aoas_check_A(A)
  if (length(L) != 1L || !is.finite(L) || L < 0 || L != floor(L)) {
    stop("L must be a non-negative integer number of grid steps.")
  }
  n <- as.numeric(L)
  out <- diag(nrow(A))
  base <- A
  while (n > 0) {
    if (n %% 2 == 1) out <- out %*% base
    n <- floor(n / 2)
    if (n > 0) base <- base %*% base
  }
  out <- pmax(out, 0)
  out / rowSums(out)
}

aoas_TV <- function(p, q) {
  0.5 * sum(abs(p - q))
}

aoas_KL <- function(p, q) {
  positive <- p > 0
  if (any(q[positive] <= 0)) return(Inf)
  sum(p[positive] * log(p[positive] / q[positive]))
}

aoas_JS <- function(p, q) {
  m <- (p + q) / 2
  0.5 * aoas_KL(p, m) + 0.5 * aoas_KL(q, m)
}

# Optimal permutation of COMPRESSED labels into GAP-AWARE label order,
# maximizing posterior-MAP contingency agreement. Labels are never used
# for the predictive NLL, which is invariant to permutations.
# Exact dynamic-programming assignment, O(K^2 2^K), no external packages.
aoas_match_labels <- function(z_gap, z_comp, K) {
  if (length(z_gap) != length(z_comp) || K > 13L || K < 1L) {
    stop("Invalid MAP sequences / K; matching supports K <= 13.")
  }
  ok <- is.finite(z_gap) & is.finite(z_comp) &
    z_gap %in% seq_len(K) & z_comp %in% seq_len(K)
  if (!any(ok)) stop("No corresponding observed labels.")
  W <- matrix(0, K, K)
  for (r in which(ok)) W[z_gap[r], z_comp[r]] <- W[z_gap[r], z_comp[r]] + 1
  N <- as.integer(2^K)
  dp <- rep(-Inf, N)
  previous <- rep(NA_integer_, N)
  last_col <- rep(NA_integer_, N)
  dp[1L] <- 0
  for (mask in 0:(N - 2L)) {
    old <- dp[mask + 1L]
    if (!is.finite(old)) next
    row_i <- sum(as.integer(intToBits(mask))[seq_len(K)]) + 1L
    for (j in seq_len(K)) {
      b <- bitwShiftL(1L, j - 1L)
      if (bitwAnd(mask, b) != 0L) next
      new_mask <- bitwOr(mask, b)
      proposal <- old + W[row_i, j]
      if (proposal > dp[new_mask + 1L]) {
        dp[new_mask + 1L] <- proposal
        previous[new_mask + 1L] <- mask
        last_col[new_mask + 1L] <- j
      }
    }
  }
  permutation <- integer(K)
  mask <- N - 1L
  for (i in K:1L) {
    permutation[i] <- last_col[mask + 1L]
    mask <- previous[mask + 1L]
  }
  list(permutation = permutation,
       optimal_MAP_agreement = dp[N] / sum(ok),
       raw_MAP_agreement = mean(z_gap[ok] == z_comp[ok]),
       contingency = W)
}

# Optional permutation-invariant comparison of MAP partitions.
aoas_adjusted_rand <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  tt <- table(x, y)
  c2 <- function(v) sum(v * (v - 1) / 2)
  N2 <- n * (n - 1) / 2
  same <- c2(tt)
  ai <- c2(rowSums(tt)); bj <- c2(colSums(tt))
  expected <- ai * bj / N2
  den <- (ai + bj) / 2 - expected
  if (abs(den) < 1e-12) return(if (same == ai && same == bj) 1 else NA_real_)
  (same - expected) / den
}

# Compare actual effective transitions across each missing-data block.
# The 'elapsed_steps' L between observed bins includes the endpoint step:
# e.g. observed grid positions 2 and 7 -> L=5, missing bins=4.
# A_comp is already aligned to gap-aware state labels, if necessary.
# No alpha / stationary / empirical-occupancy weights are introduced.
aoas_transition_distances <- function(A_gap, A_comp_aligned,
                                      observed_grid_indices, dt_seconds = 1) {
  aoas_check_A(A_gap); aoas_check_A(A_comp_aligned)
  if (nrow(A_gap) != nrow(A_comp_aligned)) stop("Different K.")
  obs <- as.integer(observed_grid_indices)
  if (length(obs) < 2L || any(!is.finite(obs)) ||
      any(diff(obs) <= 0L)) stop("Observed indices must increase.")
  delta <- diff(obs)
  r <- which(delta > 1L)
  empty_state <- data.frame(gap_id = integer(), pre_grid = integer(),
                            post_grid = integer(), missing_bins = integer(), elapsed_steps = integer(),
                            elapsed_seconds = numeric(), state = integer(), TV = numeric(),
                            KL_gap_to_comp = numeric(), JS = numeric())
  empty_gap <- data.frame(gap_id = integer(), pre_grid = integer(),
                          post_grid = integer(), missing_bins = integer(), elapsed_steps = integer(),
                          elapsed_seconds = numeric(), mean_TV_equal_states = numeric(),
                          max_TV = numeric(), mean_JS_equal_states = numeric())
  if (!length(r)) return(list(by_state = empty_state, by_gap = empty_gap))
  rows <- vector("list", length(r))
  for (ii in seq_along(r)) {
    g <- r[ii]; L <- delta[g]
    Pg <- aoas_A_power(A_gap, L)
    rows[[ii]] <- data.frame(gap_id = ii,
                             pre_grid = obs[g], post_grid = obs[g + 1L],
                             missing_bins = L - 1L, elapsed_steps = L,
                             elapsed_seconds = L * dt_seconds,
                             state = seq_len(nrow(A_gap)),
                             TV = vapply(seq_len(nrow(A_gap)), function(i)
                               aoas_TV(Pg[i, ], A_comp_aligned[i, ]), numeric(1)),
                             KL_gap_to_comp = vapply(seq_len(nrow(A_gap)), function(i)
                               aoas_KL(Pg[i, ], A_comp_aligned[i, ]), numeric(1)),
                             JS = vapply(seq_len(nrow(A_gap)), function(i)
                               aoas_JS(Pg[i, ], A_comp_aligned[i, ]), numeric(1)))
  }
  by_state <- do.call(rbind, rows)
  by_gap <- do.call(rbind, lapply(rows, function(x)
    data.frame(gap_id = x$gap_id[1], pre_grid = x$pre_grid[1],
               post_grid = x$post_grid[1], missing_bins = x$missing_bins[1],
               elapsed_steps = x$elapsed_steps[1],
               elapsed_seconds = x$elapsed_seconds[1],
               mean_TV_equal_states = mean(x$TV), max_TV = max(x$TV),
               mean_JS_equal_states = mean(x$JS))))
  rownames(by_state) <- NULL; rownames(by_gap) <- NULL
  list(by_state = by_state, by_gap = by_gap)
}

# Log emission probabilities for OBSERVED bins only. The R project
# already supplies emission_loglik_vec_skellam, _pois2 and _bivpois.
aoas_log_emission <- function(fit, model, H_obs, S_obs, expo_obs) {
  n <- length(H_obs)
  if (length(S_obs) != n || length(expo_obs) != n ||
      any(!is.finite(H_obs)) || any(!is.finite(S_obs)) ||
      any(!is.finite(expo_obs)) || any(expo_obs <= 0)) {
    stop("H/S/exposure must be complete on observed bins.")
  }
  K <- length(fit$pi1)
  out <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    out[i, ] <- switch(model,
                       skellam = emission_loglik_vec_skellam(H_obs[i] - S_obs[i],
                                                             expo_obs[i], fit$mu1, fit$mu2),
                       pois2 = emission_loglik_vec_pois2(H_obs[i], S_obs[i],
                                                         expo_obs[i], fit$lamH, fit$lamS),
                       bivpois = emission_loglik_vec_bivpois(H_obs[i], S_obs[i],
                                                             expo_obs[i], fit$lamH, fit$lamS, fit$lamC),
                       stop("Unknown emission family: ", model))
  }
  if (anyNA(out) || any(out == Inf)) stop("Invalid log emissions.")
  out
}

# True one-step-ahead FILTERING, no smoothing / Viterbi information.
# The current observation contributes to the NLL before the alpha update.
# For compressed=TRUE, each pair of retained observed bins uses A (L=1).
# Model parameters are held FIXED; this is in-sample by construction.
aoas_prequential_nll <- function(A, pi1, log_emission,
                                 observed_grid_indices, gap_aware = TRUE) {
  aoas_check_A(A)
  K <- nrow(A); n <- nrow(log_emission)
  obs <- as.integer(observed_grid_indices)
  if (ncol(log_emission) != K || n != length(obs) || n < 2L ||
      any(diff(obs) <= 0L) || length(pi1) != K ||
      any(!is.finite(pi1)) || any(pi1 < 0) || sum(pi1) <= 0) {
    stop("Dimension/index/initial probability mismatch.")
  }
  alpha <- pi1 / sum(pi1)
  nll <- numeric(n)
  elapsed <- c(NA_integer_, diff(obs))
  for (t in seq_len(n)) {
    if (t == 1L) {
      prior <- alpha
    } else {
      L <- if (gap_aware) elapsed[t] else 1L
      P <- if (L == 1L) A else aoas_A_power(A, L)
      prior <- drop(alpha %*% P)
      prior <- pmax(prior, 0); prior <- prior / sum(prior)
    }
    positive <- prior > 0
    log_terms <- rep(-Inf, K)
    log_terms[positive] <- log(prior[positive]) + log_emission[t, positive]
    m <- max(log_terms)
    if (!is.finite(m)) stop("Predictive density is zero at observed bin ", t)
    lp <- m + log(sum(exp(log_terms - m)))
    nll[t] <- -lp
    alpha <- exp(log_terms - lp)
    alpha <- alpha / sum(alpha)
  }
  r <- which(elapsed > 1L)
  details <- data.frame(gap_id = seq_along(r),
                        post_observed_index = r, pre_grid = obs[r - 1L],
                        post_grid = obs[r], missing_bins = elapsed[r] - 1L,
                        elapsed_steps = elapsed[r], NLL = nll[r])
  list(by_gap = details,
       mean_crossgap_NLL = if (length(r)) mean(nll[r]) else NA_real_,
       sum_crossgap_NLL = if (length(r)) sum(nll[r]) else NA_real_,
       n_gaps = length(r), NLL_all_bins = nll)
}

# Entire comparison. TV values use optimal MAP-label alignment for descriptive
# matrix comparisons; NLL never uses aligned labels or state smoothing.
# 'in-sample' means parameters learned from all data, NOT out-of-sample.
aoas_compare_gap_fits <- function(fit_gap, fit_comp, model,
                                  H_full, S_full, expo_full,
                                  observed_indices, dt_seconds = 1) {
  K <- length(fit_gap$pi1)
  if (length(fit_comp$pi1) != K) stop("Gap/compressed K must match.")
  obs <- as.integer(observed_indices)
  if (length(H_full) != length(S_full) ||
      length(S_full) != length(expo_full) ||
      any(obs < 1L | obs > length(H_full))) stop("Inconsistent full-grid arrays.")
  if (!is.matrix(fit_gap$gamma) || !is.matrix(fit_comp$gamma) ||
      nrow(fit_comp$gamma) != length(obs)) {
    stop("Expected gap full-grid gamma and compressed observed-bin gamma.")
  }
  z_gap <- max.col(fit_gap$gamma[obs, , drop = FALSE], ties.method = "first")
  z_comp <- max.col(fit_comp$gamma, ties.method = "first")
  match <- aoas_match_labels(z_gap, z_comp, K)
  A_comp_aligned <- fit_comp$A[match$permutation,
                               match$permutation, drop = FALSE]
  dist <- aoas_transition_distances(fit_gap$A, A_comp_aligned,
                                    obs, dt_seconds = dt_seconds)
  em_gap <- aoas_log_emission(fit_gap, model, H_full[obs],
                              S_full[obs], expo_full[obs])
  em_comp <- aoas_log_emission(fit_comp, model, H_full[obs],
                               S_full[obs], expo_full[obs])
  score_gap <- aoas_prequential_nll(fit_gap$A, fit_gap$pi1,
                                    em_gap, obs, gap_aware = TRUE)
  score_comp <- aoas_prequential_nll(fit_comp$A, fit_comp$pi1,
                                     em_comp, obs, gap_aware = FALSE)
  # Same fitted gap model, comparing propagation A^L against A:
  # isolates temporal propagation from refitting/emission differences.
  score_gap_no_power <- aoas_prequential_nll(fit_gap$A, fit_gap$pi1,
                                             em_gap, obs, gap_aware = FALSE)
  summary <- data.frame(model = model, K = K,
                        n_gaps = score_gap$n_gaps,
                        map_agreement_raw = match$raw_MAP_agreement,
                        map_agreement_aligned = match$optimal_MAP_agreement,
                        ARI = aoas_adjusted_rand(z_gap, z_comp),
                        mean_TV_equal_states = if (nrow(dist$by_gap)) mean(dist$by_gap$mean_TV_equal_states) else NA_real_,
                        max_TV = if (nrow(dist$by_gap)) max(dist$by_gap$max_TV) else NA_real_,
                        mean_NLL_gapaware = score_gap$mean_crossgap_NLL,
                        mean_NLL_compressed = score_comp$mean_crossgap_NLL,
                        mean_NLL_fixed_fit_no_power = score_gap_no_power$mean_crossgap_NLL,
                        delta_NLL_comp_minus_gap = score_comp$mean_crossgap_NLL - score_gap$mean_crossgap_NLL,
                        delta_NLL_no_power_minus_gap = score_gap_no_power$mean_crossgap_NLL - score_gap$mean_crossgap_NLL)
  gap_table <- dist$by_gap
  if (nrow(gap_table)) {
    gap_table$NLL_gapaware <- score_gap$by_gap$NLL
    gap_table$NLL_compressed <- score_comp$by_gap$NLL
    gap_table$NLL_fixed_fit_no_power <- score_gap_no_power$by_gap$NLL
    gap_table$delta_NLL_comp_minus_gap <- gap_table$NLL_compressed - gap_table$NLL_gapaware
    gap_table$delta_NLL_no_power_minus_gap <- gap_table$NLL_fixed_fit_no_power - gap_table$NLL_gapaware
  }
  list(summary = summary, by_gap = gap_table, by_state = dist$by_state,
       permutation_comp_to_gap = match$permutation,
       MAP_contingency = match$contingency)
}

