
# ---- 1. Deterministic seeds -------------------------------------------------

simulation_seed <- function(...) {
  label <- paste(..., collapse = "|")
  values <- utf8ToInt(label)
  number <- 0
  for (value in values) {
    number <- (131 * number + value) %% 2147483646
  }
  as.integer(number + 1)
}


# ---- 2. True Markov chains and source emissions -----------------------------
# Values are exactly those used for the K=2 and K=3 configurations retained
# in the R1 main simulation study. Their numerical values are not altered.
# "both_vary": the state-specific intensity and normalized hardness may vary.
# "fixed_shape": the H/S ratio is fixed across states; intensity varies.

true_parameters <- function(K_true) {

  if (K_true == 2L) {

    A <- matrix(c(
      0.97, 0.03,
      0.04, 0.96
    ), nrow = 2L, byrow = TRUE)

    mu_H <- c(2, 7)
    mu_S <- c(8, 9)

    fixed_H <- c(4, 6)
    fixed_S <- c(16, 24)

  } else if (K_true == 3L) {

    A <- matrix(c(
      0.96, 0.04, 0,
      0.03, 0.94, 0.03,
      0, 0.04, 0.96
    ), nrow = 3L, byrow = TRUE)

    mu_H <- c(2, 7, 13)
    mu_S <- c(8, 9, 6)

    fixed_H <- c(2, 4, 6)
    fixed_S <- c(8, 16, 24)

  } else {
    stop("K_true must be 2 or 3")
  }

  list(
    A = A,
    mu_H = mu_H,
    mu_S = mu_S,
    fixed_H = fixed_H,
    fixed_S = fixed_S
  )
}


# ---- 3. Simulate latent states ---------------------------------------------

simulate_states <- function(T, A, seed) {
  set.seed(seed)
  K <- nrow(A)
  Z <- integer(T)
  Z[1] <- sample.int(K, 1L)
  
  for (t in 2:T) {
    Z[t] <- sample.int(K, 1L, prob = A[Z[t - 1L], ])
  }
  
  Z
}


# ---- 4. Simulate hard and soft counts --------------------------------------

simulate_counts <- function(Z, mu_H, mu_S, b_H, b_S,
                            source_seed, background_seed) {
  
  T <- length(Z)
  period <- T / 2
  
  set.seed(source_seed)
  
  H_source <- integer(T)
  S_source <- integer(T)
  
  for (t in seq_len(T)) {
    H_source[t] <- rpois(1L, mu_H[Z[t]])
    S_source[t] <- rpois(1L, mu_S[Z[t]])
  }
  
  set.seed(background_seed)
  
  B_H <- integer(T)
  B_S <- integer(T)
  
  for (t in seq_len(T)) {
    shape <- 1 + sin(2 * pi * t / period)
    B_H[t] <- rpois(1L, b_H * shape)
    B_S[t] <- rpois(1L, b_S * shape)
  }
  
  H <- H_source + B_H
  S <- S_source + B_S
  
  data.frame(
    time = seq_len(T),
    Z_true = Z,
    exposure = rep(1, T),
    H = H,
    S = S,
    D = H - S,
    H_source = H_source,
    S_source = S_source,
    B_H = B_H,
    B_S = B_S
  )
}


# ---- 5. Quasi-periodic exposure gaps ---------------------------------------

make_quasi_mask <- function(T, n_gaps, gap_length, seed) {
  
  set.seed(seed)
  exposure <- rep(1, T)
  
  for (g in seq_len(n_gaps)) {
    
    cycle_start <- floor((g - 1L) * T / n_gaps) + 1L
    cycle_end <- floor(g * T / n_gaps)
    
    center_start <- floor(
      (cycle_start + cycle_end - gap_length + 1L) / 2
    )
    
    free_room <- min(
      center_start - cycle_start,
      cycle_end - (center_start + gap_length - 1L)
    )
    
    if (free_room < 1L) {
      stop("Gap does not fit its quasi-periodic cycle")
    }
    
    jitter <- min(8L, free_room - 1L)
    shift <- sample((-jitter):jitter, size = 1L)
    
    start <- center_start + shift
    end <- start + gap_length - 1L
    
    exposure[start:end] <- 0
  }
  
  stopifnot(sum(exposure == 0) == n_gaps * gap_length)
  
  exposure
}


# ---- 6. Apply observing mask -----------------------------------------------

apply_mask <- function(data, exposure) {
  
  stopifnot(nrow(data) == length(exposure))
  
  data$exposure <- exposure
  
  missing <- which(exposure == 0)
  
  data$H[missing] <- NA_integer_
  data$S[missing] <- NA_integer_
  data$D[missing] <- NA_integer_
  
  data
}


# ---- 7. Fit original HMM models --------------------------------------------

fit_one_model <- function(data, model, method, K, max_iter, tol) {
  
  if (method == "gap-aware") {
    
    H <- data$H
    S <- data$S
    D <- data$D
    expo <- data$exposure
    
  } else if (method == "compressed") {
    
    observed <- which(data$exposure > 0)
    
    H <- data$H[observed]
    S <- data$S[observed]
    D <- data$D[observed]
    expo <- data$exposure[observed]
    
  } else {
    stop("Unknown method: ", method)
  }
  
  if (model == "skellam") {
    
    return(hmm_skellam_EM_gapaware(
      y = D,
      expo = expo,
      K = K,
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    ))
  }
  
  if (model == "pois2") {
    
    return(hmm_pois2_EM_gapaware(
      H = H,
      S = S,
      expo = expo,
      K = K,
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    ))
  }
  
  if (model == "bivpois") {
    
    return(hmm_bivpois_EM_gapaware(
      H = H,
      S = S,
      expo = expo,
      K = K,
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    ))
  }
  
  stop("Unknown emission model: ", model)
}


# ---- 8. Adjusted Rand Index ------------------------------------------------

adjusted_rand_index <- function(z_true, z_hat) {
  
  tab <- table(z_true, z_hat)
  n <- sum(tab)
  
  if (n < 2L) return(NA_real_)
  
  both_pairs <- sum(tab * (tab - 1) / 2)
  
  true_pairs <- sum(
    rowSums(tab) * (rowSums(tab) - 1) / 2
  )
  
  fitted_pairs <- sum(
    colSums(tab) * (colSums(tab) - 1) / 2
  )
  
  expected <- true_pairs * fitted_pairs / (n * (n - 1) / 2)
  maximum <- (true_pairs + fitted_pairs) / 2
  
  if (abs(maximum - expected) < 1e-12) {
    
    if (length(unique(z_true)) == 1L &&
        length(unique(z_hat)) == 1L) {
      return(1)
    }
    
    return(NA_real_)
  }
  
  (both_pairs - expected) / (maximum - expected)
}


# ---- 9. Stationary distribution --------------------------------------------

stationary_distribution <- function(A) {
  
  K <- nrow(A)
  
  equations <- rbind(
    t(A) - diag(K),
    rep(1, K)
  )
  
  rhs <- c(rep(0, K), 1)
  
  pi <- tryCatch(
    as.numeric(qr.solve(equations, rhs)),
    error = function(e) rep(NA_real_, K)
  )
  
  if (any(!is.finite(pi))) return(pi)
  
  pi <- pmax(pi, 0)
  
  if (sum(pi) == 0) {
    return(rep(NA_real_, K))
  }
  
  pi / sum(pi)
}


# ---- 10. Model recovery metrics --------------------------------------------

evaluate_fit <- function(data, fit, method, A_true) {
  
  scores <- c(
    ARI = NA_real_,
    decoded_accuracy = NA_real_,
    TV = NA_real_,
    KL = NA_real_,
    persistence_MAE = NA_real_,
    stationary_TV = NA_real_
  )
  
  observed <- which(data$exposure > 0)
  
  if (method == "gap-aware") {
    gamma <- fit$gamma[observed, , drop = FALSE]
  } else {
    gamma <- fit$gamma
  }
  
  if (!is.matrix(gamma) ||
      nrow(gamma) != length(observed)) {
    stop("Posterior state probabilities have unexpected dimensions")
  }
  
  z_true <- data$Z_true[observed]
  z_hat <- max.col(gamma, ties.method = "first")
  
  scores["ARI"] <- adjusted_rand_index(z_true, z_hat)
  
  K_true <- nrow(A_true)
  
  if (ncol(gamma) != K_true) {
    return(scores)
  }
  
  confusion <- table(
    factor(z_true, levels = seq_len(K_true)),
    factor(z_hat, levels = seq_len(K_true))
  )
  
  permutation <- as.integer(
    clue::solve_LSAP(confusion, maximum = TRUE)
  )
  
  z_aligned <- match(z_hat, permutation)
  
  scores["decoded_accuracy"] <- mean(z_true == z_aligned)
  
  A_hat <- fit$A[permutation, permutation, drop = FALSE]
  
  tv_by_state <- numeric(K_true)
  kl_by_state <- numeric(K_true)
  
  for (k in seq_len(K_true)) {
    
    p <- A_true[k, ]
    q <- pmax(A_hat[k, ], 1e-12)
    q <- q / sum(q)
    
    tv_by_state[k] <- sum(abs(p - q)) / 2
    
    positive <- which(p > 0)
    
    kl_by_state[k] <- sum(
      p[positive] * log(p[positive] / q[positive])
    )
  }
  
  pi_true <- stationary_distribution(A_true)
  pi_hat <- stationary_distribution(A_hat)
  
  scores["TV"] <- mean(tv_by_state)
  scores["KL"] <- mean(kl_by_state)
  
  scores["persistence_MAE"] <- mean(
    abs(diag(A_true) - diag(A_hat))
  )
  
  if (all(is.finite(c(pi_true, pi_hat)))) {
    scores["stationary_TV"] <- sum(abs(pi_true - pi_hat)) / 2
  }
  
  scores
}


# ---- 11. Fit models and candidate K on one dataset -------------------------

fit_dataset <- function(
    data,
    study,
    rep_id,
    T,
    K_true,
    dgp,
    background,
    b_H,
    b_S,
    gap_design,
    n_gaps,
    gap_length,
    missing_target,
    A_true,
    methods,
    cfg
) {
  
  output <- list()
  output_index <- 0L
  
  if (cfg$K_RULE == "K_true") {
    K_values <- K_true
  } else {
    K_values <- seq_len(K_true + 2L)
  }
  
  n_observed <- sum(data$exposure > 0)
  
  for (method in methods) {
    
    for (model in cfg$MODELS) {
      
      fits_for_model <- vector("list", length(K_values))
      
      for (K in K_values) {
        
        # Original deterministic fitting seed.
        seed_fit <- simulation_seed(
          "fit",
          rep_id,
          T,
          K_true,
          dgp,
          background,
          gap_design,
          missing_target,
          model,
          K
        )
        
        set.seed(seed_fit)
        
        fit <- tryCatch(
          fit_one_model(
            data,
            model,
            method,
            K,
            cfg$MAX_ITER,
            cfg$TOL
          ),
          error = function(e) e
        )
        
        error_message <- ""
        
        if (inherits(fit, "error")) {
          error_message <- conditionMessage(fit)
          fit <- NULL
        }
        
        if (!is.null(fit)) {
          
          A_ok <- is.matrix(fit$A) &&
            identical(dim(fit$A), c(K, K)) &&
            all(is.finite(fit$A)) &&
            all(fit$A >= -1e-8) &&
            all(abs(rowSums(fit$A) - 1) < 1e-5)
          
          LL_ok <- length(fit$loglik) == 1L &&
            is.finite(fit$loglik)
          
          if (!A_ok || !LL_ok) {
            error_message <- "Invalid log-likelihood or transition matrix"
            fit <- NULL
          }
        }
        
        ll <- if (is.null(fit)) {
          NA_real_
        } else {
          as.numeric(fit$loglik)
        }
        
        bic <- if (is.finite(ll)) {
          bic_hmm(ll, K, n_observed, model)
        } else {
          NA_real_
        }
        
        converged <- if (is.null(fit)) {
          FALSE
        } else {
          isTRUE(fit$converged)
        }
        
        iterations <- if (
          is.null(fit) || is.null(fit$iter)
        ) {
          NA_integer_
        } else {
          as.integer(fit$iter)
        }
        
        scores <- c(
          ARI = NA_real_,
          decoded_accuracy = NA_real_,
          TV = NA_real_,
          KL = NA_real_,
          persistence_MAE = NA_real_,
          stationary_TV = NA_real_
        )
        
        A_text <- NA_character_
        
        if (!is.null(fit)) {
          
          A_text <- paste(
            signif(as.vector(t(fit$A)), 9),
            collapse = ";"
          )
          
          evaluation <- tryCatch(
            evaluate_fit(data, fit, method, A_true),
            error = function(e) e
          )
          
          if (inherits(evaluation, "error")) {
            
            error_message <- paste(
              error_message,
              "Metric error:",
              conditionMessage(evaluation)
            )
            
          } else {
            scores <- evaluation
          }
        }
        
        fits_for_model[[match(K, K_values)]] <- data.frame(
          study = study,
          replicate = rep_id,
          T = T,
          K_true = K_true,
          dgp = dgp,
          background = background,
          b_H = b_H,
          b_S = b_S,
          gap_pattern = "quasi-periodic",
          gap_design = gap_design,
          n_gaps = n_gaps,
          gap_length = gap_length,
          missing_target = missing_target,
          missing_fraction = mean(data$exposure == 0),
          n_observed = n_observed,
          method = method,
          model = model,
          K = K,
          loglik = ll,
          BIC = bic,
          converged = converged,
          iterations = iterations,
          n_errors = as.integer(nchar(error_message) > 0L),
          error_message = error_message,
          ARI = unname(scores["ARI"]),
          decoded_accuracy = unname(scores["decoded_accuracy"]),
          TV = unname(scores["TV"]),
          KL = unname(scores["KL"]),
          persistence_MAE = unname(scores["persistence_MAE"]),
          stationary_TV = unname(scores["stationary_TV"]),
          A_hat = A_text,
          stringsAsFactors = FALSE
        )
      }
      
      group <- do.call(rbind, fits_for_model)
      
      if (cfg$K_RULE == "K_true") {
        
        # SIM1 assumes K known: K selection is not affected much by gaps.
        
        group$K_selected <- NA_integer_
        group$bic_grid_complete <- NA
        group$bic_all_converged <- NA
        group$bic_at_upper_bound <- NA
        group$bic_provisional <- NA
        
      } else {
        
        # SIM2 search with full BIC grid: background affects K selection.
        
        grid_complete <- all(is.finite(group$BIC))
        all_converged <- all(group$converged)
        
        selected_K <- if (grid_complete) {
          group$K[which.min(group$BIC)]
        } else {
          NA_integer_
        }
        
        group$K_selected <- selected_K
        group$bic_grid_complete <- grid_complete
        group$bic_all_converged <- all_converged
        
        group$bic_at_upper_bound <- if (grid_complete) {
          selected_K == max(K_values)
        } else {
          NA
        }
        
        group$bic_provisional <- grid_complete && !all_converged
      }
      
      output_index <- output_index + 1L
      output[[output_index]] <- group
    }
  }
  
  do.call(rbind, output)
}


# ---- 12. SIM1: one complete Monte Carlo replication -------------------------

run_sim1_rep <- function(rep_id, cfg, result_dir) {

  filename <- file.path(result_dir, sprintf("rep_%03d.rds", rep_id))
  if (file.exists(filename)) stop("File already exists: ", filename)

  stopifnot(
    length(cfg$DGPS) >= 1L,
    !anyNA(cfg$DGPS),
    all(cfg$DGPS %in% c("both_vary", "fixed_shape")),
    !anyDuplicated(cfg$DGPS)
  )

  results <- list()
  result_index <- 0L

  for (T in cfg$T_GRID) {
    for (K_true in cfg$K_TRUE_GRID) {

      pars <- true_parameters(K_true)

      # Preserve the exact existing latent-state seed signature.
      Z_seed <- simulation_seed("states", rep_id, T, K_true)
      Z <- simulate_states(T, pars$A, Z_seed)

      for (dgp in cfg$DGPS) {

        if (dgp == "both_vary") {
          mu_H <- pars$mu_H
          mu_S <- pars$mu_S
        } else if (dgp == "fixed_shape") {
          mu_H <- pars$fixed_H
          mu_S <- pars$fixed_S
        }

        # Preserve source and background seed signatures from R1.
        source_seed <- simulation_seed(
          "source", rep_id, T, K_true, dgp
        )
        background_seed <- simulation_seed(
          "background", rep_id, T, K_true, dgp
        )

        full_data <- simulate_counts(
          Z = Z,
          mu_H = mu_H,
          mu_S = mu_S,
          b_H = 0,
          b_S = 0,
          source_seed = source_seed,
          background_seed = background_seed
        )

        for (missing_target in cfg$MISSING_GRID) {
          for (gap_design in cfg$GAP_DESIGNS) {

            n_gaps <- if (gap_design == "many_short") 12L else 3L
            gap_length <- as.integer(round(T * missing_target / n_gaps))

            # Same seed as the existing R1 mask generator.
            mask_seed <- simulation_seed(
              "mask", rep_id, T, K_true, n_gaps
            )

            exposure <- make_quasi_mask(
              T, n_gaps, gap_length, mask_seed
            )
            data <- apply_mask(full_data, exposure)

            result_index <- result_index + 1L
            results[[result_index]] <- fit_dataset(
              data = data,
              study = "SIM1",
              rep_id = rep_id,
              T = T,
              K_true = K_true,
              dgp = dgp,
              background = "none",
              b_H = 0,
              b_S = 0,
              gap_design = gap_design,
              n_gaps = n_gaps,
              gap_length = gap_length,
              missing_target = missing_target,
              A_true = pars$A,
              methods = c("gap-aware", "compressed"),
              cfg = cfg
            )
          }
        }
      }
    }
  }

  # Save one complete replication, never a partly written RDS.
  temp_file <- paste0(filename, ".tmp")
  saveRDS(do.call(rbind, results), temp_file)
  if (!file.rename(temp_file, filename)) {
    stop("Cannot save: ", filename)
  }

  cat("SIM1 finished replicate", rep_id, "\n")
  filename
}


# ---- 13. SIM2: one complete Monte Carlo replication -------------------------

run_sim2_rep <- function(rep_id, cfg, result_dir) {
  
  filename <- file.path(
    result_dir,
    sprintf("rep_%03d.rds", rep_id)
  )
  
  if (file.exists(filename)) {
    stop("File already exists: ", filename)
  }
  
  results <- list()
  result_index <- 0L
  
  for (T in cfg$T_GRID) {
    
    for (K_true in cfg$K_TRUE_GRID) {
      
      pars <- true_parameters(K_true)
      
      Z <- simulate_states(
        T,
        pars$A,
        simulation_seed("states", rep_id, T, K_true)
      )
      
      n_gaps <- 12L
      missing_target <- 0.12
      
      gap_length <- as.integer(
        round(T * missing_target / n_gaps)
      )
      
      mask_seed <- simulation_seed(
        "mask", rep_id, T, K_true, n_gaps
      )
      
      exposure <- make_quasi_mask(
        T, n_gaps, gap_length, mask_seed
      )
      
      for (dgp in cfg$DGPS) {
        
        if (dgp == "both_vary") {
          mu_H <- pars$mu_H
          mu_S <- pars$mu_S
          
        } else if (dgp == "fixed_shape") {
          mu_H <- pars$fixed_H
          mu_S <- pars$fixed_S
          
        } else {
          stop("Unknown source DGP: ", dgp)
        }
        
        source_seed <- simulation_seed(
          "source", rep_id, T, K_true, dgp
        )
        
        background_seed <- simulation_seed(
          "background", rep_id, T, K_true, dgp
        )
        
        for (b in seq_len(nrow(cfg$BACKGROUND))) {
          
          name <- cfg$BACKGROUND$name[b]
          b_H <- cfg$BACKGROUND$b_H[b]
          b_S <- cfg$BACKGROUND$b_S[b]
          
          full_data <- simulate_counts(
            Z, mu_H, mu_S,
            b_H, b_S,
            source_seed, background_seed
          )
          
          data <- apply_mask(full_data, exposure)
          
          result_index <- result_index + 1L
          
          results[[result_index]] <- fit_dataset(
            data = data,
            study = "SIM2",
            rep_id = rep_id,
            T = T,
            K_true = K_true,
            dgp = dgp,
            background = name,
            b_H = b_H,
            b_S = b_S,
            gap_design = "many_short",
            n_gaps = n_gaps,
            gap_length = gap_length,
            missing_target = missing_target,
            A_true = pars$A,
            methods = "gap-aware",
            cfg = cfg
          )
        }
      }
    }
  }
  
  temp_file <- paste0(filename, ".tmp")
  
  saveRDS(do.call(rbind, results), temp_file)
  
  if (!file.rename(temp_file, filename)) {
    stop("Cannot save: ", filename)
  }
  
  cat("SIM2 finished replicate", rep_id, "\n")
  
  filename
}

# ---- 14. Replication runner -----------------------------------------

run_simulation_study <- function(
    reps, cfg, result_dir, study, runner,
    expected_rows, n_workers, hmm_file, internal_file
) {

  stopifnot(
    study %in% c("SIM1", "SIM2"),
    length(reps) > 0L,
    all(is.finite(reps)),
    all(reps == as.integer(reps)),
    all(reps >= 1L),
    !anyDuplicated(reps),
    expected_rows > 0L,
    file.exists(hmm_file),
    file.exists(internal_file)
  )

  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  settings_file <- file.path(result_dir, "settings.rds")
  if (file.exists(settings_file)) {
    if (!identical(readRDS(settings_file), cfg)) {
      stop("Simulation settings changed. Use a different result directory.")
    }
  } else {
    saveRDS(cfg, settings_file)
  }

  # Keep the exact code version fixed when resuming a partially completed run.
  source_hash_file <- file.path(result_dir, "source_md5.rds")
  source_hash <- tools::md5sum(c(hmm_file, internal_file))
  if (file.exists(source_hash_file)) {
    if (!identical(readRDS(source_hash_file), source_hash)) {
      stop("Scientific source files changed. Verify before resuming.")
    }
  } else {
    saveRDS(source_hash, source_hash_file)
  }

  paths <- file.path(result_dir, sprintf("rep_%03d.rds", reps))
  done <- file.exists(paths)

  # A saved fit may contain NA diagnostics or failed K candidates, which are
  # retained by design. Check structural completeness, not convergence.
  for (i in which(done)) {
    path <- paths[i]
    part <- tryCatch(readRDS(path), error = function(e) NULL)

    if (!is.data.frame(part) || nrow(part) != expected_rows) {
      stop("Incomplete or invalid replication checkpoint: ", path)
    }
    if (!all(c("study", "replicate", "T", "K_true", "dgp",
               "model", "method", "K") %in% names(part))) {
      stop("Missing required columns in checkpoint: ", path)
    }
    if (!identical(unique(part$study), study) ||
        !identical(as.integer(unique(part$replicate)), as.integer(reps[i])) ||
        !setequal(unique(part$T), cfg$T_GRID) ||
        !setequal(unique(part$K_true), cfg$K_TRUE_GRID) ||
        !setequal(unique(part$dgp), cfg$DGPS) ||
        !setequal(unique(part$model), cfg$MODELS)) {
      stop("Replication checkpoint does not match the settings: ", path)
    }

    if (study == "SIM1") {
      keys <- c("replicate", "T", "K_true", "dgp", "method",
                "model", "gap_design", "missing_target", "K")
      if (!setequal(unique(part$method), c("gap-aware", "compressed")) ||
          !setequal(unique(part$gap_design), cfg$GAP_DESIGNS) ||
          !setequal(unique(part$missing_target), cfg$MISSING_GRID) ||
          any(part$K != part$K_true)) {
        stop("Invalid SIM1 checkpoint: ", path)
      }
    } else {
      keys <- c("replicate", "T", "K_true", "dgp", "background",
                "method", "model", "K")
      if (!setequal(unique(part$method), "gap-aware") ||
          !setequal(unique(part$background), cfg$BACKGROUND$name) ||
          any(part$K < 1L | part$K > part$K_true + 2L)) {
        stop("Invalid SIM2 checkpoint: ", path)
      }
    }

    if (anyDuplicated(part[keys])) {
      stop("Repeated experimental conditions in checkpoint: ", path)
    }
  }

  pending <- reps[!done]
  if (is.na(n_workers) || !is.finite(n_workers)) n_workers <- 1L
  workers <- max(1L, min(as.integer(n_workers), length(pending)))

  cat(study, "| complete:", sum(done), "| pending:", length(pending),
      "| workers:", if (length(pending)) workers else 0L, "\n")

  if (length(pending) == 0L) {
    cat(study, "| all requested replications already saved.\n")
    return(invisible(paths))
  }

  if (workers == 1L) {
    completed <- lapply(pending, function(rep_id) {
      runner(rep_id, cfg, result_dir)
    })
  } else {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      c("hmm_file", "internal_file"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      source(hmm_file)
      source(internal_file)
      NULL
    })
    completed <- parallel::parLapplyLB(
      cl, pending, runner, cfg = cfg, result_dir = result_dir
    )
  }

  cat(study, "| saved", length(completed), "new replications under",
      result_dir, "\n")
  invisible(paths)
}
