# ============================================================================
# AOAS R1 | SOURCE-SPECIFIC FITTING ENGINE (PROPOSED REFACTOR)
# IMPORTANT: all HMM estimation calls, random seeds and scientific
# settings below are taken unchanged from the supplied R1 scripts.
# Produces bins, complete fitted grids, loglik/BIC, convergence and selected K.
# No figures, SCORPEON comparisons, Bayesian Blocks or manuscript tables.
# This file is called by run_analysis.R from the repository root.
# ============================================================================

# =============================================================================
# AOAS / NICER: PSR J2214+0038
# =============================================================================

# SETTINGS =====================================================================
source_id <- "J2214"
source_label <- "PSR J2214+0038"
obsid <- "1070060101"
dataset <- "archive"
dt_grid <- c(10L)       # Fit each resolution; not every resolution is analysed.
analysis_dt <- 10L             # Reference dt for manuscript analysis; other dt values are fitted.
K_grid <- 1:10
models <- c("skellam", "pois2", "bivpois")
workers <- 8L
max_iter <- 1000L
tol <- 1e-3
seed <- 20261021L
reuse_fits <- FALSE            # TRUE reuses earlier HMMs already in RAM.
soft_band <- c(0.5, 2)
hard_band <- c(2, 10)
energy_tag <- "E05_2_E2_10"

# Required dependencies. Original estimators are NOT modified.
for (pkg in c("future", "future.apply", "extraDistr", "matrixStats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
}
if (!file.exists("hmm_functions.R")) stop("Missing hmm_functions.R in working directory")
source("hmm_functions.R")
project_dir <- getwd()
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# One source-specific container. Never overwrite other sources' results.
# Optional reuse requires compatible previous in-memory fits; disabled by default.
previous_container <- if (reuse_fits &&
  exists("container_J2214", envir = .GlobalEnv, inherits = FALSE)) {
  get("container_J2214", envir = .GlobalEnv)
} else NULL

container_J2214 <- list(source = source_label, obsid = obsid, dataset = dataset,
                        analysis_dt = analysis_dt, bins = list(), fits = list(),
                        fits_comp = list(), bic = data.frame(),
                        bic_comp = data.frame(), selected = data.frame(),
                        selected_comp = data.frame(),
                        tables = list(), figures = list(), checks = list())

# =============================================================================
# PART 1. COMPLETE ESTIMATION: EVERY dt, EVERY FAMILY, EVERY K
# =============================================================================

# 1.1 Read events and GTIs ONCE. These are the scientific data inputs.
events_file <- file.path("data", paste0(obsid, "_events.csv"))
gti_file <- file.path("data", paste0(obsid, "_gti.csv"))
if (!file.exists(events_file) || !file.exists(gti_file)) {
  stop("Missing event/GTI input: ", events_file, " | ", gti_file)
}
events <- read.csv(events_file)
gti <- read.csv(gti_file)
stopifnot(all(c("time", "pi") %in% names(events)),
          all(c("start", "stop") %in% names(gti)))
gti <- gti[order(gti$start), , drop = FALSE]
if (!nrow(gti) || any(gti$stop <= gti$start) ||
    (nrow(gti) > 1 && any(gti$start[-1] < head(gti$stop, -1) - 1e-7))) {
  stop("Invalid GTIs")
}
idx_gtis <- findInterval(events$time, gti$start)
inside <- idx_gtis >= 1L
inside[inside] <- events$time[inside] < gti$stop[idx_gtis[inside]]
if (!all(inside)) stop("Event times and GTI boundaries do not agree")

# 1.2 For each dt, build exactly the observed hard/soft counts and exposure.
for (dt in dt_grid) {
  dt_name <- paste0("dt", dt)
  cat("\n========== PART 1 |", source_id, "|", dt_name, "==========\n")

  t0 <- floor(min(gti$start))
  t1 <- ceiling(max(gti$stop))
  n_bins <- as.integer(ceiling((t1 - t0) / dt))
  if (n_bins <= max(K_grid) || n_bins > 5e6) stop("Unreasonable grid size")
  t_start <- t0 + (seq_len(n_bins) - 1L) * dt
  t_stop <- t_start + dt
  exposure <- numeric(n_bins)

  for (g in seq_len(nrow(gti))) {
    first <- max(1L, floor((gti$start[g] - t0) / dt) + 1L)
    last  <- min(n_bins, ceiling((gti$stop[g] - t0) / dt))
    if (first > last) next
    ii <- first:last
    exposure[ii] <- exposure[ii] + pmax(0, pmin(t_stop[ii], gti$stop[g]) -
                                        pmax(t_start[ii], gti$start[g]))
  }
  exposure <- pmin(exposure, dt)
  if (abs(sum(exposure) - sum(gti$stop - gti$start)) > 1e-5)
    stop("Exposure/GTI mismatch")
  expo <- exposure / dt                  # Offset used by original HMM functions.
  observed <- which(exposure > 0)

  eb <- floor((events$time - t0) / dt) + 1L
  keV <- events$pi * 0.01
  in_range <- eb >= 1L & eb <= n_bins
  H <- tabulate(as.integer(eb[in_range & keV >= hard_band[1] &
                                   keV <= hard_band[2]]), nbins = n_bins)
  S <- tabulate(as.integer(eb[in_range & keV >= soft_band[1] &
                                   keV < soft_band[2]]), nbins = n_bins)
  E12 <- tabulate(as.integer(eb[in_range & keV >= 12 & keV < 15]), nbins = n_bins)
  if (sum(H[expo == 0]) + sum(S[expo == 0]) != 0)
    stop("Counts outside exposure")
  H[expo == 0] <- NA_integer_
  S[expo == 0] <- NA_integer_
  E12[expo == 0] <- NA_integer_
  D <- H - S
  total <- H + S

  # TABLE DATA: one record per FULL-GRID bin. Preserve missing bins explicitly.
  bin_table <- data.frame(bin = seq_len(n_bins), t_start, t_stop,
    time_min = (t_start - t0) / 60, exposure_s = exposure, expo = expo,
    H = H, S = S, D = D, total = total, E12 = E12,
    H_rate = H / exposure, S_rate = S / exposure,
    D_rate = D / exposure, total_rate = total / exposure,
    HR = ifelse(total > 0, D / total, NA_real_))
  container_J2214$bins[[dt_name]] <- bin_table
  cat("GTIs:", nrow(gti), "| exposure (s):", sum(exposure),
      "| observed bins:", length(observed), "\n")

  # 1.3 FIT ALL (model, K) PAIRS IN PARALLEL. No figure/table generation here.
  cached <- NULL
  if (!is.null(previous_container)) {
    old_bins <- previous_container$bins[[dt_name]]
    if (!is.null(old_bins) && identical(old_bins$H, bin_table$H) &&
        identical(old_bins$S, bin_table$S) &&
        identical(old_bins$exposure_s, bin_table$exposure_s)) {
      cached <- previous_container$fits[[dt_name]]
    }
  }
  if (!is.null(cached) && all(vapply(models, function(m)
      all(paste0("K", K_grid) %in% names(cached[[m]])), logical(1)))) {
    all_fits <- cached
    cat("Reusing existing fits from previous RAM container.\n")
  } else {
    jobs <- expand.grid(model = models, K = K_grid, stringsAsFactors = FALSE)
    old_plan <- future::plan()
    future::plan(future::multisession, workers = workers)
    fits_parallel <- tryCatch(
      future.apply::future_lapply(seq_len(nrow(jobs)), function(j) {
        source(file.path(project_dir, "hmm_functions.R"), local = TRUE)
        model <- jobs$model[j]
        K <- jobs$K[j]
        set.seed(seed + dt * 1000L + K * 100L + match(model, models))
        error <- NA_character_
        fit <- tryCatch({
          if (model == "skellam") {
            hmm_skellam_EM_gapaware(D, expo, K, max_iter = max_iter,
                                    tol = tol, verbose = FALSE)
          } else if (model == "pois2") {
            hmm_pois2_EM_gapaware(H, S, expo, K, max_iter = max_iter,
                                  tol = tol, verbose = FALSE)
          } else {
            hmm_bivpois_EM_gapaware(H, S, expo, K, max_iter = max_iter,
                                    tol = tol, verbose = FALSE)
          }
        }, error = function(err) {
          error <<- conditionMessage(err)
          NULL
        })
        if (!is.null(fit) && K > 1L) {
          if (model == "skellam") fit <- sort_skellam_states_by_contrast(fit)
          if (model == "pois2") fit <- sort_pois2_states_by_contrast(fit)
          if (model == "bivpois") fit <- sort_bivpois_states_by_contrast(fit)
        }
        list(fit = fit, error = error)
      }, future.seed = TRUE),
      finally = future::plan(old_plan))

    all_fits <- list(skellam = list(), pois2 = list(), bivpois = list())
    for (j in seq_len(nrow(jobs))) {
      label <- jobs$model[j]
      kname <- paste0("K", jobs$K[j])
      all_fits[[label]][[kname]] <- fits_parallel[[j]]$fit
      if (!is.na(fits_parallel[[j]]$error)) {
        message(label, " ", kname, ": ", fits_parallel[[j]]$error)
      }
    }
  }
  container_J2214$fits[[dt_name]] <- all_fits

  # 1.4 BIC FOR ALL FITS. Only BIC is used for selection.
  bic_dt <- data.frame()
  for (model in models) {
    for (K in K_grid) {
      fit <- all_fits[[model]][[paste0("K", K)]]
      bic_dt <- rbind(bic_dt, data.frame(
        source = source_id, dataset = dataset, dt = dt, model = model, K = K,
        logLik = if (is.null(fit)) NA_real_ else fit$loglik,
        BIC = if (is.null(fit)) NA_real_ else
          bic_hmm(fit$loglik, K, length(observed), model),
        converged = if (is.null(fit)) FALSE else isTRUE(fit$converged),
        iterations = if (is.null(fit)) NA_integer_ else fit$iter))
    }
  }
  container_J2214$bic <- rbind(container_J2214$bic, bic_dt)
  # 1.5 Select the BIC optimum WITHIN each family and each dt.
  for (model in models) {
    x <- bic_dt[bic_dt$model == model & bic_dt$converged &
                  is.finite(bic_dt$BIC), ]
    if (!nrow(x)) {
      warning("No converged fit for ", model, " at dt=", dt)
      next
    }
    best <- x[which.min(x$BIC), ]
    container_J2214$selected <- rbind(container_J2214$selected,
      data.frame(source = source_id, dt = dt, model = model,
        K = best$K, BIC = best$BIC, boundary = best$K == max(K_grid)))
  }
}

# =============================================================================
# PART 1B. COMPRESSED MODEL: ALL dt, ALL FAMILIES, ALL K
# =============================================================================
# The compressed model is fitted on the SAME retained observations as the
# gap-aware model, in their observed order. It deliberately ignores elapsed
# time between GTIs (one Markov step per adjacent retained observation).
# No interpretation or figures are made here. ALL compressed fits remain in RAM.

for (dt in dt_grid) {
  dt_name <- paste0("dt", dt)
  counts_comp <- container_J2214$bins[[dt_name]]
  obs_comp <- which(counts_comp$exposure_s > 0)
  H_comp <- counts_comp$H[obs_comp]
  S_comp <- counts_comp$S[obs_comp]
  D_comp <- counts_comp$D[obs_comp]
  expo_comp <- counts_comp$expo[obs_comp]

  cat("\n========== PART 1B | COMPRESSED |", source_id,
      "| dt=", dt, "s ==========\n")

  # 1B.1. Reuse the COMPRESSED grid only if the underlying observations agree.
  cached_comp <- NULL
  if (!is.null(previous_container)) {
    old_bins <- previous_container$bins[[dt_name]]
    if (!is.null(old_bins) && identical(old_bins$H, counts_comp$H) &&
        identical(old_bins$S, counts_comp$S) &&
        identical(old_bins$exposure_s, counts_comp$exposure_s)) {
      cached_comp <- previous_container$fits_comp[[dt_name]]
    }
  }

  if (!is.null(cached_comp) && all(vapply(models, function(m)
      all(paste0("K", K_grid) %in% names(cached_comp[[m]])), logical(1)))) {
    compressed_fits <- cached_comp
    cat("Reusing compressed fits already available in RAM.\n")
  } else {
    # 1B.2. Fit 3 x 10 compressed HMMs in parallel, without any gap powers.
    jobs_comp <- expand.grid(model = models, K = K_grid,
                             stringsAsFactors = FALSE)
    old_plan <- future::plan()
    future::plan(future::multisession, workers = workers)
    comp_parallel <- tryCatch(
      future.apply::future_lapply(seq_len(nrow(jobs_comp)), function(j) {
        source(file.path(project_dir, "hmm_functions.R"), local = TRUE)
        model <- jobs_comp$model[j]
        K <- jobs_comp$K[j]
        set.seed(seed + 500000L + dt * 1000L + K * 100L +
                 match(model, models))
        error <- NA_character_
        fit <- tryCatch({
          if (model == "skellam") {
            hmm_skellam_EM_gapaware(D_comp, expo_comp, K,
                      max_iter = max_iter, tol = tol, verbose = FALSE)
          } else if (model == "pois2") {
            hmm_pois2_EM_gapaware(H_comp, S_comp, expo_comp, K,
                      max_iter = max_iter, tol = tol, verbose = FALSE)
          } else {
            hmm_bivpois_EM_gapaware(H_comp, S_comp, expo_comp, K,
                      max_iter = max_iter, tol = tol, verbose = FALSE)
          }
        }, error = function(err) {
          error <<- conditionMessage(err)
          NULL
        })
        if (!is.null(fit) && K > 1L) {
          if (model == "skellam") fit <- sort_skellam_states_by_contrast(fit)
          if (model == "pois2") fit <- sort_pois2_states_by_contrast(fit)
          if (model == "bivpois") fit <- sort_bivpois_states_by_contrast(fit)
        }
        list(fit = fit, error = error)
      }, future.seed = TRUE),
      finally = future::plan(old_plan))

    compressed_fits <- list(skellam = list(), pois2 = list(), bivpois = list())
    for (j in seq_len(nrow(jobs_comp))) {
      model <- jobs_comp$model[j]
      kname <- paste0("K", jobs_comp$K[j])
      compressed_fits[[model]][kname] <- list(comp_parallel[[j]]$fit)
      if (!is.na(comp_parallel[[j]]$error))
        message("COMPRESSED ", model, " ", kname, ": ", comp_parallel[[j]]$error)
    }
  }

  container_J2214$fits_comp[[dt_name]] <- compressed_fits

  # 1B.3. BIC for each compressed fit: same n and emission family as gap-aware.
  bic_comp_dt <- data.frame()
  for (model in models) {
    for (K in K_grid) {
      fit <- compressed_fits[[model]][[paste0("K", K)]]
      bic_comp_dt <- rbind(bic_comp_dt, data.frame(
        source = source_id, dataset = dataset, dt = dt, model = model, K = K,
        logLik = if (is.null(fit)) NA_real_ else fit$loglik,
        BIC = if (is.null(fit)) NA_real_ else
          bic_hmm(fit$loglik, K, length(obs_comp), model),
        converged = if (is.null(fit)) FALSE else isTRUE(fit$converged),
        iterations = if (is.null(fit)) NA_integer_ else fit$iter))
    }
  }
  container_J2214$bic_comp <- rbind(container_J2214$bic_comp, bic_comp_dt)
  # 1B.4. Independent BIC selection for the compressed HMM.
  for (model in models) {
    x <- bic_comp_dt[bic_comp_dt$model == model &
                     bic_comp_dt$converged & is.finite(bic_comp_dt$BIC), ]
    if (!nrow(x)) {
      warning("No converged compressed fit for ", model, " at dt=", dt)
      next
    }
    best <- x[which.min(x$BIC), ]
    container_J2214$selected_comp <- rbind(
      container_J2214$selected_comp,
      data.frame(source = source_id, dt = dt, model = model,
                 K = best$K, BIC = best$BIC,
                 boundary = best$K == max(K_grid)))
  }
}


# PART 2. RETAIN SCIENTIFIC SETTINGS IN THE MODEL CONTAINER ====================
container_J2214$settings <- list(
  source_id = source_id, source_label = source_label,
  soft_band = soft_band, hard_band = hard_band,
  energy_tag = energy_tag, dt_grid = dt_grid,
  analysis_dt = analysis_dt, K_grid = K_grid, models = models,
  workers = workers, max_iter = max_iter, tol = tol,
  seed = seed
)
cat("\nENGINE FINISHED |", source_id, "|", nrow(container_J2214$bic),
    "gap-aware and", nrow(container_J2214$bic_comp),
    "compressed candidate fits recorded.\n")
