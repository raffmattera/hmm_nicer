# ============================================================================
# AOAS R1 | Run SIM1, run SIM2, then analyze their saved replication checkpoints.
# Run from this directory: source("simulations_analysis.R")
# SIM1 and SIM2 launcher scripts SKIP completed rep_###.rds files.
# The HMM fitting and marginal MAP ARI calculations are unchanged.
# ============================================================================

# ----------------------------- SETTINGS -------------------------------------
STUDIES <- c("SIM1", "SIM2")
RESULT_ROOT <- "results"
ANALYSIS_DIR <- file.path(RESULT_ROOT)
WRITE_FIGURES <- TRUE

if (!all(STUDIES %in% c("SIM1", "SIM2"))) stop("Invalid study in STUDIES")
main_dir <- file.path(ANALYSIS_DIR, "main_paper")
supp_dir <- file.path(ANALYSIS_DIR, "web_supplement")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

read_replicates <- function(study) {
  study_dir <- file.path(RESULT_ROOT, study)
  files <- sort(list.files(study_dir, pattern = "^rep_[0-9]+[.]rds$", full.names = TRUE))
  if (!length(files)) {
    message(study, ": no finished replication RDS in ", study_dir)
    return(NULL)
  }
  parts <- lapply(files, function(path) {
    x <- tryCatch(readRDS(path), error = function(e)
      stop("Cannot read checkpoint ", path, ": ", conditionMessage(e)))
    if (!is.data.frame(x) || !nrow(x)) stop("Empty/invalid checkpoint: ", path)
    necessary <- c("study", "replicate", "T", "K_true", "dgp", "method", "model",
                   "K", "loglik", "BIC", "converged", "ARI", "TV", "KL", "decoded_accuracy")
    if (!all(necessary %in% names(x))) stop("Missing columns in checkpoint: ", path)
    if (length(unique(x$study)) != 1L || unique(x$study) != study)
      stop("Wrong study in ", path)
    expected_rep <- as.integer(sub("^rep_([0-9]+)[.]rds$", "\\1", basename(path)))
    if (length(unique(x$replicate)) != 1L || unique(x$replicate) != expected_rep)
      stop("Wrong replication ID in ", path)
    x
  })
  result <- do.call(rbind, parts)
  rownames(result) <- NULL
  key <- if (study == "SIM1")
    c("replicate", "T", "K_true", "dgp", "gap_design", "missing_target",
      "method", "model", "K")
  else
    c("replicate", "T", "K_true", "dgp", "background", "method", "model", "K")
  if (!all(key %in% names(result))) stop(study, ": checkpoint missing grouping columns")
  if (anyDuplicated(result[key])) stop(study, ": duplicate fit keys in checkpoints")

  settings_path <- file.path(study_dir, "settings.rds")
  if (file.exists(settings_path)) {
    cfg <- readRDS(settings_path)
    expected <- if (study == "SIM1")
      length(cfg$T_GRID) * length(cfg$K_TRUE_GRID) * length(cfg$DGPS) *
      length(cfg$MISSING_GRID) * length(cfg$GAP_DESIGNS) * 2L * length(cfg$MODELS)
    else
      length(cfg$T_GRID) * length(cfg$DGPS) * nrow(cfg$BACKGROUND) *
      length(cfg$MODELS) * sum(cfg$K_TRUE_GRID + 2L)
    counts <- table(result$replicate)
    if (any(counts != expected)) stop(study, ": incomplete checkpoint; expected ",
                                     expected, " rows per replication, observed ",
                                     paste(names(counts), counts, sep = ":", collapse = "; "))
    cat(study, ": verified ", length(files), " checkpoint(s); ",
        expected, " rows per replica.\n", sep = "")
  } else {
    warning(study, ": settings.rds unavailable; only schema and uniqueness checked")
  }
  hashes_path <- file.path(study_dir, "source_md5.rds")
  if (file.exists(hashes_path)) {
    previous <- unname(readRDS(hashes_path))
    current_paths <- c("hmm_functions.R", "internal_simulations.R")
    if (all(file.exists(current_paths))) {
      current <- unname(tools::md5sum(current_paths))
      if (!identical(as.character(previous), as.character(current))) {
        warning(study, ": checkpoint scientific-source hashes differ from current files. ",
                "Do not combine results from incompatible versions.")
      }
    }
  }
  result

}

mc_summary <- function(x, by, metrics, flag = NULL) {
  if (!all(c(by, metrics) %in% names(x))) stop("Missing columns in summary input")
  if (is.null(flag)) flag <- rep(TRUE, nrow(x))
  if (length(flag) != nrow(x)) stop("Wrong analysis filter length")
  out <- lapply(split(seq_len(nrow(x)), interaction(x[by], drop = TRUE, lex.order = TRUE)), function(idx) {
    d <- x[idx, , drop = FALSE]
    ok <- (!is.na(flag[idx]) & flag[idx])
    row <- d[1L, by, drop = FALSE]
    row$n_total <- nrow(d)
    row$n_eligible <- sum(ok)
    for (m in metrics) {
      v <- d[[m]][ok]
      v <- v[is.finite(v)]
      n <- length(v)
      mu <- if (n) mean(v) else NA_real_
      ss <- if (n >= 2L) stats::sd(v) else NA_real_
      half <- if (n >= 2L) stats::qt(0.975, n - 1L) * ss / sqrt(n) else NA_real_
      row[[paste0(m, "_n")]] <- n
      row[[paste0(m, "_mean")]] <- mu
      row[[paste0(m, "_sd")]] <- ss
      row[[paste0(m, "_ci_low")]] <- mu - half
      row[[paste0(m, "_ci_high")]] <- mu + half
    }
    row
  })
  if (!length(out)) return(data.frame())
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result

}

# ------------------------ SIM1: RUN + ANALYZE -------------------------------
sim1 <- NULL
sim2 <- NULL
if ("SIM1" %in% STUDIES) {
  cat("\n======================== SIM1 ========================\n")
  source("simu1.R")          # Existing rep_###.rds files are skipped.
  x <- read_replicates("SIM1")
  if (is.null(x)) stop("SIM1 produced no checkpoints")

  # Fit-quality monitoring, including nonconverged attempts.
  keys <- c("replicate", "model", "method")
  fit_quality <- do.call(rbind, lapply(
    split(seq_len(nrow(x)), interaction(x[keys], drop = TRUE, lex.order = TRUE)),
    function(ii) {
      d <- x[ii, , drop = FALSE]
      data.frame(study = "SIM1", replicate = d$replicate[1L], model = d$model[1L],
                 method = d$method[1L], n_fits = nrow(d),
                 n_finite_BIC = sum(is.finite(d$BIC)),
                 n_converged = sum(!is.na(d$converged) & d$converged),
                 convergence_rate = mean(!is.na(d$converged) & d$converged))
    }))
  utils::write.csv(fit_quality, file.path(supp_dir, "SIM1_fit_quality.csv"),
                   row.names = FALSE, na = "")

  stopifnot(all(c("gap_design", "missing_target", "persistence_MAE",
                  "stationary_TV") %in% names(x)))
  if (any(x$K != x$K_true)) stop("SIM1 contains non-oracle K fits")
  if (!setequal(unique(x$method), c("gap-aware", "compressed")))
    stop("SIM1 should contain gap-aware AND compressed")
  by <- c("replicate", "T", "K_true", "dgp", "missing_target",
          "gap_design", "model", "K")
  vars <- c("ARI", "decoded_accuracy", "TV", "KL", "persistence_MAE",
            "stationary_TV", "converged", "loglik", "BIC")
  left <- x[x$method == "gap-aware", c(by, vars), drop = FALSE]
  right <- x[x$method == "compressed", c(by, vars), drop = FALSE]
  pair <- merge(left, right, by = by, suffixes = c("_gap", "_comp"),
                all = TRUE, sort = TRUE)
  if (nrow(pair) != nrow(left) || nrow(pair) != nrow(right))
    stop("SIM1 pair matching failed: check both methods for missing fits")
  pair$both_converged <- (!is.na(pair$converged_gap) & pair$converged_gap) &
    (!is.na(pair$converged_comp) & pair$converged_comp)
  pair$delta_ARI_gap_minus_comp <- pair$ARI_gap - pair$ARI_comp
  pair$delta_accuracy_gap_minus_comp <- pair$decoded_accuracy_gap -
    pair$decoded_accuracy_comp
  pair$delta_TV_comp_minus_gap <- pair$TV_comp - pair$TV_gap
  pair$delta_KL_comp_minus_gap <- pair$KL_comp - pair$KL_gap
  pair$delta_persistence_MAE_comp_minus_gap <- pair$persistence_MAE_comp -
    pair$persistence_MAE_gap
  pair$delta_stationary_TV_comp_minus_gap <- pair$stationary_TV_comp -
    pair$stationary_TV_gap

  utils::write.csv(pair, file.path(supp_dir, "SIM1_paired_replications.csv"), row.names = FALSE, na = "")
  by_summary <- c("T", "K_true", "dgp", "gap_design", "missing_target", "model")
  outcome_cols <- c("delta_ARI_gap_minus_comp", "delta_accuracy_gap_minus_comp",
                    "delta_TV_comp_minus_gap", "delta_KL_comp_minus_gap",
                    "delta_persistence_MAE_comp_minus_gap",
                    "delta_stationary_TV_comp_minus_gap")
  main <- mc_summary(pair, by_summary, outcome_cols, flag = pair$both_converged)
  utils::write.csv(main, file.path(main_dir, "SIM1_paired_effects.csv"), row.names = FALSE, na = "")
  # Marginal state-recovery results in the two methods; helps distinguish
  # negative pair differences from a genuine loss of recovery accuracy.
  fit_by <- c("T", "K_true", "dgp", "gap_design", "missing_target", "model", "method")
  margins <- mc_summary(x, fit_by,
                              c("ARI", "decoded_accuracy", "TV", "KL"),
                              flag = (!is.na(x$converged) & x$converged))
  margins$convergence_rate <- margins$n_eligible / margins$n_total
  utils::write.csv(margins, file.path(supp_dir, "SIM1_model_recovery.csv"), row.names = FALSE, na = "")
  cat("SIM1: ", sum(pair$both_converged), "/", nrow(pair),
      " paired fits converge; inference uses only those pairs.\n", sep = "")
  sim1 <- list(pairs = pair, main = main, marginal = margins)

}

# ------------------------ SIM2: RUN + ANALYZE -------------------------------
if ("SIM2" %in% STUDIES) {
  cat("\n======================== SIM2 ========================\n")
  source("simu2.R")          # Existing rep_###.rds files are skipped.
  x <- read_replicates("SIM2")
  if (is.null(x)) stop("SIM2 produced no checkpoints")

  keys <- c("replicate", "model", "method")
  fit_quality <- do.call(rbind, lapply(
    split(seq_len(nrow(x)), interaction(x[keys], drop = TRUE, lex.order = TRUE)),
    function(ii) {
      d <- x[ii, , drop = FALSE]
      data.frame(study = "SIM2", replicate = d$replicate[1L], model = d$model[1L],
                 method = d$method[1L], n_fits = nrow(d),
                 n_finite_BIC = sum(is.finite(d$BIC)),
                 n_converged = sum(!is.na(d$converged) & d$converged),
                 convergence_rate = mean(!is.na(d$converged) & d$converged))
    }))
  utils::write.csv(fit_quality, file.path(supp_dir, "SIM2_fit_quality.csv"),
                   row.names = FALSE, na = "")

  stopifnot(all(c("background", "b_H", "b_S", "K_selected",
                  "bic_grid_complete", "bic_all_converged") %in% names(x)))
  if (!setequal(unique(x$method), "gap-aware"))
    stop("SIM2 should use only gap-aware fitting")
  by <- c("replicate", "T", "K_true", "dgp", "background", "model")
  selected <- lapply(split(seq_len(nrow(x)), interaction(x[by], drop = TRUE, lex.order = TRUE)), function(ii) {
    d <- x[ii, , drop = FALSE]
    expected_K <- seq_len(d$K_true[1L] + 2L)
    if (!setequal(d$K, expected_K)) stop("SIM2 incomplete K grid at replication ",
                                         d$replicate[1L], ", ", d$model[1L])
    complete <- all(is.finite(d$BIC))
    all_converged <- all((!is.na(d$converged) & d$converged))
    # Reconstruct selected K, do not blindly trust repeated checkpoint fields.
    ind <- if (complete) which.min(d$BIC) else NA_integer_
    out <- d[1L, by, drop = FALSE]
    out$b_H <- d$b_H[1L]
    out$b_S <- d$b_S[1L]
    out$background_contrast_H_minus_S <- out$b_H - out$b_S
    out$grid_complete <- complete
    out$all_converged <- all_converged
    out$primary_valid <- complete && all_converged
    out$provisional <- complete && !all_converged
    out$K_selected <- if (complete) d$K[ind] else NA_integer_
    out$K_correct <- if (complete) d$K[ind] == d$K_true[ind] else NA
    out$K_under <- if (complete) d$K[ind] < d$K_true[ind] else NA
    out$K_over <- if (complete) d$K[ind] > d$K_true[ind] else NA
    out$K_at_boundary <- if (complete) d$K[ind] == max(expected_K) else NA
    out$selected_converged <- if (complete) isTRUE(d$converged[ind]) else FALSE
    out$selected_BIC <- if (complete) d$BIC[ind] else NA_real_
    out$selected_ARI <- if (complete) d$ARI[ind] else NA_real_
    out$selected_accuracy <- if (complete) d$decoded_accuracy[ind] else NA_real_
    out$selected_TV <- if (complete) d$TV[ind] else NA_real_
    out$selected_KL <- if (complete) d$KL[ind] else NA_real_
    # Keep checkpoint's selected K as a consistency check when populated.
    saved <- unique(stats::na.omit(d$K_selected))
    if (complete && length(saved) && (length(saved) != 1L || saved != out$K_selected))
      stop("SIM2 saved and reconstructed selected K disagree")
    out
  })
  selection <- do.call(rbind, selected)
  rownames(selection) <- NULL
  utils::write.csv(selection, file.path(supp_dir, "SIM2_BIC_selection_by_replication.csv"), row.names = FALSE, na = "")

  sel_by <- c("T", "K_true", "dgp", "background", "model")
  summary <- mc_summary(selection, sel_by,
                               c("K_correct", "K_under", "K_over", "K_at_boundary",
                                 "selected_ARI", "selected_accuracy", "selected_TV"),
                               flag = selection$primary_valid)
  summary$complete_grid_rate <- vapply(split(seq_len(nrow(selection)), interaction(selection[sel_by], drop = TRUE, lex.order = TRUE)),
                                      function(ii) mean(selection$grid_complete[ii]), numeric(1))
  summary$all_converged_rate <- vapply(split(seq_len(nrow(selection)), interaction(selection[sel_by], drop = TRUE, lex.order = TRUE)),
                                     function(ii) mean(selection$all_converged[ii]), numeric(1))
  summary$provisional_grid_rate <- vapply(split(seq_len(nrow(selection)), interaction(selection[sel_by], drop = TRUE, lex.order = TRUE)),
                                         function(ii) mean(selection$provisional[ii]), numeric(1))
  utils::write.csv(summary, file.path(main_dir, "SIM2_background_and_K_selection.csv"), row.names = FALSE, na = "")

  # Supplementary sensitivity: all finite BIC grids, even when other
  # candidate fits did not converge. Never mislabel as definitive selection.
  # The selected fit itself must have converged to enter this sensitivity.
  loose_flag <- selection$grid_complete & selection$selected_converged
  sensitivity <- mc_summary(selection, sel_by,
                                 c("K_correct", "K_under", "K_over",
                                   "K_at_boundary", "selected_ARI"),
                                 flag = loose_flag)
  utils::write.csv(sensitivity, file.path(supp_dir, "SIM2_BIC_finite_grid_sensitivity.csv"), row.names = FALSE, na = "")

  # Selection frequencies for exact K, conditional on strict primary validity.
  good <- selection[selection$primary_valid & !is.na(selection$K_selected), , drop = FALSE]
  if (nrow(good)) {
    freq_by <- c(sel_by, "K_selected")
    freq <- do.call(rbind, lapply(split(seq_len(nrow(good)), interaction(good[freq_by], drop = TRUE, lex.order = TRUE)), function(ii) {
      d <- good[ii, , drop = FALSE]
      out <- d[1L, freq_by, drop = FALSE]
      out$n_selected <- nrow(d)
      out
    }))
    totals <- aggregate(n_selected ~ T + K_true + dgp + background + model,
                        data = freq, FUN = sum)
    names(totals)[names(totals) == "n_selected"] <- "n_valid_grids"
    freq <- merge(freq, totals, by = sel_by, all.x = TRUE, sort = TRUE)
    freq$frequency <- freq$n_selected / freq$n_valid_grids
    utils::write.csv(freq, file.path(supp_dir, "SIM2_selected_K_frequencies.csv"), row.names = FALSE, na = "")
  }

  base_by <- c("replicate", "T", "K_true", "dgp", "model")
  baseline <- selection[selection$background == "none",
                        c(base_by, "primary_valid", "K_selected", "K_correct",
                          "selected_ARI", "selected_TV"), drop = FALSE]
  alternative <- selection[selection$background != "none", , drop = FALSE]
  paired <- merge(alternative, baseline, by = base_by,
                  suffixes = c("_bg", "_none"), all.x = TRUE, sort = TRUE)
  if (nrow(paired) != nrow(alternative) || anyNA(paired$primary_valid_none))
    stop("SIM2: background/non-background pairing failed")
  paired$both_valid <- paired$primary_valid_bg & paired$primary_valid_none
  paired$delta_correct_bg_minus_none <- as.numeric(paired$K_correct_bg) -
    as.numeric(paired$K_correct_none)
  paired$delta_K_bg_minus_none <- paired$K_selected_bg - paired$K_selected_none
  paired$delta_ARI_bg_minus_none <- paired$selected_ARI_bg - paired$selected_ARI_none
  paired$delta_TV_bg_minus_none <- paired$selected_TV_bg - paired$selected_TV_none
  utils::write.csv(paired, file.path(supp_dir, "SIM2_paired_background_replications.csv"), row.names = FALSE, na = "")
  bg_by <- c("T", "K_true", "dgp", "background", "model")
  effects <- mc_summary(paired, bg_by,
                              c("delta_correct_bg_minus_none", "delta_K_bg_minus_none",
                                "delta_ARI_bg_minus_none", "delta_TV_bg_minus_none"),
                              flag = paired$both_valid)
  utils::write.csv(effects, file.path(main_dir, "SIM2_paired_background_effects.csv"), row.names = FALSE, na = "")
  cat("SIM2: ", sum(selection$primary_valid), "/", nrow(selection),
      " grids complete and fully converged. Provisional grids: ",
      sum(selection$provisional), ".\n", sep = "")
  sim2 <- list(selection = selection, main = summary, paired = paired, effects = effects)

}

# ------------------------------- FIGURES ------------------------------------
if (WRITE_FIGURES && requireNamespace("ggplot2", quietly = TRUE)) {
  gg <- asNamespace("ggplot2")
  if (!is.null(sim1) && nrow(sim1$main)) {
    a <- sim1$main
    a <- a[a$delta_ARI_gap_minus_comp_n > 0, , drop = FALSE]
    a$model <- factor(a$model, levels = c("skellam", "pois2", "bivpois"),
                      labels = c("Skellam", "Independent Poisson", "Bivariate Poisson"))
    a$case <- paste0("K=", a$K_true, " | ", a$dgp)
    p <- gg$ggplot(a, gg$aes(x = missing_target,
      y = delta_ARI_gap_minus_comp_mean, group = model, linetype = model,
      shape = model)) +
      gg$geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
      gg$geom_line(color = "grey25", linewidth = 0.4) +
      gg$geom_point(color = "grey10", size = 1.9) +
      gg$geom_errorbar(gg$aes(ymin = delta_ARI_gap_minus_comp_ci_low,
                              ymax = delta_ARI_gap_minus_comp_ci_high),
                       width = 0.008, color = "grey40", linewidth = 0.3) +
      gg$facet_grid(case ~ gap_design + T, scales = "free_y") +
      gg$scale_x_continuous(breaks = sort(unique(a$missing_target)),
                            labels = function(v) paste0(round(100*v), "%")) +
      gg$labs(x = "Fraction of missing observations",
              y = "Paired difference in ARI (gap-aware minus compressed)",
              linetype = "Emission", shape = "Emission") +
      gg$theme_bw(base_size = 9) +
      gg$theme(panel.grid.minor = gg$element_blank(),
               legend.position = "bottom", strip.background = gg$element_rect(fill = "grey94"))
    gg$ggsave(file.path(main_dir, "SIM1_gap_ARI_paired.pdf"), p,
              width = 12, height = 8, units = "in", limitsize = FALSE)
  }
  if (!is.null(sim2) && nrow(sim2$main)) {
    a <- sim2$main
    a <- a[a$K_correct_n > 0, , drop = FALSE]
    a$background <- factor(a$background,
       levels = c("none", "symmetric", "near_symmetric", "moderate_asymmetric"))
    a$model <- factor(a$model, levels = c("skellam", "pois2", "bivpois"),
                      labels = c("Skellam", "Independent Poisson", "Bivariate Poisson"))
    a$case <- paste0("K=", a$K_true, " | ", a$dgp)
    p <- gg$ggplot(a, gg$aes(x = background,
                             y = K_correct_mean, group = model,
                             linetype = model, shape = model)) +
      gg$geom_line(color = "grey25", linewidth = 0.4) +
      gg$geom_point(color = "grey10", size = 1.9) +
      gg$geom_errorbar(gg$aes(ymin = pmax(0, K_correct_ci_low),
                              ymax = pmin(1, K_correct_ci_high)),
                       width = 0.12, color = "grey40", linewidth = 0.3) +
      gg$facet_grid(case ~ T) +
      gg$coord_cartesian(ylim = c(0, 1)) +
      gg$labs(x = "Instrumental background scenario",
              y = "Proportion of replications with correct BIC-selected K",
              linetype = "Emission", shape = "Emission") +
      gg$theme_bw(base_size = 9) +
      gg$theme(panel.grid.minor = gg$element_blank(),
               axis.text.x = gg$element_text(angle = 30, hjust = 1),
               legend.position = "bottom", strip.background = gg$element_rect(fill = "grey94"))
    gg$ggsave(file.path(main_dir, "SIM2_BIC_K_recovery.pdf"), p,
              width = 11.5, height = 8, units = "in", limitsize = FALSE)
  }

} else if (WRITE_FIGURES) {
  message("ggplot2 not installed: all tables saved, figures skipped")
}

writeLines(c(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "SIM1/SIM2 launched from this script (existing complete checkpoints skipped).",
  "Decoder: marginal posterior MAP, as in the unchanged source files.",
  "SIM1: oracle K, paired gap-aware/compressed comparisons on the same data.",
  "SIM2: K chosen with BIC over 1:(K_true+2).",
  "Main BIC selection: finite complete grid and all candidates converged.",
  "Monte Carlo CI: 95% two-sided Student-t intervals across replications."
), file.path(ANALYSIS_DIR, "analysis_metadata.txt"))
writeLines(capture.output(sessionInfo()), file.path(ANALYSIS_DIR, "R_sessionInfo.txt"))
cat("\nSimulation analysis complete. Tables in: ", ANALYSIS_DIR, "\n", sep = "")
