# =============================================================================
# Bootstrap_2596012703.R
# Run AFTER Application_2596012703.R, from the project root.
#   source("Bootstrap_2596012703.R")
# =============================================================================

# ---- 1. Settings -------------------------------------------------------------
BOOT_B <- as.integer(Sys.getenv("BOOT_B", "200"))
BOOT_TOL <- as.numeric(Sys.getenv("BOOT_TOL", "1e-4"))
BOOT_MAXIT <- as.integer(Sys.getenv("BOOT_MAXIT", "1000"))
BOOT_KMEANS <- Sys.getenv("BOOT_KMEANS", "0") == "1"
LEVEL <- 0.95
MIN_RECOVERY <- 0.95
BASE_SEED <- 20261112L
SEED_OFFSET <- 300000L # original GRS_skellam position in the nine-model script
DT <- 1L
EXPECTED_K <- 7L # the model used in the revised manuscript; fail if changed
PAIR <- c(4L, 5L)
FREE_CORES <- 2L
n_phys <- parallel::detectCores(logical = FALSE)
if (is.na(n_phys)) {
  n_log <- parallel::detectCores(logical = TRUE)
  n_phys <- if (is.na(n_log)) 1L else max(1L, n_log %/% 2L)
}
BOOT_CORES <- as.integer(Sys.getenv("BOOT_CORES",
                                    as.character(max(1L, n_phys - FREE_CORES))))
OUT <- Sys.getenv("BOOT_OUT", file.path("paper_outputs", "analysis_R1",
                                        "web_supplement", "bootstrap"))
FIT_FILE <- file.path("paper_outputs", "GRS.rds")
HMM_FILE <- "hmm_functions.R"
# Also accept the uploaded filename if it has not yet been renamed.
if (!file.exists(HMM_FILE) && file.exists("hmm_functions(8).R"))
  HMM_FILE <- "hmm_functions(8).R"

stopifnot(is.finite(BOOT_B), BOOT_B >= 1L, is.finite(BOOT_CORES),
          BOOT_CORES >= 1L, is.finite(BOOT_TOL), BOOT_TOL > 0,
          is.finite(BOOT_MAXIT), BOOT_MAXIT >= 1L)
if (!file.exists(FIT_FILE) || !file.exists(HMM_FILE))
  stop("Run from the project root after Application_2596012703.R.")
for (pkg in c("extraDistr", "matrixStats"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
HMM_FILE <- normalizePath(HMM_FILE)
source(HMM_FILE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(OUT, "BOOT_GRS_skellam")

# ---- 2. Read the actual BIC-selected gap-aware model --------------------------
x <- readRDS(FIT_FILE)
stopifnot(as.character(x$obsid) == "2596012703")
sel <- x$selected
K <- sel$K[sel$dt == DT & sel$model == "skellam"]
if (length(K) != 1L || is.na(K) || K != EXPECTED_K)
  stop("Expected the paper's selected gap-aware Skellam model with K = 7.")
fit <- x$fits$dt1$skellam[[paste0("K", K)]]
if (is.null(fit) || !isTRUE(fit$converged)) stop("Selected fit missing/not converged.")
fit <- sort_skellam_states_by_contrast(fit)
expo <- x$bins$dt1$expo
stopifnot(length(expo) > 1L, all(is.finite(expo)), all(expo >= 0 & expo <= 1),
          any(expo > 0), all(dim(fit$A) == c(K, K)),
          length(fit$pi1) == K, length(fit$mu1) == K, length(fit$mu2) == K,
          all(is.finite(c(fit$pi1, fit$A, fit$mu1, fit$mu2))),
          all(fit$A >= 0), all(abs(rowSums(fit$A) - 1) < 1e-8),
          all(fit$pi1 >= 0), abs(sum(fit$pi1) - 1) < 1e-8,
          all(fit$mu1 > 0), all(fit$mu2 > 0))
fit$gamma <- NULL
rm(x)
cfg <- list(K = K, dt = DT, expo = expo, fit = fit,
            B = BOOT_B, tol = BOOT_TOL, maxit = BOOT_MAXIT,
            kmeans = BOOT_KMEANS, seed = BASE_SEED + SEED_OFFSET)

# ---- 3. Parametric bootstrap (original estimator and random draws) ------------
targets <- function(m, dt) {
  K <- length(m$mu1)
  a <- diag(m$A)
  c(setNames((m$mu1 - m$mu2) / dt, paste0("contrast_", seq_len(K))),
    setNames(a, paste0("a_", seq_len(K))),
    setNames(dt / pmax(1 - a, 1e-15), paste0("tau_", seq_len(K))))
}

refit <- function(cfg, D, from_estimate) {
  f <- cfg$fit
  tryCatch({
    if (from_estimate)
      hmm_skellam_EM_gapaware(D, cfg$expo, cfg$K,
                              max_iter = cfg$maxit, tol = cfg$tol,
                              mu1_init = f$mu1, mu2_init = f$mu2, A_init = f$A,
                              pi1_init = f$pi1, verbose = FALSE)
    else
      hmm_skellam_EM_gapaware(D, cfg$expo, cfg$K,
                              max_iter = cfg$maxit, tol = cfg$tol, verbose = FALSE)
  }, error = function(e) structure(list(msg = conditionMessage(e)),
                                   class = "refit_error"))
}

label_agree <- function(m, Z, obs, K) {
  if (!requireNamespace("clue", quietly = TRUE) ||
      is.null(m$gamma) || !is.matrix(m$gamma) ||
      nrow(m$gamma) != length(Z)) return(NA)
  g <- m$gamma[obs, , drop = FALSE]
  z <- Z[obs]
  M <- t(sapply(seq_len(K), function(l) colSums(g[z == l, , drop = FALSE])))
  tryCatch(all(as.integer(clue::solve_LSAP(M, maximum = TRUE)) == seq_len(K)),
           error = function(e) NA)
}

one_rep <- function(r, cfg) {
  tryCatch({
    set.seed(cfg$seed + r)
    f <- cfg$fit
    K <- cfg$K
    TT <- length(cfg$expo)
    obs <- cfg$expo > 0
    e <- cfg$expo[obs]
    Z <- integer(TT)
    Z[1] <- sample.int(K, 1, prob = f$pi1)
    for (t in 2:TT) Z[t] <- sample.int(K, 1, prob = f$A[Z[t - 1], ])
    zo <- Z[obs]
    n <- sum(obs)
    # Same draw order as boot_GRS_skellam.R (preserves random stream).
    h <- rpois(n, e * f$mu1[zo])
    s <- rpois(n, e * f$mu2[zo])
    D <- rep(NA_real_, TT)
    D[obs] <- h - s
    adj <- which(obs[-TT] & obs[-1])
    exits <- sapply(seq_len(K), function(k) sum(Z[adj] == k & Z[adj + 1] != k))
    fits <- list(refit(cfg, D, TRUE))
    if (cfg$kmeans) fits[[2]] <- refit(cfg, D, FALSE)
    ok <- !vapply(fits, inherits, logical(1), what = "refit_error")
    ll <- vapply(seq_along(fits), function(i)
      if (ok[i]) fits[[i]]$loglik else -Inf, numeric(1))
    ll[!is.finite(ll)] <- -Inf
    if (!any(is.finite(ll))) {
      errors <- vapply(fits[!ok], function(z) z$msg, character(1))
      stop("No finite refit likelihood. ", paste(errors, collapse = " | "))
    }
    best <- which.max(ll)
    m <- sort_skellam_states_by_contrast(fits[[best]])
    data.frame(rep = r, t(targets(m, cfg$dt)), loglik = m$loglik,
               converged = isTRUE(m$converged), iter = m$iter,
               start_used = c("estimate", "kmeans")[best],
               loglik_gap_starts = if (length(ll) == 2L && all(is.finite(ll)))
                 abs(diff(ll)) else NA_real_,
               label_agree = label_agree(m, Z, obs, K),
               t(setNames(exits, paste0("true_exits_", seq_len(K)))),
               check.names = FALSE)
  }, error = function(e) list(rep = r, error = conditionMessage(e)))
}

# ---- 4. Resume safely, without mixing different reference fits ----------------
ck <- paste0(prefix, "_reps.rds")
sig_file <- paste0(prefix, "_config.rds")
signature <- list(obsid = "2596012703", dt = DT, K = K, expo = expo,
                  fit = fit[c("mu1", "mu2", "A", "pi1")],
                  seed = cfg$seed, tol = cfg$tol, maxit = cfg$maxit,
                  kmeans = cfg$kmeans, hmm_md5 = unname(tools::md5sum(HMM_FILE)))
done <- if (file.exists(ck)) readRDS(ck) else NULL
if (!is.null(done)) {
  req <- c("rep", names(targets(fit, DT)), "loglik", "converged", "iter",
           "start_used", "loglik_gap_starts", "label_agree",
           paste0("true_exits_", seq_len(K)))
  stopifnot(is.data.frame(done), all(req %in% names(done)),
            !anyNA(done$rep), !anyDuplicated(done$rep),
            all(done$rep >= 1 & done$rep == floor(done$rep)))
  if (file.exists(sig_file)) {
    if (!identical(readRDS(sig_file), signature))
      stop("Checkpoint configuration differs. Choose a new BOOT_OUT; do not mix fits.")
  } else {
    if (Sys.getenv("BOOT_TRUST_LEGACY", "0") != "1")
      stop("Legacy checkpoint has no reference-fit signature. Verify provenance ",
           "then set BOOT_TRUST_LEGACY=1, or choose a new BOOT_OUT.")
    # Check the legacy metadata and point estimates when available.
    old_meta <- paste0(prefix, "_meta.csv")
    if (file.exists(old_meta)) {
      m <- read.csv(old_meta)
      stopifnot(nrow(m) == 1L, m$K == K, m$dt == DT,
                isTRUE(all.equal(m$tol, cfg$tol)), m$max_iter == cfg$maxit,
                as.logical(m$kmeans_start) == cfg$kmeans)
    }
    old_summary <- paste0(prefix, "_matched_summary.csv")
    if (!file.exists(old_summary)) old_summary <- paste0(prefix, "_summary.csv")
    if (file.exists(old_summary)) {
      os <- read.csv(old_summary)
      nm <- paste0(os$quantity, "_", os$state)
      est <- targets(fit, DT)
      stopifnot(all(nm %in% names(est)),
                isTRUE(all.equal(as.numeric(os$estimate),
                                 unname(est[nm]), tolerance = 1e-7)))
    }
    warning("Using a legacy checkpoint at your request; full provenance cannot be verified.")
  }
}
saveRDS(signature, sig_file)

run_missing <- function(done, cfg) {
  todo <- setdiff(seq_len(cfg$B), if (is.null(done)) integer(0) else done$rep)
  if (!length(todo)) {
    message("All requested replicates already exist; rebuilding tables only.")
    return(done)
  }
  if (Sys.getenv("BOOT_SUMMARY_ONLY", "0") == "1") {
    warning("Summary-only mode: ", length(todo), " requested replicates are missing.")
    return(done)
  }
  cores <- min(BOOT_CORES, length(todo))
  cl <- NULL
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(f) { source(f); NULL }, HMM_FILE)
    parallel::clusterExport(cl, c("targets", "refit", "label_agree", "one_rep"),
                            envir = environment(one_rep))
  }
  chunk_size <- max(3L * cores, 2L)
  consecutive_empty <- 0L
  while (length(todo)) {
    ids <- head(todo, chunk_size)
    todo <- todo[-seq_along(ids)]
    results <- if (is.null(cl)) lapply(ids, one_rep, cfg = cfg) else
      parallel::parLapplyLB(cl, ids, one_rep, cfg = cfg, chunk.size = 1L)
    good <- Filter(is.data.frame, results)
    bad <- Filter(function(z) !is.data.frame(z), results)
    if (length(bad)) {
      errors <- do.call(rbind, lapply(bad, as.data.frame))
      ef <- paste0(prefix, "_errors.csv")
      write.table(errors, ef, sep = ",", row.names = FALSE,
                  col.names = !file.exists(ef), append = file.exists(ef))
      message(length(bad), " failed replicate(s); see ", ef)
    }
    if (length(good)) {
      done <- rbind(done, do.call(rbind, good))
      # Save after each chunk; retain the old file until replacement is written.
      tmp <- paste0(ck, ".tmp")
      saveRDS(done, tmp)
      if (!file.rename(tmp, ck)) stop("Could not replace checkpoint: ", ck)
    }
    consecutive_empty <- if (length(good)) 0L else consecutive_empty + 1L
    if (consecutive_empty >= 2L) stop("Two consecutive chunks failed; inspect error log.")
    message("Completed ", if (is.null(done)) 0L else sum(done$rep <= cfg$B),
            "/", cfg$B, " replicates; converged: ",
            if (is.null(done)) 0L else sum(done$rep <= cfg$B & done$converged %in% TRUE))
  }
  done
}
done <- run_missing(done, cfg)
if (is.null(done) || !nrow(done)) stop("No bootstrap replicates available.")
d_all <- done[done$rep <= BOOT_B, , drop = FALSE]
d_all <- d_all[order(d_all$rep), , drop = FALSE]
d <- d_all[d_all$converged %in% TRUE, , drop = FALSE]
if (!nrow(d)) stop("No converged bootstrap replicate available.")
if (nrow(d) != BOOT_B)
  warning("Only ", nrow(d), " of ", BOOT_B,
          " requested replicates converged. Do not claim all converged.")

# ---- 5. Match regimes (same nearest-contrast rule as 09_bootstrap_tables.R) ---
est <- targets(fit, DT)
qs <- c("contrast", "a", "tau")
est_con <- unname(est[paste0("contrast_", seq_len(K))])
X <- as.matrix(d[, paste0("contrast_", seq_len(K)), drop = FALSE])
M <- setNames(lapply(qs, function(q)
  as.matrix(d[, paste0(q, "_", seq_len(K)), drop = FALSE])), qs)
matched <- setNames(lapply(qs, function(q) matrix(NA_real_, nrow(d), K)), qs)
sc <- sd(est_con)
if (!is.finite(sc) || sc <= 0) sc <- 1
for (i in seq_len(nrow(d))) {
  D <- outer(X[i, ], est_con, "-")^2 / sc^2
  nearest <- apply(D, 1, which.min)
  for (k in seq_len(K)) {
    j <- which(nearest == k)
    if (length(j)) {
      jj <- j[which.min(D[j, k])]
      for (q in qs) matched[[q]][i, k] <- M[[q]][i, jj]
    }
  }
}
matched_df <- data.frame(rep = d$rep)
for (q in qs)
  for (k in seq_len(K)) matched_df[[paste0(q, "_", k)]] <- matched[[q]][, k]
write.csv(matched_df, paste0(prefix, "_matched_replicates.csv"), row.names = FALSE)

# ---- 6. Conditional percentile intervals and recovery diagnostics ------------
al <- (1 - LEVEL) / 2
ci <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(v, c(al, 1 - al), names = FALSE))
}
rows <- list()
for (k in seq_len(K)) {
  rec_n <- sum(!is.na(matched$contrast[, k]))
  rec <- rec_n / nrow(d) # same denominator as original matched tables
  for (q in qs) {
    v <- matched[[q]][, k]
    interval <- ci(v)
    rows[[length(rows) + 1L]] <- data.frame(
      source = "GRS 1915+105", dt = DT, model = "skellam", K = K, state = k,
      quantity = q, estimate = unname(est[paste0(q, "_", k)]),
      boot_pct_lo = interval[1], boot_pct_hi = interval[2],
      boot_median = if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE),
      recovered = rec, n_recovered = rec_n, B_converged = nrow(d),
      median_true_exits = median(d[[paste0("true_exits_", k)]]),
      report_interval = rec_n > 0L && (q == "contrast" || rec >= MIN_RECOVERY))
  }
}
S <- do.call(rbind, rows)
write.csv(S, paste0(prefix, "_matched_summary.csv"), row.names = FALSE)
meta <- data.frame(source = "GRS 1915+105", dt = DT, model = "skellam", K = K,
                   tag = "GRS_skellam", B_requested = BOOT_B, B_done = nrow(d_all),
                   B_converged = nrow(d), median_iter = median(d$iter), tol = cfg$tol,
                   max_iter = cfg$maxit, kmeans_start = cfg$kmeans,
                   seed_first = cfg$seed + 1L,
                   share_label_agree = if (all(is.na(d$label_agree))) NA_real_
                   else mean(d$label_agree, na.rm = TRUE))
write.csv(meta, paste0(prefix, "_meta.csv"), row.names = FALSE)

# ---- 7. Paired contrast difference for the two low-intensity regimes ----------
k1 <- PAIR[1]; k2 <- PAIR[2]
diff_boot <- matched$contrast[, k2] - matched$contrast[, k1]
ok <- !is.na(diff_boot)
pc <- ci(diff_boot)
pair <- data.frame(state1 = k1, state2 = k2, contrast_difference = est_con[k2] - est_con[k1],
                   boot_pct_lo = pc[1], boot_pct_hi = pc[2], n_both_recovered = sum(ok),
                   B_converged = nrow(d), share_positive = if (any(ok)) mean(diff_boot[ok] > 0) else NA_real_)
write.csv(pair, paste0(prefix, "_contrast_difference_5_minus_4.csv"), row.names = FALSE)
