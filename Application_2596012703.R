# =============================================================================
# AOAS R1 | GRS 1915+105 | COMPLETE REAL-DATA ANALYSIS
# =============================================================================

# =============================================================================
# 0. SCIENTIFIC SETTINGS AND PROJECT DIRECTORIES
# =============================================================================
source_id <- "GRS"
source_label <- "GRS 1915+105"
obsid <- "2596012703"
dataset <- "prepared"

dt_grid <- c(10L, 5L, 1L)
analysis_dt <- 1L
K_grid <- 1:15
models <- c("skellam", "pois2", "bivpois")

workers <- 8L
max_iter <- 1000L
tol <- 1e-3
seed <- 20261021L

soft_band <- c(2, 5)
hard_band <- c(5, 12)
energy_tag <- "E2_5_E5_12"

# Optional, precomputed Bayesian Blocks JSONs; skip this benchmark if absent.
run_bayesian_blocks <- TRUE

output_dir <- file.path("paper_outputs_v2", "analysis_R1")
main_dir <- file.path(output_dir, "main_paper")
supp_dir <- file.path(output_dir, "web_supplement")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("future", "future.apply", "extraDistr", "matrixStats", "ggplot2", "cowplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
}
suppressPackageStartupMessages(library(ggplot2))

if (!file.exists("hmm_functions.R")) {
  stop("Missing original hmm_functions.R; run from the repository root")
}
source("hmm_functions.R")

project_dir <- getwd()
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# The files in data/background keep their original filenames and location.
main_sources <- "GRS"
background_sources <- "GRS"
source_specs <- data.frame(
  source = "GRS", obsid = obsid, label = source_label,
  dataset = dataset, dt = 1L,
  soft_lo = 2, soft_hi = 5, hard_lo = 5, hard_hi = 12,
  energy_tag = energy_tag, stringsAsFactors = FALSE
)

container_GRS <- list(
  source = source_label, obsid = obsid, dataset = dataset,
  analysis_dt = analysis_dt, bins = list(), fits = list(),
  fits_comp = list(), bic = data.frame(), bic_comp = data.frame(),
  selected = data.frame(), selected_comp = data.frame(),
  tables = list(), figures = list(), checks = list()
)

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
  container_GRS$bins[[dt_name]] <- bin_table
  cat("GTIs:", nrow(gti), "| exposure (s):", sum(exposure),
      "| observed bins:", length(observed), "\n")

  # 1.3 FIT ALL (model, K) PAIRS IN PARALLEL.
  # Estimation calls, seeds and sorting routines are unchanged.
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
  container_GRS$fits[[dt_name]] <- all_fits

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
  container_GRS$bic <- rbind(container_GRS$bic, bic_dt)
  # 1.5 Select the BIC optimum WITHIN each family and each dt.
  for (model in models) {
    x <- bic_dt[bic_dt$model == model & bic_dt$converged &
                  is.finite(bic_dt$BIC), ]
    if (!nrow(x)) {
      warning("No converged fit for ", model, " at dt=", dt)
      next
    }
    best <- x[which.min(x$BIC), ]
    container_GRS$selected <- rbind(container_GRS$selected,
      data.frame(source = source_id, dt = dt, model = model,
        K = best$K, BIC = best$BIC, boundary = best$K == max(K_grid)))
  }
}

# =============================================================================
# PART 1B. COMPRESSED MODEL (ORIGINAL ESTIMATORS, ORIGINAL SEEDS)
# =============================================================================
# The compressed time series ignores elapsed time at observation gaps.
# It uses the SAME retained bins as its gap-aware counterpart.

for (dt in dt_grid) {
  dt_name <- paste0("dt", dt)
  counts_comp <- container_GRS$bins[[dt_name]]
  obs_comp <- which(counts_comp$exposure_s > 0)
  H_comp <- counts_comp$H[obs_comp]
  S_comp <- counts_comp$S[obs_comp]
  D_comp <- counts_comp$D[obs_comp]
  expo_comp <- counts_comp$expo[obs_comp]

  cat("\n========== PART 1B | COMPRESSED |", source_id,
      "| dt=", dt, "s ==========\n")

  # 1B.1. Fit all compressed (model, K) pairs in parallel.
  # Same observed bins, original estimators and original random seeds.
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
    compressed_fits[[model]][[kname]] <- comp_parallel[[j]]$fit
    if (!is.na(comp_parallel[[j]]$error)) {
      message("COMPRESSED ", model, " ", kname, ": ", comp_parallel[[j]]$error)
    }
  }
  container_GRS$fits_comp[[dt_name]] <- compressed_fits

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
  container_GRS$bic_comp <- rbind(container_GRS$bic_comp, bic_comp_dt)
  # 1B.4. Independent BIC selection for the compressed HMM.
  for (model in models) {
    x <- bic_comp_dt[bic_comp_dt$model == model &
                     bic_comp_dt$converged & is.finite(bic_comp_dt$BIC), ]
    if (!nrow(x)) {
      warning("No converged compressed fit for ", model, " at dt=", dt)
      next
    }
    best <- x[which.min(x$BIC), ]
    container_GRS$selected_comp <- rbind(
      container_GRS$selected_comp,
      data.frame(source = source_id, dt = dt, model = model,
                 K = best$K, BIC = best$BIC,
                 boundary = best$K == max(K_grid)))
  }
}

# =============================================================================
# FIT CONTAINER: FINISH, VALIDATE AND SAVE
# =============================================================================
container_GRS$settings <- list(
  source_id = source_id, source_label = source_label,
  soft_band = soft_band, hard_band = hard_band,
  energy_tag = energy_tag, dt_grid = dt_grid,
  analysis_dt = analysis_dt, K_grid = K_grid, models = models,
  workers = workers, max_iter = max_iter, tol = tol, seed = seed
)

# Save the complete fitted grid from this run, including all K=1:15 fits.
# The file is a reproducible output, not a required input.
dir.create("paper_outputs", recursive = TRUE, showWarnings = FALSE)
saveRDS(container_GRS, file.path("paper_outputs", "GRS.rds"))
cat("Saved complete GRS fit container: paper_outputs/GRS.rds\n")
cat("COMPLETE GRS MODEL GRID:", nrow(container_GRS$bic), "gap-aware and",
    nrow(container_GRS$bic_comp), "compressed candidates.\n")

GRS <- container_GRS
containers <- list(GRS = GRS)
observations <- list(GRS = list(
  spec = source_specs, events = events, gti = gti,
  bins = GRS$bins$dt1, dt = 1L
))

# Application-level indexing uses the source label rather than a scalar dt.
analysis_dt <- c(GRS = 1L)

# =============================================================================
# 1. OBSERVATION GAP AND SCORPEON INSTRUMENTAL OVERVIEW
# =============================================================================

gap_geometry <- data.frame()
for (src in main_sources) {
  b <- observations[[src]]$bins
  dt <- analysis_dt[[src]]
  retained <- which(b$exposure_s>0)
  gaps <- pmax(0L,diff(retained)-1L)
  positive <- gaps[gaps>0]
  gap_geometry <- rbind(gap_geometry,data.frame(source=src,dt=dt,
    observed_bins=length(retained),total_bins=nrow(b),
    missing_bins_pct=100*(1-length(retained)/nrow(b)),
    n_internal_gaps=length(positive),
    median_gap_s=if(length(positive))median(positive)*dt else NA_real_,
    max_gap_s=if(length(positive))max(positive)*dt else NA_real_))
}

cat("\nOBSERVATION GAPS (GRS):\n")
print(gap_geometry, row.names = FALSE, digits = 4)

read_instrumental_on_grid <- function(src, bins, gti, sp, dt) {
  base <- file.path("data", "background", as.character(sp$obsid))
  sc_file <- file.path(base, "scorpeon_timeseries.csv")
  if (!file.exists(sc_file) && dt != 10L) {
    fallback <- file.path(base,paste0("scorpeon_dt10_",sp$energy_tag),
                          "scorpeon_timeseries.csv")
    if (file.exists(fallback)) {
      sc_file <- fallback
      message(src, ": using dt10 SCORPEON aligned in real time to dt", dt)
    }
  }
  if (!file.exists(sc_file)) stop("Missing SCORPEON CSV: ",sc_file)
  sc <- read.csv(sc_file)
  needed <- c("t_start_event_frame","t_stop_event_frame","valid_estimate",
              "hard_background_rate","soft_background_rate",
              "hard_fracexp","soft_fracexp")
  if (!all(needed %in% names(sc))) stop("Missing SCORPEON columns: ",sc_file)
  sc$valid_estimate <- tolower(trimws(as.character(sc$valid_estimate))) %in%
    c("true","t","1")
  sc <- sc[order(sc$t_start_event_frame), , drop=FALSE]

  coverage <- numeric(nrow(bins))
  bgH_raw <- numeric(nrow(bins))
  bgS_raw <- numeric(nrow(bins))
  for (j in seq_len(nrow(sc))) {
    if (!sc$valid_estimate[j]) next
    vals <- as.numeric(unlist(sc[j,c("t_start_event_frame","t_stop_event_frame",
                    "hard_background_rate","soft_background_rate",
                    "hard_fracexp","soft_fracexp")],use.names=FALSE))
    if (!all(is.finite(vals))) next
    left_sc <- sc$t_start_event_frame[j]
    right_sc <- sc$t_stop_event_frame[j]
    usable_fraction <- min(1,sc$hard_fracexp[j],sc$soft_fracexp[j])
    if (right_sc <= left_sc || usable_fraction <= 0) next
    first <- max(1L,floor((left_sc-min(bins$t_start))/dt)+1L)
    last <- min(nrow(bins),ceiling((right_sc-min(bins$t_start))/dt))
    if (first>last) next

    for (i in first:last) {
      if (bins$exposure_s[i] <= 0) next
      left <- max(left_sc,bins$t_start[i])
      right <- min(right_sc,bins$t_stop[i])
      if (right<=left) next
      seconds <- sum(pmax(0,pmin(right,gti$stop)-pmax(left,gti$start))) *
                 usable_fraction
      if (seconds<=0) next
      coverage[i] <- coverage[i]+seconds
      bgH_raw[i] <- bgH_raw[i]+seconds*sc$hard_background_rate[j]
      bgS_raw[i] <- bgS_raw[i]+seconds*sc$soft_background_rate[j]
    }
  }
  if (any(coverage>bins$exposure_s+1e-4))
    stop("SCORPEON intervals overlap excessively for ",src)
  valid <- bins$exposure_s>0 & coverage>=0.9*bins$exposure_s
  if (!any(valid)) stop("No sufficiently covered SCORPEON bins for ",src)
  bgH <- bgS <- rep(NA_real_,nrow(bins))
  bgH[valid] <- bgH_raw[valid]/coverage[valid]*bins$exposure_s[valid]
  bgS[valid] <- bgS_raw[valid]/coverage[valid]*bins$exposure_s[valid]

  # Independent instrumental diagnostics, not corrections to the HMM.
  n_det <- dead <- NA_real_
  for (kind in c("proxy","deadtime")) {
    f <- if (kind == "proxy") {
      file.path(base, "mkf_proxy_timeseries.csv")
    } else {
      file.path(base, "deadtime_mkf_timeseries.csv")
    }
    if (!file.exists(f)) stop("Missing instrument diagnostic: ",f)
    hk <- read.csv(f)
    if (!"time_event_frame" %in% names(hk))
      stop("Missing time_event_frame: ",f)
    ti <- hk$time_event_frame
    bi <- findInterval(ti,bins$t_start)
    gi <- findInterval(ti,gti$start)
    keep <- which(is.finite(ti)&bi>0L&bi<=nrow(bins)&gi>0L&gi<=nrow(gti))
    keep <- keep[ti[keep]<bins$t_stop[bi[keep]] &
                 ti[keep]<gti$stop[gi[keep]] &
                 bins$exposure_s[bi[keep]]>0]
    if (!length(keep)) next
    if ("num_fpm_on" %in% names(hk))
      n_det <- median(as.numeric(hk$num_fpm_on[keep]),na.rm=TRUE)
    if (kind=="deadtime" && "fpm_deadtime_fraction_mean" %in% names(hk))
      dead <- mean(as.numeric(hk$fpm_deadtime_fraction_mean[keep]),na.rm=TRUE)
  }
  hh <- sum(bins$H[valid]); ss <- sum(bins$S[valid])
  summary <- data.frame(source=src,obsid=as.character(sp$obsid),dt=dt,
    exposure_s=sum(bins$exposure_s),observed_bins=sum(bins$exposure_s>0),
    detectors_median=n_det,deadtime_fraction=dead,
    background_H_pct=if(hh>0)100*sum(bgH[valid])/hh else NA_real_,
    background_S_pct=if(ss>0)100*sum(bgS[valid])/ss else NA_real_,
    sc_coverage_pct=100*sum(bins$exposure_s[valid])/sum(bins$exposure_s))
  list(valid=valid,bgH=bgH,bgS=bgS,coverage=coverage,
       scorpeon_file=sc_file,summary=summary)
}

background_data <- list()
instrument_summary_main <- data.frame()
for (src in main_sources) {
  dat <- observations[[src]]
  bg <- read_instrumental_on_grid(src,dat$bins,dat$gti,
                                  dat$spec,dat$dt)
  background_data[[src]] <- bg
  instrument_summary_main <- rbind(instrument_summary_main,bg$summary)
}
cat("\nINSTRUMENTAL OVERVIEW: MAIN SOURCES (deadtime_fraction NOT percent)\n")
print(instrument_summary_main,row.names=FALSE,digits=5)

# =============================================================================
# 2. BIC TABLES, SELECTED MODELS AND VITERBI DECODING
# =============================================================================

K_sensitivity <- data.frame()
BIC_grid <- data.frame()
for (src in main_sources) {
  x <- containers[[src]]
  for (method in c("gap-aware","compressed")) {
    z <- if (method=="gap-aware") x$selected else x$selected_comp
    q <- if (method=="gap-aware") x$bic else x$bic_comp
    if (!nrow(z) || !nrow(q)) stop("Missing selected/BIC tables for ",src)
    z$source <- src;z$method <- method
    q$source <- src;q$method <- method
    K_sensitivity <- rbind(K_sensitivity,z)
    BIC_grid <- rbind(BIC_grid,q)
  }
}
cat("\nSUPPLEMENT: GRS K SELECTION ACROSS dt, EMISSIONS AND METHODS\n")
print(K_sensitivity,row.names=FALSE)
write.csv(K_sensitivity,file.path(supp_dir,"K_sensitivity_GRS.csv"),row.names=FALSE)
cat("\nCOMPLETE BIC GRID (all dt and both methods; numbers within each emission only):\n")
print(BIC_grid,row.names=FALSE,digits=5)
write.csv(BIC_grid,file.path(supp_dir,"BIC_complete_grid_GRS.csv"),row.names=FALSE)

# 2B. Select ONLY the final dt for scientific interpretation.
# Draw one BIC curve and three Viterbi-coloured count curves per emission model.
# For figures use Viterbi; for observed regime statistics use posterior gamma.
Viterbi_paths <- list()
state_summary_all <- data.frame()
for (src in main_sources) {
  x <- containers[[src]]
  dt <- analysis_dt[[src]]
  b <- x$bins[[paste0("dt",dt)]]
  obs <- which(b$exposure_s>0)
  lc <- b[obs,,drop=FALSE]
  Viterbi_paths[[src]] <- list()

  for (model in models) {
    sel <- x$selected[x$selected$dt==dt & x$selected$model==model,,drop=FALSE]
    if (nrow(sel)!=1L) {warning("No selected ",model," for ",src);next}
    K <- sel$K
    fit <- x$fits[[paste0("dt",dt)]][[model]][[paste0("K",K)]]
    if (is.null(fit)) stop("Selected fit missing: ",src," / ",model)
    # BIC curves show gap-aware candidates ONLY. Compressed BIC is in the
    # complete tables and in the gap-aware/compressed comparison below.
    bic <- BIC_grid[BIC_grid$source==src & BIC_grid$dt==dt &
                    BIC_grid$model==model & BIC_grid$method=="gap-aware" &
                    is.finite(BIC_grid$BIC),,drop=FALSE]
    p_BIC <- ggplot(bic,aes(K,BIC))+
      geom_line(colour="grey45",linewidth=0.6)+
      geom_point(colour="black",size=1.8)+
      scale_x_continuous(breaks=sort(unique(bic$K)))+
      theme_bw(base_size=11)+
      theme(panel.grid.minor=element_blank())+
      labs(x="Number of latent regimes K",y="BIC (gap-aware)",
           title=paste(src,model,"| dt",dt,"s"))
    print(p_BIC)
    bic_dir <- if(model=="skellam")main_dir else supp_dir
    ggsave(file.path(bic_dir,paste0(src,"_",model,"_BIC.pdf")),
           p_BIC,width=6.4,height=4)

    gamma <- fit$gamma[obs,,drop=FALSE]
    if (nrow(gamma)!=length(obs)||ncol(gamma)!=K) stop("Gamma mismatch for ",src,model)
    # Decode the FULL time grid using the original HMM Viterbi routines.
    path <- if(model=="skellam") {
      viterbi_skellam_maxplus(b$D,b$expo,fit$pi1,fit$A,fit$mu1,fit$mu2,fill_gaps="NA")
    } else if(model=="pois2") {
      viterbi_pois2_maxplus(b$H,b$S,b$expo,fit$pi1,fit$A,fit$lamH,fit$lamS,
                           fill_gaps="NA")
    } else {
      viterbi_bivpois_maxplus(b$H,b$S,b$expo,fit$pi1,fit$A,
                             fit$lamH,fit$lamS,fit$lamC,fill_gaps="NA")
    }
    if(length(path)!=nrow(b)||anyNA(path[obs])) stop("Viterbi failed for ",src,model)
    Viterbi_paths[[src]][[model]] <- path

    plot_data <- lc
    plot_data$state <- factor(path[obs],levels=seq_len(K))
    # Muted, consistent discrete colours across H/S/D for this fit.
    state_palette <- grDevices::hcl.colors(K, palette="Dark 3")
    names(state_palette) <- as.character(seq_len(K))
    plot_dir <- if(model=="skellam")main_dir else supp_dir
    p_H <- ggplot(plot_data,aes(time_min,H_rate,colour=state))+
      geom_point(size=0.58)+scale_colour_manual(values=state_palette)+theme_bw()+
      labs(x="Elapsed minutes",y="Hard counts/s",colour="Viterbi regime",
           title=paste(src,model,"| H(t)"))
    p_S <- ggplot(plot_data,aes(time_min,S_rate,colour=state))+
      geom_point(size=0.58)+scale_colour_manual(values=state_palette)+theme_bw()+
      labs(x="Elapsed minutes",y="Soft counts/s",colour="Viterbi regime",
           title=paste(src,model,"| S(t)"))
    p_D <- ggplot(plot_data,aes(time_min,D_rate,colour=state))+
      geom_point(size=0.58)+scale_colour_manual(values=state_palette)+theme_bw()+
      labs(x="Elapsed minutes",y="H - S counts/s",colour="Viterbi regime",
           title=paste(src,model,"| D(t)"))
    print(p_H);print(p_S);print(p_D)
    ggsave(file.path(plot_dir,paste0(src,"_",model,"_H_viterbi.pdf")),p_H,width=9,height=3.5)
    ggsave(file.path(plot_dir,paste0(src,"_",model,"_S_viterbi.pdf")),p_S,width=9,height=3.5)
    ggsave(file.path(plot_dir,paste0(src,"_",model,"_D_viterbi.pdf")),p_D,width=9,height=3.5)

    # Empirical observed count intensity and contrast within each inferred state.
    # Posterior weights gamma, not hard classification and not model lambda alone.
    by_state <- data.frame()
    for (k in seq_len(K)) {
      w <- gamma[,k]
      seconds <- sum(w*lc$exposure_s)
      hh <- sum(w*lc$H);ss <- sum(w*lc$S)
      full <- which(lc$expo>=0.999)
      wf <- w[full];df <- lc$D[full]
      mean_D <- if(sum(wf)>0)sum(wf*df)/sum(wf) else NA_real_
      var_D <- if(sum(wf)>0)sum(wf*(df-mean_D)^2)/sum(wf) else NA_real_
      if(model=="skellam") {
        muH <- fit$mu1[k];muS <- fit$mu2[k]
      } else if(model=="pois2") {
        muH <- fit$lamH[k];muS <- fit$lamS[k]
      } else {
        muH <- fit$lamH[k]+fit$lamC[k]
        muS <- fit$lamS[k]+fit$lamC[k]
      }
      by_state <- rbind(by_state,data.frame(source=src,dt=dt,model=model,K=K,
        state=k,occupancy_pct=100*seconds/sum(lc$exposure_s),
        n_Viterbi=sum(path[obs]==k),
        observed_H_rate=if(seconds>0)hh/seconds else NA_real_,
        observed_S_rate=if(seconds>0)ss/seconds else NA_real_,
        observed_D_rate=if(seconds>0)(hh-ss)/seconds else NA_real_,
        observed_total_rate=if(seconds>0)(hh+ss)/seconds else NA_real_,
        observed_hardness=if(hh+ss>0)(hh-ss)/(hh+ss) else NA_real_,
        observed_D_mean_full_bins=mean_D,observed_D_variance_full_bins=var_D,
        fitted_H_per_full_bin=muH,fitted_S_per_full_bin=muS,
        fitted_D_per_second=(muH-muS)/dt,
        A_kk=fit$A[k,k],mean_sojourn_s=if(fit$A[k,k]<1)dt/(1-fit$A[k,k]) else Inf))
    }
    cat("\nSELECTED STATE SUMMARY:",src,model,"K=",K,"\n")
    print(by_state,row.names=FALSE,digits=4)
    cat("TRANSITION MATRIX A:",src,model,"\n")
    print(round(fit$A,4))
    write.csv(by_state,file.path(plot_dir,paste0(src,"_",model,"_states.csv")),row.names=FALSE)
    write.csv(fit$A,file.path(plot_dir,paste0(src,"_",model,"_transition_matrix.csv")),
              row.names=FALSE)
    state_summary_all <- rbind(state_summary_all,by_state)
  }
}

# =============================================================================
# 2C. REGIME CHARACTERISTICS, GTI ALLOCATION AND NEILSEN (2020) FLARE PHASES
# =============================================================================
# This section analyses the selected GAP-AWARE models only: Skellam, independent
# bivariate Poisson and correlated bivariate Poisson, all at dt=1 second.
#
# Two distinct quantities must NOT be confused:
#   fitted lambda_H, lambda_S: parameters of the emission distribution;
#   observed H/S: posterior-weighted summaries of the actual two count series.
# For Skellam, lambda_H+lambda_S is the dispersion of D and need NOT equal
# the observed total H+S rate under emission misspecification.
#
# Phase boundaries are operationally reconstructed from durations published by
# Neilsen et al. (2020). They are relative to the START of our flaring GTI.
# Their screening and exact time origin need not coincide with our extraction.
# The original phase 4 and phase 5 overlap; 4a denotes the nonoverlapping part.
# =============================================================================

b <- GRS$bins$dt1
obs <- which(b$exposure_s > 0)
lc <- b[obs, , drop = FALSE]
mid <- (lc$t_start + lc$t_stop) / 2
lc$gti_id <- findInterval(mid, gti$start)

stopifnot(all(lc$gti_id > 0L))
stopifnot(all(mid < gti$stop[lc$gti_id]))

flare_gti <- lc$gti_id[which.max(lc$H + lc$S)]
flare_start <- gti$start[flare_gti]
flare_end <- gti$stop[flare_gti]
flare_duration <- flare_end - flare_start

# Specific observation used in the AOAS R1 analysis.
stopifnot(flare_gti == 6L, nrow(gti) == 8L)
stopifnot(abs(flare_duration - 756) < 1e-6)
stopifnot(abs(sum(lc$exposure_s) - 5218) < 1e-5)

lc$flare_seconds <- mid - flare_start

# First four published intervals: 333, 146, 158, 99 s.
# Fifth interval: final 43 s of the observed GTI, overlapping phase 4 by 23 s.
neilsen_phases <- data.frame(
  phase = c("1", "2", "3", "4", "4a", "5"),
  start_s = c(0, 333, 479, 637, 637, flare_duration - 43),
  stop_s = c(333, 479, 637, 736, flare_duration - 43, flare_duration),
  is_disjoint_partition = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

write.csv(neilsen_phases,
          file.path(supp_dir, "GRS_Neilsen_operational_phase_boundaries.csv"),
          row.names = FALSE)

# The common GTI-based observed summaries do NOT depend on the HMM emission.
GTI_observed <- data.frame()
for (g in seq_len(nrow(gti))) {
  ii <- which(lc$gti_id == g)
  E <- sum(lc$exposure_s[ii])
  hh <- sum(lc$H[ii]); ss <- sum(lc$S[ii])

  GTI_observed <- rbind(GTI_observed, data.frame(
    GTI = g, exposure_s = E, observed_bins = length(ii),
    H_rate = hh / E, S_rate = ss / E, D_rate = (hh - ss) / E,
    total_rate = (hh + ss) / E,
    H_over_S = if (ss > 0) hh / ss else NA_real_,
    gap_after_GTI_s = if (g < nrow(gti))
      gti$start[g + 1L] - gti$stop[g] else NA_real_
  ))
}
write.csv(GTI_observed,
          file.path(supp_dir, "GRS_observed_counts_by_GTI.csv"),
          row.names = FALSE)

regime_parameters <- data.frame()
regime_empirical <- data.frame()
regime_runs <- data.frame()
regime_GTI <- data.frame()
regime_phases <- data.frame()
phase_counts <- data.frame()
flare_plot_data <- data.frame()
max_K <- max(GRS$selected$K[GRS$selected$dt == 1L & GRS$selected$model %in% models])

stopifnot(max_K == 14L)

for (model in models) {

  selected <- GRS$selected[
    GRS$selected$dt == 1L & GRS$selected$model == model, , drop = FALSE
  ]
  stopifnot(nrow(selected) == 1L)

  K <- selected$K
  fit <- GRS$fits$dt1[[model]][[paste0("K", K)]]
  stopifnot(!is.null(fit), !is.null(fit$gamma))

  gamma <- fit$gamma
  if (nrow(gamma) == nrow(b)) {
    gamma <- gamma[obs, , drop = FALSE]
  }
  if (nrow(gamma) != nrow(lc) || ncol(gamma) != K) {
    stop("Posterior dimensions disagree with observed bins: ", model)
  }

  path_full <- Viterbi_paths[["GRS"]][[model]]
  stopifnot(length(path_full) == nrow(b))
  path <- as.integer(path_full[obs])
  stopifnot(!anyNA(path))

  # ----------------------------------------------------------
  # A. State-dependent emission parameters and their moments
  # ----------------------------------------------------------
  if (model == "skellam") {
    lambdaH <- fit$mu1
    lambdaS <- fit$mu2
    shared <- rep(NA_real_, K)
    varianceD <- fit$mu1 + fit$mu2
  } else if (model == "pois2") {
    lambdaH <- fit$lamH
    lambdaS <- fit$lamS
    shared <- rep(0, K)
    varianceD <- fit$lamH + fit$lamS
  } else {
    lambdaH <- fit$lamH + fit$lamC
    lambdaS <- fit$lamS + fit$lamC
    shared <- fit$lamC
    varianceD <- fit$lamH + fit$lamS
  }

  # Stationary distribution of the fitted one-second transition matrix.
  system <- t(fit$A) - diag(K)
  system[K, ] <- rep(1, K)
  rhs <- rep(0, K)
  rhs[K] <- 1
  stationary <- tryCatch(
    as.numeric(solve(system, rhs)),
    error = function(err) rep(NA_real_, K)
  )

  # ----------------------------------------------------------
  # B. Consecutive Viterbi runs. NEVER join segments across a gap.
  # ----------------------------------------------------------
  begin_run <- c(
    TRUE,
    diff(lc$bin) != 1L |
      diff(lc$gti_id) != 0L |
      diff(path) != 0L
  )
  run_number <- cumsum(begin_run)
  run_groups <- split(seq_len(nrow(lc)), run_number)

  for (run in run_groups) {
    regime_runs <- rbind(regime_runs, data.frame(
      model = model, K = K, state = path[run[1]],
      GTI = lc$gti_id[run[1]],
      first_bin = lc$bin[run[1]],
      last_bin = lc$bin[run[length(run)]],
      start_time_s = lc$t_start[run[1]],
      stop_time_s = lc$t_stop[run[length(run)]],
      observed_bins = length(run),
      observed_duration_s = sum(lc$exposure_s[run])
    ))
  }

  # ----------------------------------------------------------
  # C. Tables of model parameters and empirical state summaries
  # ----------------------------------------------------------
  for (k in seq_len(K)) {
    w <- gamma[, k]
    E <- sum(w * lc$exposure_s)

    Hsum <- sum(w * lc$H)
    Ssum <- sum(w * lc$S)
    Hrate <- Hsum / E
    Srate <- Ssum / E

    full <- which(lc$expo >= 0.999)
    w_full <- w[full]
    D_full <- lc$D[full]
    D_mean <- sum(w_full * D_full) / sum(w_full)
    D_variance <- sum(w_full * (D_full - D_mean)^2) / sum(w_full)

    lengths <- regime_runs$observed_duration_s[
      regime_runs$model == model & regime_runs$state == k
    ]

    sojourn <- if (fit$A[k, k] < 1) {
      1 / (1 - fit$A[k, k])
    } else {
      Inf
    }

    regime_parameters <- rbind(regime_parameters, data.frame(
      model = model, K = K, state = k,
      lambda_H = lambdaH[k], lambda_S = lambdaS[k],
      lambda_shared = shared[k],
      delta = lambdaH[k] - lambdaS[k],
      intensity = lambdaH[k] + lambdaS[k],
      q = (lambdaH[k] - lambdaS[k]) / (lambdaH[k] + lambdaS[k]),
      variance_D = varianceD[k],
      A_kk = fit$A[k, k],
      expected_sojourn_s = sojourn,
      stationary_prob = stationary[k]
    ))

    regime_empirical <- rbind(regime_empirical, data.frame(
      model = model, K = K, state = k,
      effective_exposure_s = E,
      posterior_occupancy_pct = 100 * E / sum(lc$exposure_s),
      Viterbi_bins = sum(path == k),
      observed_H_rate = Hrate,
      observed_S_rate = Srate,
      observed_D_rate = Hrate - Srate,
      observed_total_rate = Hrate + Srate,
      observed_H_over_S = if (Srate > 0) Hrate / Srate else NA_real_,
      observed_q = if (Hrate + Srate > 0)
        (Hrate - Srate) / (Hrate + Srate) else NA_real_,
      observed_D_variance_full_bins = D_variance,
      variance_to_observed_total = D_variance / (Hrate + Srate),
      fitted_varD_to_observed_total = varianceD[k] / (Hrate + Srate),
      n_Viterbi_runs = length(lengths),
      median_observed_run_s = if (length(lengths)) median(lengths) else NA_real_,
      max_observed_run_s = if (length(lengths)) max(lengths) else NA_real_
    ))
  }

  # Export the ORIGINAL full fitted transition matrix as well.
  write.csv(fit$A,
            file.path(supp_dir, paste0("GRS_", model, "_transition_matrix.csv")),
            row.names = FALSE)

  # ----------------------------------------------------------
  # D. Occupancy by each of the eight NICER observing GTIs
  # ----------------------------------------------------------
  GTI_table <- data.frame()
  for (g in seq_len(nrow(gti))) {
    ii <- which(lc$gti_id == g)
    E <- sum(lc$exposure_s[ii])

    for (k in seq_len(K)) {
      posterior_seconds <- sum(gamma[ii, k] * lc$exposure_s[ii])
      viterbi_seconds <- sum(lc$exposure_s[ii][path[ii] == k])

      GTI_table <- rbind(GTI_table, data.frame(
        model = model, K = K, GTI = g, state = k,
        GTI_exposure_s = E,
        state_posterior_seconds = posterior_seconds,
        posterior_pct_within_GTI = 100 * posterior_seconds / E,
        Viterbi_pct_within_GTI = 100 * viterbi_seconds / E
      ))
    }
  }

  # Among the observed seconds assigned to state k, which GTIs contain them?
  
  GTI_table$pct_state_exposure_in_GTI <- NA_real_
  
  for (k in seq_len(K)) {
    
    rows <- which(GTI_table$state == k)
    
    state_total <- sum(
      GTI_table$state_posterior_seconds[rows]
    )
    
    GTI_table$pct_state_exposure_in_GTI[rows] <-
      100 * GTI_table$state_posterior_seconds[rows] /
      state_total
  }
  
  regime_GTI <- rbind(regime_GTI, GTI_table)
  
  # ----------------------------------------------------------
  # E. Occupancy in Neilsen intervals and observational summaries
  # The original interval 4 and final interval 5 overlap by 23 seconds.
  # These rows must NOT be summed. The 4a row gives a disjoint alternative.
  # ----------------------------------------------------------
  for (r in seq_len(nrow(neilsen_phases))) {
    phase <- neilsen_phases$phase[r]
    ii <- which(
      lc$gti_id == flare_gti &
      lc$flare_seconds >= neilsen_phases$start_s[r] &
      lc$flare_seconds < neilsen_phases$stop_s[r]
    )
    E <- sum(lc$exposure_s[ii])
    Hphase <- sum(lc$H[ii]); Sphase <- sum(lc$S[ii])

    if (model == "skellam") {
      phase_counts <- rbind(phase_counts, data.frame(
        phase = phase, exposure_s = E,
        H_rate = Hphase / E, S_rate = Sphase / E,
        total_rate = (Hphase + Sphase) / E,
        D_rate = (Hphase - Sphase) / E,
        H_over_S = if (Sphase > 0) Hphase / Sphase else NA_real_
      ))
    }

    for (k in seq_len(K)) {
      posterior_s <- sum(gamma[ii, k] * lc$exposure_s[ii])
      viterbi_s <- sum(lc$exposure_s[ii][path[ii] == k])
      regime_phases <- rbind(regime_phases, data.frame(
        model = model, K = K, phase = phase,
        state = k, phase_exposure_s = E,
        posterior_pct = 100 * posterior_s / E,
        Viterbi_pct = 100 * viterbi_s / E
      ))
    }
  }

  # ----------------------------------------------------------
  # F. Common flare figure: observed D coloured by model-specific Viterbi state
  # ----------------------------------------------------------
  ii <- which(lc$gti_id == flare_gti)
  flare_plot_data <- rbind(flare_plot_data, data.frame(
    model = model,
    seconds = lc$flare_seconds[ii],
    D = lc$D[ii],
    state = factor(path[ii], levels = seq_len(max_K))
  ))

  cat("\nREGIME OCCUPANCY BY GTI:", model, "K=", K, "\n")
  print(round(xtabs(posterior_pct_within_GTI ~ GTI + state,
                     data = GTI_table), 1))
}

# ------------------------------------------------------------
# G. Export MAIN-PAPER tables (Skellam) and SUPPLEMENT tables (all models)
# ------------------------------------------------------------
write.csv(regime_parameters,
          file.path(supp_dir, "GRS_regime_fitted_parameters_all.csv"),
          row.names = FALSE)
write.csv(regime_empirical,
          file.path(supp_dir, "GRS_regime_observed_statistics_all.csv"),
          row.names = FALSE)
write.csv(regime_runs,
          file.path(supp_dir, "GRS_Viterbi_observed_runs_all.csv"),
          row.names = FALSE)
write.csv(regime_GTI,
          file.path(supp_dir, "GRS_state_occupancy_by_GTI_all.csv"),
          row.names = FALSE)
write.csv(regime_phases,
          file.path(supp_dir, "GRS_Neilsen_phase_occupancy_all.csv"),
          row.names = FALSE)
write.csv(phase_counts,
          file.path(supp_dir, "GRS_Neilsen_phase_observed_counts.csv"),
          row.names = FALSE)

sk_params <- regime_parameters[regime_parameters$model == "skellam", ]
sk_empirical <- regime_empirical[regime_empirical$model == "skellam", ]

write.csv(sk_params,
          file.path(main_dir, "GRS_skellam_fitted_parameters.csv"),
          row.names = FALSE)
write.csv(sk_empirical,
          file.path(main_dir, "GRS_skellam_regime_description.csv"),
          row.names = FALSE)
write.csv(data.frame(
  state = sk_params$state,
  fitted_varD = sk_params$variance_D,
  observed_varD = sk_empirical$observed_D_variance_full_bins,
  observed_total = sk_empirical$observed_total_rate,
  ratio = sk_empirical$variance_to_observed_total
), file.path(supp_dir, "GRS_skellam_dispersion_diagnostic.csv"),
row.names = FALSE)

# ------------------------------------------------------------
# H. Publishable flare figure: three emission models on the SAME observed data
# ------------------------------------------------------------
selected_one_second <- GRS$selected[GRS$selected$dt == 1L, ]
selected_one_second <- selected_one_second[
  match(models, selected_one_second$model), , drop = FALSE
]
flare_plot_data$model <- factor(
  flare_plot_data$model, levels = models,
  labels = paste0(
    c("Skellam", "Independent bivariate Poisson", "Correlated bivariate Poisson"),
    " (K = ", selected_one_second$K, ")"
  )
)

state_colors <- grDevices::hcl.colors(max_K, palette = "Dark 3")
names(state_colors) <- as.character(seq_len(max_K))

p_three_flare <- ggplot(
  flare_plot_data,
  aes(x = seconds, y = D, colour = state)
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.25) +
  geom_point(size = 0.65, alpha = 0.85) +
  geom_vline(xintercept = c(333, 479, 637, 713),
             linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  facet_grid(model ~ .) +
  scale_colour_manual(values = state_colors, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 333, 479, 637, 713, 756)) +
  coord_cartesian(xlim = c(0, flare_duration)) +
  labs(x = "Seconds since the beginning of flaring GTI 6",
       y = "Observed hard-soft contrast (counts per 1-s bin)",
       colour = "Regime (model-specific)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

print(p_three_flare)
ggsave(file.path(supp_dir, "GRS_three_emissions_Neilsen_flare.pdf"),
       p_three_flare, width = 11.3, height = 9.2)

# Model-specific state numbers have no one-to-one correspondence.
# The extra gray vertical lines mark RECONSTRUCTED spectral intervals.


# ------------------------------------------------------------
# I. BIC curves for the three GAP-AWARE emission models
# ------------------------------------------------------------

bic_plot <- GRS$bic[
  GRS$bic$dt == 1L &
    GRS$bic$model %in% models &
    is.finite(GRS$bic$BIC),
  c("model", "K", "BIC"),
  drop = FALSE
]

bic_plot$model <- factor(
  bic_plot$model,
  levels = models,
  labels = c(
    "Skellam",
    "Independent BP",
    "Correlated BP"
  )
)

p_bic <- ggplot(
  bic_plot,
  aes(x = K, y = BIC)
) +
  geom_line(
    colour = "grey40",
    linewidth = 0.65
  ) +
  geom_point(
    colour = "black",
    size = 1.8
  ) +
  facet_wrap(
    ~ model,
    scales = "free_y",
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = 1:15
  ) +
  labs(
    x = "Number of latent regimes",
    y = "BIC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey93")
  )

print(p_bic)

ggsave(
  file.path(main_dir, "GRS_BIC_three_emissions_K15.pdf"),
  p_bic,
  width = 11.5,
  height = 4.2
)


# =============================================================================
# 3. GAP-AWARE VERSUS COMPRESSED: SAME DATA, SAME EMISSION
# =============================================================================

gap_comparison_all <- data.frame()
for (src in main_sources) {
  x <- containers[[src]]
  dt <- analysis_dt[[src]]
  b <- x$bins[[paste0("dt",dt)]]
  obs <- which(b$exposure_s>0)
  for (model in models) {
    sg <- x$selected[x$selected$dt==dt & x$selected$model==model,,drop=FALSE]
    sc <- x$selected_comp[x$selected_comp$dt==dt & x$selected_comp$model==model,,drop=FALSE]
    if(nrow(sg)!=1L||nrow(sc)!=1L) {warning("Missing selection for ",src,model);next}
    Kg <- sg$K;Kc <- sc$K
    gap <- x$fits[[paste0("dt",dt)]][[model]][[paste0("K",Kg)]]
    comp_ind <- x$fits_comp[[paste0("dt",dt)]][[model]][[paste0("K",Kc)]]
    comp_same <- x$fits_comp[[paste0("dt",dt)]][[model]][[paste0("K",Kg)]]
    if (is.null(gap) || is.null(comp_ind) || is.null(comp_same)) {
      stop("Comparison fit missing for ", src, " / ", model)
    }
    # Independent-K comparison: only the same emission family/data, not A rows.
    z_gap <- max.col(gap$gamma[obs,,drop=FALSE],ties.method="first")
    z_comp <- max.col(comp_ind$gamma,ties.method="first")
    ari_agreement <- aoas_adjusted_rand(z_gap,z_comp)
    e_gap <- aoas_log_emission(gap,model,b$H[obs],b$S[obs],b$expo[obs])
    e_comp_ind <- aoas_log_emission(comp_ind,model,b$H[obs],b$S[obs],b$expo[obs])
    score_gap <- aoas_prequential_nll(gap$A,gap$pi1,e_gap,obs,gap_aware=TRUE)
    score_comp_ind <- aoas_prequential_nll(comp_ind$A,comp_ind$pi1,
                                           e_comp_ind,obs,gap_aware=FALSE)
    # Fixed-K comparison: same K as gap-aware, state matching by the original
    # helper. TV is computed on transition behaviour across observation gaps.
    e_comp_same <- aoas_log_emission(comp_same,model,b$H[obs],b$S[obs],b$expo[obs])
    score_comp_same <- aoas_prequential_nll(comp_same$A,comp_same$pi1,
                                            e_comp_same,obs,gap_aware=FALSE)
    diag <- aoas_compare_gap_fits(gap,comp_same,model,b$H,b$S,b$expo,obs,dt)
    tv <- diag$summary$mean_TV_equal_states[1]
    nll_gain <- score_comp_same$mean_crossgap_NLL-score_gap$mean_crossgap_NLL
    bic_comp_same <- bic_hmm(comp_same$loglik,Kg,length(obs),model)
    row <- data.frame(source=src,dt=dt,model=model,
      K_gap=Kg,K_compressed_independent=Kc,
      BIC_gap=sg$BIC,BIC_compressed_independent=sc$BIC,
      BIC_compressed_same_K=bic_comp_same,
      delta_BIC_same_K=bic_comp_same-sg$BIC,
      ARI_gap_vs_compressed_independent_K=ari_agreement,
      n_evaluated_gaps=score_gap$n_gaps,
      postgap_NLL_gap=score_gap$mean_crossgap_NLL,
      postgap_NLL_compressed_independent_K=score_comp_ind$mean_crossgap_NLL,
      postgap_NLL_compressed_same_K=score_comp_same$mean_crossgap_NLL,
      delta_postgap_NLL_comp_minus_gap_same_K=nll_gain,
      postgap_NLL_relative_improvement_pct=if(is.finite(nll_gain)&&
        is.finite(score_comp_same$mean_crossgap_NLL)&&
        score_comp_same$mean_crossgap_NLL>0)
        100*nll_gain/score_comp_same$mean_crossgap_NLL else NA_real_,
      transition_TV_at_gaps=tv,
      transition_TV_discrepancy_pct=100*tv,
      mean_NLL_all_bins_gap=mean(score_gap$NLL_all_bins,na.rm=TRUE),
      mean_NLL_all_bins_comp=mean(score_comp_same$NLL_all_bins,na.rm=TRUE))
    gap_comparison_all <- rbind(gap_comparison_all,row)
    write.csv(diag$by_gap,file.path(supp_dir,paste0(src,"_",model,"_gap_details.csv")),
              row.names=FALSE)
    write.csv(diag$by_state,file.path(supp_dir,paste0(src,"_",model,"_transition_TV_by_state.csv")),
              row.names=FALSE)
  }
}
cat("\nGAP-AWARE VERSUS COMPRESSED, SELECTED dt:\n")
print(gap_comparison_all,row.names=FALSE,digits=4)
write.csv(gap_comparison_all,file.path(supp_dir,"gapaware_vs_compressed_all_models.csv"),
          row.names=FALSE)
write.csv(gap_comparison_all[gap_comparison_all$model=="skellam",,drop=FALSE],
          file.path(main_dir,"gapaware_vs_compressed_skellam.csv"),row.names=FALSE)

# =============================================================================
# 5. BAYESIAN BLOCKS: PRECOMPUTED GRS JSON (OPTIONAL)
# =============================================================================
# To reproduce BB benchmarks supply the existing two GRS JSON files under
# results/bayesian_blocks/2596012703/. Their generation is independent of R.

bb_tolerance_bins <- 1L

# Maximum-cardinality, one-to-one matching of two sorted boundary-position sets.
# Both inputs are ORIGINAL time-grid bin indices, never compressed indices.
match_bb_switches <- function(hmm_boundaries, bb_boundaries, tolerance) {
  i <- j <- 1L
  matches <- 0L
  while (i <= length(hmm_boundaries) && j <= length(bb_boundaries)) {
    d <- hmm_boundaries[i] - bb_boundaries[j]
    if (abs(d) <= tolerance) {
      matches <- matches + 1L
      i <- i + 1L
      j <- j + 1L
    } else if (d < -tolerance) {
      i <- i + 1L
    } else {
      j <- j + 1L
    }
  }
  matches
}

bb_paths <- character(0)
for (src in main_sources) {
  sp <- source_specs[source_specs$source == src, , drop = FALSE]
  dt <- analysis_dt[[src]]
  for (signal in c("contrast", "total")) {
    bb_paths <- c(bb_paths, file.path(
      "results", "bayesian_blocks", as.character(sp$obsid),
      paste0("bayesian_blocks_compressed_", signal, ".json")
    ))
  }
}

missing_bb <- bb_paths[!file.exists(bb_paths)]
if (!run_bayesian_blocks) {
  message("BAYESIAN BLOCKS: disabled by user setting.")
} else if (length(missing_bb)) {
  message("BAYESIAN BLOCKS NOT ANALYSED: missing ", length(missing_bb), " JSON file(s):\n",
          paste(missing_bb, collapse = "\n"),
          "\nGenerate them in Python, place them at the exact paths and rerun run_analysis.R.")
} else {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Install R package jsonlite to read Bayesian Blocks JSONs")

  bb_summary <- data.frame()
  for (src in main_sources) {
    x <- containers[[src]]
    sp <- source_specs[source_specs$source == src, , drop = FALSE]
    dt <- analysis_dt[[src]]
    b <- x$bins[[paste0("dt", dt)]]
    gti <- observations[[src]]$gti

    for (signal in c("contrast", "total")) {
      json_file <- file.path(
        "results", "bayesian_blocks", as.character(sp$obsid),
        paste0("bayesian_blocks_compressed_", signal, ".json")
      )
      js <- jsonlite::fromJSON(json_file)
      if (!isTRUE(all.equal(as.numeric(js$p0), 0.05)) ||
          as.character(js$source) != src ||
          as.character(js$obsid) != as.character(sp$obsid) ||
          as.character(js$variant) != as.character(sp$dataset) ||
          as.integer(js$dt) != dt || as.character(js$signal) != signal)
        stop("Stale or incompatible Bayesian Blocks metadata: ", json_file)
      if (!isTRUE(all.equal(as.numeric(unlist(js$bands_keV$soft)),
                            as.numeric(c(sp$soft_lo, sp$soft_hi)))) ||
          !isTRUE(all.equal(as.numeric(unlist(js$bands_keV$hard)),
                            as.numeric(c(sp$hard_lo, sp$hard_hi)))))
        stop("Bayesian Blocks energy bands disagree with HMM: ", json_file)

      rows <- match(round(as.numeric(js$bin_left), 6), round(b$t_start, 6))
      if (!length(rows) || anyNA(rows) || anyDuplicated(rows) ||
          is.unsorted(rows, strictly = TRUE) ||
          any(b$exposure_s[rows] < dt - 1e-6))
        stop("BB full-exposure grid does not match HMM: ", json_file)
      block <- as.integer(js$bin_block)
      if (length(block) != length(rows) || anyNA(block) ||
          any(diff(block) < 0L) || any(block < 1L) ||
          as.integer(js$n_bins) != length(rows))
        stop("Invalid BB block labels / bin count: ", json_file)
      if (!is.data.frame(js$blocks) || !nrow(js$blocks) ||
          length(unique(block)) != nrow(js$blocks))
        stop("Inconsistent number of Bayesian Blocks: ", json_file)

      # Validate the actual HMM-observed count series against every BB segment,
      # guarding against stale JSON generated from different processed events.
      value <- if (signal == "contrast") b$D[rows] else b$total[rows]
      for (q in seq_len(nrow(js$blocks))) {
        start_q <- as.integer(js$blocks$first_bin[q]) + 1L
        end_q <- start_q + as.integer(js$blocks$n_bins[q]) - 1L
        if (!is.finite(start_q) || !is.finite(end_q) ||
            start_q < 1L || end_q > length(rows) || end_q < start_q ||
            !all(block[start_q:end_q] == q) ||
            abs(sum(value[start_q:end_q]) - as.numeric(js$blocks$sum[q])) > 1e-6)
          stop("BB counts/segments differ from HMM observed counts: ", json_file)
      }

      # Match only switch boundaries within continuous, fully exposed GTIs.
      gti_id <- findInterval(b$t_start[rows] + dt / 2, gti$start)
      eligible <- which(diff(rows) == 1L &
                        diff(gti_id) == 0L &
                        abs(diff(b$t_start[rows]) - dt) < 1e-6) + 1L
      bb_changepoints <- eligible[block[eligible] != block[eligible - 1L]]

      for (model in models) {
        zfull <- Viterbi_paths[[src]][[model]]
        if (is.null(zfull)) {
          warning("No Viterbi path: ", src, " / ", model)
          next
        }
        z <- as.integer(zfull[rows])
        if (anyNA(z)) stop("Missing Viterbi label: ", src, " / ", model)
        hmm_changepoints <- eligible[z[eligible] != z[eligible - 1L]]
        # Match separately within each GTI using ORIGINAL-grid indices.
        # This prevents an apparent one-bin match at two ends of a large gap.
        n_matched <- 0L
        for (gi in sort(unique(gti_id[eligible]))) {
          hmm_g <- rows[hmm_changepoints[gti_id[hmm_changepoints] == gi]]
          bb_g <- rows[bb_changepoints[gti_id[bb_changepoints] == gi]]
          n_matched <- n_matched + match_bb_switches(
            hmm_g, bb_g, bb_tolerance_bins)
        }
        n_hmm <- length(hmm_changepoints)
        n_bb <- length(bb_changepoints)
        precision <- if (n_bb) n_matched / n_bb else NA_real_
        recall <- if (n_hmm) n_matched / n_hmm else NA_real_
        f1 <- if (n_bb + n_hmm) 2 * n_matched / (n_bb + n_hmm) else NA_real_

        bb_row <- data.frame(
          source = src, dt = dt, signal = signal, model = model,
          bb_p0 = as.numeric(js$p0), matching_tolerance_bins = bb_tolerance_bins,
          matching_tolerance_seconds = bb_tolerance_bins * dt,
          n_full_exposure_bins = length(rows),
          n_eligible_within_GTI_boundaries = length(eligible),
          n_blocks_BB = length(unique(block)),
          n_BB_changepoints_within_GTI = n_bb,
          n_HMM_changepoints_within_GTI = n_hmm,
          n_matched_changepoints = n_matched,
          BB_change_precision_vs_HMM = precision,
          BB_change_recall_of_HMM = recall,
          BB_HMM_change_F1 = f1,
          n_BB_boundaries_excluded_at_discontinuities =
            sum(diff(block) != 0L) - n_bb,
          n_HMM_switches_excluded_at_discontinuities =
            sum(diff(z) != 0L) - n_hmm,
          n_distinct_HMM_states = length(unique(z)),
          n_BB_blocks_spanning_gaps = as.integer(js$n_blocks_spanning_gaps),
          n_partial_bin_photons_excluded = as.integer(js$excluded_partial_bin_counts))
        bb_summary <- rbind(bb_summary, bb_row)
        cat("BAYESIAN BLOCKS:", src, signal, model,
            "| BB within-GTI changes=", n_bb,
            "| HMM within-GTI changes=", n_hmm,
            "| matched=", n_matched, "\n")
      }

      # Show all segmented bins, but mark only BB boundaries inside continuous
      # GTIs. Boundaries at observation discontinuities are not meaningful times.
      plot_data <- b[rows, , drop = FALSE]
      plot_data$bb_value <- if (signal == "contrast")
        plot_data$D_rate else plot_data$total_rate
      p_bb <- ggplot(plot_data, aes(time_min, bb_value)) +
        geom_point(size = 0.5, colour = "grey30") +
        geom_vline(xintercept = plot_data$time_min[bb_changepoints],
                   colour = "grey65", linewidth = 0.3, linetype = "dotted") +
        theme_bw() +
        labs(x = "Elapsed minutes",
             y = if (signal == "contrast") "H - S counts/s" else "H + S counts/s",
             title = paste(src, signal, "| Bayesian Blocks (within-GTI changes)"))
      print(p_bb)
      ggsave(file.path(supp_dir, paste0(src, "_BB_", signal, ".pdf")),
             p_bb, width = 9, height = 3.5)
    }
  }
  if (nrow(bb_summary) != length(main_sources) * 2L * length(models))
    warning("Bayesian Blocks comparison is incomplete; inspect missing fits")
  cat("\nBAYESIAN BLOCKS CHANGEPOINT COMPARISON:\n")
  print(bb_summary, row.names = FALSE, digits = 4)
  write.csv(bb_summary, file.path(supp_dir, "Bayesian_Blocks_comparison.csv"),
            row.names = FALSE)
}

# =============================================================================
# 6. SCORPEON BACKGROUND DIAGNOSTICS FOR GRS ONLY
# =============================================================================
# Background estimates remain external diagnostics: no count subtraction or refit.

background_summary_all <- instrument_summary_main
background_by_state <- data.frame()
background_pairwise <- data.frame()

for (src in background_sources) {
  sp <- source_specs[source_specs$source==src,,drop=FALSE]
  dt <- analysis_dt[[src]]
  x <- containers[[src]]
  b <- x$bins[[paste0("dt",dt)]]

  bg <- background_data[[src]]
  valid <- which(bg$valid)
  observed <- which(b$exposure_s>0)
  loc <- match(valid,observed)
  if (anyNA(loc)) stop("Background/HMM observation alignment failed: ",src)

  cat("\nSCORPEON ROBUSTNESS | ",src," | dt=",dt," s\n",sep="")
  for (model in models) {
    sel <- x$selected[x$selected$dt==dt & x$selected$model==model,,drop=FALSE]
    if (nrow(sel)!=1L) stop("No selected fit: ",src," / ",model)
    fit <- x$fits[[paste0("dt",dt)]][[model]][[paste0("K",sel$K)]]
    if (is.null(fit)) stop("Selected fit missing: ",src," / ",model)
    gamma <- fit$gamma[observed,,drop=FALSE]
    if (ncol(gamma)!=sel$K || nrow(gamma)!=length(observed))
      stop("Wrong posterior dimensions: ",src," / ",model)

    state_table <- data.frame()
    for (k in seq_len(sel$K)) {
      w <- gamma[loc,k]                   # posterior Pr(Z_t=k | data)
      seconds <- sum(w*b$exposure_s[valid])
      if (seconds<=0) next
      Hobs <- sum(w*b$H[valid])           # actually observed photons
      Sobs <- sum(w*b$S[valid])
      Hbg <- sum(w*bg$bgH[valid])         # SCORPEON estimated photons
      Sbg <- sum(w*bg$bgS[valid])

      state_table <- rbind(state_table,data.frame(
        source=src,dt=dt,model=model,K=sel$K,state=k,
        exposure_s=seconds,
        H_observed_rate=Hobs/seconds,S_observed_rate=Sobs/seconds,
        D_observed_rate=(Hobs-Sobs)/seconds,
        D_background_rate=(Hbg-Sbg)/seconds,
        D_residual_rate=(Hobs-Sobs-Hbg+Sbg)/seconds,
        background_total_rate=(Hbg+Sbg)/seconds,
        background_H_pct=if(Hobs>0)100*Hbg/Hobs else NA_real_,
        background_S_pct=if(Sobs>0)100*Sbg/Sobs else NA_real_,
        background_sd_D_rate=sqrt(((Hbg+Sbg)/seconds)/dt)))
    }
    print(state_table,row.names=FALSE,digits=5)
    background_by_state <- rbind(background_by_state,state_table)

    # Compare differences between every pair of inferred regimes. This is
    # more interpretable than plotting hundreds of background points over time.
    if (nrow(state_table)>=2L) {
      pairs <- combn(seq_len(nrow(state_table)),2L)
      for (j in seq_len(ncol(pairs))) {
        a <- pairs[1,j]; c <- pairs[2,j]
        delta_obs <- state_table$D_observed_rate[c]-state_table$D_observed_rate[a]
        delta_bg <- state_table$D_background_rate[c]-state_table$D_background_rate[a]
        background_pairwise <- rbind(background_pairwise,data.frame(
          source=src,dt=dt,model=model,
          state_from=state_table$state[a],state_to=state_table$state[c],
          delta_observed_D=delta_obs,
          delta_background_D=delta_bg,
          delta_residual_D=delta_obs-delta_bg,
          delta_background_total=state_table$background_total_rate[c]-
                                 state_table$background_total_rate[a],
          abs_background_to_observed_ratio=if(abs(delta_obs)>1e-10)
            abs(delta_bg)/abs(delta_obs) else NA_real_))
      }
    }
  }
}

cat("\nINSTRUMENTAL SUMMARY: GRS (deadtime is a fraction, not %)\n")
print(background_summary_all,row.names=FALSE,digits=5)
cat("\nBACKGROUND PAIRWISE STATE CONTRASTS (all emission families)\n")
print(background_pairwise,row.names=FALSE,digits=4)

# Only final interpretive results are exported; no redundant graphics or
# time-series CSV dumps. The per-state tables record K explicitly.
write.csv(background_summary_all,
          file.path(supp_dir,"background_instrument_GRS.csv"),row.names=FALSE)
write.csv(background_by_state,
          file.path(supp_dir,"background_by_state_GRS.csv"),row.names=FALSE)
write.csv(background_pairwise,
          file.path(supp_dir,"background_pairwise_GRS.csv"),row.names=FALSE)

cat("\nAOAS R1 APPLICATION ANALYSIS COMPLETE\n")
cat("Main paper directory: ",main_dir,"\n",sep="")
cat("Web supplement directory: ",supp_dir,"\n",sep="")

# =============================================================================
# 6B. SCORPEON BACKGROUND LIGHT CURVES (ORIGINAL GRS BAND DEFINITIONS)
# =============================================================================
# Rates are aligned with the HMM grid using the original SCORPEON integration
# above. Lines are broken when the observed grid contains an interruption.

bg <- background_data[["GRS"]]
b <- GRS$bins$dt1
valid <- which(bg$valid)

background_plot <- data.frame(
  bin = b$bin[valid],
  time_min = b$time_min[valid],
  Soft = bg$bgS[valid] / b$exposure_s[valid],
  Hard = bg$bgH[valid] / b$exposure_s[valid]
)
background_plot$Contrast <- background_plot$Hard - background_plot$Soft
background_plot$segment <- cumsum(c(
  1L, as.integer(diff(background_plot$bin) > 1L)
))

background_long <- rbind(
  data.frame(time_min = background_plot$time_min,
             rate = background_plot$Soft, series = "Soft",
             segment = background_plot$segment,
             panel = "(a) Hard and soft background"),
  data.frame(time_min = background_plot$time_min,
             rate = background_plot$Hard, series = "Hard",
             segment = background_plot$segment,
             panel = "(a) Hard and soft background"),
  data.frame(time_min = background_plot$time_min,
             rate = background_plot$Contrast, series = "Contrast",
             segment = background_plot$segment,
             panel = "(b) Hard-soft background contrast")
)
background_long$series <- factor(background_long$series,
                                 levels = c("Soft", "Hard", "Contrast"))
background_long$panel <- factor(background_long$panel,
  levels = c("(a) Hard and soft background",
             "(b) Hard-soft background contrast"))

p_background <- ggplot(background_long, aes(
  x = time_min, y = rate, colour = series, shape = series,
  group = interaction(series, segment)
)) +
  geom_line(linewidth = 0.35, alpha = 0.65) +
  geom_point(size = 0.40, alpha = 0.75) +
  facet_grid(panel ~ ., scales = "free_y") +
  scale_colour_manual(values = c(
    Soft = "#2878B5", Hard = "#C35D35", Contrast = "#303030")) +
  scale_shape_manual(values = c(Soft = 16, Hard = 17, Contrast = 16)) +
  labs(title = "GRS 1915+105: SCORPEON background",
       subtitle = "NICER ObsID 2596012703 | bin width 1 s",
       x = "Elapsed time (minutes)",
       y = "Estimated background (counts/s)",
       colour = "Series", shape = "Series") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.spacing.y = grid::unit(0.7, "lines"))

print(p_background)
ggsave(file.path(supp_dir, "GRS_SCORPEON_HS_D.pdf"),
       p_background, width = 9.5, height = 6)


# =============================================================================
# 7. ORIGINAL GRS OBSERVATIONAL FIGURES
# =============================================================================
# The following plotting code is from the supplied figures_main_paper.R.
# Keep its original three-bin photon figure, H/S contrast figures, quiet/flare
# zooms and final main-paper three-panel plot.

# ============================================================
# AOAS | GRS 1915+105 | ORIGINAL OBSERVATION
# ObsID 2596012703, Delta = 1 second
# Original figure code, now run within the single-source application script.
# ============================================================

library(ggplot2)

grs_obsid <- "2596012703"
grs_events <- read.csv(
  file.path("data", paste0(grs_obsid, "_events.csv"))
)
grs_gti <- read.csv(
  file.path("data", paste0(grs_obsid, "_gti.csv"))
)
grs_gti <- grs_gti[order(grs_gti$start), ]

grs_t0 <- floor(min(grs_gti$start))
grs_soft <- c(2, 5)
grs_hard <- c(5, 12)

grs_events$energy <- grs_events$pi * 0.01
grs_events$band <- ifelse(
  grs_events$energy >= grs_soft[1] &
    grs_events$energy < grs_soft[2], "Soft",
  ifelse(grs_events$energy >= grs_hard[1] &
           grs_events$energy <= grs_hard[2],
         "Hard", NA_character_)
)

grs_photons <- grs_events[
  !is.na(grs_events$band),
  c("time", "energy", "band")
]

# Reconstruct the exact original grid and GTI exposure.
make_grs_bins <- function(dt = 1L) {
  
  n <- as.integer(ceiling(
    (ceiling(max(grs_gti$stop)) - grs_t0) / dt
  ))
  
  left <- grs_t0 + (seq_len(n) - 1L) * dt
  right <- left + dt
  exposed <- numeric(n)
  
  for (g in seq_len(nrow(grs_gti))) {
    first <- max(
      1L, floor((grs_gti$start[g] - grs_t0) / dt) + 1L
    )
    last <- min(
      n, ceiling((grs_gti$stop[g] - grs_t0) / dt)
    )
    if (first > last) next
    
    ii <- first:last
    exposed[ii] <- exposed[ii] +
      pmax(0, pmin(right[ii], grs_gti$stop[g]) -
             pmax(left[ii], grs_gti$start[g]))
  }
  
  exposed <- pmin(exposed, dt)
  
  bin_id <- floor((grs_photons$time - grs_t0) / dt) + 1L
  valid <- bin_id >= 1L & bin_id <= n
  
  H <- tabulate(
    as.integer(bin_id[valid &
                        grs_photons$band == "Hard"]),
    nbins = n
  )
  S <- tabulate(
    as.integer(bin_id[valid &
                        grs_photons$band == "Soft"]),
    nbins = n
  )
  
  stopifnot(sum(H[exposed == 0]) == 0)
  stopifnot(sum(S[exposed == 0]) == 0)
  
  H[exposed == 0] <- NA_integer_
  S[exposed == 0] <- NA_integer_
  
  data.frame(
    bin = seq_len(n),
    t_start = left,
    t_stop = right,
    exposure_s = exposed,
    expo = exposed / dt,
    H = H,
    S = S,
    D = H - S,
    time_min = (left - grs_t0) / 60
  )
}

grs_b1 <- make_grs_bins(1L)

# ============================================================
# FIGURE S1: PHOTON ARRIVALS WITH THREE 1-S BINS
# ============================================================

candidate <- seq_len(nrow(grs_b1) - 2L)

# Prefer a partial-exposure bin followed by two full bins.
good <- candidate[
  grs_b1$expo[candidate] >= 0.2 &
    grs_b1$expo[candidate] < 0.95 &
    grs_b1$expo[candidate + 1L] > 0.999 &
    grs_b1$expo[candidate + 2L] > 0.999
]

if (!length(good)) {
  good <- candidate[
    grs_b1$expo[candidate] > 0.999 &
      grs_b1$expo[candidate + 1L] > 0.999 &
      grs_b1$expo[candidate + 2L] > 0.999
  ]
}
stopifnot(length(good) > 0)

# Find an illustrative window with both photon bands.
n_soft <- vapply(good, function(i)
  sum(grs_b1$S[i:(i + 2L)]), numeric(1))
n_hard <- vapply(good, function(i)
  sum(grs_b1$H[i:(i + 2L)]), numeric(1))

eligible <- which(n_soft >= 5 & n_hard >= 5)
if (!length(eligible)) eligible <- seq_along(good)

# Avoid selecting an extremely crowded 3-s interval.
target <- 120
score <- abs(
  n_soft[eligible] + n_hard[eligible] - target
)
i <- good[eligible[which.min(score)]]

grs_selected <- grs_b1[i:(i + 2L), ]
start <- grs_selected$t_start[1]
end <- tail(grs_selected$t_stop, 1)

grs_ph <- grs_photons[
  grs_photons$time >= start &
    grs_photons$time < end, ,
  drop = FALSE
]
grs_ph$x <- grs_ph$time - start
grs_ph$band <- factor(
  grs_ph$band, levels = c("Soft", "Hard")
)

# Shade only the actual exposed portions.
grs_rectangles <- do.call(
  rbind,
  lapply(i:(i + 2L), function(j) {
    a <- pmax(grs_b1$t_start[j], grs_gti$start)
    b <- pmin(grs_b1$t_stop[j], grs_gti$stop)
    ok <- b > a
    data.frame(
      xmin = a[ok] - start,
      xmax = b[ok] - start
    )
  })
)

grs_labels <- data.frame(
  x = c(0.5, 1.5, 2.5),
  y = 13.5,
  label = sprintf(
    "S=%d H=%d\ne=%.2f",
    grs_selected$S,
    grs_selected$H,
    grs_selected$expo
  )
)

p_GRS_photons <- ggplot() +
  geom_rect(
    data = grs_rectangles,
    aes(xmin = xmin, xmax = xmax,
        ymin = 0, ymax = 14.5),
    fill = "grey83", alpha = 0.65
  ) +
  geom_vline(
    xintercept = 0:3,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 5,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  geom_segment(
    data = grs_ph,
    aes(x = x, xend = x,
        y = 0, yend = energy),
    colour = "grey55",
    linewidth = 0.18,
    alpha = 0.5
  ) +
  geom_point(
    data = grs_ph,
    aes(x = x, y = energy, shape = band),
    colour = "black",
    size = 1,
    alpha = 0.65
  ) +
  geom_text(
    data = grs_labels,
    aes(x = x, y = y, label = label),
    size = 3,
    lineheight = 1
  ) +
  scale_shape_manual(
    values = c(Soft = 16, Hard = 17)
  ) +
  scale_x_continuous(breaks = 0:3) +
  coord_cartesian(
    xlim = c(0, 3), ylim = c(0, 14.5)
  ) +
  labs(
    x = "Elapsed time (s)",
    y = "Photon energy (keV)",
    shape = "Energy band",
    title = "GRS 1915+105",
    subtitle = "Photon arrivals, 1-s bins and GTI exposure"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(p_GRS_photons)
print(grs_selected[
  , c("bin", "exposure_s", "expo", "H", "S", "D")
])

grs_output <- "paper_outputs/analysis_R1/web_supplement"
dir.create(
  grs_output, recursive = TRUE, showWarnings = FALSE
)

ggsave(
  file.path(grs_output, "GRS_photons_binning_1s.pdf"),
  p_GRS_photons,
  width = 9.5,
  height = 4.4
)

# ============================================================
# FIGURE S2: HARD + SOFT ABOVE, CONTRAST BELOW
# ============================================================

grs_lc <- grs_b1[
  grs_b1$expo > 0, ,
  drop = FALSE
]

grs_lc$segment <- cumsum(c(
  1L, as.integer(diff(grs_lc$bin) > 1L)
))

grs_hs <- rbind(
  data.frame(
    time_min = grs_lc$time_min,
    counts = grs_lc$S,
    band = "Soft",
    segment = grs_lc$segment
  ),
  data.frame(
    time_min = grs_lc$time_min,
    counts = grs_lc$H,
    band = "Hard",
    segment = grs_lc$segment
  )
)

grs_hs$panel <- "(a) Hard and soft counts"

grs_d <- data.frame(
  time_min = grs_lc$time_min,
  counts = grs_lc$D,
  band = "Contrast",
  segment = grs_lc$segment,
  panel = "(b) Hard-soft contrast"
)

grs_plot <- rbind(grs_hs, grs_d)

grs_plot$panel <- factor(
  grs_plot$panel,
  levels = c(
    "(a) Hard and soft counts",
    "(b) Hard-soft contrast"
  )
)
grs_plot$band <- factor(
  grs_plot$band,
  levels = c("Soft", "Hard", "Contrast")
)

p_GRS_HSD <- ggplot(
  grs_plot,
  aes(
    x = time_min,
    y = counts,
    colour = band,
    shape = band,
    group = interaction(band, segment)
  )
) +
  geom_line(linewidth = 0.25, alpha = 0.4) +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_grid(panel ~ ., scales = "free_y") +
  scale_colour_manual(values = c(
    Soft = "#2878B5",
    Hard = "#C35D35",
    Contrast = "#303030"
  )) +
  scale_shape_manual(values = c(
    Soft = 16,
    Hard = 17,
    Contrast = 16
  )) +
  labs(
    x = "Elapsed time (minutes)",
    y = "Photon counts per 1-s bin",
    colour = "Series",
    shape = "Series"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey93"),
    panel.grid.minor = element_blank(),
    panel.spacing.y = grid::unit(0.7, "lines")
  )

print(p_GRS_HSD)

ggsave(
  file.path(grs_output, "GRS_HS_D_two_panels.pdf"),
  p_GRS_HSD,
  width = 10,
  height = 6
)



# ============================================================
# GRS 1915+105 | QUIET AND FLARING ZOOMS
# Original ObsID 2596012703 | Delta = 1 s
# ============================================================

library(ggplot2)

zoom_seconds <- 120L
min_block_bins <- 60L

# ------------------------------------------------------------
# 1. Identify fully exposed bins belonging to each GTI
# ------------------------------------------------------------

b <- grs_b1
b$total <- b$H + b$S

midpoint <- (b$t_start + b$t_stop) / 2
b$gti_id <- findInterval(midpoint, grs_gti$start)

valid <- which(
  b$expo >= 0.999 &
    is.finite(b$total) &
    b$gti_id >= 1L
)

valid <- valid[
  midpoint[valid] < grs_gti$stop[b$gti_id[valid]]
]

blocks <- split(valid, b$gti_id[valid])

# Only sufficiently long, continuous observing blocks
blocks <- blocks[vapply(blocks, function(ii) {
  length(ii) >= min_block_bins &&
    all(diff(ii) == 1L)
}, logical(1))]

stopifnot(length(blocks) >= 2L)

# ------------------------------------------------------------
# 2. Summarise variability in each observing block
# ------------------------------------------------------------

block_stats <- do.call(rbind, lapply(names(blocks), function(id) {
  
  ii <- blocks[[id]]
  y <- b$total[ii]
  
  data.frame(
    GTI = as.integer(id),
    duration_s = length(ii),
    mean_total = mean(y),
    median_total = median(y),
    sd_total = sd(y),
    Fano = var(y) / max(mean(y), 1e-8),
    q99_total = as.numeric(quantile(y, 0.99)),
    max_total = max(y)
  )
}))

block_stats <- block_stats[
  order(block_stats$GTI), ,
  drop = FALSE
]

cat("\nGRS: OBSERVING BLOCK DIAGNOSTICS\n")
print(block_stats, row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# 3. Select two distinct GTIs
#
# Quiet: lowest count overdispersion
# Flaring candidate: highest count overdispersion
#
# Override manually after inspecting block_stats.
# ------------------------------------------------------------

quiet_gti_manual <- NA_integer_
flaring_gti_manual <- NA_integer_

quiet_gti <- if (is.na(quiet_gti_manual)) {
  block_stats$GTI[which.min(block_stats$Fano)]
} else {
  quiet_gti_manual
}

flaring_gti <- if (is.na(flaring_gti_manual)) {
  candidates <- block_stats[
    block_stats$GTI != quiet_gti, ,
    drop = FALSE
  ]
  candidates$GTI[which.max(candidates$Fano)]
} else {
  flaring_gti_manual
}

stopifnot(
  quiet_gti != flaring_gti,
  as.character(quiet_gti) %in% names(blocks),
  as.character(flaring_gti) %in% names(blocks)
)

cat("\nSelected QUIET GTI:", quiet_gti, "\n")
cat("Selected FLARING-CANDIDATE GTI:", flaring_gti, "\n")

# ------------------------------------------------------------
# 4. Choose a representative window within each GTI
# ------------------------------------------------------------

choose_zoom <- function(gti_id, mode = c("quiet", "flaring")) {
  
  mode <- match.arg(mode)
  ii <- blocks[[as.character(gti_id)]]
  y <- b$total[ii]
  
  n <- length(ii)
  k <- min(zoom_seconds, n)
  
  starts <- seq_len(n - k + 1L)
  
  if (mode == "quiet") {
    
    # Rolling variance / mean: find a locally stable window.
    cs <- c(0, cumsum(y))
    cs2 <- c(0, cumsum(y^2))
    
    roll_sum <- cs[starts + k] - cs[starts]
    roll_sum2 <- cs2[starts + k] - cs2[starts]
    
    roll_mean <- roll_sum / k
    roll_var <- pmax(
      0,
      (roll_sum2 - roll_sum^2 / k) / (k - 1)
    )
    
    roll_fano <- roll_var / pmax(roll_mean, 1e-8)
    start_pos <- starts[which.min(roll_fano)]
    
  } else {
    
    # Include the highest observed count excursion
    # in the flaring-candidate GTI.
    peak_pos <- which.max(y)
    
    start_pos <- max(
      1L,
      min(n - k + 1L, peak_pos - floor(k / 2))
    )
  }
  
  ii[start_pos:(start_pos + k - 1L)]
}

quiet_idx <- choose_zoom(quiet_gti, "quiet")
flare_idx <- choose_zoom(flaring_gti, "flaring")

# ------------------------------------------------------------
# 5. Common plotting function
# ------------------------------------------------------------

plot_grs_zoom <- function(idx, label, gti_id) {
  
  z <- b[idx, , drop = FALSE]
  
  # Time relative to start of the selected window
  z$seconds <- z$t_start - z$t_start[1]
  
  hs <- rbind(
    data.frame(
      seconds = z$seconds,
      counts = z$S,
      series = "Soft",
      panel = "(a) Hard and soft counts"
    ),
    data.frame(
      seconds = z$seconds,
      counts = z$H,
      series = "Hard",
      panel = "(a) Hard and soft counts"
    )
  )
  
  dd <- data.frame(
    seconds = z$seconds,
    counts = z$D,
    series = "Contrast",
    panel = "(b) Hard-soft contrast"
  )
  
  plot_data <- rbind(hs, dd)
  
  plot_data$series <- factor(
    plot_data$series,
    levels = c("Soft", "Hard", "Contrast")
  )
  
  plot_data$panel <- factor(
    plot_data$panel,
    levels = c(
      "(a) Hard and soft counts",
      "(b) Hard-soft contrast"
    )
  )
  
  p <- ggplot(
    plot_data,
    aes(
      x = seconds,
      y = counts,
      colour = series,
      shape = series,
      group = series
    )
  ) +
    geom_line(linewidth = 0.35, alpha = 0.65) +
    geom_point(size = 0.9, alpha = 0.85) +
    facet_grid(panel ~ ., scales = "free_y") +
    scale_colour_manual(values = c(
      Soft = "#2878B5",
      Hard = "#C35D35",
      Contrast = "#303030"
    )) +
    scale_shape_manual(values = c(
      Soft = 16,
      Hard = 17,
      Contrast = 16
    )) +
    labs(
      title = paste0("GRS 1915+105: ", label),
      subtitle = paste0(
        "GTI ", gti_id,
        " | ", nrow(z), " s | ",
        "starting at elapsed minute ",
        round(z$time_min[1], 2)
      ),
      x = "Time within selected window (s)",
      y = "Photon counts per 1-s bin",
      colour = "Series",
      shape = "Series"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.background = element_rect(fill = "grey93"),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.6, "lines")
    )
  
  p
}

p_quiet <- plot_grs_zoom(
  quiet_idx, "Quiet interval", quiet_gti
)

p_flaring <- plot_grs_zoom(
  flare_idx, "High-variability interval", flaring_gti
)

print(p_quiet)
print(p_flaring)

# ------------------------------------------------------------
# 6. Save the two figures for the Web Supplement
# ------------------------------------------------------------

dir.create(
  grs_output,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(grs_output, "GRS_quiet_zoom.pdf"),
  p_quiet,
  width = 9,
  height = 5.5
)

ggsave(
  file.path(grs_output, "GRS_flaring_zoom.pdf"),
  p_flaring,
  width = 9,
  height = 5.5
)

# Report the selected windows for reproducibility
zoom_summary <- data.frame(
  window = c("Quiet", "Flaring candidate"),
  GTI = c(quiet_gti, flaring_gti),
  start_elapsed_min = c(
    b$time_min[quiet_idx[1]],
    b$time_min[flare_idx[1]]
  ),
  duration_s = c(length(quiet_idx), length(flare_idx)),
  mean_total_counts = c(
    mean(b$total[quiet_idx]),
    mean(b$total[flare_idx])
  ),
  Fano = c(
    var(b$total[quiet_idx]) / mean(b$total[quiet_idx]),
    var(b$total[flare_idx]) / mean(b$total[flare_idx])
  )
)

print(zoom_summary, row.names = FALSE)

# ============================================================
# GRS 1915+105 | FINAL SECTION 2 FIGURE
# (a) First four observing blocks, excluding >60 total counts
# (b) Flaring observing block
# (c) First post-flare observing block
# Only H_t and S_t | Delta = 1 s
# ============================================================

library(ggplot2)
library(cowplot)

out_dir <- "paper_outputs/analysis_R1/main_paper"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Adjustable graphical settings
threshold_total <- 60
major_gap_seconds <- 15 * 60
n_early_blocks <- 4L

# ------------------------------------------------------------
# 1. Observed bins and major observing blocks
# ------------------------------------------------------------

z <- grs_b1
z$total <- z$H + z$S

z <- z[
  z$expo > 0 &
    is.finite(z$total),
  , drop = FALSE
]

z <- z[order(z$bin), , drop = FALSE]

# Group nearby GTIs into major observing episodes.
# Short internal gaps remain visible, but a separation
# exceeding 15 minutes starts a new major block.

z$major_block <- cumsum(c(
  1L,
  as.integer(diff(z$t_start) > major_gap_seconds)
))

block_ids <- sort(unique(z$major_block))

block_info <- do.call(rbind, lapply(block_ids, function(k) {
  
  x <- z[z$major_block == k, , drop = FALSE]
  
  data.frame(
    block = k,
    start_min = min(x$time_min),
    end_min = max(x$time_min),
    observed_bins = nrow(x),
    max_total_counts = max(x$total),
    bins_above_60 = sum(x$total > threshold_total)
  )
}))

cat("\nMAJOR OBSERVING BLOCKS\n")
print(block_info, row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# 2. Flaring block and the following block
# ------------------------------------------------------------

peak_row <- which.max(z$total)
flare_block <- z$major_block[peak_row]

before_flare <- block_ids[block_ids < flare_block]
after_flare <- block_ids[block_ids > flare_block]

stopifnot(
  length(before_flare) >= n_early_blocks,
  length(after_flare) >= 1L
)

early_blocks <- head(before_flare, n_early_blocks)
post_block <- after_flare[1L]

cat("\nSelected early blocks:",
    paste(early_blocks, collapse = ", "), "\n")
cat("Flaring block:", flare_block, "\n")
cat("Post-flare block:", post_block, "\n")

# ------------------------------------------------------------
# 3. Extract the three portions
# ------------------------------------------------------------

# FIRST FOUR BLOCKS
# Exclude high-total-count bins only from the
# descriptive main panel, not from the analysis.

early <- z[
  z$major_block %in% early_blocks &
    z$total <= threshold_total,
  , drop = FALSE
]

# FULL FLARING BLOCK, including rise and decay
flare <- z[
  z$major_block == flare_block,
  , drop = FALSE
]

# FIRST BLOCK AFTER THE FLARE
post <- z[
  z$major_block == post_block,
  , drop = FALSE
]

stopifnot(
  nrow(early) > 1L,
  nrow(flare) > 1L,
  nrow(post) > 1L
)

# Preserve absolute elapsed time in the main panel.
early$x <- early$time_min

# Local time for the two zooms.
flare$x <- (flare$t_start - min(flare$t_start)) / 60
post$x <- (post$t_start - min(post$t_start)) / 60

# ------------------------------------------------------------
# 4. Convert H/S into long plotting format
# ------------------------------------------------------------

prepare_hs <- function(dat) {
  
  # Restart lines after gaps OR removed bins.
  dat$segment <- cumsum(c(
    1L,
    as.integer(diff(dat$bin) > 1L)
  ))
  
  ans <- rbind(
    data.frame(
      x = dat$x,
      count = dat$S,
      band = "Soft",
      segment = dat$segment
    ),
    data.frame(
      x = dat$x,
      count = dat$H,
      band = "Hard",
      segment = dat$segment
    )
  )
  
  ans$band <- factor(ans$band, levels = c("Soft", "Hard"))
  ans
}

early_hs <- prepare_hs(early)
flare_hs <- prepare_hs(flare)
post_hs <- prepare_hs(post)

# ------------------------------------------------------------
# 5. Common plotting style
# ------------------------------------------------------------

colours <- c(
  Soft = "#2878B5",
  Hard = "#C35D35"
)

shapes <- c(
  Soft = 16,
  Hard = 17
)

make_hs_plot <- function(dat, title, xlabel,
                         point_size = 0.55) {
  
  ggplot(
    dat,
    aes(
      x = x,
      y = count,
      colour = band,
      shape = band,
      group = interaction(band, segment)
    )
  ) +
    geom_line(
      linewidth = 0.30,
      alpha = 0.55
    ) +
    geom_point(
      size = point_size,
      alpha = 0.80
    ) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = shapes) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.07))
    ) +
    labs(
      title = title,
      x = xlabel,
      y = "Counts per 1-s bin",
      colour = "Energy band",
      shape = "Energy band"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 6, 8)
    )
}

# ------------------------------------------------------------
# 6. MAIN: first four observing blocks
# ------------------------------------------------------------

p_early <- make_hs_plot(
  early_hs,
  "(a) First four observing blocks",
  "Elapsed time (minutes)",
  point_size = 0.55
)

# Vertical scale follows ONLY the first four blocks
# after excluding bins with H+S > 60.
early_ymax <- max(early$H, early$S)

p_early <- p_early +
  coord_cartesian(
    ylim = c(0, max(5, early_ymax * 1.08))
  )

# ------------------------------------------------------------
# 7. ZOOM: flare and post-flare
# ------------------------------------------------------------

p_flare <- make_hs_plot(
  flare_hs,
  "(b) Flaring block",
  "Time within block (minutes)",
  point_size = 0.43
)

p_post <- make_hs_plot(
  post_hs,
  "(c) Post-flare block",
  "Time within block (minutes)",
  point_size = 0.60
)

# Each zoom has its own vertical scale.

# ------------------------------------------------------------
# 8. Shared legend
# ------------------------------------------------------------

legend_plot <- make_hs_plot(
  early_hs,
  "",
  "",
  point_size = 1
) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

legend <- cowplot::get_legend(legend_plot)

# ------------------------------------------------------------
# 9. Compose: one wide main + two zooms below
# ------------------------------------------------------------

zoom_row <- cowplot::plot_grid(
  p_flare,
  p_post,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "hv",
  axis = "tblr"
)

figure_body <- cowplot::plot_grid(
  p_early,
  zoom_row,
  ncol = 1,
  rel_heights = c(1.05, 1),
  align = "v",
  axis = "lr"
)

p_GRS_final <- cowplot::plot_grid(
  legend,
  figure_body,
  ncol = 1,
  rel_heights = c(0.10, 1)
)

print(p_GRS_final)

# ------------------------------------------------------------
# 10. Save final figure
# ------------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    "GRS_HS_first4_flare_postflare.pdf"
  ),
  p_GRS_final,
  width = 10.5,
  height = 7.2
)

cat("\nFINAL FIGURE SUMMARY\n")
cat("Early blocks:", paste(early_blocks, collapse = ", "), "\n")
cat("Flaring block:", flare_block, "\n")
cat("Post-flare block:", post_block, "\n")
cat("Early bins displayed:", nrow(early), "\n")
cat("Early bins excluded above 60:",
    sum(z$major_block %in% early_blocks &
          z$total > threshold_total), "\n")
cat("Flare bins:", nrow(flare), "\n")
cat("Post-flare bins:", nrow(post), "\n")
cat("Maximum H/S in main:", early_ymax, "\n")


# =============================================================================
# END OF COMPLETE GRS 1915+105 REPRODUCTION
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "GRS_R_session_info.txt"))
cat("\nReproduction finished. Fit container: paper_outputs/GRS.rds\n")
cat("Main outputs:", main_dir, "\n")
cat("Supplement outputs:", supp_dir, "\n")
