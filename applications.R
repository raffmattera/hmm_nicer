# =============================================================================
# AOAS R1 | REAL-DATA APPLICATIONS: OBSERVATIONS -> HMM -> SCIENTIFIC ANALYSES
# =============================================================================
# RUN FROM THE REPOSITORY ROOT:
#   source("run_analysis.R")
#
# Section 1. Descriptive observations: GRS (1 s), J1023 (10 s) ONLY.
# Section 2. Fit missing source containers here; interpret selected regimes for
#            two applications. K sensitivity uses every available dt and family.
# Section 3. Gap-aware versus compressed: GRS and J1023 ONLY.
# Section 4. Bayesian Blocks: GRS and J1023 ONLY.
# Section 5. Instrumental-background robustness: GRS, J1023, J2214 (10 s).
#
# This is the main public R entry point. It calls the three fit-only source
# scripts directly, but only if a corresponding fitted RDS is missing.
# No fitting is performed if the required containers already exist.
# =============================================================================

# PART 0. SCIENTIFIC SETTINGS ==================================================
main_sources <- c("GRS", "J1023")
background_sources <- c("GRS", "J1023", "J2214")
analysis_dt <- c(GRS = 1L, J1023 = 10L, J2214 = 10L)
models <- c("skellam", "pois2", "bivpois")

# Bayesian Blocks JSON is generated separately using 05_bayesian_blocks.py.
# This R script only reads the four precomputed JSON files (no Python required).

# Source-specific observational settings (matching the source R1 scripts).
source_specs <- data.frame(
  source = c("GRS", "J1023", "J2214"),
  obsid = c("2596012703", "4034060103", "1070060101"),
  label = c("GRS 1915+105", "PSR J1023+0038", "J2214 (verify target designation)"),
  dataset = c("prepared", "prepared", "archive"),
  dt = c(1L, 10L, 10L),
  soft_lo = c(2, 0.5, 0.5), soft_hi = c(5, 2, 2),
  hard_lo = c(5, 2, 2), hard_hi = c(12, 5, 10),
  energy_tag = c("E2_5_E5_12", "E05_2_E2_5", "E05_2_E2_10"),
  stringsAsFactors = FALSE
)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install R package ggplot2")
suppressPackageStartupMessages(library(ggplot2))
if (!file.exists("hmm_functions.R")) stop("Run this script from the repository root")

output_dir <- file.path("paper_outputs", "analysis_R1")
main_dir <- file.path(output_dir, "main_paper")
supp_dir <- file.path(output_dir, "web_supplement")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. DESCRIPTIVE ANALYSIS: ONLY GRS AND J1023
# =============================================================================
# Input files are the processed events and GTIs. This is intentionally placed
# BEFORE model fitting to show what is actually observed.
observations <- list()

# 1A. A window of 30 actual photons, from ONE GTI.
# The vertical sticks show photon energy; they are not intensity curves.
# This is the only dedicated plotting function: it reproduces the user-selected
# event-window rule exactly and records which photons were displayed.
plot_R1_photons <- function(events, gtis, cfg, n_photons=30L) {
  energy <- events$pi * 0.01
  soft <- energy>=cfg$soft[1] & energy<cfg$soft[2]
  hard <- energy>=cfg$hard[1] & energy<=cfg$hard[2]
  keep <- soft | hard
  if (sum(keep)<n_photons) stop("Too few photons: ",cfg$label)

  ph <- data.frame(time=events$time[keep],energy=energy[keep],
                   band=ifelse(soft[keep],"Soft","Hard"))
  ph <- ph[order(ph$time), , drop=FALSE]
  ph$gti <- findInterval(ph$time,gtis$start)
  n <- nrow(ph)
  starts <- seq_len(n-n_photons+1L)
  ends <- starts+n_photons-1L
  cs <- c(0L,cumsum(ph$band=="Soft"))
  nsoft <- cs[ends+1L]-cs[starts]
  candidates <- starts[
    ph$gti[starts]==ph$gti[ends] &
    nsoft>=3L & (n_photons-nsoft)>=3L
  ]
  if (!length(candidates))
    stop("No single-GTI photon window with both bands: ",cfg$label)

  i <- candidates[which.min(abs(candidates-(n-n_photons+2L)/2))]
  z <- ph[i:(i+n_photons-1L), , drop=FALSE]
  z$time_s <- z$time-z$time[1L]
  z$band <- factor(z$band,levels=c("Soft","Hard"))

  p <- ggplot2::ggplot(z,ggplot2::aes(x=time_s,y=energy))+
    ggplot2::geom_segment(
      ggplot2::aes(xend=time_s,y=0,yend=energy),
      colour="grey77",linewidth=0.3
    )+
    ggplot2::geom_point(ggplot2::aes(shape=band),size=2.3,
                        colour="black")+
    ggplot2::geom_hline(yintercept=cfg$hard[1],linetype="dashed",
                        linewidth=0.4)+
    ggplot2::scale_shape_manual(values=c(Soft=16,Hard=17))+
    ggplot2::coord_cartesian(ylim=c(0,max(cfg$hard[2],z$energy)*1.08))+
    ggplot2::labs(
      title=cfg$label,
      subtitle=paste0("Observed NICER photons; soft ",cfg$soft[1],"–",
                      cfg$soft[2]," keV, hard ",cfg$hard[1],"–",
                      cfg$hard[2]," keV"),
      x="Time since first displayed photon (s)",
      y="Photon energy (keV)",shape="Band"
    )+
    ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(legend.position="top",
                   panel.grid.minor=ggplot2::element_blank())

  list(plot=p,info=data.frame(
    source=cfg$label, obsid=cfg$obsid,
    gti_index=z$gti[1L], photons=nrow(z),
    soft_photons=sum(z$band=="Soft"),
    hard_photons=sum(z$band=="Hard"),
    first_event_time=z$time[1L],last_event_time=z$time[nrow(z)],
    span_s=max(z$time_s)
  ))
}


for (src in main_sources) {
  spec <- source_specs[source_specs$source == src, , drop = FALSE]
  dt <- analysis_dt[[src]]
  event_file <- file.path("data", paste0(spec$obsid, "_events.csv"))
  gti_file <- file.path("data", paste0(spec$obsid, "_gti.csv"))
  if (!file.exists(event_file) || !file.exists(gti_file)) {
    stop("Missing processed input for ", src, ": ", event_file, " / ", gti_file)
  }
  events <- read.csv(event_file)
  gti <- read.csv(gti_file)
  stopifnot(all(c("time", "pi") %in% names(events)),
            all(c("start", "stop") %in% names(gti)))
  gti <- gti[order(gti$start), , drop = FALSE]
  if (!nrow(gti) || any(gti$stop <= gti$start) ||
      (nrow(gti) > 1L && any(gti$start[-1] < head(gti$stop, -1) - 1e-7))) {
    stop("Invalid GTIs for ", src)
  }
  gi <- findInterval(events$time, gti$start)
  inside <- gi >= 1L
  inside[inside] <- events$time[inside] < gti$stop[gi[inside]]
  if (!all(inside)) stop("Events outside GTI for ", src)

  # Replicate the exact bin boundaries, exposure fractions and band endpoints
  # used by the original R1 fitting scripts, including partially exposed bins.
  t0 <- floor(min(gti$start))
  t1 <- ceiling(max(gti$stop))
  n_bins <- as.integer(ceiling((t1-t0) / dt))
  starts <- t0 + (seq_len(n_bins)-1L)*dt
  stops <- starts + dt
  exposure <- numeric(n_bins)
  for (g in seq_len(nrow(gti))) {
    first <- max(1L, floor((gti$start[g]-t0)/dt)+1L)
    last <- min(n_bins, ceiling((gti$stop[g]-t0)/dt))
    if (first > last) next
    ii <- first:last
    exposure[ii] <- exposure[ii] + pmax(0, pmin(stops[ii],gti$stop[g]) -
                                       pmax(starts[ii],gti$start[g]))
  }
  exposure <- pmin(exposure,dt)
  if (abs(sum(exposure)-sum(gti$stop-gti$start)) > 1e-5)
    stop("Exposure mismatch for ", src)
  bin_id <- floor((events$time-t0)/dt)+1L
  keV <- events$pi*0.01
  in_grid <- bin_id>=1L & bin_id<=n_bins
  H <- tabulate(as.integer(bin_id[in_grid & keV>=spec$hard_lo &
                                            keV<=spec$hard_hi]),nbins=n_bins)
  S <- tabulate(as.integer(bin_id[in_grid & keV>=spec$soft_lo &
                                            keV<spec$soft_hi]),nbins=n_bins)
  if (sum(H[exposure==0])+sum(S[exposure==0]) != 0)
    stop("Counts outside exposure for ", src)
  H[exposure==0] <- NA_integer_
  S[exposure==0] <- NA_integer_
  bins <- data.frame(bin=seq_len(n_bins), t_start=starts,t_stop=stops,
    time_min=(starts-t0)/60,exposure_s=exposure,expo=exposure/dt,
    H=H,S=S,D=H-S,total=H+S,H_rate=H/exposure,S_rate=S/exposure,
    D_rate=(H-S)/exposure,total_rate=(H+S)/exposure)
  observations[[src]] <- list(spec=spec,events=events,gti=gti,bins=bins,dt=dt)
  cat("\nDESCRIPTIVE SOURCE:",src,"| dt=",dt,"| exposure=",sum(exposure),"s\n")

  # 1A. Real photon arrivals: soft circles and hard triangles, all black.
  photon_cfg <- list(
    label=as.character(spec$label), obsid=as.character(spec$obsid),
    soft=c(spec$soft_lo,spec$soft_hi), hard=c(spec$hard_lo,spec$hard_hi)
  )
  photon_result <- plot_R1_photons(events, gti, photon_cfg,
                                   n_photons=30L)
  p_events <- photon_result$plot
  print(p_events)
  print(photon_result$info,row.names=FALSE)
  ggplot2::ggsave(file.path(main_dir,paste0(src,"_photon_arrivals.pdf")),
                  p_events,width=8.5,height=4.6)

  # 1B. Three independent observed series: H_t, S_t, D_t = H_t - S_t.
  lc <- bins[bins$exposure_s>0,,drop=FALSE]
  lc$segment <- cumsum(c(1L,as.integer(diff(lc$bin)>1L)))
  p_H <- ggplot(lc,aes(time_min,H_rate,group=segment))+
    geom_line(linewidth=0.35)+geom_point(size=0.35)+theme_bw()+
    labs(x="Elapsed minutes",y="Hard counts/s",title=paste(src,"H(t)"))
  p_S <- ggplot(lc,aes(time_min,S_rate,group=segment))+
    geom_line(linewidth=0.35)+geom_point(size=0.35)+theme_bw()+
    labs(x="Elapsed minutes",y="Soft counts/s",title=paste(src,"S(t)"))
  p_D <- ggplot(lc,aes(time_min,D_rate,group=segment))+
    geom_line(linewidth=0.35)+geom_point(size=0.35)+theme_bw()+
    labs(x="Elapsed minutes",y="H - S counts/s",title=paste(src,"D(t)"))
  print(p_H);print(p_S);print(p_D)
  ggsave(file.path(main_dir,paste0(src,"_H_observed.pdf")),p_H,width=9,height=3)
  ggsave(file.path(main_dir,paste0(src,"_S_observed.pdf")),p_S,width=9,height=3)
  ggsave(file.path(main_dir,paste0(src,"_D_observed.pdf")),p_D,width=9,height=3)
}

# 1C. Gap pattern at the chosen resolution. It is NOT assumed that J1023 has
#     exactly one gap: the observed number/lengths are calculated from the GTIs.
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
cat("\nOBSERVATION GAPS (GRS/J1023):\n")
print(gap_geometry,row.names=FALSE,digits=4)

# 1D. Background and deadtime overview. This step is descriptive: it DOES NOT
#     subtract background counts or alter any HMM likelihood.
#
#     SCORPEON reports rates over its own time intervals. For each HMM bin:
#       (1) intersect the SCORPEON interval with the bin and the NICER GTIs;
#       (2) integrate hard/soft background rates over the shared exposure;
#       (3) keep bins with >=90% coverage and extrapolate <=10% uncovered time.
#     This is a single technical helper because the same measurement alignment
#     must be used for GRS, J1023 and J2214, without triplicating 100 lines.
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
# 2. FITTING WHEN NEEDED + SELECTED REGIME INTERPRETATION
# =============================================================================
# Fitting starts HERE, after the descriptive section. There is no separate
# run_all_R1.R. Each source-specific R script contains ONLY count construction,
# gap-aware/compressed fitting, BIC diagnostics, and parameter settings.
#
# The source files must run from the repository root. They are called with
# ordinary source() exactly as in the original run_all.R. Scientific settings
# are restored after the source calls to avoid variable-name conflicts.
# Fit algorithm/seed formulas are unchanged; existing RDS are not overwritten.
# The repository root must contain hmm_functions.R and the source scripts.
# These are the SAME source calls used by the original run_all.R.
for (f in c("hmm_functions.R", "01_J1023.R", "02_GRS1915.R", "03_J2214.R")) {
  if (!file.exists(f)) stop("Missing file in working directory ", getwd(), ": ", f)
}

dir.create("paper_outputs", recursive=TRUE, showWarnings=FALSE)

if (!file.exists("paper_outputs/J1023.rds")) {
  cat("\nFITTING J1023 (all dt, families and K)\n")
  source("01_J1023.R")
  saveRDS(container_J1023, "paper_outputs/J1023.rds")
}

if (!file.exists("paper_outputs/GRS.rds")) {
  cat("\nFITTING GRS (all dt, families and K)\n")
  source("02_GRS1915.R")
  saveRDS(container_GRS, "paper_outputs/GRS.rds")
}

if (!file.exists("paper_outputs/J2214.rds")) {
  cat("\nFITTING J2214 (10 s only, all families and K)\n")
  source("03_J2214.R")
  saveRDS(container_J2214, "paper_outputs/J2214.rds")
}

# Source-specific scripts set their own global analysis_dt and models.
# Restore the main-analysis settings before reading the saved results.
analysis_dt <- c(GRS=1L, J1023=10L, J2214=10L)
models <- c("skellam", "pois2", "bivpois")

GRS <- readRDS("paper_outputs/GRS.rds")
J1023 <- readRDS("paper_outputs/J1023.rds")
J2214 <- readRDS("paper_outputs/J2214.rds")
containers <- list(GRS=GRS,J1023=J1023,J2214=J2214)

# Minimal integrity check, without recomputing all bands or comparing tables.
# If processed inputs or energy bands change, archive and regenerate the RDS.
for (src in background_sources) {
  x <- containers[[src]]
  sp <- source_specs[source_specs$source==src,,drop=FALSE]
  if (as.character(x$obsid)!=as.character(sp$obsid) ||
      is.null(x$bins[[paste0("dt",analysis_dt[[src]])]]))
    stop("Wrong observation or missing selected dt in ",src," RDS")
}
source("hmm_functions.R") # Original likelihood, Viterbi, ARI and gap diagnostics.

# 2A. ONE sensitivity table across ALL dt for the two main sources.
# All emission families and gap-handling methods are included. This is the only
# section of 1-4 that reads fit grids at dt other than the selected dt.
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
cat("\nSUPPLEMENT: K SELECTION ACROSS dt, EMISSIONS, METHODS (GRS/J1023)\n")
print(K_sensitivity,row.names=FALSE)
write.csv(K_sensitivity,file.path(supp_dir,"K_sensitivity_GRS_J1023.csv"),row.names=FALSE)
cat("\nCOMPLETE BIC GRID (all dt and both methods; numbers within each emission only):\n")
print(BIC_grid,row.names=FALSE,digits=5)
write.csv(BIC_grid,file.path(supp_dir,"BIC_complete_grid_GRS_J1023.csv"),row.names=FALSE)

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
    if (is.null(fit)||!isTRUE(fit$converged)) {
      warning("Nonconverged selected fit: ",src," / ",model);next
    }
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
# 3. GAP-AWARE VERSUS COMPRESSED: ONLY GRS AND J1023 AT CHOSEN dt
# =============================================================================
# Real-world TV measures DISAGREEMENT between estimated transition behaviour,
# not improvement against an unknown true A as in Simulation 1.
# NLL below is IN-SAMPLE prequential; never label it held-out prediction.
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
    if(is.null(gap)||is.null(comp_ind)||is.null(comp_same)||
       !isTRUE(gap$converged)||!isTRUE(comp_ind$converged)||!isTRUE(comp_same$converged)) {
      warning("Comparison fit missing/nonconverged for ",src," / ",model);next
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
# 4. BAYESIAN BLOCKS: ONLY GRS AND J1023 AT CHOSEN dt
# =============================================================================
# The four lightweight JSONs are PRECOMPUTED with 05_bayesian_blocks.py:
#   python3 05_bayesian_blocks.py --source GRS --dt 1 --signal contrast --p0 0.05 --compressed
#   python3 05_bayesian_blocks.py --source GRS --dt 1 --signal total --p0 0.05 --compressed
#   python3 05_bayesian_blocks.py --source J1023 --dt 10 --signal contrast --p0 0.05 --compressed
#   python3 05_bayesian_blocks.py --source J1023 --dt 10 --signal total --p0 0.05 --compressed
# These commands are NOT invoked from R. No reticulate, Python, or Astropy
# dependency is required for the R analysis when the four JSONs are supplied.
#
# HMM state labels recur; BB block IDs identify unique temporal segments.
# Thus ARI between them is NOT an appropriate measure of changepoint agreement.
# We compare switch POSITIONS only, with one-to-one matching of changepoints.
# Only adjacent FULL-EXPOSURE bins belonging to the SAME GTI are eligible;
# compressed-series adjacency across observation gaps is never scored.
# Matching tolerance is +/- one ORIGINAL bin (1 s GRS, 10 s J1023).
# Rates and temporal precision/recall quantify descriptive agreement, not
# comparative accuracy against ground truth or a formally calibrated test.

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
if (length(missing_bb)) {
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
# 5. SCORPEON ROBUSTNESS: GRS 1 s, J1023 10 s, J2214 10 s
# =============================================================================
# HMMs have already been fitted to the observed counts, NOT background-corrected
# counts. SCORPEON is an external estimate. For each inferred state k, calculate
# posterior-weighted observed and background rates over the SAME well-covered
# bins. The nominal residual D_obs - D_bg is a diagnostic, NOT a refitted state.
#
# Useful scientific questions:
#   1. Is the observed difference in contrast between states larger than the
#      corresponding SCORPEON background difference?
#   2. Does any pair of states appear mostly background-associated?
#   3. How does this depend on emission family and source?
# These descriptive comparisons do NOT prove that a state is astrophysical.
background_summary_all <- instrument_summary_main
background_by_state <- data.frame()
background_pairwise <- data.frame()

for (src in background_sources) {
  sp <- source_specs[source_specs$source==src,,drop=FALSE]
  dt <- analysis_dt[[src]]
  x <- containers[[src]]
  b <- x$bins[[paste0("dt",dt)]]

  # The two main sources were aligned in Section 1; J2214 enters ONLY here.
  if (src=="J2214") {
    gti <- read.csv(file.path("data",paste0(sp$obsid,"_gti.csv")))
    gti <- gti[order(gti$start),,drop=FALSE]
    bg <- read_instrumental_on_grid(src,b,gti,sp,dt)
    background_data[[src]] <- bg
    background_summary_all <- rbind(background_summary_all,bg$summary)
  }
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
    if (is.null(fit)||!isTRUE(fit$converged))
      stop("Non-converged selected fit: ",src," / ",model)
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

cat("\nINSTRUMENTAL SUMMARY: GRS/J1023/J2214 (deadtime is a fraction, not %)\n")
print(background_summary_all,row.names=FALSE,digits=5)
cat("\nBACKGROUND PAIRWISE STATE CONTRASTS (all emission families)\n")
print(background_pairwise,row.names=FALSE,digits=4)

# Only final interpretive results are exported; no redundant graphics or
# time-series CSV dumps. In particular, J2214 K is already explicit in the
# state-level table, so no separate J2214 K-selection CSV is needed.
write.csv(background_summary_all,
          file.path(supp_dir,"background_instrument_all_sources.csv"),row.names=FALSE)
write.csv(background_by_state,
          file.path(supp_dir,"background_by_state_all_sources.csv"),row.names=FALSE)
write.csv(background_pairwise,
          file.path(supp_dir,"background_pairwise_all_sources.csv"),row.names=FALSE)

cat("\nAOAS R1 APPLICATION ANALYSIS COMPLETE\n")
cat("Main paper directory: ",main_dir,"\n",sep="")
cat("Web supplement directory: ",supp_dir,"\n",sep="")
