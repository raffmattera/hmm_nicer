# ============================================================================
# Simulation 2: emission structure and time-varying background
#
# Run from the replication directory:
#   Rscript simu2.R
# or from an interactive R session:
#   source("simu2.R")
#
# Every replication is saved immediately after completion as rep_###.rds.
# Re-running this launcher skips existing, structurally complete replications.
# To split the run across machines, use non-overlapping REPS ranges.
# ============================================================================

# -------------------------- COMPUTING OPTIONS -------------------------------
REPS <- 1:100                 # E.g. 1:30, 31:60, 61:100
MAX_ITER <- 100L
TOL <- 0.01
cores <- parallel::detectCores()
N_WORKERS <- if (is.na(cores)) 1L else max(1L, as.integer(cores) - 2L)

# -------------------------- SCIENTIFIC DESIGN -------------------------------
T_GRID <- c(1000L, 2000L)
K_TRUE_GRID <- c(2L, 3L)
DGPS <- c("both_vary", "fixed_shape")
MODELS <- c("skellam", "pois2", "bivpois")

BACKGROUND <- data.frame(
  name = c("none", "symmetric", "near_symmetric", "moderate_asymmetric"),
  b_H = c(0, 3, 3.25, 4),
  b_S = c(0, 3, 2.75, 2),
  stringsAsFactors = FALSE
)

HMM_FILE <- "hmm_functions.R"
INTERNAL_FILE <- "internal_simulations.R"
RESULT_DIR <- file.path("results", "SIM2")

# -------------------------- LOAD FUNCTIONS ----------------------------------
if (!file.exists(HMM_FILE) || !file.exists(INTERNAL_FILE)) {
  stop("Place hmm_functions.R and internal_simulations.R in the working directory.")
}
if (!requireNamespace("clue", quietly = TRUE)) {
  stop("Install the R package 'clue' before running the simulations.")
}
source(HMM_FILE)
source(INTERNAL_FILE)
stopifnot(exists("bic_hmm", mode = "function"))

CFG <- list(
  design = "AOAS_R1_SIM2",
  T_GRID = T_GRID,
  K_TRUE_GRID = K_TRUE_GRID,
  MODELS = MODELS,
  DGPS = DGPS,
  BACKGROUND = BACKGROUND,
  MAX_ITER = MAX_ITER,
  TOL = TOL,
  N_STARTS = 1L,
  K_RULE = "1:(K_true+2)",
  gap_design = "many_short",
  n_gaps = 12L,
  missing_target = 0.12,
  transition_matrices = "R0",
  emissions = "R1"
)

# Candidate state counts are 1:(K_true + 2), for every K_true.
expected_rows <- (
  length(T_GRID) * length(DGPS) * nrow(BACKGROUND) *
  length(MODELS) * sum(K_TRUE_GRID + 2L)
)

run_simulation_study(
  reps = REPS,
  cfg = CFG,
  result_dir = RESULT_DIR,
  study = "SIM2",
  runner = run_sim2_rep,
  expected_rows = expected_rows,
  n_workers = N_WORKERS,
  hmm_file = HMM_FILE,
  internal_file = INTERNAL_FILE
)
