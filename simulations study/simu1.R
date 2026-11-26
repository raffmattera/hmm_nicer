# ============================================================================
# Run from the replication directory:
#   Rscript simu1.R
# or from an interactive R session:
#   source("simu1.R")
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
# These settings reproduce the selected R1 simulation conditions.
T_GRID <- c(1000L, 2000L)
K_TRUE_GRID <- c(2L, 3L)
DGPS <- c("both_vary", "fixed_shape")
MODELS <- c("skellam", "pois2", "bivpois")
MISSING_GRID <- c(0.12, 0.24, 0.36)
GAP_DESIGNS <- c("many_short", "few_long")

# Keep these paths relative to the working (replication) directory.
HMM_FILE <- "hmm_functions.R"
INTERNAL_FILE <- "internal_simulations.R"
RESULT_DIR <- file.path("results", "SIM1")

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
  design = "AOAS_R1_SIM1",
  T_GRID = T_GRID,
  K_TRUE_GRID = K_TRUE_GRID,
  MODELS = MODELS,
  DGPS = DGPS,
  MISSING_GRID = MISSING_GRID,
  GAP_DESIGNS = GAP_DESIGNS,
  MAX_ITER = MAX_ITER,
  TOL = TOL,
  N_STARTS = 1L,
  K_RULE = "K_true",
  transition_matrices = "R0",
  emissions = "R1"
)

expected_rows <- (
  length(T_GRID) * length(K_TRUE_GRID) * length(DGPS) *
  length(MISSING_GRID) * length(GAP_DESIGNS) *
  2L * length(MODELS)  # 2 temporal treatments
)

run_simulation_study(
  reps = REPS,
  cfg = CFG,
  result_dir = RESULT_DIR,
  study = "SIM1",
  runner = run_sim1_rep,
  expected_rows = expected_rows,
  n_workers = N_WORKERS,
  hmm_file = HMM_FILE,
  internal_file = INTERNAL_FILE
)
