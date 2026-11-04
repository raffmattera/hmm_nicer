# AOAS R1 | One-step Monte Carlo replication and analysis

This package includes the **four uploaded R scripts unchanged**:
`hmm_functions.R`, `internal_simulations.R`, `simu1.R` and `simu2.R`.
Only `simulations_analysis.R` is newly revised. The original EM, DGP, seed
signatures, simulation design, model grids and marginal posterior-MAP ARI
remain unchanged.

## Run everything from one entry point

Open R in the folder containing these five scripts and execute:

```r
source("simulations_analysis.R")
```

This **runs** `source("simu1.R")` followed by its analysis, then
`source("simu2.R")` followed by its analysis, and finally creates the
figures. The fitting launchers skip structurally complete `rep_###.rds`
checkpoints. New replications are saved after each completed replication.
**IMPORTANT:** with `REPS <- 1:100` in both launchers, this command will
compute all *missing* replications up to 100. For a test or splitting across
machines, edit `REPS` in **both** `simu1.R` and `simu2.R` before executing.
For example, `REPS <- 1L` is a one-replication smoke test. Use non-overlapping
ranges on different computers, and do not merge different source versions.

Requirements:

```r
install.packages(c("extraDistr", "matrixStats", "clue", "ggplot2"))
```

`ggplot2` is optional; missing it skips only the two PDF figures.
Set `STUDIES <- "SIM1"` or `STUDIES <- "SIM2"` at the top of
`simulations_analysis.R` to run just one study. Keep the working directory
at the package root so that `source()` resolves the correct files.

## Scientific design

- **SIM1**: temporal fragmentation, fixed oracle `K=K_true`, both gap-aware
  and compressed fits, paired on the same simulated dataset. Compare
  posterior-MAP ARI, state accuracy and transition-matrix TV/KL. Positive
  ARI difference (gap-aware minus compressed) favors gap-aware; positive
  transition-TV difference (compressed minus gap-aware) favors gap-aware.
- **SIM2**: three emission models, four background settings, gap-aware fits
  and candidate `K` grid `1:(K_true+2)`. Main results count correct, under-
  and over-selection, selected-state recovery and paired changes from
  the zero-background baseline.
- Only jointly converged paired fits enter the primary SIM1 comparison.
  In SIM2 the primary BIC results require a **finite, fully converged K grid**;
  finite yet not fully converged grids are retained only in the labelled
  supplementary sensitivity output. Neither filtering choice deletes the
  underlying raw checkpoints.
- Summary `_mean`, `_sd`, `_n`, `_ci_low` and `_ci_high` columns give the
  sample mean, empirical SD, number of finite eligible observations and a
  two-sided **95% Student-t Monte Carlo confidence interval** across
  replications. This is not a confidence interval for an individual HMM
  parameter or an exact binomial proportion interval.
- The EM algorithm uses posterior smoothing probabilities. The simulation
  *decoding* metrics use marginal MAP (`max.col(gamma)`), as agreed. Viterbi
  remains a distinct trajectory-visualisation procedure in applications.

## Outputs

The scripts create `results/SIM1/rep_###.rds` and
`results/SIM2/rep_###.rds`, along with their settings and source checksums.
The analysis writes:

`results/analysis_R1/main_paper/`
- `SIM1_paired_effects.csv`
- `SIM2_background_and_K_selection.csv`
- `SIM2_paired_background_effects.csv`
- `SIM1_gap_ARI_paired.pdf` (if ggplot2 installed)
- `SIM2_BIC_K_recovery.pdf` (if ggplot2 installed)

`results/analysis_R1/web_supplement/`
- `SIM1_fit_quality.csv`
- `SIM1_paired_replications.csv`
- `SIM1_model_recovery.csv`
- `SIM2_fit_quality.csv`
- `SIM2_BIC_selection_by_replication.csv`
- `SIM2_BIC_finite_grid_sensitivity.csv`
- `SIM2_selected_K_frequencies.csv` (if at least one strictly valid grid)
- `SIM2_paired_background_replications.csv`

`results/analysis_R1/` also contains `analysis_metadata.txt` and
`R_sessionInfo.txt`.

## Reproducibility and limitations

Each complete replication is saved atomically, and the launchers refuse to
reuse incompatible settings/source checksum files. Do not edit the four
scientific scripts after starting an experiment. Running the consolidated
entry point with a source-hash mismatch stops rather than silently reusing
incompatible checkpoints. Matching random seeds does not guarantee
bitwise-identical fitted solutions across different R/package versions.

The analysis was statically checked, but no full R execution test was
possible in the authoring runtime because R is not installed there.
