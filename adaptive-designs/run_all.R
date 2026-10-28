# =============================================================================
# Run all simulations for the adaptive-designs post.
# Everything is cached in sims/, so re-running only computes what changed.
# Afterwards, render adaptive-designs.qmd (it only reads results/).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})
stopifnot(requireNamespace("adaptr", quietly = TRUE))

# All paths are relative to this topic's folder (found via here::here(),
# which uses the project root of the repository)
topic_dir   <- here::here("adaptive-designs")
dir_sims    <- file.path(topic_dir, "sims/")     # raw simulations and calibrations (large; not in git)
dir_results <- file.path(topic_dir, "results/")  # small summary objects read by the post (in git)
walk(c(dir_sims, dir_results), dir.create, showWarnings = FALSE, recursive = TRUE)

source(file.path(topic_dir, "R", "helpers.R"))



# =============================================================================
# Settings
# Draft vs. final: for final results use n_rep = 10000 and n_rep_calib = 5000
# (this takes hours the first time; everything is cached in sims/, so later
# runs only re-load).
# =============================================================================

settings <- list(
  n_rep         = 5000,   # simulated trials per design and scenario 
  n_rep_calib   = 5000,   # simulated trials per calibration step  
  n_rep_history = 500,    # full-detail trials for allocation-over-time plots
  n_boot        = 1000,   # bootstrap samples for Monte Carlo 95% CIs used in adaptr::check_performance
  target        = 0.10,   # false-go rate the designs are calibrated to
  n_cores       = 4,      # cores for parallel simulation
  seed          = 4131,
  include_two_arm = FALSE # also run the two-arm example (control + 1 arm, own calibration)
)

# Scenarios: id, title, design family, true response rates (first arm = control)
scenarios <- list(
  pick_the_winner = list(
    title = "Pick the winner (E = 40%, others 15%)", family = "six_arm",
    ys = c(Control = 0.15, A = 0.15, B = 0.15, C = 0.15, D = 0.15, E = 0.40)),
  futility_screen = list(
    title = "Futility screen (all arms 15%)", family = "six_arm",
    ys = c(Control = 0.15, A = 0.15, B = 0.15, C = 0.15, D = 0.15, E = 0.15)),
  close_effects = list(
    title = "Close effects (E = 25%, D = 20%, others 15%)", family = "six_arm",
    ys = c(Control = 0.15, A = 0.15, B = 0.15, C = 0.15, D = 0.20, E = 0.25)),
  two_arm = list(
    title = "Two arms (E = 40% vs. control 20%)", family = "two_arm",
    ys = c(Control = 0.20, E = 0.40))
)
# The two-arm example is optional; its family (and calibration) is only used if included
if (!settings$include_two_arm) scenarios$two_arm <- NULL


adaptr::setup_cluster(settings$n_cores)


# =============================================================================
# 01: Calibrate the superiority threshold of every design in every family
#     to the target false-go rate under the null scenario.
#     Output: results/calibration.rds
# =============================================================================

families <- unique(map_chr(scenarios, "family"))

# Each family: the three main designs (futility vs. the original control only)
# plus the two adaptive designs with futility also vs. a new control (FNC)
calibration <- map(set_names(families), \(fam) {
  message("Calibrating family: ", fam)
  cals   <- calibrate_family(fam, target = settings$target,
                             n_rep = settings$n_rep_calib, seed = settings$seed)
  cals_v <- calibrate_family(fam, target = settings$target,
                             n_rep = settings$n_rep_calib, seed = settings$seed,
                             designs = design_names[-1], futility_only_first = FALSE)
  all <- c(cals, cals_v)
  list(thresholds = get_thresholds(all),
       table      = threshold_table(all),
       steps      = calibration_data(all))
})

calibration$settings <- settings
saveRDS(calibration, paste0(dir_results, "calibration.rds"))
message("Saved results/calibration.rds")



# =============================================================================
# 02: Simulate all scenarios with the calibrated thresholds and save small
#     summary objects (no simulations inside).
#     Output: results/scenario_<id>.rds (five designs each)
# =============================================================================


#calibration <- readRDS(paste0(dir_results, "calibration.rds"))

iwalk(scenarios, \(sc, id) {
  message("Scenario: ", id)
  thr  <- calibration[[sc$family]]$thresholds
  sims <- c(run_designs(id, sc$ys, sc$family, thr, n_rep = settings$n_rep, seed = settings$seed),
            run_designs(id, sc$ys, sc$family, thr, n_rep = settings$n_rep, seed = settings$seed,
                        designs = design_names[-1], futility_only_first = FALSE))
  
  # Monte Carlo uncertainty (bootstrap CIs) for the main example only
  S <- summarise_designs(sims, sc$ys, sc$family, id = id, thresholds = thr,
                         uncertainty = id == "pick_the_winner",
                         n_boot = settings$n_boot, seed = settings$seed)
  S$title <- sc$title
  S$extra <- list()
  
  # Main example and null scenario: convergence of key metrics
  if (id == "pick_the_winner") {
    S$extra$convergence <- map(sims, convergence_data,
                               metrics = c("size_mean", "prob_select_arm_E", "ratio_ys_mean"))
    S$extra$status <- list(
      `Stop/drop` = status_data(sims[["Stop/drop"]], arms = c("Control", "A", "E")),
      `Stop/drop + RAR`                  = status_data(sims[["Stop/drop + RAR"]], arms = c("Control", "A", "E")))
  }
  if (id == "futility_screen") {
    S$extra$convergence <- map(sims, convergence_data,
                               metrics = c("size_mean", "prob_select_arm_Control"))
  }
  
  saveRDS(S, paste0(dir_results, "scenario_", id, ".rds"))
})

message("Saved scenario results")

# =============================================================================
# 03: Single simulated trials, allocation over time and allocation tables
#     (pick-the-winner scenario).
#     Output: results/single_trials.rds
# =============================================================================

#calibration <- readRDS(paste0(dir_results, "calibration.rds"))
sc  <- scenarios$pick_the_winner
s   <- family_settings[[sc$family]]
thr <- calibration[[sc$family]]$thresholds

spec_b   <- make_design("Stop/drop", sc$ys, s,
                        superiority = thr[["Stop/drop"]])
spec_c   <- make_design("Stop/drop + RAR", sc$ys, s, superiority = thr[["Stop/drop + RAR"]])
spec_c_v <- make_design("Stop/drop + RAR", sc$ys, s, superiority = thr[["Stop/drop + RAR (FNC)"]],
                        futility_only_first = FALSE)
useless  <- c("A", "B", "C", "D")

# Conditions for choosing illustrative trials (first seed that meets each)
cond_typical <- \(h) any(h$status == "superior" & h$arm == "E") && max(h$n) < s$max_n
cond_promotion <- \(h) {
  p <- promotions(h)
  nrow(p) > 0 && p$arm[1] == "E" && p$look[1] <= max(h$look) - 2 &&
    sum(h$look == p$look[1] & h$arm %in% useless & h$status == "active") >= 2
}
cond_false <- \(h) { p <- promotions(h); nrow(p) > 0 && p$arm[1] %in% useless }

seeds <- c(typical   = find_seed(spec_c, cond_typical),
           promotion = find_seed(spec_c, cond_promotion),
           false_win = find_seed(spec_c, cond_false))

lab_b  <- "(b) Stop/drop"
lab_c  <- "(c) Stop/drop + RAR"
lab_cv <- "(c) Stop/drop + RAR (FNC)"

single <- list(
  seeds = seeds,
  typical = list(
    hists = list(trial_history(spec_b, seeds[["typical"]]),
                 trial_history(spec_c, seeds[["typical"]])) |> set_names(c(lab_b, lab_c)),
    thr = c(thr[["Stop/drop"]], thr[["Stop/drop + RAR"]]) |> set_names(c(lab_b, lab_c))),
  promotion = list(
    hists = list(trial_history(spec_c, seeds[["promotion"]]),
                 trial_history(spec_c_v, seeds[["promotion"]])) |> set_names(c(lab_c, lab_cv)),
    thr = c(thr[["Stop/drop + RAR"]], thr[["Stop/drop + RAR (FNC)"]]) |> set_names(c(lab_c, lab_cv))),
  false_win = list(
    hists = list(trial_history(spec_c, seeds[["false_win"]])) |> set_names(lab_c),
    thr = thr[["Stop/drop + RAR"]] |> set_names(lab_c)),
  max_n = s$max_n,
  alloc_fixed = alloc_table_fixed(),
  alloc_rar = alloc_table_rar()
)

# Allocation over time across many trials (needs full-detail simulations)
sims_hist <- run_designs("history", sc$ys, sc$family, thr, n_rep = settings$n_rep_history,
                         seed = settings$seed,
                         designs = c("Stop/drop", "Stop/drop + RAR"),
                         sparse = FALSE)
single$history <- map(sims_hist, history_data)

saveRDS(single, paste0(dir_results, "single_trials.rds"))
message("Saved results/single_trials.rds")

