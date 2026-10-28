# =============================================================================
# Helper functions: adaptive phase II designs with adaptr
#
# Conventions
# - Every function from the adaptr package is called as adaptr::<function>(),
#   so it is always visible what adaptr does and what is our own code.
# - Each function's documentation says whether it WRAPS an adaptr function
#   or is OWN CODE (possibly working on adaptr output).
# - Code that reads adaptr's internal object structure (rather than its
#   exported functions) is marked "relies on adaptr 1.5.0 object structure";
#   these are the places most likely to need changes for a new adaptr version.
#
# Help: the roxygen comments can be viewed as help pages with
#   docstring::docstring(make_design)
#


# =============================================================================
# 1. Constants
# =============================================================================

design_names <- c("Fixed", "Stop/drop", "Stop/drop + RAR")

design_cols <- c(
  "Fixed"                                                    = "#2a78d6",
  "Stop/drop"                                                = "#eb6834",
  "Stop/drop + RAR"                                          = "#1baf7a",
  "Stop/drop (FNC)"                                          = "#a8431a",
  "Stop/drop + RAR (FNC)"                                    = "#0f6b4a"
)

arm_cols <- c(Control = "#3d3d3d", A = "#2a78d6", B = "#eda100",
              C = "#e87ba4", D = "#8a63b8", E = "#008300")

outcome_cols <- c(
  "Winner declared early"                        = "#008300",
  "Winner declared at final analysis"            = "#86c286",
  "Winner kept; other arms dropped for futility" = "#2a78d6",
  "No winner by the maximum"                     = "#a6a6a6",
  "All experimental arms dropped (no-go)"        = "#eda100",
  "Stopped for equivalence"                      = "#6d6d6d"
)

status_cols <- c(Recruiting = "grey75", Inferiority = "#d6604d",
                 Futility = "#c2a35a", Equivalence = "#4393c3",
                 Superiority = "#4dac26")

#' Label of a design, with a suffix for the futility-vs-new-control variant
#'
#' OWN CODE.
#'
#' @param d design name, one of `design_names`.
#' @param futility_only_first logical; `FALSE` marks the variant in which
#'   futility is also checked against an arm promoted to control.
#' @return A character string.
design_label <- function(d, futility_only_first = TRUE) {
  if (futility_only_first) d else paste0(d, " (FNC)")
}

#' File-system-safe version of a design label ("/" is not allowed in file names)
#'
#' OWN CODE.
#' @param x character vector.
#' @return A character vector.
file_safe <- function(x) gsub("/", "-", x, fixed = TRUE)

#' Design families
#'
#' OWN CODE. Each family fixes the arms, maximum sample size, analysis
#' schedule and allocation rules. Scenarios (true response rates) are supplied
#' separately, so a family and its calibrated thresholds can be reused across
#' scenarios. The first arm is always the common control.
family_settings <- list(
  six_arm = list(
    label         = "Control + 5 experimental arms, max 300",
    arms          = c("Control", "A", "B", "C", "D", "E"),
    null_rate     = 0.15,
    max_n         = 300,
    data_looks    = seq(from = 60, to = 300, by = 30),
    lag           = 0,      # participants randomized but awaiting outcome (adaptr supports > 0)
    futility_diff = 0.05,
    fixed_alloc   = list(control_prob_fixed = "sqrt-based fixed"),
    rar           = list(control_prob_fixed = "sqrt-based",
                         min_probs = c(NA, rep(0.05, 5)),
                         rescale_probs = "limits",
                         soften_power = 0.7),
    search_range  = c(0.9, 0.9999)
  ),
  # optional example, used only if run_all.R sets include_two_arm = TRUE
  two_arm = list(
    label         = "Control + 1 experimental arm, max 160",
    arms          = c("Control", "E"),
    null_rate     = 0.20,
    max_n         = 160,
    data_looks    = seq(from = 40, to = 160, by = 10),
    lag           = 0,
    futility_diff = 0.05,
    fixed_alloc   = list(fixed_probs = c(0.5, 0.5)),
    rar           = list(min_probs = c(0.1, 0.1), soften_power = 0.5),
    search_range  = c(0.8, 0.9999)
  )
)

# =============================================================================
# 2. Trial specification
# =============================================================================

#' Trial specification for one design, family and scenario
#'
#' WRAPS `adaptr::setup_trial_binom()`. Translates our three designs into
#' adaptr arguments:
#' * Fixed: a single analysis at the maximum sample size, fixed allocation.
#' * Stop/drop: interim analyses (stopping and arm dropping), fixed allocation.
#' * Stop/drop + RAR: interim analyses and response-adaptive randomization.
#'
#' All designs compare each experimental arm with the current control:
#' * superiority, P(better than control) > x: the arm becomes the new control
#'   and the old control is dropped;
#' * inferiority, P(better than control) < 1 - x: the arm is dropped;
#' * futility, P(not better than control by `futility_diff`) > 0.90: the arm
#'   is dropped (checked against the original control only, unless
#'   `futility_only_first = FALSE`).
#' The trial stops when one arm is left or at the maximum sample size.
#'
#' @param design one of `design_names`.
#' @param ys true response rates, in the order of the family's arms.
#' @param settings one element of `family_settings`.
#' @param superiority superiority threshold x; inferiority is 1 - x.
#' @param futility_only_first logical; `TRUE` (default) checks futility only
#'   against the original control, `FALSE` also against promoted controls.
#' @return An adaptr `trial_spec` object.
make_design <- function(design, ys, settings, superiority = 0.99,
                        futility_only_first = TRUE) {
  stopifnot(design %in% design_names, length(ys) == length(settings$arms))
  max_n <- settings$max_n
  looks <- settings$data_looks
  common <- list(
    arms                = settings$arms,
    true_ys             = unname(ys),
    control             = settings$arms[1],
    highest_is_best     = TRUE,                  # response: higher is better
    superiority         = superiority,
    inferiority         = 1 - superiority,
    futility_prob       = 0.9,
    futility_diff       = settings$futility_diff,
    futility_only_first = futility_only_first,
    n_draws             = 5000
  )
  timing <- if (design == "Fixed") {
    list(data_looks = max_n)                       # one final analysis only
  } else {
    list(data_looks = looks,
         randomised_at_looks = pmin(looks + settings$lag, max_n))
  }
  alloc <- if (design == "Stop/drop + RAR") settings$rar else settings$fixed_alloc
  do.call(adaptr::setup_trial_binom, c(common, timing, alloc))
}

# =============================================================================
# 3. Calibration
# =============================================================================

#' Calibration function: false-go rate as a function of the threshold
#'
#' OWN CODE (not an adaptr function), written in the format that the `fun`
#' argument of `adaptr::calibrate_trial()` expects: a function of `x` (the
#' candidate threshold) and `trial_spec` that returns a list with `sims`,
#' `trial_spec` and `y` (the value to calibrate). `calibrate_trial()` only
#' passes `x` and `trial_spec`, so `n_rep`, `base_seed` and `cores` are fixed
#' here when the function is created (a "closure"), instead of being looked up
#' as global variables as in the template in `?adaptr::calibrate_trial`.
#'
#' Target (y): the false-go rate, P(an arm other than the original control is
#' selected), evaluated in the null scenario (all arms equal). With
#' `select_strategy = "control or best"`, a trial selects the original
#' control if it is still active at the end (no-go); otherwise the declared
#' winner, or at the maximum sample size the active arm most likely to be
#' best (go). So false go = 1 - P(original control selected).
#' This differs from adaptr's default target, P(stop for superiority), which
#' misses trials in which the original control was beaten but no single
#' winner emerged by the maximum sample size.
#' Assumes the control arm is called "Control" (true for all families).
#'
#' @param n_rep,base_seed,cores passed to `adaptr::run_trials()`.
#' @return A function of `x` and `trial_spec`.
make_false_go_fun <- function(n_rep, base_seed, cores = NULL) {
  function(x, trial_spec) {
    trial_spec$superiority <- x
    trial_spec$inferiority <- 1 - x
    sims <- adaptr::run_trials(trial_spec, n_rep = n_rep, cores = cores,
                               base_seed = base_seed)
    perf <- adaptr::check_performance(sims, select_strategy = "control or best")
    list(sims = sims, trial_spec = trial_spec,
         y = 1 - perf$est[perf$metric == "prob_select_arm_Control"])
  }
}

#' Calibrate the superiority threshold of each design in a family
#'
#' WRAPS `adaptr::calibrate_trial()`. Each design is calibrated under the null
#' scenario (all arms at the family's `null_rate`) so that the false-go rate
#' is at most `target`. Results are cached in `sims/`.
#'
#' @param family name of a family in `family_settings`.
#' @param target target false-go rate.
#' @param n_rep simulations per calibration step.
#' @param seed base random seed.
#' @param designs designs to calibrate.
#' @param futility_only_first see `make_design()`.
#' @return A named list of adaptr `trial_calibration` objects.
calibrate_family <- function(family, target = 0.10, n_rep = 5000, seed = 4131,
                             designs = design_names, futility_only_first = TRUE) {
  s <- family_settings[[family]]
  null_ys <- rep(s$null_rate, length(s$arms))
  suffix <- if (futility_only_first) "" else " - futility vs new control"
  designs |>
    set_names(map_chr(designs, design_label, futility_only_first = futility_only_first)) |>
    map(\(d) adaptr::calibrate_trial(
      make_design(d, null_ys, s, futility_only_first = futility_only_first),
      n_rep        = n_rep,
      base_seed    = seed,
      fun          = make_false_go_fun(n_rep, seed),
      target       = target,
      search_range = s$search_range,
      tol          = target / 20,
      dir          = -1,          # only accept false-go rates at or below target
      iter_max     = 30,
      path         = sprintf("%scalib - %s - %s%s - target %.3f - n%d.RDS",
                             dir_sims, family, file_safe(d), suffix, target, n_rep)
    ))
}

#' Calibrated thresholds, rounded up to 4 decimals
#'
#' OWN CODE (reads `best_x` from adaptr calibration objects). Rounding up
#' makes the threshold slightly stricter, so the false-go rate stays at or
#' below the target.
#'
#' @param cals output of `calibrate_family()`.
#' @return A named numeric vector of superiority thresholds.
get_thresholds <- function(cals) {
  map_dbl(cals, \(cl) ceiling(cl$best_x * 1e4) / 1e4)
}

#' Summary table of a calibration
#'
#' OWN CODE (reads adaptr calibration objects).
#'
#' @param cals output of `calibrate_family()`.
#' @return A tibble with one row per design.
threshold_table <- function(cals) {
  tibble(
    Design                             = names(cals),
    `Superiority threshold`            = get_thresholds(cals),
    `Inferiority threshold`            = 1 - get_thresholds(cals),
    `False-go rate (calibration sims)` = map_dbl(cals, \(cl) cl$best_y),
    `Calibration successful`           = map_lgl(cals, \(cl) cl$success),
    Evaluations                        = map_int(cals, \(cl) nrow(cl$evaluations))
  )
}

#' All calibration steps, for plotting
#'
#' OWN CODE (reads `evaluations` from adaptr calibration objects).
#'
#' @param cals output of `calibrate_family()`.
#' @return A tibble with design, x (threshold), y (false-go rate), best.
calibration_data <- function(cals) {
  imap(cals, \(cl, d) tibble(design = d, x = cl$evaluations$x,
                             y = cl$evaluations$y,
                             best = cl$evaluations$x == cl$best_x)) |>
    list_rbind()
}

#' Exact false-go rate of the fixed design (no simulation)
#'
#' OWN CODE (no adaptr). For the fixed design with one analysis, `n0`
#' participants on control and `nk` on each of `K` experimental arms, and all
#' arms with true response rate `p`:
#'
#'   false go(x) = sum_{s0} Bin(s0; n0, p) * [1 - F(x | s0)^K],
#'   F(x | s0)   = sum_{s: g(s, s0) <= x} Bin(s; nk, p),
#'
#' where g(s, s0) = P(p_k > p_0 | data) under Beta(1, 1) priors, computed by
#' numerical integration. Given the control count s0, the K arms are
#' independent, which gives the power K. Arm sizes are fixed at their
#' expected values (adaptr randomizes, so in the simulations they vary a
#' little).
#'
#' @param n0,nk participants on control and on each experimental arm.
#' @param K number of experimental arms.
#' @param p common true response rate.
#' @param target false-go rate to solve for.
#' @param n_grid grid points for the integral.
#' @return A list: `f` (false-go rate as a function of x) and `x_target`
#'   (smallest x with f(x) <= target).
fixed_false_go_exact <- function(n0, nk, K = 5, p = 0.15, target = 0.10,
                                 n_grid = 2000) {
  t  <- (seq_len(n_grid) - 0.5) / n_grid
  fk <- sapply(0:nk, \(s) dbeta(t, 1 + s, 1 + nk - s))   # density, arm k
  c0 <- sapply(0:n0, \(s) pbeta(t, 1 + s, 1 + n0 - s))   # cdf, control
  g  <- crossprod(fk, c0) / n_grid                      # P(p_k > p_0 | s_k, s_0)
  wk <- dbinom(0:nk, nk, p)
  w0 <- dbinom(0:n0, n0, p)
  f  <- \(x) sum(w0 * (1 - colSums(wk * (g <= x))^K))
  xs <- sort(unique(as.vector(g)))
  xs <- xs[xs > 0.5]
  fx <- vapply(xs, f, numeric(1))
  list(f = f, x_target = xs[which(fx <= target)[1]])
}

# =============================================================================
# 4. Simulation
# =============================================================================

#' Simulate designs for one scenario
#'
#' WRAPS `adaptr::run_trials()`. File names include the threshold and
#' `n_rep`, so changing either never clashes with saved simulations, and
#' existing files are re-loaded instead of re-simulated.
#'
#' @param id scenario identifier (used in file names).
#' @param ys true response rates.
#' @param family name of a family in `family_settings`.
#' @param thresholds named vector of superiority thresholds (names = design
#'   labels, as returned by `get_thresholds()`).
#' @param n_rep number of simulated trials per design.
#' @param seed base random seed.
#' @param designs,futility_only_first see `calibrate_family()`.
#' @param sparse passed to `adaptr::run_trials()`; `FALSE` keeps the
#'   look-by-look details needed by `adaptr::plot_history()`.
#' @return A named list of adaptr `trial_results` objects.
run_designs <- function(id, ys, family, thresholds, n_rep, seed = 4131,
                        designs = design_names, futility_only_first = TRUE,
                        sparse = TRUE) {
  s <- family_settings[[family]]
  labels <- map_chr(designs, design_label, futility_only_first = futility_only_first)
  map2(designs, labels, \(d, lab) adaptr::run_trials(
    make_design(d, ys, s, superiority = thresholds[[lab]],
                futility_only_first = futility_only_first),
    n_rep     = n_rep,
    base_seed = seed,
    sparse    = sparse,
    path      = sprintf("%s%s - %s - sup %.4f - n%d%s.RDS", dir_sims, id, file_safe(lab),
                        thresholds[[lab]], n_rep, if (sparse) "" else " - full")
  )) |>
    set_names(labels)
}

# =============================================================================
# 5. Summary objects (no simulations inside, so they are small)
# =============================================================================

#' Summarize simulated designs into a small object
#'
#' OWN CODE built on `adaptr::check_performance()` and
#' `adaptr::extract_results()`. The per-arm sample sizes are read from each
#' simulated trial (relies on adaptr 1.5.0 object structure:
#' `trial_results[[i]]$trial_res$ns_all`).
#'
#' The arm taken forward uses `select_strategy = "control or best"`: the
#' declared winner if the trial stopped with one arm left; otherwise the
#' original control if still active (a no-go), else the active arm most
#' likely to be best at the last analysis (a go).
#'
#' @param sims named list of adaptr `trial_results` (from `run_designs()`).
#' @param ys true response rates.
#' @param family name of a family in `family_settings`.
#' @param id scenario identifier.
#' @param thresholds named vector of thresholds used.
#' @param uncertainty logical; also compute bootstrap 95% CIs
#'   (`adaptr::check_performance(uncertainty = TRUE)`).
#' @param n_boot,seed bootstrap settings.
#' @return A list (summary object) with performance metrics, per-trial
#'   results, per-arm sample sizes and scenario information.
summarise_designs <- function(sims, ys, family, id = NA_character_,
                              thresholds = NULL, uncertainty = FALSE,
                              n_boot = 1000, seed = 4131) {
  s     <- family_settings[[family]]
  arms  <- s$arms
  ys    <- set_names(ys, arms)
  ctrl  <- arms[1]

  perf <- imap(sims, \(x, lab) {
    p <- adaptr::check_performance(x, select_strategy = "control or best",
                                   uncertainty = uncertainty, n_boot = n_boot,
                                   boot_seed = seed)
    p <- as_tibble(as.data.frame(p))
    if (!"lo_ci" %in% names(p)) p <- mutate(p, lo_ci = NA_real_, hi_ci = NA_real_)
    p |> select(metric, est, lo_ci, hi_ci) |> mutate(design = lab, .before = 1)
  }) |> list_rbind()

  res <- imap(sims, \(x, lab)
    adaptr::extract_results(x, select_strategy = "control or best") |>
      as_tibble() |>
      mutate(design = lab, .before = 1)
  ) |> list_rbind()

  # relies on adaptr 1.5.0 object structure: trial_res$ns_all per simulated trial
  arm_n <- imap(sims, \(x, lab) {
    n_all <- map(x$trial_results, \(tr) as.numeric(tr$trial_res$ns_all))
    tibble(design = lab,
           sim    = rep(seq_along(n_all), each = length(arms)),
           arm    = rep(arms, times = length(n_all)),
           n      = unlist(n_all))
  }) |> list_rbind()

  mean_sizes <- arm_n |>
    group_by(design, arm) |>
    summarise(mean_n = mean(n), .groups = "drop")

  list(
    id         = id,
    family     = family,
    designs    = names(sims),
    ys         = ys,
    arms       = arms,
    control    = ctrl,
    max_n      = s$max_n,
    best       = if (n_distinct(ys) == 1) NULL else names(ys)[which.max(ys)],
    useless    = setdiff(names(ys)[ys <= ys[[ctrl]]], ctrl),
    thresholds = thresholds,
    n_rep      = sims[[1]]$n_rep,
    perf       = perf,
    res        = res,
    arm_n      = arm_n,
    mean_sizes = mean_sizes,
    provenance = list(date = Sys.time(),
                      adaptr = as.character(utils::packageVersion("adaptr")),
                      R = R.version.string)
  )
}

#' One performance metric for one design
#'
#' OWN CODE.
#' @param S summary object.
#' @param d design label.
#' @param m metric
#' name as in `adaptr::check_performance()`.
#' @param what `"est"`, `"lo_ci"`
#' or `"hi_ci"`.
#' @return A number (`NA` if not available).
get_metric <- function(S, d, m, what = "est") {
  v <- S$perf[[what]][S$perf$design == d & S$perf$metric == m]
  if (length(v) == 0) NA_real_ else v
}

#' Key metrics for one design (named vector; NA = not applicable)
#'
#' OWN CODE, combining adaptr metrics; "false go" sums the selection
#' probabilities of arms that are no better than control.
#'
#' @param S summary object.
#' @param d design label.
#' @param what see
#' `get_metric()`.
#' @return A named numeric vector.
metric_values <- function(S, d, what = "est") {
  g <- \(m) get_metric(S, d, m, what)
  has_best <- !is.null(S$best)
  c(
    "Expected sample size"                     = g("size_mean"),
    "Sample size, 25th percentile"             = g("size_p25"),
    "Sample size, 75th percentile"             = g("size_p75"),
    "P(stop for superiority)"                  = g("prob_superior"),
    "P(stop for futility)"                     = g("prob_futility"),
    "P(reach max without decision)"            = g("prob_max"),
    "Power: P(best arm taken forward)"                   = if (has_best) g(paste0("prob_select_arm_", S$best)) else NA,
    "P(false go: a useless arm taken forward)" = if (length(S$useless))
      sum(map_dbl(S$useless, \(a) g(paste0("prob_select_arm_", a)))) else NA,
    "P(no-go: control taken forward)"          = g(paste0("prob_select_arm_", S$control)),
    "Response rate among trial participants"   = g("ratio_ys_mean"),
    "Total responders in trial"                = g("sum_ys_mean"),
    "RMSE of treatment effect (selected arm)"  = g("rmse_te")
  )
}

#' Table of key metrics: one row per metric, one column per design
#'
#' OWN CODE. Rows that are not applicable for any design are dropped.
#'
#' @param S summary object.
#' @param digits decimals for display.
#' @param digits_n decimals for sample-size rows.
#' @param ci logical; add bootstrap 95% CIs (if computed) as "est (lo, hi)".
#' @param bold_best logical; mark the best design per row in bold (markdown),
#'   for rows where "best" is well defined.
#' @return A tibble.
metric_table <- function(S, digits = 2, digits_n = 0, ci = FALSE, bold_best = FALSE) {
  est <- map(set_names(S$designs), \(d) metric_values(S, d))
  keep <- reduce(map(est, \(v) !is.na(v)), `|`)
  # sample-size rows get digits_n decimals, all other rows digits decimals
  row_digits <- set_names(if_else(grepl("sample size", names(est[[1]]), ignore.case = TRUE),
                                  digits_n, digits), names(est[[1]]))
  fmt_num <- \(x, dg) set_names(map2_chr(x, dg, \(v, k) formatC(v, format = "f", digits = k)), names(x))
  higher_better <- c("Power: P(best arm taken forward)" = TRUE,
                     "P(false go: a useless arm taken forward)" = FALSE,
                     "Expected sample size" = FALSE,
                     "Response rate among trial participants" = TRUE,
                     "RMSE of treatment effect (selected arm)" = FALSE)
  fmt <- imap(est, \(v, d) {
    out <- fmt_num(v, row_digits)
    if (ci) {
      lo <- metric_values(S, d, "lo_ci"); hi <- metric_values(S, d, "hi_ci")
      has_ci <- !is.na(lo) & !is.na(hi)
      out[has_ci] <- paste0(out[has_ci], " (", fmt_num(lo[has_ci], row_digits[has_ci]),
                            ", ", fmt_num(hi[has_ci], row_digits[has_ci]), ")")
    }
    out[is.na(v)] <- ""
    out
  })
  if (bold_best && length(S$designs) > 1) {
    for (m in intersect(names(higher_better), names(keep)[keep])) {
      vals <- map_dbl(est, \(v) round(v[[m]], row_digits[[m]]))  # compare as displayed
      if (all(is.na(vals)) || diff(range(vals, na.rm = TRUE)) < 1e-12) next  # nothing to compare
      best <- if (higher_better[[m]]) max(vals, na.rm = TRUE) else min(vals, na.rm = TRUE)
      for (d in names(est)) if (isTRUE(abs(vals[[d]] - best) < 1e-12))
        fmt[[d]][[m]] <- paste0("**", fmt[[d]][[m]], "**")
    }
  }
  as_tibble(fmt) |>
    mutate(Metric = names(est[[1]]), .before = 1) |>
    filter(keep)
}

#' Key metrics in long format (numeric), for combining scenarios
#'
#' OWN CODE.
#' @param S summary object.
#' @return A tibble (design, metric, value).
metrics_long <- function(S) {
  map(S$designs, \(d) {
    v <- metric_values(S, d)
    tibble(design = d, metric = names(v), value = unname(v))
  }) |>
    list_rbind() |>
    filter(!is.na(value))
}

#' "Bigger is better" metrics on a 0-1 scale
#'
#' OWN CODE. Smaller-is-better metrics are converted: 1 - P(false go),
#' 1 - expected n / max n, and relative accuracy = smallest RMSE / RMSE.
#'
#' @param S summary object.
#' @return A tibble (design, group, metric, value).
glance_data <- function(S) {
  has_best <- !is.null(S$best)
  rmse <- map_dbl(set_names(S$designs), \(d) get_metric(S, d, "rmse_te"))
  map(S$designs, \(d) {
    g  <- \(m) get_metric(S, d, m)
    mv <- metric_values(S, d)
    tibble(
      design = d,
      group  = c("Pick the right arm", "Pick the right arm",
                 "Efficiency", "Efficiency",
                 "Participant benefit", "Participant benefit",
                 "Accuracy"),
      metric = c("Power: P(best arm taken forward)",
                 "Avoid taking a useless arm forward",
                 "P(conclusive decision)", "Sample size saved",
                 "Share of participants on the best arm", "Response rate in trial",
                 "Estimation accuracy (relative)"),
      value  = c(
        if (has_best) mv[["Power: P(best arm taken forward)"]] else NA,
        1 - mv[["P(false go: a useless arm taken forward)"]],
        g("prob_conclusive"),
        1 - g("size_mean") / S$max_n,
        if (has_best) {
          ms <- S$mean_sizes[S$mean_sizes$design == d, ]
          ms$mean_n[ms$arm == S$best] / sum(ms$mean_n)
        } else NA,
        if (has_best) g("ratio_ys_mean") else NA,
        if (has_best && all(is.finite(rmse))) min(rmse) / rmse[[d]] else NA
      )
    )
  }) |>
    list_rbind() |>
    filter(!is.na(value))
}

#' Probability that each arm is taken forward
#'
#' OWN CODE (reads `prob_select_arm_*` from adaptr performance metrics).
#' @param S summary object.
#' @return A tibble (design, arm, prob, type).
selection_data <- function(S) {
  S$perf |>
    filter(startsWith(metric, "prob_select_arm_")) |>
    mutate(arm = sub("prob_select_arm_", "", metric, fixed = TRUE),
           type = case_when(
             arm == S$control                    ~ "Control (no-go)",
             !is.null(S$best) & arm %in% S$best  ~ "Best arm",
             arm %in% S$useless                  ~ "Useless arm",
             TRUE                                ~ "Better than control, not best"
           )) |>
    select(design, arm, prob = est, type)
}

#' Bias, SD and RMSE of the estimated treatment effect
#'
#' OWN CODE on `adaptr::extract_results()` output (`err_te` = estimated minus
#' true effect of the selected arm vs. the original control). Only trials
#' taking an experimental arm forward contribute.
#'
#' @param S summary object.
#' @return A tibble, one row per design.
bias_table <- function(S) {
  ctrl_n <- S$mean_sizes |> filter(arm == S$control) |> select(design, `Mean n, control` = mean_n)
  best_n <- if (is.null(S$best)) NULL else
    S$mean_sizes |> filter(arm == S$best) |> select(design, `Mean n, best arm` = mean_n)
  out <- S$res |>
    filter(!is.na(err_te)) |>
    mutate(true_te = unname(S$ys[selected_arm] - S$ys[[S$control]])) |>
    group_by(design) |>
    summarise(`Trials included` = n(),
              Bias = mean(err_te),
              `Bias (% of true effect)` = if (mean(true_te) > 0) 100 * mean(err_te) / mean(true_te) else NA_real_,
              SD = sd(err_te),
              RMSE = sqrt(mean(err_te^2)),
              .groups = "drop") |>
    left_join(ctrl_n, by = "design")
  if (!is.null(best_n)) out <- left_join(out, best_n, by = "design")
  out |>
    select(Design = design, `Trials included`, any_of(c("Mean n, control", "Mean n, best arm")),
           Bias, any_of("Bias (% of true effect)"), SD, RMSE) |>
    arrange(match(Design, S$designs))
}

#' How each simulated trial ended
#'
#' OWN CODE. Endings are classified by adaptr's `final_status` AND by whether
#' an arm other than the original control was selected (a go):
#' * adaptr records "futility" when the last comparators are dropped for
#'   futility. With FNC this can happen after a winner has replaced the
#'   control, so those trials are a go, not a no-go.
#' * adaptr records "superiority" whenever one arm is left, even if that arm
#'   is the original control (all others dropped, the last for inferiority);
#'   that is a no-go.
#' * "max" trials can be a go or a no-go (see `select_strategy` in
#'   `summarise_designs()`); they are kept as one category.
#'
#' @param S summary object.
#' @return A factor (one value per simulated trial).
classify_endings <- function(S) {
  go <- !is.na(S$res$selected_arm) & S$res$selected_arm != S$control
  out <- case_when(
    S$res$final_status == "superiority" & go & S$res$final_n < S$max_n ~ "Winner declared early",
    S$res$final_status == "superiority" & go                    ~ "Winner declared at final analysis",
    S$res$final_status == "futility" & go                       ~ "Winner kept; other arms dropped for futility",
    S$res$final_status %in% c("superiority", "futility")        ~ "All experimental arms dropped (no-go)",
    S$res$final_status == "equivalence"                         ~ "Stopped for equivalence",
    TRUE                                                        ~ "No winner by the maximum"
  )
  factor(out, levels = names(outcome_cols))
}

#' Share of trials per ending, with label positions
#'
#' OWN CODE.
#' @param S summary object.
#' @return A tibble.
endings_data <- function(S) {
  S$res |>
    mutate(outcome = classify_endings(S)) |>
    count(design, outcome) |>
    group_by(design) |>
    mutate(share = n / sum(n)) |>
    arrange(design, outcome) |>
    mutate(x_mid = cumsum(share) - share / 2) |>
    ungroup()
}

#' Share of trials still running after n participants
#'
#' OWN CODE.
#' @param S summary object.
#' @return A tibble (design, n, still_running).
running_data <- function(S) {
  grid <- sort(unique(c(0, S$res$final_n)))
  map(S$designs, \(d) {
    n_d <- S$res$final_n[S$res$design == d]
    tibble(design = d, n = grid, still_running = map_dbl(grid, \(x) mean(n_d > x)))
  }) |>
    list_rbind()
}

#' Remove per-trial details to make a summary object even smaller
#'
#' OWN CODE. Keeps what the scenario overview needs.
#' @param S summary object.
#' @return A list.
compact_summary <- function(S) {
  S[c("id", "family", "designs", "ys", "arms", "control", "max_n", "best",
      "useless", "thresholds", "n_rep", "perf", "mean_sizes", "provenance")]
}

# =============================================================================
# 6. Data behind adaptr's own plots (so the post does not need the sims)
# =============================================================================

#' Arm statuses over time
#'
#' WRAPS `adaptr::plot_status()` and keeps only the plotted data
#' (columns x, p, status and arm_facet for arm-level plots).
#'
#' @param sims an adaptr `trial_results` object.
#' @param arms arms to include (`NULL` = overall trial status).
#' @return A tibble.
status_data <- function(sims, arms = NULL) {
  as_tibble(adaptr::plot_status(sims, x_value = "total n", arm = arms)$data)
}

#' Allocation probabilities over time (median and spread across trials)
#'
#' WRAPS `adaptr::plot_history()` and keeps only the plotted data
#' (columns x, lo, mid, hi, arm). Needs simulations run with `sparse = FALSE`.
#'
#' @param sims an adaptr `trial_results` object (sparse = FALSE).
#' @return A tibble.
history_data <- function(sims) {
  as_tibble(adaptr::plot_history(sims, x_value = "total n", y_value = "prob")$data)
}

#' Convergence of performance metrics with the number of simulations
#'
#' WRAPS `adaptr::plot_convergence()` and keeps only the plotted data
#' (columns res_points = number of simulations, y, split, metric, labels).
#'
#' @param sims an adaptr `trial_results` object.
#' @param metrics metrics to track (see `adaptr::plot_convergence()`).
#' @param n_split number of consecutive batches (1 = cumulative curve).
#' @return A tibble.
convergence_data <- function(sims, metrics, n_split = 1) {
  as_tibble(adaptr::plot_convergence(sims, metrics = metrics,
                                     select_strategy = "control or best",
                                     n_split = n_split)$data)
}

# =============================================================================
# 7. Single simulated trials
# =============================================================================

#' Look-by-look history of one simulated trial
#'
#' WRAPS `adaptr::run_trial(sparse = FALSE)`; reads the look-by-look results
#' (relies on adaptr 1.5.0 object structure: `all_looks`, with fields
#' `arms`, `old_status`, `new_status`, `old_alloc`, `new_alloc`,
#' `probs_better_first`, `ns_all`, `sum_ys_all`). `probs_better_first` holds
#' the probabilities computed at the start of each analysis, i.e. the values
#' that triggered that analysis' decisions.
#'
#' @param spec an adaptr `trial_spec`.
#' @param seed random seed.
#' @return A tibble with one row per arm and analysis.
trial_history <- function(spec, seed) {
  tr <- adaptr::run_trial(spec, seed = seed, sparse = FALSE)
  imap(tr$all_looks, \(x, i) tibble(
    look       = i,
    n          = sum(x$ns_all),
    arm        = x$arms,
    old_status = x$old_status,
    status     = x$new_status,
    old_alloc  = as.numeric(x$old_alloc),
    alloc      = as.numeric(x$new_alloc),
    p_better   = if (is.null(x$probs_better_first)) NA_real_ else as.numeric(x$probs_better_first),
    responders = as.numeric(x$sum_ys_all)
  )) |>
    list_rbind()
}

#' Status changes in a trial history (drops, promotions, winner)
#'
#' OWN CODE.
#' @param h output of `trial_history()`.
#' @return A tibble.
trial_events <- function(h) {
  h |>
    filter(status != old_status) |>
    mutate(event = case_when(
      status == "control"                               ~ "promoted to control",
      status == "inferior" & old_status == "control"    ~ "replaced as control",
      status == "inferior"                              ~ "dropped: inferior",
      status == "futile"                                ~ "dropped: futility",
      status == "superior"                              ~ "declared winner",
      TRUE                                              ~ status
    ))
}

#' Promotions of an experimental arm to control
#'
#' OWN CODE.
#' @param h output of `trial_history()`.
#' @return A tibble (look, n, arm).
promotions <- function(h) {
  h |> filter(status == "control", old_status == "active") |> select(look, n, arm)
}

#' First seed whose simulated trial meets a condition
#'
#' OWN CODE (uses `trial_history()`).
#' @param spec trial specification.
#' @param cond function of a trial history
#' returning `TRUE`/`FALSE`.
#' @param seeds seeds to try.
#' @return An integer seed (or `NA`).
find_seed <- function(spec, cond, seeds = 1:500) {
  detect(seeds, \(sd) isTRUE(cond(trial_history(spec, sd)))) %||% NA_integer_
}

# =============================================================================
# 8. Allocation rules in tables
# =============================================================================

#' Allocation shares under fixed square-root allocation (design b)
#'
#' OWN CODE reproducing adaptr 1.5.0 behavior with
#' `control_prob_fixed = "sqrt-based fixed"`: the control's share is
#' recalculated by the square-root rule for the k active experimental arms,
#' but the experimental arms keep their initial weight 1 / (sqrt(K) + K);
#' all shares are then rescaled to sum to 1. The pure square-root rule is
#' shown for contrast.
#'
#' @param K initial number of experimental arms.
#' @return A tibble.
alloc_table_fixed <- function(K = 5) {
  tibble(k = K:1) |>
    mutate(ctrl_raw = sqrt(k) / (sqrt(k) + k),
           exp_raw  = 1 / (sqrt(K) + K),
           tot      = ctrl_raw + k * exp_raw) |>
    transmute(`Active experimental arms`                      = k,
              `Control share (adaptr)`                        = ctrl_raw / tot,
              `Each experimental arm (adaptr)`                = exp_raw / tot,
              `Control share (pure square-root rule)`         = sqrt(k) / (sqrt(k) + k),
              `Each experimental arm (pure square-root rule)` = 1 / (sqrt(k) + k))
}

#' Control share and minimum shares with RAR (design c)
#'
#' OWN CODE reproducing adaptr 1.5.0 behavior with
#' `control_prob_fixed = "sqrt-based"`, `min_probs` and
#' `rescale_probs = "limits"`.
#'
#' @param K initial number of experimental arms.
#' @param min_prob initial
#' minimum share per experimental arm.
#' @return A tibble.
alloc_table_rar <- function(K = 5, min_prob = 0.05) {
  tibble(k = K:1) |>
    transmute(`Active experimental arms`              = k,
              `Control share (fixed)`                 = sqrt(k) / (sqrt(k) + k),
              `Shared by RAR among experimental arms` = 1 - sqrt(k) / (sqrt(k) + k),
              `Minimum share per experimental arm`    = if_else(k > 1, min_prob * (K + 1) / (k + 1), NA_real_))
}

# =============================================================================
# 9. Plots (own ggplot2 code; adaptr plots are re-drawn from their data)
# =============================================================================

theme_post <- function() {
  theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.position = "top")
}

#' One interim analysis as a flowchart (designs b and c)
#'
#' OWN CODE (ggplot2 only, so it renders identically to HTML and PDF).
#' @return A ggplot.
plot_flowchart <- function() {
  H <- 0.72; MW <- 4.6; RW <- 3.3; RX <- 5.2
  nodes <- tibble(
    x = c(rep(0, 8), rep(RX, 5)),
    y = c(8:1, 6, 5, 4, 3, 2),
    w = c(rep(MW, 8), rep(RW, 5)),
    type = c("action", "action", "decision", "decision", "decision", "decision",
             "decision", "action", "effect", "effect", "effect", "stop", "stop"),
    label = c(
      "Interim analysis\n(every 30 participants with outcome data)",
      "Update posteriors for each active arm:\nBeta(1 + responders, 1 + non-responders)",
      "Inferiority:\nP(arm better than control) < 1 - x ?",
      "Superiority:\nP(arm better than control) > x ?",
      "Futility (original control only):\nP(not better by >= 5 points) > 0.90 ?",
      "Only one arm left?",
      "Maximum sample size reached?",
      "Update allocation for the next participants",
      "Drop the arm",
      "Arm becomes the NEW control;\nold control dropped",
      "Drop the arm",
      "STOP: superiority\n(futility if all others were futile)",
      "STOP: maximum reached\n(selection: see Step 3)"
    )
  )
  fills <- c(action = "#e8eef7", decision = "#fff3d6", effect = "#eeeeee", stop = "#f6d5d0")
  ys_dec <- c(6, 5, 4, 3, 2)
  down  <- tibble(x = 0, y = 8:2 - H / 2, xend = 0, yend = 7:1 + H / 2)
  right <- tibble(x = MW / 2, y = ys_dec, xend = RX - RW / 2, yend = ys_dec,
                  lt = c("dashed", "dashed", "dashed", "solid", "solid"))
  no_lab  <- tibble(x = 0.15, y = ys_dec - 0.5, label = "no")
  yes_lab <- tibble(x = (MW / 2 + RX - RW / 2) / 2, y = ys_dec + 0.12, label = "yes")
  XR <- RX + RW / 2 + 0.35
  recheck <- tibble(x = c(RX + RW / 2, XR, XR, 0.8, 0.8), y = c(5, 5, 6.5, 6.5, 6 + H / 2))
  LX <- -MW / 2 - 0.4
  loop <- tibble(x = c(-MW / 2, LX, LX, -MW / 2), y = c(1, 1, 8, 8))
  arr <- grid::arrow(length = grid::unit(2, "mm"), type = "closed")
  ggplot() +
    geom_rect(data = nodes, aes(xmin = x - w / 2, xmax = x + w / 2,
                                ymin = y - H / 2, ymax = y + H / 2, fill = type),
              colour = "grey55", linewidth = 0.3) +
    geom_text(data = nodes, aes(x = x, y = y, label = label), size = 3.1, lineheight = 0.95) +
    geom_segment(data = down, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arr, colour = "grey35") +
    geom_segment(data = right, aes(x = x, y = y, xend = xend, yend = yend, linetype = lt),
                 arrow = arr, colour = "grey35") +
    geom_path(data = recheck, aes(x = x, y = y), linetype = "dotted", arrow = arr, colour = "grey35") +
    geom_path(data = loop, aes(x = x, y = y), arrow = arr, colour = "grey35") +
    geom_text(data = no_lab, aes(x = x, y = y, label = label), hjust = 0, size = 2.8, colour = "grey35") +
    geom_text(data = yes_lab, aes(x = x, y = y, label = label), vjust = 0, size = 2.8, colour = "grey35") +
    annotate("text", x = XR + 0.1, y = 5.75, label = "re-check against\nthe new control",
             hjust = 0, size = 2.8, colour = "grey35") +
    annotate("text", x = LX - 0.12, y = 4.5, label = "next analysis",
             angle = 90, vjust = 0, size = 2.8, colour = "grey35") +
    scale_fill_manual(values = fills, guide = "none") +
    scale_linetype_identity() +
    coord_cartesian(xlim = c(LX - 0.45, XR + 1.6), ylim = c(0.5, 8.5), expand = FALSE) +
    theme_void()
}

#' Single simulated trials, look by look
#'
#' OWN CODE. Three panels per trial (columns):
#' * Events: what happened at which analysis (labels).
#' * Allocation share: one stacked bar per interval between analyses; the
#'   segments are the shares of the active arms (control at the bottom),
#'   labeled when >= 10%.
#' * P(better than the current control), with the superiority and
#'   inferiority thresholds.
#' Events are also vertical lines across all panels (color = arm, line type
#' = event); simultaneous events are offset slightly so none is hidden.
#'
#' @param hists named list of trial histories (`trial_history()`).
#' @param thr named vector of superiority thresholds (same names).
#' @param max_n maximum sample size.
#' @return A ggplot.
plot_anatomy <- function(hists, thr, max_n) {
  p_ev <- "Events"; p_al <- "Allocation share"; p_pb <- "P(better than current control)"
  panels <- c(p_ev, p_al, p_pb)
  lev <- names(hists)
  dx <- max_n / 120                                  # offset for simultaneous events

  # allocation as steps: initial share from 0 to the first analysis, then the
  # share set at each analysis until the next analysis
  steps <- imap(hists, \(h, nm) {
    looks <- h |> distinct(look, n) |> arrange(look)
    first <- h |> filter(look == 1) |>
      transmute(arm, x = 0, xend = looks$n[1], y = old_alloc)
    # only arms still active after the analysis: adaptr 1.5.0 can report a
    # share for an arm dropped in the re-check after a promotion
    later <- h |> filter(!is.na(alloc), status %in% c("active", "control")) |>
      left_join(looks |> transmute(look = look - 1, xend = n), by = "look") |>
      filter(!is.na(xend)) |>
      transmute(arm, x = n, xend, y = alloc)
    bind_rows(first, later) |> filter(!is.na(y)) |> mutate(trial = nm)
  }) |> list_rbind()

  # stack the shares of the active arms within each interval (they sum to 1),
  # in a fixed arm order (control at the bottom)
  bars <- steps |>
    mutate(arm = factor(arm, levels = names(arm_cols))) |>
    arrange(trial, x, arm) |>
    group_by(trial, x) |>
    mutate(ymax = cumsum(y), ymin = ymax - y) |>
    ungroup() |>
    mutate(arm = as.character(arm),
           label = if_else(y >= 0.10, scales::percent(y, accuracy = 1), ""))

  probs <- imap(hists, \(h, nm) h |> filter(!is.na(p_better)) |>
                  transmute(trial = nm, arm, n, y = p_better)) |> list_rbind()

  # "replaced as control" is implied by the promotion at the same analysis
  events <- imap(hists, \(h, nm) trial_events(h) |>
                   filter(event != "replaced as control") |>
                   group_by(n) |>
                   mutate(k = row_number(), n_k = n()) |>
                   ungroup() |>
                   mutate(trial = nm,
                          x = n + (k - (n_k + 1) / 2) * dx,
                          label = paste0(arm, ": ", event))) |>
    list_rbind()

  thr_lines <- tibble(trial = rep(lev, each = 2),
                      y = c(rbind(thr[lev], 1 - thr[lev])),
                      panel = p_pb)

  look_n <- sort(unique(unlist(map(hists, \(h) unique(h$n)))))
  fix <- \(df, pnl) mutate(df, trial = factor(trial, levels = lev),
                           panel = factor(pnl, levels = panels))
  all_panels <- \(df) bind_rows(fix(df, p_ev), fix(df, p_al), fix(df, p_pb))
  winners <- filter(events, event == "declared winner")
  winner_layers <- if (nrow(winners) == 0) NULL else list(
    geom_rect(data = all_panels(winners),
              aes(xmin = x - 1.5 * dx, xmax = x + 1.5 * dx, ymin = -Inf, ymax = Inf, fill = arm),
              alpha = 0.2, inherit.aes = FALSE, show.legend = FALSE),
    geom_vline(data = all_panels(winners),
               aes(xintercept = x, colour = arm, linetype = event), linewidth = 1.6)
  )
  ev_lt <- c("promoted to control" = "solid", "declared winner" = "solid",
             "dropped: inferior" = "dashed",
             "dropped: futility" = "dotted")

  ggplot() +
    # allocation: one stacked bar per interval between analyses
    geom_rect(data = fix(bars, p_al),
              aes(xmin = x, xmax = xend, ymin = ymin, ymax = ymax, fill = arm),
              colour = "white", linewidth = 0.3, alpha = 0.85, show.legend = FALSE) +
    geom_text(data = fix(bars, p_al),
              aes(x = (x + xend) / 2, y = (ymin + ymax) / 2, label = label),
              colour = "white", size = 2.1) +
    # the winner (if any): a shaded band and a thick line across all panels
    winner_layers +
    # other events: thin lines across all panels
    geom_vline(data = all_panels(filter(events, event != "declared winner")),
               aes(xintercept = x, colour = arm, linetype = event), linewidth = 0.6) +
    # fixed y ranges: 0-100% for allocation and probabilities; the (short)
    # events panel uses negative values, which never get axis labels
    geom_blank(data = fix(tibble(trial = rep(lev, 2), y = rep(c(-0.6, 0), each = length(lev))), p_ev),
               aes(y = y)) +
    geom_blank(data = fix(tibble(trial = rep(lev, 2), y = rep(c(0, 1), each = length(lev))), p_al),
               aes(y = y)) +
    geom_blank(data = fix(tibble(trial = rep(lev, 2), y = rep(c(0, 1), each = length(lev))), p_pb),
               aes(y = y)) +
    # event labels, stacked in the events panel
    geom_text(data = fix(events, p_ev) |> mutate(y = -0.08 - 0.11 * (k - 1)),
              aes(x = x, y = y, label = label, colour = arm),
              hjust = 0, nudge_x = dx / 2, size = 2.5, show.legend = FALSE) +
    # P(better than current control)
    geom_hline(data = fix(thr_lines, p_pb), aes(yintercept = y), linetype = "dashed", colour = "grey50") +
    geom_line(data = fix(probs, p_pb), aes(x = n, y = y, colour = arm), linewidth = 0.6) +
    geom_point(data = fix(probs, p_pb), aes(x = n, y = y, colour = arm), size = 1.4) +
    facet_grid(panel ~ trial, scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = arm_cols, name = "Arm") +
    scale_linetype_manual(values = ev_lt, name = "Event") +
    scale_fill_manual(values = arm_cols, guide = "none") +
    scale_x_continuous(limits = c(0, max_n + 4 * dx), breaks = c(0, look_n)) +
    scale_y_continuous(labels = scales::percent,
                       breaks = \(lims) if (lims[1] < -0.1) numeric(0) else seq(0, 1, 0.25)) +
    labs(x = "Participants randomized", y = NULL) +
    theme_post() +
    theme(legend.position = "bottom", legend.box = "vertical",
          strip.placement = "outside",
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
}

#' "Bigger is better" chart for one scenario
#'
#' OWN CODE.
#' @param S summary object.
#' @return A ggplot.
plot_glance <- function(S) {
  gd <- glance_data(S)
  gd |>
    mutate(group  = factor(group, levels = unique(gd$group)),
           metric = factor(metric, levels = rev(unique(gd$metric))),
           design = factor(design, levels = rev(S$designs))) |>
    ggplot(aes(x = value, y = metric, fill = design)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75,
             colour = "white", linewidth = 0.4) +
    geom_text(aes(label = scales::percent(value, accuracy = 1)),
              position = position_dodge(width = 0.8), hjust = -0.15, size = 3, colour = "grey20") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                       labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = design_cols, breaks = S$designs) +
    labs(x = "Bigger is better", y = NULL, fill = NULL) +
    theme_post() +
    theme(panel.grid.major.y = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0))
}

#' Which arm is taken forward, per design
#'
#' OWN CODE.
#' @param S summary object.
#' @return A ggplot.
plot_selection <- function(S) {
  type_cols <- c("Best arm" = "#1b7837", "Better than control, not best" = "#86c286",
                 "Useless arm" = "#d6604d", "Control (no-go)" = "grey60")
  sd_ <- selection_data(S) |>
    mutate(design = factor(design, levels = S$designs), arm = factor(arm, levels = S$arms))
  ggplot(sd_, aes(x = arm, y = prob, fill = type)) +
    geom_col() +
    geom_text(aes(label = if_else(prob >= 0.01, scales::percent(prob, accuracy = 1), "")),
              vjust = -0.3, size = 3) +
    facet_wrap(~ design) +
    scale_fill_manual(values = type_cols, breaks = intersect(names(type_cols), sd_$type)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.08), expand = c(0, 0)) +
    labs(x = NULL, y = "P(arm taken forward)", fill = NULL) +
    theme_post()
}

#' Mean participants per arm, per design
#'
#' OWN CODE.
#' @param S summary object.
#' @return A ggplot.
plot_arm_sizes <- function(S) {
  S$mean_sizes |>
    mutate(design  = factor(design, levels = S$designs),
           arm     = factor(arm, levels = S$arms),
           is_best = !is.null(S$best) & arm %in% S$best) |>
    ggplot(aes(x = arm, y = mean_n, fill = is_best)) +
    geom_col() +
    facet_wrap(~ design) +
    scale_fill_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "grey65"), guide = "none") +
    labs(x = NULL, y = "Mean participants randomized",
         subtitle = if (is.null(S$best)) "All arms have the same true response rate"
                    else paste0("Green = the truly best arm (", S$best, ")")) +
    theme_post()
}

#' How trials end, per design
#'
#' OWN CODE.
#' @param S summary object.
#' @return A ggplot.
plot_endings <- function(S) {
  ends <- endings_data(S) |>
    mutate(design = factor(design, levels = rev(S$designs)),
           txt_col = if_else(outcome %in% c("Winner declared early",
                                            "Winner kept; other arms dropped for futility"),
                             "white", "grey10"))
  ggplot(ends, aes(x = share, y = design, fill = outcome, group = outcome)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.4,
             position = position_stack(reverse = TRUE)) +
    geom_text(aes(x = x_mid, y = design,
                  label = if_else(share >= 0.04, scales::percent(share, accuracy = 1), ""),
                  colour = txt_col),
              size = 3.2, inherit.aes = FALSE) +
    scale_colour_identity() +
    scale_fill_manual(values = outcome_cols,
                      breaks = intersect(names(outcome_cols), as.character(ends$outcome))) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = "% of simulated trials", y = NULL, fill = NULL) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_post() +
    theme(panel.grid = element_blank())
}

#' Share of trials still running, per design
#'
#' OWN CODE.
#' @param S summary object.
#' @return A ggplot.
plot_running <- function(S) {
  rd <- running_data(S) |> mutate(design = factor(design, levels = S$designs))
  ggplot(rd, aes(x = n, y = still_running, colour = design)) +
    geom_step(linewidth = 1, direction = "hv") +
    scale_colour_manual(values = design_cols, breaks = S$designs) +
    scale_x_continuous(breaks = unique(rd$n)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Participants randomized", y = "% of trials still running", colour = NULL) +
    theme_post() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
}

#' Calibration steps per design
#'
#' OWN CODE.
#' @param cal_data output of `calibration_data()`.
#' @param target
#' target false-go rate.
#' @return A ggplot.
plot_calibration <- function(cal_data, target) {
  cal_data |>
    mutate(design = factor(design, levels = unique(cal_data$design))) |>
    ggplot(aes(x = x, y = y, colour = design)) +
    geom_hline(yintercept = target, linetype = "dashed") +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_point(data = \(d) filter(d, best), shape = 21, size = 5, stroke = 1.2, fill = NA) +
    facet_wrap(~ design, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = design_cols, guide = "none") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Superiority threshold x (inferiority = 1 - x)",
         y = "False-go rate under the null",
         subtitle = "Each point is one calibration step; the circled point is the chosen threshold") +
    theme_post()
}

#' Arm statuses over time (re-drawn from `status_data()`)
#'
#' OWN CODE re-drawing `adaptr::plot_status()` data.
#' @param sd output of
#' `status_data()`.
#' @return A ggplot.
plot_status_data <- function(sd) {
  sd <- sd |> mutate(status = factor(status, levels = intersect(names(status_cols), unique(status))))
  p <- ggplot(sd, aes(x = x, y = p, fill = status)) +
    geom_area(alpha = 0.8) +
    scale_fill_manual(values = status_cols) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Participants randomized", y = "% of simulated trials", fill = NULL) +
    theme_post()
  if ("arm_facet" %in% names(sd)) p <- p + facet_wrap(~ arm_facet, nrow = 1)
  p
}

#' Allocation probabilities over time (re-drawn from `history_data()`)
#'
#' OWN CODE re-drawing `adaptr::plot_history()` data.
#' @param hd named list of `history_data()` outputs (one per design).
#' @return A ggplot.
plot_history_data <- function(hd) {
  imap(hd, \(d, nm) mutate(d, design = nm)) |>
    list_rbind() |>
    mutate(design = factor(design, levels = names(hd))) |>
    ggplot(aes(x = x)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = arm), alpha = 0.15) +
    geom_line(aes(y = mid, colour = arm), linewidth = 0.9) +
    facet_wrap(~ design) +
    scale_colour_manual(values = arm_cols) +
    scale_fill_manual(values = arm_cols) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Participants randomized", y = "Allocation probability",
         colour = "Arm", fill = "Arm") +
    theme_post()
}

#' Convergence of metrics (re-drawn from `convergence_data()`)
#'
#' OWN CODE re-drawing `adaptr::plot_convergence()` data.
#' @param cd named list of `convergence_data()` outputs (one per design).
#' @return A ggplot.
plot_convergence_data <- function(cd) {
  imap(cd, \(d, nm) mutate(d, design = nm)) |>
    list_rbind() |>
    mutate(design = factor(design, levels = names(cd))) |>
    ggplot(aes(x = res_points, y = y, colour = design,
               group = interaction(design, split))) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ labels, scales = "free_y") +
    scale_colour_manual(values = design_cols) +
    labs(x = "Number of simulated trials", y = NULL, colour = NULL) +
    theme_post()
}

#' Compare designs across scenarios (one panel per metric)
#'
#' OWN CODE. Scenarios are rows; the designs are points on the same row, so
#' designs are compared left to right on a common scale.
#'
#' @param gl tibble with columns scenario, design, metric, value
#'   (from `glance_data()` of several scenarios).
#' @param scenarios scenario order.
#' @param designs design order.
#' @return A ggplot.
plot_scenarios <- function(gl, scenarios, designs = design_names) {
  gl <- gl |>
    filter(design %in% designs) |>
    mutate(scenario = factor(scenario, levels = rev(scenarios)),
           design   = factor(design, levels = designs),
           metric   = factor(metric, levels = unique(gl$metric)))
  rng <- gl |> group_by(scenario, metric) |>
    summarise(lo = min(value), hi = max(value), .groups = "drop")
  ggplot(gl, aes(y = scenario)) +
    geom_linerange(data = rng, aes(xmin = lo, xmax = hi), colour = "grey70", linewidth = 1.2) +
    geom_point(aes(x = value, colour = design, shape = design), size = 2.8) +
    facet_wrap(~ metric, ncol = 2) +
    scale_colour_manual(values = design_cols, breaks = designs) +
    scale_shape_manual(values = c(16, 17, 15, 1, 2)[seq_along(designs)] |> set_names(designs)) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Bigger is better", y = NULL, colour = NULL, shape = NULL) +
    theme_post()
}

