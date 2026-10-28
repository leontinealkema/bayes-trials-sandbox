# Adaptive designs

How Bayesian adaptive designs work: a phase II pick-the-winner trial with five experimental arms and a control, simulated with the R package [adaptr](https://inceptdk.github.io/adaptr/index.html).

**Read the post:** [Bayesian adaptive designs: a worked example](adaptive-design101.md) ([PDF](adaptive-design101.pdf))

## Files

| File | What it is |
|---|---|
| `adaptive-design101.qmd` | the post (Quarto source) |
| `adaptive-design101.md` | the post as GitHub markdown |
| `run_all.R` | runs all simulations and writes `results/` |
| `R/helpers.R` | helper functions (roxygen documentation) |
| `references.bib` | references cited in the post |
| `results/` | saved results that the post reads |

## Reproduce

Open `bayes-trials-sandbox.Rproj` (repository root), run `adaptive-designs/run_all.R`, then render `adaptive-design101.qmd`. Requires adaptr 1.5.0. The first run of `run_all.R` takes a while; simulations are cached in `sims/` (not in git), so later runs only compute what changed.
