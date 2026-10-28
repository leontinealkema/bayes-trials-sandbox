# Bayesian methods in clinical trials

This is [Leontine Alkema](https://leontinealkema.github.io/)'s sandbox for learning and trying out Bayesian methods for clinical trial design and analysis. Each topic folder contains a worked example built around an open-source R package: a post that explains the method, the R code that produces the results, and the saved results the post reads. The examples are meant to be didactic: they use simple settings to show how a method works.

## Topics

| Folder | Topic | Post |
|----|----|----|
| [`adaptive-designs/`](adaptive-designs/) | How adaptive designs work: a phase II pick-the-winner example (stopping, arm dropping and response-adaptive randomization, simulated with the adaptr R package) | [Bayesian adaptive designs: a worked example](adaptive-designs/adaptive-design101.md) ([PDF](adaptive-designs/adaptive-design101.pdf)) |

## Structure of a topic folder

```         
<topic>/
├── <post>.qmd      a post (renders to HTML, PDF and GitHub markdown)
├── <post>.md       the post as GitHub markdown (figures in <post>_files/)
├── run_all.R       runs all simulations; writes results/
├── R/helpers.R     helper functions (with roxygen documentation)
├── references.bib  references cited in the post(s)
├── results/        small saved results read by the post(s) (in git)
├── sims/           cached raw simulations (not in git)
└── private/        local material such as papers (not in git)
```

## About this sandbox

**Credits**

These examples stand on the work of the authors of the R packages and tutorials they use. Thank you for making your work openly available. Each post cites the package (with the version used) and the papers it relies on.

**Use of Claude**

I developed these examples with the help of Claude, to draft R code and figures, check settings against the package source code and documentation, and draft text. I direct the work and review the code, results and text, and I take responsibility for the content.

**Questions and input welcome**

I’m a Bayesian statistician who is doing a deep dive into clinical trial design, and I’m sharing my explorations as I learn. Questions, corrections and suggestions are very welcome. Please email me at [lalkema\@umass.edu](mailto:lalkema@umass.edu).

**License**

Text and figures are licensed under [CC BY 4.0](LICENSE-CC-BY-4.0) ([summary](https://creativecommons.org/licenses/by/4.0/)); code is licensed under the [MIT License](LICENSE). © 2026 Leontine Alkema.
