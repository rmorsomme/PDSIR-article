
<!-- README.md is generated from README.Rmd. Please edit that file -->

## PDSIR-article – a project for reproducibility

The `R` project `PDSIR-article` provides `R` scripts to reproduce the
analyses present in the article Exact Bayesian Inference for Stochastic
Epidemic Models via Uniformly Ergodic Block Sampling by R. Morsomme and
J. Xu. The `R` package `PDSIR` contains the functions to run the
sampling algorithm introduced in the article, `PDSIR-article` is a
lightweight project for reproducibility.

The `R` scripts are located in the folder `R`. The five scripts
`analysis-*.R` run Markov chain Monte Carlo samplers whose draws are
saved in the folder `mcmc-draws`. The script `make_figures.R` makes the
various figures and computes the summary statistics found in the
article. The figures are saved in the folder `figures`.

For convenience, `.sh` scripts to run the `R` scripts on `SLURM` are
provided in the folder `slurm`. Running the sampling algorithm on
`SLURM` is convenient for computation-heavy tasks, such as computing the
posterior coverage rate across 1000 simulations.
