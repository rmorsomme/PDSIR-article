
<!-- README.md is generated from README.Rmd. Please edit that file -->

## PDSIR-article – a repository for reproducibility

This repository provides `R` scripts for reproducing the result 
of the article Exact Bayesian Inference for Partially Observed Stochastic Epidemic Models via Uniformly Ergodic Block Sampling by R. Morsomme and
J. Xu. The `R` package [PDSIR](https://github.com/rmorsomme/PDSIR) 
contains the functions for running the sampling algorithm introduced
in the article.

The scripts are located in the folder `R`. The five scripts
`analysis-*.R` run Markov chain Monte Carlo samplers whose draws are
saved in  `mcmc-draws`. The script `make_figures.R` makes the
various figures, which are saved in `figures`, and computes
the summary statistics found in the article.

For convenience, `.sh` scripts for running the simulations on a cluster via
SLURM are also provided in `slurm`. This is particularly convenient for 
computing the posterior coverage rate across 1,000 simulations, which is
computationally expensive.
