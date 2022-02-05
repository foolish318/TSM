# TSM
Computational requirements
 The code was last run with R version 3.5.0.
 The following R packages are required: KFAS, MCMCpack, mvtnorm, Rcpp, RcppArmadillo,
numDeriv, sandwich, xtable, urca, dynlm, dplyr, VAR.etp.
{ The le setup.R lists the dependencies and installs all missing packages (latest version).
{ To use the package Rcpp a C++ compiler is needed. On Windows it is necessary to
install Rtools. See here for the Rcpp documentation and here for unocial Rtools
Windows installation instructions.
 The computational costs of the estimation of the ESE model are relatively high, but they
can be substantially reduced if the MCMC sampling is run in parallel.
{ Running the MCMC sampler for a single chain (100,000 iterations) takes about 10-20
hours on a typical modern computer, depending on its processing power. The results
in the paper are based on the combination of ve separate chains.
{ To keep the run time manageable, the sampling the MCMC chains can be done with
parallel computation on multiple cores (see instructions below).
{ The parallel MCMC chains were run on a Linux high-performance computing (HPC)
cluster, specically on a 16-core Intel server with 300 GB of RAM.
