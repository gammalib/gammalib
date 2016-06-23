CTA module unit tests
=====================
This folder contains scripts and executables that are used for CTA
module unit testing.

benchmark_irf_computing.py
  This script performs a benchmark of the response computation for all
  different analysis methods (unbinned, binned, stacked). Specifically it
  compares the total number of modelled counts to the expected number from
  Monte-Carlo simulations. This script is used to validate the absolute
  normalization of the CTA response computations.

benchmark_ml_fitting.py
  This script performs a benchmark for maximum likelihood fitting of CTA
  data.

test_CTA.cpp
  Unit test C++ file for CTA specific GammaLib components.

test_CTA.py
  Unit test Python script for CTA specific GammaLib components.
