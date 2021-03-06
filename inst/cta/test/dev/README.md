Developer test scripts
======================
This folder contains scripts that are used for various testing purposes
during code development. None of the scripts actually needs to work, they
may actually be outdated, but could still be useful for further testing
if needed.

`benchmark_irf_computing.py`

>  This script performs a benchmark of the response computation for all
>  different analysis methods (unbinned, binned, stacked). Specifically it
>  compares the total number of modelled counts to the expected number from
>  Monte-Carlo simulations. This script is used to validate the absolute
>  normalization of the CTA response computations.

`benchmark_ml_fitting.py`
>  This script performs a benchmark for maximum likelihood fitting of CTA
>  data.

`example_binned_ml_fit.py`
>  Illustrates a binned maximum likelihood fit of CTA data.
>  This script implements a full analysis chain, and computes
>  also the TS statistics value that estimates the significance
>  of the source detection.

`example_make_model.py`
>  Illustrates how to make a model in counts space for a counts
>  map given a specific model in XML format.

`example_sim_photons.py`
>  Illustrates how photons are simulated from a source model. If
>  matplotlib is installed, the simulated photons are shown as
>  a spectrum.

`show_model_mc.py`
>  For a given model, compares Monte Carlo simulation to model
>  ecaluation.

`test_background_integration.py`
>  Test the event simulation using GCTAModelBackground.

`test_cta_bkgcube.py`
>  Test energy integration of GCTABackground3D.

`test_gauss.py`
>  Test Gaussian.

`test_irf_elliptical.py`
>  Test the response computation for elliptical models.

`test_irf_offset.py`
>  Creates images of source models for various offset angles. This
>  tests the IRF variation with offset angle.

`test_irf_radial.py`
>  Test the response computation for radial models.

`test_irf_trafo.py`
>  Test coordinate transformations for IRF computation.

`test_likelihood_profile.py`
>  Performs likelihood profile computation.

`test_ml_fitting.py`
>  Performs a binned maximum likelihood analysis.

`test_model.py`
>  Creates images of sources models convolved with the CTA instrument
>  response function. Also allows the creation of gradient images.

`test_npred_computation.py`
>  Test the Npred computation of the CTA response.

`test_npred_elliptical.py`
>  Test the Npred computation of the CTA response for elliptical models.

`test_npred_integration.py`
>  Compare numerical to analytical Npred computation for CTA response.

`test_psf_radial.py`
>  Tests the radial response computation.

`test_psf_shape.py`
>  Simulates the PSF distribution.

`test_radial_acceptance.py`
>  Displays the various radial acceptance models.

`test_response_table.py`
>  Display CTA response table.

`test_shell.py`
>  Shows the shell model radial distribution.

`test_sim_edisp.py`
>  Compares energy dispersion simulation to energy dispersion
>  model.

`test_sim_psf.py`
>  Compares Point Spread function simulation to Point Spread function
>  model.
