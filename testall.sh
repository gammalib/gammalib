#!/bin/sh
#
# Run all test and example scripts that come with the package
# ===========================================================
base=$PWD


# test
# ====
echo
echo "=====> test"
cd test
./example_model_xml_io.py
./example_radial_models.py
cd $base


# inst/cta/test
# =============
echo
echo "=====> inst/cta/test"
cd inst/cta/test
./test_model.py
./test_gauss.py
./example_binned_ml_fit.py
./example_sim_photons.py
cd $base


# inst/lat/test
# =============
echo
echo "=====> inst/lat/test"
cd inst/lat/test
./test_python.py
cd $base


# inst/mwl/test
# =============
echo
echo "=====> inst/mwl/test"
cd inst/mwl/test
./test_python.py
cd $base


# Signal completion
# =================
echo 
echo "All scripts completed."
