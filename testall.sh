#!/bin/sh
#
# Run all test and example scripts that come with the package
#
# Copyright (C) 2011-2016 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =====================================================================
base=$PWD

# test/dev
# ========
echo
echo "=====> test/dev"
cd test/dev
./example_model_xml_io.py
./example_radial_models.py
cd $base


# examples/python
# ===============
echo
echo "=====> examples/python"
cd examples/python
./matrix_howto.py
./models_howto.py
./xml_howto.py
./xml_html_create.py
cd $base


# inst/cta/test/dev
# =================
echo
echo "=====> inst/cta/test/dev"
cd inst/cta/test/dev
./example_binned_ml_fit.py
./example_make_model.py
./example_sim_photons.py
./test_gauss.py
./test_model.py
./test_sim_psf.py
./test_sim_edisp.py
cd $base


# inst/lat/test/dev
# =================
echo
echo "=====> inst/lat/test/dev"
cd inst/lat/test/dev
#./test_python.py
cd $base


# inst/mwl/test/dev
# =================
echo
echo "=====> inst/mwl/test/dev"
cd inst/mwl/test/dev
./test_python.py
cd $base


# Signal completion
# =================
echo 
echo "All scripts completed."
