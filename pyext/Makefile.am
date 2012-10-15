#############################################################################
# Makefile for gammalib Python extension module                             #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2012 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
# ------------------------------------------------------------------------- #
# This Makefile will build the Python interface files (*.py) and the Python #
# wrapper files using SWIG (if available; otherwise we have to make sure    #
# the the files exist already). The Python distutils are used to build and  #
# install the Python modules (see file setup.py.in).                        #
#                                                                           #
# Process this file with automake to produce Makefile.in.                   #
#############################################################################

# Works only if Python is available
if WITH_PYTHON

# Set instrument specific SWIG targets
if WITH_INST_MWL
  MWL_SWIG_TARGETS = gammalib/mwl_wrap.cpp gammalib/mwl.py
endif
if WITH_INST_CTA
  CTA_SWIG_TARGETS = gammalib/cta_wrap.cpp gammalib/cta.py
endif
if WITH_INST_LAT
  LAT_SWIG_TARGETS = gammalib/lat_wrap.cpp gammalib/lat.py
endif
INST_SWIG_TARGETS = $(MWL_SWIG_TARGETS) $(CTA_SWIG_TARGETS) $(LAT_SWIG_TARGETS)

# Sources that are required for the build. All targets in the list will be built
# before anything. These are honored by 'make all', 'make check' and 'make install'.
# We use this here to make sure that all these targets are build by swig before we
# start compiling the Python modules.
BUILT_SOURCES = gammalib/app_wrap.cpp  gammalib/app.py \
                gammalib/base_wrap.cpp gammalib/base.py \
                gammalib/fits_wrap.cpp gammalib/fits.py \
                gammalib/linalg_wrap.cpp gammalib/linalg.py \
                gammalib/model_wrap.cpp gammalib/model.py \
                gammalib/numerics_wrap.cpp gammalib/numerics.py \
                gammalib/obs_wrap.cpp gammalib/obs.py \
                gammalib/opt_wrap.cpp gammalib/opt.py \
                gammalib/sky_wrap.cpp gammalib/sky.py \
                gammalib/support_wrap.cpp gammalib/support.py \
                gammalib/xml_wrap.cpp gammalib/xml.py \
                gammalib/test_wrap.cpp gammalib/test.py \
                $(INST_SWIG_TARGETS)

# Dependencies. Obtain the list of dependencies by executing ls.
DEPFILES = $(shell ls gammalib/*.Pi 2> /dev/null)
include $(DEPFILES)

# Set files to be removed in clean. We remove only the SWIG interface files if SWIG
# is available to rebuild them. 
if HAVE_SWIG
  CLEAN_FILES = $(BUILT_SOURCES)
endif
CLEANFILES = $(CLEAN_FILES) $(DEPFILES)

# Set SWIG compile command
SWIGCOMPILE = $(SWIG) -c++ -python -Wall

# Things we want to do when 'make' is invoked.
all-local: build

# Build the gammalib extension module
build: $(BUILT_SOURCES)
	@PYTHON_BUILD_PREFIX@ $(PYTHON) setup.py build_ext

# Do nothing for 'make check'
#check-local:

# Install the gammalib extension module
install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

# Remove build directory and test libraries
clean-local:
	rm -rf build
	rm -rf gammalib/__init__.py? gammalib/_*.so gammalib/*.pyc

# Rules for building Python wrappers using SWIG. The rules are only active if
# SWIG is avialable. The is a general rule for core modules and specific rules
# for each instrument module.
if HAVE_SWIG

# Rule to build SWIG wrappers
swig: $(BUILT_SOURCES)

# General rule for core modules
%_wrap.cpp %.py: %.i
	if $(SWIGCOMPILE) -MMD -MF "$*.Tpi" -o $@ -outdir gammalib $<; \
	then mv -f "$*.Tpi" "$*.Pi"; else rm -f "$*.Tpi"; exit 1; fi

# Rule for MWL module
if WITH_INST_MWL
gammalib/mwl.py: gammalib/mwl_wrap.cpp
gammalib/mwl_wrap.cpp: $(top_srcdir)/inst/mwl/pyext/mwl.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/mwl.Tpi" -o gammalib/mwl_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/mwl.Tpi" "gammalib/mwl.Pi"; else rm -f "gammalib/mwl.Tpi"; exit 1; fi
endif

# Rule for CTA module
if WITH_INST_CTA
gammalib/cta.py: gammalib/cta_wrap.cpp
gammalib/cta_wrap.cpp: $(top_srcdir)/inst/cta/pyext/cta.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/cta.Tpi" -o gammalib/cta_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/cta.Tpi" "gammalib/cta.Pi"; else rm -f "gammalib/cta.Tpi"; exit 1; fi
endif

# Rule for LAT module
if WITH_INST_LAT
gammalib/lat.py: gammalib/lat_wrap.cpp
gammalib/lat_wrap.cpp: $(top_srcdir)/inst/lat/pyext/lat.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/lat.Tpi" -o gammalib/lat_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/lat.Tpi" "gammalib/lat.Pi"; else rm -f "gammalib/lat.Tpi"; exit 1; fi
endif

endif # HAVE_SWIG

endif # WITH_PYTHON
