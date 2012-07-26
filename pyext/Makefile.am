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
if WITH_INST_COM
  COM_SWIG_TARGETS = gammalib/com_wrap.cpp gammalib/com.py
endif
INST_SWIG_TARGETS = $(MWL_SWIG_TARGETS) \
                    $(CTA_SWIG_TARGETS) \
                    $(COM_SWIG_TARGETS)

# Sources that are required for the build. All targets in the list will be built
# before anything. These are honored by 'make all', 'make check' and 'make install'.
# We use this here to make sure that all these targets are build by swig before we
# start compiling the Python modules.
BUILT_SOURCES = gammalib/app_wrap.cpp  gammalib/app.py \
                gammalib/fits_wrap.cpp gammalib/fits.py \
                gammalib/linalg_wrap.cpp gammalib/linalg.py \
                gammalib/model_wrap.cpp gammalib/model.py \
                gammalib/numerics_wrap.cpp gammalib/numerics.py \
                gammalib/obs_wrap.cpp gammalib/obs.py \
                gammalib/opt_wrap.cpp gammalib/opt.py \
                gammalib/sky_wrap.cpp gammalib/sky.py \
                gammalib/support_wrap.cpp gammalib/support.py \
                gammalib/xml_wrap.cpp gammalib/xml.py \
                $(INST_SWIG_TARGETS)

# Set files to be removed in clean. We remove only the SWIG interface files if SWIG
# is available to rebuild them. 
if HAVE_SWIG
  CLEAN_FILES = $(BUILT_SOURCES)
endif
CLEANFILES = $(CLEAN_FILES)

# Things we want to do when 'make' is invoked.
all-local: build

# Build the gammalib extension module
build: $(BUILT_SOURCES)
	$(PYTHON) setup.py build_ext

# If 'make check' is invoked we build the extension modules inplace, so that
# putting the 'pyext' directory in the PYTHONPATH allows importing the gammalib
# module. We also catch this special build in the setup script to force static
# linking with libgamma, which will allow a successful 'make check' on systems
# (such as Mac OS X) that have fixed shared library paths.
check-local:
	$(PYTHON) setup.py build_ext --inplace

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
	$(SWIG) -c++ -python -Wall -o $@ -outdir gammalib $<

# Rule for MWL module
if WITH_INST_MWL
gammalib/mwl_wrap.cpp gammalib/mwl.py: $(top_srcdir)/inst/mwl/pyext/mwl.i
	$(SWIG) -c++ -python -Wall -o gammalib/mwl_wrap.cpp -outdir gammalib $<
endif

# Rule for CTA module
if WITH_INST_CTA
gammalib/cta_wrap.cpp gammalib/cta.py: $(top_srcdir)/inst/cta/pyext/cta.i
	$(SWIG) -c++ -python -Wall -o gammalib/cta_wrap.cpp -outdir gammalib $<
endif

# Rule for LAT module
if WITH_INST_LAT
gammalib/lat_wrap.cpp gammalib/lat.py: $(top_srcdir)/inst/lat/pyext/lat.i
	$(SWIG) -c++ -python -Wall -o gammalib/lat_wrap.cpp -outdir gammalib $<
endif

# Rule for COM module
if WITH_INST_COM
gammalib/com_wrap.cpp gammalib/com.py: $(top_srcdir)/inst/com/pyext/com.i
	$(SWIG) -c++ -python -Wall -o gammalib/com_wrap.cpp -outdir gammalib $<
endif

endif # HAVE_SWIG

endif # WITH_PYTHON
