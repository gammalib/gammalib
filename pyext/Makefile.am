#############################################################################
# Makefile for gammalib Python extension module                             #
# ------------------------------------------------------------------------- #
# Copyright (C) 2011-2017 Juergen Knoedlseder                               #
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

# Set wrapper directory depending on the whether the wrappers were actually
# shipped with the distribution of whether the wrappers have to be build
# using swig
if BUILD_WRAPPERS
  wrapperdir = $(top_builddir)/pyext
else
  wrapperdir = $(srcdir)
endif

# Set instrument specific SWIG targets
if WITH_INST_MWL
  MWL_SWIG_TARGETS = $(wrapperdir)/gammalib/mwl_wrap.cpp \
                     $(wrapperdir)/gammalib/mwl.py
  MWL_TESTS        = $(top_srcdir)/inst/mwl/test/test_MWL.py
endif
if WITH_INST_CTA
  CTA_SWIG_TARGETS = $(wrapperdir)/gammalib/cta_wrap.cpp \
                     $(wrapperdir)/gammalib/cta.py
  CTA_TESTS        = $(top_srcdir)/inst/cta/test/test_CTA.py
endif
if WITH_INST_LAT
  LAT_SWIG_TARGETS = $(wrapperdir)/gammalib/lat_wrap.cpp \
                     $(wrapperdir)/gammalib/lat.py
  LAT_TESTS        = $(top_srcdir)/inst/lat/test/test_LAT.py
endif
if WITH_INST_COM
  COM_SWIG_TARGETS = $(wrapperdir)/gammalib/com_wrap.cpp \
                     $(wrapperdir)/gammalib/com.py
  COM_TESTS        = $(top_srcdir)/inst/com/test/test_COM.py
endif
INST_SWIG_TARGETS = $(MWL_SWIG_TARGETS) \
                    $(CTA_SWIG_TARGETS) \
                    $(LAT_SWIG_TARGETS) \
                    $(COM_SWIG_TARGETS)

# Set SWIG targets
PYTHON_WRAPPERS = $(wrapperdir)/gammalib/app_wrap.cpp \
                  $(wrapperdir)/gammalib/app.py \
                  $(wrapperdir)/gammalib/base_wrap.cpp \
                  $(wrapperdir)/gammalib/base.py \
                  $(wrapperdir)/gammalib/fits_wrap.cpp \
                  $(wrapperdir)/gammalib/fits.py \
                  $(wrapperdir)/gammalib/linalg_wrap.cpp \
                  $(wrapperdir)/gammalib/linalg.py \
                  $(wrapperdir)/gammalib/model_wrap.cpp \
                  $(wrapperdir)/gammalib/model.py \
                  $(wrapperdir)/gammalib/numerics_wrap.cpp \
                  $(wrapperdir)/gammalib/numerics.py \
                  $(wrapperdir)/gammalib/obs_wrap.cpp \
                  $(wrapperdir)/gammalib/obs.py \
                  $(wrapperdir)/gammalib/opt_wrap.cpp \
                  $(wrapperdir)/gammalib/opt.py \
                  $(wrapperdir)/gammalib/sky_wrap.cpp \
                  $(wrapperdir)/gammalib/sky.py \
                  $(wrapperdir)/gammalib/support_wrap.cpp \
                  $(wrapperdir)/gammalib/support.py \
                  $(wrapperdir)/gammalib/test_wrap.cpp \
                  $(wrapperdir)/gammalib/test.py \
                  $(wrapperdir)/gammalib/xml_wrap.cpp \
                  $(wrapperdir)/gammalib/xml.py \
                  $(wrapperdir)/gammalib/xspec_wrap.cpp \
                  $(wrapperdir)/gammalib/xspec.py \
                  $(wrapperdir)/gammalib/vo_wrap.cpp \
                  $(wrapperdir)/gammalib/vo.py \
                  $(INST_SWIG_TARGETS)

# Test scripts to be distributed in $(pyexecdir)/gammalib/tests
testscripts = gammalib/tests/__init__.py \
              $(top_srcdir)/test/test_python.py \
              $(top_srcdir)/test/test_GApplication.py \
              $(top_srcdir)/test/test_GFits.py \
              $(top_srcdir)/test/test_GLinalg.py \
              $(top_srcdir)/test/test_GModel.py \
              $(top_srcdir)/test/test_GNumerics.py \
              $(top_srcdir)/test/test_GObservation.py \
              $(top_srcdir)/test/test_GOptimizer.py \
              $(top_srcdir)/test/test_GSky.py \
              $(top_srcdir)/test/test_GSupport.py \
              $(top_srcdir)/test/test_GTest.py \
              $(top_srcdir)/test/test_GVO.py \
              $(top_srcdir)/test/test_GXml.py \
              $(top_srcdir)/test/test_GXspec.py \
              $(MWL_TESTS) $(CTA_TESTS) $(LAT_TESTS) $(COM_TESTS)

# Test data to be distributed in $(pyexecdir)/gammalib/tests/data
testdata = $(top_srcdir)/test/data/initial_coords.txt \
           $(top_srcdir)/test/data/fk5_to_galactic.txt \
           $(top_srcdir)/test/data/file.fits.gz \
           $(top_srcdir)/test/data/test_GApplication.par \
           $(top_srcdir)/test/data/model_point_plaw.xml

# Distribute scripts
testscriptsdir   = $(pyexecdir)/gammalib/tests
testscripts_DATA = $(testscripts)
testdatadir      = $(pyexecdir)/gammalib/tests/data
testdata_DATA    = $(testdata)

# Distribute CTA test data (optionally)
if WITH_INST_CTA
  testdatacta       = $(top_srcdir)/inst/cta/test/data/irf_1dc.xml \
                      $(top_srcdir)/inst/cta/test/data/irf_unbinned.xml \
                      $(top_srcdir)/inst/cta/test/data/crab_irf.xml \
                      $(top_srcdir)/inst/cta/test/data/onoff_model.xml \
                      $(top_srcdir)/inst/cta/test/data/crab_events.fits.gz \
                      $(top_srcdir)/inst/cta/test/data/crab_cntmap.fits.gz
  testdatactadir    = $(pyexecdir)/gammalib/tests/cta/data
  testdatacta_DATA  = $(testdatacta)
  testcaldbcta      = $(top_srcdir)/inst/cta/test/caldb/cta_dummy_irf.dat \
                      $(top_srcdir)/inst/cta/test/caldb/prod1_gauss.fits.gz \
                      $(top_srcdir)/inst/cta/test/caldb/prod1_king.fits.gz
  testcaldbctadir   = $(pyexecdir)/gammalib/tests/cta/caldb
  testcaldbcta_DATA = $(testcaldbcta)
  testdc1cta        = $(top_srcdir)/inst/cta/test/caldb/dc1/arf.fits.gz \
                      $(top_srcdir)/inst/cta/test/caldb/dc1/psf.fits.gz \
                      $(top_srcdir)/inst/cta/test/caldb/dc1/psf_magic.fits.gz \
                      $(top_srcdir)/inst/cta/test/caldb/dc1/rmf.fits.gz
  testdc1ctadir     = $(pyexecdir)/gammalib/tests/cta/caldb/dc1
  testdc1cta_DATA   = $(testdc1cta)
endif

# If building of Python wrappers was requested then add them to the build
# source
if BUILD_WRAPPERS

# Set SWIG compile command
SWIGCOMPILE = $(SWIG) -c++ -python -Wall

# Sources that are required for the build. All targets in the list will be
# built before anything. These are honored by 'make all', 'make check' and
# 'make install'. We use this here to make sure that all these targets are
# build by swig before we start compiling the Python modules.
BUILT_SOURCES = $(PYTHON_WRAPPERS)

# Add the wrappers to the files that will be deleted in a 'make clean'
CLEAN_FILES = $(BUILT_SOURCES)

# Rule to build SWIG wrappers
swig: $(BUILT_SOURCES)

# General rule for core modules
%_wrap.cpp %.py: %.i
	if $(SWIGCOMPILE) -MMD -MF "$*.Tpi" -I$(top_srcdir)/pyext -o $@ -outdir gammalib $<; \
	then mv -f "$*.Tpi" "$*.Pi"; else rm -f "$*.Tpi"; exit 1; fi

# Rule for MWL module
if WITH_INST_MWL
$(wrapperdir)/gammalib/mwl.py: $(wrapperdir)/gammalib/mwl_wrap.cpp
$(wrapperdir)/gammalib/mwl_wrap.cpp: $(top_srcdir)/inst/mwl/pyext/mwl.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/mwl.Tpi" -I$(top_srcdir)/pyext -o gammalib/mwl_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/mwl.Tpi" "gammalib/mwl.Pi"; else rm -f "gammalib/mwl.Tpi"; exit 1; fi
endif

# Rule for CTA module
if WITH_INST_CTA
$(wrapperdir)/gammalib/cta.py: $(wrapperdir)/gammalib/cta_wrap.cpp
$(wrapperdir)/gammalib/cta_wrap.cpp: $(top_srcdir)/inst/cta/pyext/cta.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/cta.Tpi" -I$(top_srcdir)/pyext -o gammalib/cta_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/cta.Tpi" "gammalib/cta.Pi"; else rm -f "gammalib/cta.Tpi"; exit 1; fi
endif

# Rule for LAT module
if WITH_INST_LAT
$(wrapperdir)/gammalib/lat.py: $(wrapperdir)/gammalib/lat_wrap.cpp
$(wrapperdir)/gammalib/lat_wrap.cpp: $(top_srcdir)/inst/lat/pyext/lat.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/lat.Tpi" -I$(top_srcdir)/pyext -o gammalib/lat_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/lat.Tpi" "gammalib/lat.Pi"; else rm -f "gammalib/lat.Tpi"; exit 1; fi
endif

# Rule for COM module
if WITH_INST_COM
$(wrapperdir)/gammalib/com.py: $(wrapperdir)/gammalib/com_wrap.cpp
$(wrapperdir)/gammalib/com_wrap.cpp: $(top_srcdir)/inst/com/pyext/com.i
	if $(SWIGCOMPILE) -MMD -MF "gammalib/com.Tpi" -I$(top_srcdir)/pyext -o gammalib/com_wrap.cpp -outdir gammalib $<; \
	then mv -f "gammalib/com.Tpi" "gammalib/com.Pi"; else rm -f "gammalib/com.Tpi"; exit 1; fi
endif

# ... otherwise, if building of Python wrappers has not been required then
# they are located in the source directory and we don't have to build them
# using swig
else
CLEAN_FILES =
BUILT_SOURCES =
endif # BUILD_WRAPPERS

# Extra files to be distributed in the tar file
EXTRA_DIST = gammalib/tests $(PYTHON_WRAPPERS)
dist-hook:
	cp -r $(srcdir)/*.i $(distdir)
	cp -r $(srcdir)/gammalib/*.i $(distdir)/gammalib

# Dependencies. Obtain the list of dependencies by executing ls.
DEPFILES = $(shell ls gammalib/*.Pi 2> /dev/null)
include $(DEPFILES)

# Things we want to do when 'make' is invoked.
all-local: build

# Build the gammalib extension module
build: $(BUILT_SOURCES)
	@PYTHON_BUILD_PREFIX@ $(PYTHON) setup.py build_ext
if IS_ELCAPITAN
	-for f in build/gammalib/_*.so; do \
		echo "Build post process Python module "$$f; \
		libgamma=`otool -L $$f | grep -o libgamma.*.dylib`; \
		install_name_tool -change $(libdir)/$$libgamma @rpath/$$libgamma $$f; \
		install_name_tool -add_rpath $(abs_top_builddir)/src/.libs $$f; \
	done
endif

# Do nothing for 'make check'
#check-local:

# Install the gammalib Python module
install-exec-local:
if IS_ELCAPITAN
	-for f in build/gammalib/_*.so; do \
		echo "Install pre process Python module "$$f; \
		install_name_tool -rpath $(abs_top_builddir)/src/.libs $(libdir) $$f; \
	done
endif
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

# Uninstall the gammalib Python module
uninstall-local:
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/*.py
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/*.so
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/tests/*.pyc
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib/tests/__pycache__
	rm -rf $(DESTDIR)$(pyexecdir)/gammalib-@VERSION@*.egg-info
	rm -rf $(DESTDIR)$(datarootdir)/gammalib

# Set files to be removed in 'make clean'
CLEANFILES = $(CLEAN_FILES) $(DEPFILES)

# Files to be removed in build directory. Only remove wrapper files if the
# have been built or if the were copied in the build directory
clean-local:
	rm -rf build
if BUILD_WRAPPERS
	rm -rf $(top_builddir)/pyext/gammalib/*.py
	rm -rf $(top_builddir)/pyext/gammalib/*_wrap.cpp
else
	if test "x$(top_builddir)" != "x$(top_srcdir)"; then \
		rm -rf $(top_builddir)/pyext/gammalib/*.py; \
		rm -rf $(top_builddir)/pyext/gammalib/*_wrap.cpp; \
	fi
endif
	rm -rf $(top_builddir)/pyext/gammalib/*.pyc
	rm -rf $(top_builddir)/pyext/gammalib/_*.so
	rm -rf $(top_builddir)/pyext/gammalib/__pycache__

endif # WITH_PYTHON
