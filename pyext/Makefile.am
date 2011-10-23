## Process this file with automake to produce Makefile.in
#
# The actual implemented logic is to use swig if swig is
# available on the system

if WITH_PYTHON

# Set instrument specific SWIG targets
if WITH_INST_MWL
  MWL_SWIG_TARGETS = $(top_srcdir)/inst/mwl/pyext/mwl_wrap.cpp \
					 $(top_srcdir)/inst/mwl/pyext/mwl.py
endif
if WITH_INST_CTA
  CTA_SWIG_TARGETS = $(top_srcdir)/inst/cta/pyext/cta_wrap.cpp \
					 $(top_srcdir)/inst/cta/pyext/cta.py
endif
if WITH_INST_LAT
  LAT_SWIG_TARGETS = $(top_srcdir)/inst/lat/pyext/lat_wrap.cpp \
					 $(top_srcdir)/inst/lat/pyext/lat.py
endif
INST_SWIG_TARGETS = $(MWL_SWIG_TARGETS) $(CTA_SWIG_TARGETS) $(LAT_SWIG_TARGETS)

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

# Set files to be removed in clean. We remove only the SWIG interface files if swig
# is available to rebuild them. 
if HAVE_SWIG
  CLEAN_FILES = $(BUILT_SOURCES)
endif
CLEANFILES = $(CLEAN_FILES)

all-local: build

# We build the python extension twice since the install step needs the result
# in a specific directory and the test make check has no access to the name
# of this specific directory. So we build one instance for checking, the other
# for installation. If we could get the name built by distutils we could avoid
# this.
#$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build_ext --build-platlib build
build: $(BUILT_SOURCES)
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build_ext

# make check installs the gammalib extension for checking in the install directory
check-local:
	$(PYTHON) setup.py install --home=install

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build install
	
if HAVE_SWIG
%_wrap.cpp %.py: %.i
	$(SWIG) -c++ -python -Wall -o $@ -outdir gammalib $<
endif # HAVE_SWIG

endif # WITH_PYTHON
