## Process this file with automake to produce Makefile.in
#
# The actual implemented logic is to use swig if swig is
# available on the system

if ENABLE_PYTHON_BINDING
if HAVE_PYTHON
if WITH_PYTHON_BINDING

# Set SWIG interface filename
SWIG_INTERFACE = gammalib.i

# Define the targets that should be built before anything else.
# These are honored by 'make all', 'make check' and 'make install'.
if HAVE_SWIG
  BUILT_SOURCES = gammalib_wrap.cpp gammalib.py
endif

# Set files to be removed in clean
if HAVE_SWIG
  CLEANFILES = gammalib_wrap.cpp gammalib.py
endif

# Set extra files to distribute
#EXTRA_DIST = runme.py

# Set instrument specific SWIG options (those will be catched in gammalib.i)
if WITH_INST_MWL
  OPT_MWL = -DWITH_INST_MWL
  INC_MWL = -I$(top_srcdir)/inst/mwl/pyext
endif
if WITH_INST_CTA
  OPT_CTA = -DWITH_INST_CTA
  INC_CTA = -I$(top_srcdir)/inst/cta/pyext
endif
if WITH_INST_LAT
  OPT_LAT = -DWITH_INST_LAT
  INC_LAT = -I$(top_srcdir)/inst/lat/pyext
endif
INC_INST = $(INC_MWL) $(INC_CTA) $(INC_LAT)
OPT_INST = $(OPT_MWL) $(OPT_CTA) $(OPT_LAT)

all-local: build

# We build the python extension twice since the install step needs the result
# in a specific directory and the test make check has no access to the name
# of this specific directory. So we build one instance for checking, the other
# for installation. If we could get the name built by distutils we could avoid
# this.
build: gammalib_wrap.cpp gammalib.py
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build --build-platlib build
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
if HAVE_SWIG
gammalib_wrap.cpp gammalib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -I$(top_srcdir)/src $(INC_INST) $(OPT_INST) -o gammalib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # WITH_PYTHON_BINDING
endif # HAVE_PYTHON
endif # ENABLE_PYTHON_BINDING
