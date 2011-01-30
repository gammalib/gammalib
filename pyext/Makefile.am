## Process this file with automake to produce Makefile.in

if ENABLE_PYTHON_BINDING
if HAVE_PYTHON

# Set SWIG interface filename
SWIG_INTERFACE = gammalib.i

# Set files to be built
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

all-local: gammalib_wrap.cpp gammalib.py
	$(PYTHON_BUILD_PREFIX) $(PYTHON) setup.py build --build-platlib build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
if HAVE_SWIG
gammalib_wrap.cpp gammalib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I$(top_srcdir)/src $(INC_INST) $(OPT_INST) -o gammalib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # HAVE_PYTHON
endif # ENABLE_PYTHON_BINDING
