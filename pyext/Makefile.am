## Process this file with automake to produce Makefile.in
#
# Note that we could also ship gammalib_wrap.cpp with the code
# so that the call to SWIG is not required. We should add some
# test to the master configure.ac script that will test if
# gammalib_wrap.cpp is shipped with the distribution. In this
# case the SWIG call may be omitted and the gammalib_wrap.cpp
# file will be simply compiled.

if ENABLE_PYTHON_BINDING
if HAVE_PYTHON

SWIG_INTERFACE = gammalib.i
BUILT_SOURCES = gammalib_wrap.cpp gammalib.py
#CLEANFILES = gammalib_wrap.cpp gammalib.py
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
	@if test "x`uname`" = "xSunOS"; then \
		echo "WARNING: This is Solaris"; \
		CC=CC; \
	fi
	@if test "x`uname`" = "xDarwin"; then \
		echo "WARNING: This is Mac OS X"; \
	fi
	$(PYTHON) setup.py build --build-platlib build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
if HAVE_SWIG
gammalib_wrap.cpp gammalib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I$(top_srcdir)/src $(INC_INST) $(OPT_INST) -o gammalib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
#	@if test "x`uname`" = "xCygwin"; then \
#		echo "WARNING: working around bug in SWIG on Windows/Cygwin platforms ONLY"; \
#		sed -i.bak -f $(srcdir)/fix_wrapper.sed gammalib_wrap.cpp; \
#	fi
else
gammalib_wrap.cpp gammalib.py: $(SWIG_INTERFACE)
endif # HAVE_SWIG

endif # HAVE_PYTHON
endif # ENABLE_PYTHON_BINDING
