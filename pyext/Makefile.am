## Process this file with automake to produce Makefile.in

if HAVE_SWIG
if HAVE_PYTHON
if ENABLE_PYTHON_BINDING

SWIG_INTERFACE = gammalib.i
BUILT_SOURCES = gammalib_wrap.cpp gammalib.py
CLEANFILES = gammalib_wrap.cpp gammalib.py
#EXTRA_DIST =	runme.py

all-local: gammalib_wrap.cpp gammalib.py
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
gammalib_wrap.cpp gammalib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I$(top_srcdir)/src -o gammalib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
	@if test "x`uname`" = "xCygwin"; then \
		echo "WARNING: working around bug in SWIG on Windows/Cygwin platforms ONLY"; \
		sed -i.bak -f $(srcdir)/fix_wrapper.sed gammalib_wrap.cpp; \
	fi

endif # ENABLE_PYTHON_BINDING
endif # HAVE_PYTHON
endif # HAVE_SWIG
