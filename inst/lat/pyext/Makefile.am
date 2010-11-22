## Process this file with automake to produce Makefile.in

if ENABLE_PYTHON_BINDING

SWIG_INTERFACE = latlib.i
BUILT_SOURCES = latlib_wrap.cpp latlib.py
CLEANFILES = latlib_wrap.cpp latlib.py

all-local: latlib_wrap.cpp latlib.py
	$(PYTHON) setup.py build --build-platlib build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
latlib_wrap.cpp latlib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I../../../pyext -o latlib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
	@if test "x`uname`" = "xCygwin"; then \
		echo "WARNING: working around bug in SWIG on Windows/Cygwin platforms ONLY"; \
		sed -i.bak -f $(srcdir)/fix_wrapper.sed latlib_wrap.cpp; \
	fi

endif # ENABLE_PYTHON_BINDING
