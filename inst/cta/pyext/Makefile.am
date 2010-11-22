## Process this file with automake to produce Makefile.in

if ENABLE_PYTHON_BINDING

SWIG_INTERFACE = ctalib.i
BUILT_SOURCES = ctalib_wrap.cpp ctalib.py
CLEANFILES = ctalib_wrap.cpp ctalib.py

all-local: ctalib_wrap.cpp ctalib.py
	$(PYTHON) setup.py build --build-platlib build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(prefix)

clean-local:
	rm -rf build
	
ctalib_wrap.cpp ctalib.py: $(SWIG_INTERFACE)
	$(SWIG) -c++ -python -Wall -includeall -I../../../pyext -o ctalib_wrap.cpp -outdir . $(srcdir)/$(SWIG_INTERFACE)
	@if test "x`uname`" = "xCygwin"; then \
		echo "WARNING: working around bug in SWIG on Windows/Cygwin platforms ONLY"; \
		sed -i.bak -f $(srcdir)/fix_wrapper.sed ctalib_wrap.cpp; \
	fi

endif # ENABLE_PYTHON_BINDING
