## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/cta/include

# Define sources for this directory
sources = src/GCTAException.cpp \
          src/GCTAObservation.cpp \
          src/GCTAEventList.cpp \
          src/GCTAEventAtom.cpp \
	  src/GCTAResponse.cpp \
	  src/GCTAInstDir.cpp \
	  src/GCTAPointing.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCTAException.hpp \
                     include/GCTAObservation.hpp \
		     include/GCTAEventList.hpp \
		     include/GCTAEventAtom.hpp \
		     include/GCTAPointing.hpp \
		     include/GCTAInstDir.hpp \
		     include/GCTAResponse.hpp

# Build libtool library
noinst_LTLIBRARIES = libcta.la
libcta_la_SOURCES = $(sources)
libcta_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Programs to run for make check
#check_PROGRAMS = test_CTA

# Tests
#TESTS = test_CTA

#test_CTA_SOURCES = test/test_CTA.cpp
#test_CTA_LDADD = ../../src/libgamma.la
