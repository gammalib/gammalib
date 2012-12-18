## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Define sources for this directory
sources = GIntegral.cpp \
          GDerivative.cpp \
          GFunction.cpp \
          GNumerics.cpp \
          GException_numerics.cpp

# Build libtool library
noinst_LTLIBRARIES = libnumerics.la
libnumerics_la_SOURCES = $(sources)
libnumerics_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
