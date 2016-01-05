## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/mwl/include

# Define sources for this directory
sources = src/GMWLException.cpp \
          src/GMWLObservation.cpp \
          src/GMWLSpectrum.cpp \
          src/GMWLDatum.cpp \
          src/GMWLInstDir.cpp \
          src/GMWLResponse.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GMWLException.hpp \
                     include/GMWLObservation.hpp \
                     include/GMWLSpectrum.hpp \
                     include/GMWLDatum.hpp \
                     include/GMWLInstDir.hpp \
                     include/GMWLResponse.hpp

# Build libtool library
noinst_LTLIBRARIES = libmwl.la
libmwl_la_SOURCES = $(sources)
libmwl_la_LDFLAGS =
