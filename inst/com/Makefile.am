## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/com/include \
           -I$(top_srcdir)/inst/com/src

# Define sources for this directory
sources = src/GCOMSupport.cpp \
          src/GCOMException.cpp \
          src/GCOMEventCube.cpp \
          src/GCOMEventBin.cpp \
          src/GCOMInstDir.cpp \
          src/GCOMPointing.cpp \
          src/GCOMResponse.cpp \
          src/GCOMObservation.cpp \
          src/GCOMModelDRBFitting.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCOMException.hpp \
                     include/GCOMEventCube.hpp \
                     include/GCOMEventBin.hpp \
                     include/GCOMInstDir.hpp \
                     include/GCOMPointing.hpp \
                     include/GCOMResponse.hpp \
                     include/GCOMObservation.hpp \
                     include/GCOMModelDRBFitting.hpp \
                     include/GCOMLib.hpp

# Calibration files to distribute
irfs = caldb/u47512_iaq.fits.gz \
       caldb/u47569_iaq.fits.gz \
       caldb/u47741_iaq.fits.gz \
       caldb/u48091_iaq.fits.gz

# Build libtool library
noinst_LTLIBRARIES = libcom.la
libcom_la_SOURCES = $(sources)
libcom_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs)
caldbdir   = $(datadir)/caldb/com
caldb_DATA = $(irfs)
