## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/com/include \
           -I$(top_srcdir)/inst/com/src

# Define sources for this directory
sources = src/GCOMException.cpp \
          src/GCOMObservation.cpp \
          src/GCOMEventCube.cpp \
          src/GCOMEventBin.cpp \
          src/GCOMResponse.cpp \
          src/GCOMInstDir.cpp \
          src/GCOMPointing.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCOMException.hpp \
                     include/GCOMObservation.hpp \
                     include/GCOMEventCube.hpp \
                     include/GCOMEventBin.hpp \
                     include/GCOMPointing.hpp \
                     include/GCOMInstDir.hpp \
                     include/GCOMResponse.hpp \
                     include/GCOMLib.hpp

# Calibration files to distribute
irfs = 

# Build libtool library
noinst_LTLIBRARIES = libcom.la
libcom_la_SOURCES = $(sources)
libcom_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs)
caldbdir   = $(datadir)/caldb/com
caldb_DATA = $(irfs)
