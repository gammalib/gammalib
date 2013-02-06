## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/xxx/include \
           -I$(top_srcdir)/inst/xxx/src

# Define sources for this directory
sources = src/GXXXException.cpp \
          src/GXXXEventCube.cpp \
          src/GXXXEventBin.cpp \
          src/GXXXInstDir.cpp \
          src/GXXXPointing.cpp \
          src/GXXXResponse.cpp \
          src/GXXXObservation.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GXXXException.hpp \
                     include/GXXXEventCube.hpp \
                     include/GXXXEventBin.hpp \
                     include/GXXXInstDir.hpp \
                     include/GXXXPointing.hpp \
                     include/GXXXResponse.hpp \
                     include/GXXXObservation.hpp \
                     include/GXXXLib.hpp

# Calibration files to distribute
irfs = 

# Build libtool library
noinst_LTLIBRARIES = libxxx.la
libxxx_la_SOURCES = $(sources)
libxxx_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs)
caldbdir   = $(datadir)/caldb/xxx
caldb_DATA = $(irfs)
