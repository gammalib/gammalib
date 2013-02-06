## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/spi/include \
           -I$(top_srcdir)/inst/spi/src

# Define sources for this directory
sources = src/GSPIException.cpp \
          src/GSPIEventCube.cpp \
          src/GSPIEventBin.cpp \
          src/GSPIInstDir.cpp \
          src/GSPIPointing.cpp \
          src/GSPIResponse.cpp \
          src/GSPIObservation.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GSPIException.hpp \
                     include/GSPIEventCube.hpp \
                     include/GSPIEventBin.hpp \
                     include/GSPIInstDir.hpp \
                     include/GSPIPointing.hpp \
                     include/GSPIResponse.hpp \
                     include/GSPIObservation.hpp \
                     include/GSPILib.hpp

# Calibration files to distribute
irfs = 

# Build libtool library
noinst_LTLIBRARIES = libspi.la
libspi_la_SOURCES = $(sources)
libspi_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs)
caldbdir   = $(datadir)/caldb/spi
caldb_DATA = $(irfs)
