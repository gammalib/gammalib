## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/com/include \
           -I$(top_srcdir)/inst/com/src

# Define sources for this directory
sources = src/GCOMSupport.cpp \
          src/GCOMEventCube.cpp \
          src/GCOMEventBin.cpp \
          src/GCOMInstDir.cpp \
          src/GCOMResponse.cpp \
          src/GCOMObservation.cpp \
          src/GCOMModelDRBFitting.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCOMEventCube.hpp \
                     include/GCOMEventBin.hpp \
                     include/GCOMInstDir.hpp \
                     include/GCOMResponse.hpp \
                     include/GCOMObservation.hpp \
                     include/GCOMModelDRBFitting.hpp \
                     include/GCOMLib.hpp

# Calibration files to distribute
caldb_indx = caldb/data/cgro/comptel/caldb.indx
caldb_irfs = caldb/data/cgro/comptel/bcf/u47512_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47569_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47741_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48091_iaq.fits.gz

# Build libtool library
noinst_LTLIBRARIES = libcom.la
libcom_la_SOURCES = $(sources)
libcom_la_LDFLAGS =

# Distribute calibration database files
EXTRA_DIST     = $(caldb_indx) $(caldb_irfs)
caldbindxdir   = $(datadir)/caldb/data/cgro/comptel
caldbindx_DATA = $(caldb_indx)
caldbirfsdir   = $(datadir)/caldb/data/cgro/comptel/bcf
caldbirfs_DATA = $(caldb_irfs)
