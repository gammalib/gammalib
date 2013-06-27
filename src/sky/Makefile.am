## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GSkyDir.cpp \
          GSkyPixel.cpp \
          GSkymap.cpp \
          GSkyRegions.cpp \
          GSkyRegion.cpp \
          GSkyRegionCircle.cpp \
          GWcs.cpp \
          GWcsRegistry.cpp \
          GWcslib.cpp \
          GWcsAZP.cpp \
          GWcsCAR.cpp \
          GWcsMER.cpp \
          GWcsTAN.cpp \
          GWcsSTG.cpp \
          GWcsHPX.cpp \
          GException_sky.cpp

# Build libtool library
noinst_LTLIBRARIES = libsky.la
libsky_la_SOURCES = $(sources)
libsky_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
