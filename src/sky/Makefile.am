## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GSkyDir.cpp \
	  GHorizDir.cpp \
          GSkyPixel.cpp \
          GSkyMap.cpp \
          GSkyRegions.cpp \
          GSkyRegion.cpp \
          GSkyRegionCircle.cpp \
          GSkyProjection.cpp \
          GHealpix.cpp \
          GWcsRegistry.cpp \
          GWcs.cpp \
          GWcsAIT.cpp \
          GWcsAZP.cpp \
          GWcsCAR.cpp \
          GWcsMER.cpp \
          GWcsMOL.cpp \
          GWcsSIN.cpp \
          GWcsSTG.cpp \
          GWcsTAN.cpp \
          GException_sky.cpp

# Build libtool library
noinst_LTLIBRARIES = libsky.la
libsky_la_SOURCES = $(sources)
libsky_la_LDFLAGS =
