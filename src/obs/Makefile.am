## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GData.cpp \
	  GData_iterator.cpp \
	  GData_optimizer.cpp \
	  GObservation.cpp \
	  GEvents.cpp \
	  GEventList.cpp \
	  GEventCube.cpp \
	  GEvent.cpp \
	  GEventAtom.cpp \
	  GEventBin.cpp \
	  GGti.cpp \
	  GEbounds.cpp \
	  GResponse.cpp

# Build libtool library
noinst_LTLIBRARIES = libobs.la
libobs_la_SOURCES = $(sources)
libobs_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
