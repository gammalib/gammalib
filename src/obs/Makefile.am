## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GEnergy.cpp \
          GTimes.cpp \
          GTime.cpp \
          GTimeReference.cpp \
          GCaldb.cpp \
          GObservations.cpp \
          GObservations_iterator.cpp \
          GObservations_optimizer.cpp \
          GObservation.cpp \
          GObservationRegistry.cpp \
          GEvents.cpp \
          GEventList.cpp \
          GEventCube.cpp \
          GEvent.cpp \
          GEventAtom.cpp \
          GEventBin.cpp \
          GGti.cpp \
          GRoi.cpp \
          GEbounds.cpp \
          GResponse.cpp \
          GInstDir.cpp \
          GPointing.cpp \
          GPhotons.cpp \
          GPhoton.cpp \
          GSource.cpp \
          GException_obs.cpp

# Build libtool library
noinst_LTLIBRARIES = libobs.la
libobs_la_SOURCES = $(sources)
libobs_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
