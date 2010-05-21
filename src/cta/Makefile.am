## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GCTAEventAtom.cpp
#sources = GCTAObservation.cpp \
#          GCTAEventList.cpp \
#          GCTAEventCube.cpp \
#          GCTAEventAtom.cpp \
#          GCTAEventBin.cpp \
#          GCTAResponse.cpp \
#          GCTAResponseTable.cpp \
#          GCTAAeff.cpp \
#          GCTAPsf.cpp \
#          GCTAEdisp.cpp

# Build libtool library
noinst_LTLIBRARIES = libcta.la
libcta_la_SOURCES = $(sources)
libcta_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
