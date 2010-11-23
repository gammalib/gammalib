## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/lat/include

# Define sources for this directory
sources = src/GLATObservation.cpp \
          src/GLATEventList.cpp \
          src/GLATEventCube.cpp \
          src/GLATEventAtom.cpp \
          src/GLATEventBin.cpp \
          src/GLATResponse.cpp \
          src/GLATInstDir.cpp \
          src/GLATRoi.cpp \
          src/GLATPointing.cpp \
          src/GLATResponseTable.cpp \
          src/GLATAeff.cpp \
          src/GLATPsf.cpp \
          src/GLATEdisp.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GLATObservation.hpp \
                     include/GLATEventList.hpp \
                     include/GLATEventCube.hpp \
                     include/GLATEventAtom.hpp \
                     include/GLATEventBin.hpp \
                     include/GLATPointing.hpp \
                     include/GLATInstDir.hpp \
                     include/GLATRoi.hpp \
                     include/GLATResponse.hpp \
                     include/GLATResponseTable.hpp \
                     include/GLATLib.hpp

# Build libtool library
noinst_LTLIBRARIES = liblat.la
liblat_la_SOURCES = $(sources)
liblat_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
