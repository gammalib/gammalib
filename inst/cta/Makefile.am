## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_srcdir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/cta/include \
	   -I$(top_srcdir)/inst/cta/src

# Define sources for this directory
sources = src/GCTASupport.cpp \
          src/GCTAException.cpp \
          src/GCTAObservation.cpp \
          src/GCTAEventList.cpp \
          src/GCTAEventAtom.cpp \
          src/GCTAEventCube.cpp \
          src/GCTAEventBin.cpp \
          src/GCTAResponse.cpp \
          src/GCTAInstDir.cpp \
          src/GCTARoi.cpp \
          src/GCTAPointing.cpp \
	  src/GCTAModelRadialRegistry.cpp \
          src/GCTAModelRadial.cpp \
          src/GCTAModelRadialGauss.cpp \
          src/GCTAModelRadialAcceptance.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCTAException.hpp \
                     include/GCTAObservation.hpp \
                     include/GCTAEventList.hpp \
                     include/GCTAEventAtom.hpp \
                     include/GCTAEventCube.hpp \
                     include/GCTAEventBin.hpp \
                     include/GCTAPointing.hpp \
                     include/GCTAInstDir.hpp \
                     include/GCTARoi.hpp \
                     include/GCTAResponse.hpp \
		     include/GCTAModelRadialRegistry.hpp \
                     include/GCTAModelRadial.hpp \
                     include/GCTAModelRadialGauss.hpp \
                     include/GCTAModelRadialAcceptance.hpp \
                     include/GCTALib.hpp

# Calibration files to distribute
irfs = caldb/kb_A_50h_v3.dat \
       caldb/kb_B_50h_v3.dat \
       caldb/kb_C_50h_v3.dat \
       caldb/kb_D_50h_v3.dat \
       caldb/kb_E_50h_v3.dat \
       caldb/kb_F_50h_v3.dat \
       caldb/kb_G_50h_v3.dat \
       caldb/kb_H_50h_v3.dat \
       caldb/kb_I_50h_v3.dat \
       caldb/kb_J_50h_v3.dat \
       caldb/kb_K_50h_v3.dat \
       caldb/kb_NA_50h_v3.dat \
       caldb/kb_NB_50h_v3.dat

# Build libtool library
noinst_LTLIBRARIES = libcta.la
libcta_la_SOURCES = $(sources)
libcta_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs)
caldbdir = $(datadir)/caldb/cta
caldb_DATA = $(irfs)
