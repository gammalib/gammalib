## Process this file with automake to produce Makefile.in

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

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
          src/GCTAOnOffObservation.cpp \
          src/GCTAOnOffObservations.cpp \
          src/GCTAEventList.cpp \
          src/GCTAEventAtom.cpp \
          src/GCTAEventCube.cpp \
          src/GCTAEventBin.cpp \
          src/GCTAResponse.cpp \
          src/GCTAResponseIrf.cpp \
          src/GCTAResponseCube.cpp \
          src/GCTAResponse_helpers.cpp \
          src/GCTAResponseTable.cpp \
          src/GCTAAeff.cpp \
          src/GCTAAeffPerfTable.cpp \
          src/GCTAAeffArf.cpp \
          src/GCTAAeff2D.cpp \
          src/GCTAPsf.cpp \
          src/GCTAPsfPerfTable.cpp \
          src/GCTAPsfVector.cpp \
          src/GCTAPsf2D.cpp \
          src/GCTAPsfKing.cpp \
          src/GCTAEdisp.cpp \
          src/GCTAEdispPerfTable.cpp \
          src/GCTAEdispRmf.cpp \
          src/GCTAEdisp2D.cpp \
          src/GCTABackground.cpp \
          src/GCTABackgroundPerfTable.cpp \
          src/GCTABackground3D.cpp \
          src/GCTACubeExposure.cpp \
          src/GCTACubeBackground.cpp \
          src/GCTACubePsf.cpp \
          src/GCTACubeSource.cpp \
          src/GCTACubeSourcePoint.cpp \
          src/GCTACubeSourceDiffuse.cpp \
          src/GCTAInstDir.cpp \
          src/GCTARoi.cpp \
          src/GCTAPointing.cpp \
          src/GCTAModelCubeBackground.cpp \
          src/GCTAModelIrfBackground.cpp \
          src/GCTAModelAeffBackground.cpp \
          src/GCTAModelRadialRegistry.cpp \
          src/GCTAModelRadial.cpp \
          src/GCTAModelRadialGauss.cpp \
          src/GCTAModelRadialPolynom.cpp \
          src/GCTAModelRadialProfile.cpp \
          src/GCTAModelRadialAcceptance.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCTATypemaps.hpp \
                     include/GCTAException.hpp \
                     include/GCTAObservation.hpp \
                     include/GCTAOnOffObservation.hpp \
                     include/GCTAOnOffObservations.hpp \
                     include/GCTAEventList.hpp \
                     include/GCTAEventAtom.hpp \
                     include/GCTAEventCube.hpp \
                     include/GCTAEventBin.hpp \
                     include/GCTAPointing.hpp \
                     include/GCTAInstDir.hpp \
                     include/GCTARoi.hpp \
                     include/GCTAResponse.hpp \
                     include/GCTAResponseIrf.hpp \
                     include/GCTAResponseCube.hpp \
                     include/GCTAResponseTable.hpp \
                     include/GCTAAeff.hpp \
                     include/GCTAAeffPerfTable.hpp \
                     include/GCTAAeffArf.hpp \
                     include/GCTAAeff2D.hpp \
                     include/GCTAPsf.hpp \
                     include/GCTAPsfPerfTable.hpp \
                     include/GCTAPsfVector.hpp \
                     include/GCTAPsf2D.hpp \
                     include/GCTAPsfKing.hpp \
                     include/GCTAEdisp.hpp \
                     include/GCTAEdispPerfTable.hpp \
                     include/GCTAEdispRmf.hpp \
                     include/GCTAEdisp2D.hpp \
                     include/GCTABackground.hpp \
                     include/GCTABackgroundPerfTable.hpp \
                     include/GCTABackground3D.hpp \
                     include/GCTACubeExposure.hpp \
                     include/GCTACubePsf.hpp \
                     include/GCTACubeBackground.hpp \
                     include/GCTACubeSource.hpp \
                     include/GCTACubeSourcePoint.hpp \
                     include/GCTACubeSourceDiffuse.hpp \
                     include/GCTAModelCubeBackground.hpp \
                     include/GCTAModelIrfBackground.hpp \
                     include/GCTAModelAeffBackground.hpp \
                     include/GCTAModelRadialRegistry.hpp \
                     include/GCTAModelRadial.hpp \
                     include/GCTAModelRadialGauss.hpp \
                     include/GCTAModelRadialPolynom.hpp \
                     include/GCTAModelRadialProfile.hpp \
                     include/GCTAModelRadialAcceptance.hpp \
                     include/GCTALib.hpp

# Build libtool library
noinst_LTLIBRARIES = libcta.la
libcta_la_SOURCES = $(sources)
libcta_la_LDFLAGS =

# Distribute calibration database files
caldb_prod2_indxdir      = $(datadir)/caldb/data/cta/prod2
caldb_prod2_indx_DATA    = caldb/data/cta/prod2/caldb.indx
caldb_prod2_north1dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_0.5h
caldb_prod2_north1_DATA  = caldb/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz
caldb_prod2_north2dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_5h
caldb_prod2_north2_DATA  = caldb/data/cta/prod2/bcf/North_5h/irf_file.fits.gz
caldb_prod2_north3dir    = $(datadir)/caldb/data/cta/prod2/bcf/North_50h
caldb_prod2_north3_DATA  = caldb/data/cta/prod2/bcf/North_50h/irf_file.fits.gz
#
caldb_prod2_south1dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_0.5h
caldb_prod2_south1_DATA  = caldb/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz
caldb_prod2_south2dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_5h
caldb_prod2_south2_DATA  = caldb/data/cta/prod2/bcf/South_5h/irf_file.fits.gz
caldb_prod2_south3dir    = $(datadir)/caldb/data/cta/prod2/bcf/South_50h
caldb_prod2_south3_DATA  = caldb/data/cta/prod2/bcf/South_50h/irf_file.fits.gz
