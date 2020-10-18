#############################################################################
# Makefile.am for cta module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2020 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

# Set sub directories
SUBDIRS = test

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
sources = src/cta_helpers_response_stacked_vector.cpp \
          src/GCTASupport.cpp \
          src/GCTAException.cpp \
          src/GCTAObservation.cpp \
          src/GCTAOnOffObservation.cpp \
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
          src/GCTAPsfTable.cpp \
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
          src/GCTACubeEdisp.cpp \
          src/GCTACubeSource.cpp \
          src/GCTAInstDir.cpp \
          src/GCTARoi.cpp \
          src/GCTAPointing.cpp \
          src/GCTAModelBackground.cpp \
          src/GCTAModelSpatial.cpp \
          src/GCTAModelSpatialLookup.cpp \
          src/GCTAModelSpatialGaussSpectrum.cpp \
          src/GCTAModelSpatialGradient.cpp \
          src/GCTAModelSpatialMultiplicative.cpp \
          src/GCTAModelSpatialRegistry.cpp \
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
                     include/GCTAPsfTable.hpp \
                     include/GCTAEdisp.hpp \
                     include/GCTAEdispPerfTable.hpp \
                     include/GCTAEdispRmf.hpp \
                     include/GCTAEdisp2D.hpp \
                     include/GCTABackground.hpp \
                     include/GCTABackgroundPerfTable.hpp \
                     include/GCTABackground3D.hpp \
                     include/GCTACubeExposure.hpp \
                     include/GCTACubePsf.hpp \
                     include/GCTACubeEdisp.hpp \
                     include/GCTACubeBackground.hpp \
                     include/GCTACubeSource.hpp \
                     include/GCTAModelBackground.hpp \
                     include/GCTAModelSpatial.hpp \
                     include/GCTAModelSpatialLookup.hpp \
                     include/GCTAModelSpatialGaussSpectrum.hpp \
                     include/GCTAModelSpatialGradient.hpp \
                     include/GCTAModelSpatialMultiplicative.hpp \
                     include/GCTAModelSpatialRegistry.hpp \
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

# Test data to distribute
testdatacta       = test/data/irf_1dc.xml \
                    test/data/irf_unbinned.xml \
                    test/data/irf_hess_aeff.fits.gz \
                    test/data/irf_hess_bkg.fits.gz \
                    test/data/irf_hess_edisp.fits.gz \
                    test/data/irf_hess_psf.fits.gz \
                    test/data/crab_irf.xml \
                    test/data/onoff_model.xml \
                    test/data/crab_events.fits.gz \
                    test/data/crab_cntmap.fits.gz
testdatactadir    = $(pyexecdir)/gammalib/tests/cta/data
testdatacta_DATA  = $(testdatacta)
testcaldbcta      = test/caldb/cta_dummy_irf.dat \
                    test/caldb/prod1_gauss.fits.gz \
                    test/caldb/prod1_king.fits.gz
testcaldbctadir   = $(pyexecdir)/gammalib/tests/cta/caldb
testcaldbcta_DATA = $(testcaldbcta)
testdc1cta        = test/caldb/dc1/arf.fits.gz \
                    test/caldb/dc1/psf.fits.gz \
                    test/caldb/dc1/psf_magic.fits.gz \
                    test/caldb/dc1/rmf.fits.gz
testdc1ctadir     = $(pyexecdir)/gammalib/tests/cta/caldb/dc1
testdc1cta_DATA   = $(testdc1cta)


# Extra files to be distributed in the tar file
EXTRA_DIST = README.md caldb pyext \
             include/GCTALib.hpp \
             src/cta_helpers_response_stacked_vector.hpp \
             src/GCTAResponse_helpers.hpp src/GCTASupport.hpp \
             test/test_CTA.hpp test/test_CTA.py test/caldb test/data
