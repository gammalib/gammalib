#############################################################################
# Makefile.am for spi module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2020 by Juergen Knoedlseder                                 #
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
           -I$(top_srcdir)/inst/spi/include \
           -I$(top_srcdir)/inst/spi/src

# Define sources for this directory
sources = src/GSPITools.cpp \
          src/GSPIObservation.cpp \
          src/GSPIEventBin.cpp \
          src/GSPIEventCube.cpp \
          src/GSPIInstDir.cpp \
          src/GSPIResponse.cpp \
          src/GSPIModelDataSpace.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GSPITools.hpp \
                     include/GSPIObservation.hpp \
                     include/GSPIEventBin.hpp \
                     include/GSPIEventCube.hpp \
                     include/GSPIInstDir.hpp \
                     include/GSPIResponse.hpp \
                     include/GSPIModelDataSpace.hpp

# Build libtool library
noinst_LTLIBRARIES = libspi.la
libspi_la_SOURCES = $(sources)
libspi_la_LDFLAGS =

# Calibration files to distribute
# TODO: Follow the example below to distribution your calibration data
#caldb_indx = caldb/data/cgro/comptel/caldb.indx
#caldb_irfs = caldb/data/cgro/comptel/bcf/m10371_sda.fits.gz \
#             caldb/data/cgro/comptel/bcf/m10564_sdb.fits.gz
#
# Distribute calibration database files
#caldbindxdir   = $(datadir)/caldb/data/cgro/comptel
#caldbindx_DATA = $(caldb_indx)
#caldbirfsdir   = $(datadir)/caldb/data/cgro/comptel/bcf
#caldbirfs_DATA = $(caldb_irfs)

# Test data to distribute
testdata            = test/data/obs.xml \
                      test/data/models.xml \
                      test/data/models_noinit.xml
testdatadir         = $(pyexecdir)/gammalib/tests/spi/data
testdata_DATA       = $(testdata)
#
testdataobs         = test/data/obs/og_spi.fits.gz \
                      test/data/obs/spi_bvar_gedsat_se.fits.gz \
                      test/data/obs/spi_pnt_filter.fits.gz
testdataobsdir      = $(pyexecdir)/gammalib/tests/spi/data/obs
testdataobs_DATA    = $(testdataobs)
#
testdataobsspi      = test/data/obs/spi/back_model.fits.gz \
                      test/data/obs/spi/back_model_idx.fits.gz \
                      test/data/obs/spi/dead_time.fits.gz \
                      test/data/obs/spi/energy_boundaries.fits.gz \
                      test/data/obs/spi/evts_det_spec.fits.gz \
                      test/data/obs/spi/gti.fits.gz \
                      test/data/obs/spi/pointing.fits.gz
testdataobsspidir   = $(pyexecdir)/gammalib/tests/spi/data/obs/spi
testdataobsspi_DATA = $(testdataobsspi)

# Extra files to be distributed in the tar file
#EXTRA_DIST = README.md caldb pyext \
#             include/GSPILib.hpp \
#             test/test_SPI.hpp test/test_SPI.py test/data
EXTRA_DIST = README.md pyext \
             include/GSPILib.hpp \
             test/test_SPI.hpp test/test_SPI.py test/data
