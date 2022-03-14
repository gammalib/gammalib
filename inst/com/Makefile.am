#############################################################################
# Makefile.am for com module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2022 Juergen Knoedlseder                               #
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
           -I$(top_srcdir)/inst/com/include \
           -I$(top_srcdir)/inst/com/src

# Define sources for this directory
sources = src/com_helpers_response_vector.cpp \
          src/GCOMTools.cpp \
          src/GCOMSupport.cpp \
          src/GCOMSelection.cpp \
          src/GCOMStatus.cpp \
          src/GCOMDri.cpp \
          src/GCOMOads.cpp \
          src/GCOMOad.cpp \
          src/GCOMBvcs.cpp \
          src/GCOMBvc.cpp \
          src/GCOMTim.cpp \
          src/GCOMRoi.cpp \
          src/GCOMEventAtom.cpp \
          src/GCOMEventList.cpp \
          src/GCOMEventCube.cpp \
          src/GCOMEventBin.cpp \
          src/GCOMInstDir.cpp \
          src/GCOMResponse.cpp \
          src/GCOMD1Response.cpp \
          src/GCOMD2Response.cpp \
          src/GCOMInstChars.cpp \
          src/GCOMIaq.cpp \
          src/GCOMObservation.cpp \
          src/GCOMModelDRBPhibarNodes.cpp \
          src/GCOMModelDRBPhibarBins.cpp \
          src/GCOMModelDRM.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCOMTools.hpp \
                     include/GCOMSelection.hpp \
                     include/GCOMStatus.hpp \
                     include/GCOMDri.hpp \
                     include/GCOMOads.hpp \
                     include/GCOMOad.hpp \
                     include/GCOMBvcs.hpp \
                     include/GCOMBvc.hpp \
                     include/GCOMTim.hpp \
                     include/GCOMRoi.hpp \
                     include/GCOMEventAtom.hpp \
                     include/GCOMEventList.hpp \
                     include/GCOMEventCube.hpp \
                     include/GCOMEventBin.hpp \
                     include/GCOMInstDir.hpp \
                     include/GCOMResponse.hpp \
                     include/GCOMD1Response.hpp \
                     include/GCOMD2Response.hpp \
                     include/GCOMInstChars.hpp \
                     include/GCOMIaq.hpp \
                     include/GCOMObservation.hpp \
                     include/GCOMModelDRBPhibarNodes.hpp \
                     include/GCOMModelDRBPhibarBins.hpp \
                     include/GCOMModelDRM.hpp \
                     include/GCOMLib.hpp

# Calibration files to distribute
caldb_indx = caldb/data/cgro/comptel/caldb.indx
caldb_irfs = caldb/data/cgro/comptel/bcf/module_status.fits.gz \
             caldb/data/cgro/comptel/bcf/m10371_sda.fits.gz \
             caldb/data/cgro/comptel/bcf/m10564_sdb.fits.gz \
             caldb/data/cgro/comptel/bcf/r00004_ict.fits.gz \
             caldb/data/cgro/comptel/bcf/u09517_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47512_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47569_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47741_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47857_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47979_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48091_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48173_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48174_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48175_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48176_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48177_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48178_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48179_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48180_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48181_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48182_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48183_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48184_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48185_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48186_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48187_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48188_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48189_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48190_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48191_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48192_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48193_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48194_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48195_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48196_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48197_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48198_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48199_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/m16845_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/plaw2_0.75-1MeV_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/plaw2_1-3MeV_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/plaw2_10-30MeV_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/plaw2_3-10MeV_iaq.fits.gz

# Distribute calibration database files
caldbindxdir   = $(datadir)/caldb/data/cgro/comptel
caldbindx_DATA = $(caldb_indx)
caldbirfsdir   = $(datadir)/caldb/data/cgro/comptel/bcf
caldbirfs_DATA = $(caldb_irfs)

# Test data to distribute
testdatacom       = test/data/obs.xml \
                    test/data/crab.xml \
                    test/data/m32171_drx.fits.gz \
                    test/data/m34997_drg.fits.gz \
                    test/data/m50438_dre.fits.gz \
                    test/data/m50439_dre.fits.gz \
                    test/data/m50440_dre.fits.gz \
                    test/data/m50441_dre.fits.gz \
                    test/data/m20039_oad.fits.gz \
                    test/data/m20041_oad.fits.gz \
                    test/data/s10150_10000rows_bvc.fits.gz
testdatacomdir    = $(pyexecdir)/gammalib/tests/com/data
testdatacom_DATA  = $(testdatacom)
testcaldbcom      = test/caldb/data/cgro/comptel/caldb.indx
testcaldbcomdir   = $(pyexecdir)/gammalib/tests/com/caldb/data/cgro/comptel
testcaldbcom_DATA = $(testcaldbcom)
testirfcom        = test/caldb/data/cgro/comptel/bcf/m10371_sda.fits.gz \
                    test/caldb/data/cgro/comptel/bcf/m10564_sdb.fits.gz \
                    test/caldb/data/cgro/comptel/bcf/u47569_iaq.fits.gz
testirfcomdir     = $(pyexecdir)/gammalib/tests/com/caldb/data/cgro/comptel/bcf
testirfcom_DATA   = $(testirfcom)

# Build libtool library
noinst_LTLIBRARIES = libcom.la
libcom_la_SOURCES = $(sources)
libcom_la_LDFLAGS =

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md caldb pyext \
             include/GCOMLib.hpp \
             src/GCOMSupport.hpp \
             src/com_helpers_response_vector.hpp \
             test/test_COM.hpp test/test_COM.py test/caldb test/data
