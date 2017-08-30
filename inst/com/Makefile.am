#############################################################################
# Makefile.am for com module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2017 Juergen Knoedlseder                               #
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
sources = src/GCOMSupport.cpp \
          src/GCOMEventCube.cpp \
          src/GCOMEventBin.cpp \
          src/GCOMInstDir.cpp \
          src/GCOMResponse.cpp \
          src/GCOMD1Response.cpp \
          src/GCOMD2Response.cpp \
          src/GCOMInstChars.cpp \
          src/GCOMIaq.cpp \
          src/GCOMObservation.cpp \
          src/GCOMModelDRBFitting.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GCOMEventCube.hpp \
                     include/GCOMEventBin.hpp \
                     include/GCOMInstDir.hpp \
                     include/GCOMResponse.hpp \
                     include/GCOMD1Response.hpp \
                     include/GCOMD2Response.hpp \
                     include/GCOMInstChars.hpp \
                     include/GCOMIaq.hpp \
                     include/GCOMObservation.hpp \
                     include/GCOMModelDRBFitting.hpp \
                     include/GCOMLib.hpp

# Calibration files to distribute
caldb_indx = caldb/data/cgro/comptel/caldb.indx
caldb_irfs = caldb/data/cgro/comptel/bcf/m10371_sda.fits.gz \
             caldb/data/cgro/comptel/bcf/m10564_sdb.fits.gz \
             caldb/data/cgro/comptel/bcf/r00004_ict.fits.gz \
             caldb/data/cgro/comptel/bcf/u47512_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47569_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u47741_iaq.fits.gz \
             caldb/data/cgro/comptel/bcf/u48091_iaq.fits.gz

# Build libtool library
noinst_LTLIBRARIES = libcom.la
libcom_la_SOURCES = $(sources)
libcom_la_LDFLAGS =

# Distribute calibration database files
caldbindxdir   = $(datadir)/caldb/data/cgro/comptel
caldbindx_DATA = $(caldb_indx)
caldbirfsdir   = $(datadir)/caldb/data/cgro/comptel/bcf
caldbirfs_DATA = $(caldb_irfs)

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md caldb pyext \
             include/GCOMLib.hpp \
             src/GCOMSupport.hpp \
             test/test_COM.hpp test/test_COM.py test/data
