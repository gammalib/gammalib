#############################################################################
# Makefile.am for mwl module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2021 Juergen Knoedlseder                               #
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
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/mwl/include

# Define sources for this directory
sources = src/GMWLObservation.cpp \
          src/GMWLSpectrum.cpp \
          src/GMWLDatum.cpp \
          src/GMWLInstDir.cpp \
          src/GMWLResponse.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GMWLObservation.hpp \
                     include/GMWLSpectrum.hpp \
                     include/GMWLDatum.hpp \
                     include/GMWLInstDir.hpp \
                     include/GMWLResponse.hpp

# Build libtool library
noinst_LTLIBRARIES = libmwl.la
libmwl_la_SOURCES = $(sources)
libmwl_la_LDFLAGS =

# Test data to distribute
testdatamwl       = test/data/crab_mwl.fits
testdatamwldir    = $(pyexecdir)/gammalib/tests/mwl/data
testdatamwl_DATA  = $(testdatamwl)

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md pyext \
             include/GMWLLib.hpp \
             test/test_MWL.hpp test/test_MWL.py test/data
