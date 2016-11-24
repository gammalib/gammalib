#############################################################################
# Makefile.am for fits module                                               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016 Juergen Knoedlseder                                    #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GFitsHeaderCard.cpp \
          GFitsHeader.cpp \
          GFitsImage.cpp \
          GFitsImageByte.cpp \
          GFitsImageSByte.cpp \
          GFitsImageUShort.cpp \
          GFitsImageShort.cpp \
          GFitsImageULong.cpp \
          GFitsImageLong.cpp \
          GFitsImageLongLong.cpp \
          GFitsImageFloat.cpp \
          GFitsImageDouble.cpp \
          GFitsTable.cpp \
          GFitsAsciiTable.cpp \
          GFitsBinTable.cpp \
          GFitsTableCol.cpp \
          GFitsTableBitCol.cpp \
          GFitsTableBoolCol.cpp \
          GFitsTableByteCol.cpp \
          GFitsTableStringCol.cpp \
          GFitsTableUShortCol.cpp \
          GFitsTableShortCol.cpp \
          GFitsTableULongCol.cpp \
          GFitsTableLongCol.cpp \
          GFitsTableLongLongCol.cpp \
          GFitsTableFloatCol.cpp \
          GFitsTableDoubleCol.cpp \
          GFitsTableCFloatCol.cpp \
          GFitsTableCDoubleCol.cpp \
          GFitsHDU.cpp \
          GFits.cpp \
          GException_fits.cpp

# Build libtool library
noinst_LTLIBRARIES = libfits.la
libfits_la_SOURCES = $(sources)
libfits_la_LDFLAGS =

# Extra files to be distributed in the tar file
EXTRA_DIST = GFitsCfitsio.hpp
