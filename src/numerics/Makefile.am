#############################################################################
# Makefile.am for numerics module                                           #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Define sources for this directory
sources = GIntegral.cpp \
          GIntegrals.cpp \
          GDerivative.cpp \
          GFunction.cpp \
          GFunctions.cpp \
          GMath.cpp \
          GNdarray.cpp \
          GFft.cpp \
          GFftWavetable.cpp

# Build libtool library
noinst_LTLIBRARIES = libnumerics.la
libnumerics_la_SOURCES = $(sources)
libnumerics_la_LDFLAGS =
