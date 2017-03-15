#############################################################################
# Makefile.am for sky module                                                #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GSkyDir.cpp \
          GHorizDir.cpp \
          GSkyPixel.cpp \
          GSkyMap.cpp \
          GSkyRegions.cpp \
          GSkyRegion.cpp \
          GSkyRegionCircle.cpp \
          GSkyProjection.cpp \
          GHealpix.cpp \
          GWcsRegistry.cpp \
          GWcs.cpp \
          GWcsAIT.cpp \
          GWcsAZP.cpp \
          GWcsCAR.cpp \
          GWcsGLS.cpp \
          GWcsMER.cpp \
          GWcsMOL.cpp \
          GWcsSIN.cpp \
          GWcsSFL.cpp \
          GWcsSTG.cpp \
          GWcsTAN.cpp \
          GException_sky.cpp

# Build libtool library
noinst_LTLIBRARIES = libsky.la
libsky_la_SOURCES = $(sources)
libsky_la_LDFLAGS =
