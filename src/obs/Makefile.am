#############################################################################
# Makefile.am for obs module                                                #
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

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GEnergies.cpp \
          GEnergy.cpp \
          GTimes.cpp \
          GTime.cpp \
          GTimeReference.cpp \
          GCaldb.cpp \
          GObservations.cpp \
          GObservations_likelihood.cpp \
          GObservation.cpp \
          GObservationRegistry.cpp \
          GEvents.cpp \
          GEventList.cpp \
          GEventCube.cpp \
          GEvent.cpp \
          GEventAtom.cpp \
          GEventBin.cpp \
          GGti.cpp \
          GRoi.cpp \
          GEbounds.cpp \
          GPhases.cpp \
          GResponse.cpp \
          GResponseCache.cpp \
          GResponseVectorCache.cpp \
          GInstDir.cpp \
          GPhotons.cpp \
          GPhoton.cpp \
          GSource.cpp \
          GPulsar.cpp \
          GPulsarEphemerides.cpp

# Build libtool library
noinst_LTLIBRARIES = libobs.la
libobs_la_SOURCES = $(sources)
libobs_la_LDFLAGS =
