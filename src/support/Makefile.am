#############################################################################
# Makefile.am for support module                                            #
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
           -I$(top_srcdir)/src \
           -I$(top_builddir)/src/fits \
           -I$(top_srcdir)/src/fits

# Define sources for this directory
sources = GException.cpp \
          GTools.cpp \
          GNodeArray.cpp \
          GBilinear.cpp \
          GCsv.cpp \
          GRan.cpp \
          GUrl.cpp \
          GUrlFile.cpp \
          GUrlString.cpp \
          GFilename.cpp \
          GDaemon.cpp

# Build libtool library
noinst_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = $(sources)
libsupport_la_LDFLAGS =

