#############################################################################
# Makefile.am for template module                                           #
# ------------------------------------------------------------------------- #
# Copyright (C) 2017 Juergen Knoedlseder                                    #
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
sources = GTPLClass.cpp \
          GXXXObservation.cpp \
          GXXXResponse.cpp \
          GXXXEventList.cpp \
          GXXXEventAtom.cpp \
          GXXXEventCube.cpp \
          GXXXEventBin.cpp \
          GXXXInstDir.cpp \
          GXXXRoi.cpp

# Build libtool library (note that this library is not linked included into
# GammaLib, but we make sure that the library building step works correctly
# for the code)
noinst_LTLIBRARIES = libtemplate.la
libtemplate_la_SOURCES = $(sources)
libtemplate_la_LDFLAGS =

# Additional files to distribute
EXTRA_DIST = GTPLClass.hpp \
             GTPLClass..i \
             GXXXLib.hpp \
             GXXXObservation.hpp \
             GXXXObservation.i \
             GXXXResponse.hpp \
             GXXXResponse.i \
             GXXXEventList.hpp \
             GXXXEventList.i \
             GXXXEventAtom.hpp \
             GXXXEventAtom.i \
             GXXXEventCube.hpp \
             GXXXEventCube.i \
             GXXXEventBin.hpp \
             GXXXEventBin.i \
             GXXXInstDir.hpp \
             GXXXInstDir.i \
             GXXXRoi.hpp \
             GXXXRoi.i \
             inst_Makefile.am \
             inst_README.md \
             inst_test_Makefile.am \
             inst_test_XXX.cpp \
             inst_test_XXX.hpp \
             inst_test_XXX.py \
             xxx.i
