#############################################################################
# Makefile.am for example/cpp module                                        #
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

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Include files
INCLUDES = -I$(top_srcdir)/include

# Program name
bin_PROGRAMS = createfits \
               createskymap \
               createxml \
               interpolate \
               numerics \
               optimize \
               readmodel

# Test programs (needed so that files are compiled for 'make check')
check_PROGRAMS = createfits \
                 createskymap \
                 createxml \
                 interpolate \
                 numerics \
                 optimize \
                 readmodel

# Program definition
createfits_SOURCES = createfits.cpp
createfits_LDFLAGS = @LDFLAGS@
createfits_CPPFLAGS = @CPPFLAGS@
createfits_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
createskymap_SOURCES = createskymap.cpp
createskymap_LDFLAGS = @LDFLAGS@
createskymap_CPPFLAGS = @CPPFLAGS@
createskymap_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
createxml_SOURCES = createxml.cpp
createxml_LDFLAGS = @LDFLAGS@
createxml_CPPFLAGS = @CPPFLAGS@
createxml_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
interpolate_SOURCES = interpolate.cpp
interpolate_LDFLAGS = @LDFLAGS@
interpolate_CPPFLAGS = @CPPFLAGS@
interpolate_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
numerics_SOURCES = numerics.cpp
numerics_LDFLAGS = @LDFLAGS@
numerics_CPPFLAGS = @CPPFLAGS@
numerics_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
optimize_SOURCES = optimize.cpp
optimize_LDFLAGS = @LDFLAGS@
optimize_CPPFLAGS = @CPPFLAGS@
optimize_LDADD = $(top_builddir)/src/libgamma.la

# Program definition
readmodel_SOURCES = readmodel.cpp
readmodel_LDFLAGS = @LDFLAGS@
readmodel_CPPFLAGS = @CPPFLAGS@
readmodel_LDADD = $(top_builddir)/src/libgamma.la
