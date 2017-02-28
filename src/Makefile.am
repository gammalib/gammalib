#############################################################################
# Makefile.am for source module                                             #
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

INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Set sub directories
SUBDIRS = base \
          support \
          linalg \
          numerics \
          fits \
          xml \
          vo \
          xspec \
          sky \
          opt \
          obs \
          model \
          app \
          test

# Define sources for this directory
sources =

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional instrument specific libraries
if WITH_INST_MWL
  INST_MWL = $(top_builddir)/inst/mwl/libmwl.la
endif
if WITH_INST_CTA
  INST_CTA = $(top_builddir)/inst/cta/libcta.la
endif
if WITH_INST_LAT
  INST_LAT = $(top_builddir)/inst/lat/liblat.la
endif
if WITH_INST_COM
  INST_COM = $(top_builddir)/inst/com/libcom.la
endif

# Scripts to be distributed
scripts = gammalib-setup \
          gammalib-init.sh \
          gammalib-init.csh

# Build libtool library
lib_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = $(sources)
libgamma_la_LIBADD = $(top_builddir)/src/base/libbase.la \
                     $(top_builddir)/src/support/libsupport.la \
                     $(top_builddir)/src/linalg/liblinalg.la \
                     $(top_builddir)/src/numerics/libnumerics.la \
                     $(top_builddir)/src/fits/libfits.la \
                     $(top_builddir)/src/xml/libxml.la \
                     $(top_builddir)/src/vo/libvo.la \
                     $(top_builddir)/src/xspec/libxspec.la \
                     $(top_builddir)/src/sky/libsky.la \
                     $(top_builddir)/src/opt/libopt.la \
                     $(top_builddir)/src/obs/libobs.la \
                     $(top_builddir)/src/model/libmodel.la \
                     $(top_builddir)/src/app/libapp.la \
                     $(top_builddir)/src/test/libtest.la \
                     $(INST_MWL) \
                     $(INST_CTA) \
                     $(INST_LAT) \
                     $(INST_COM)
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION) @OPENMP_LIB_LDFLAGS@

# Distribute scripts
EXTRA_DIST   = $(scripts)
scriptsdir   = $(bindir)
scripts_DATA = $(scripts)

# Add hook rule to make scripts executable
install-data-hook:
	@list='$(scripts_DATA)';\
	for p in $$list; do \
	  if test -f "$$p"; then chmod +x "$(DESTDIR)$(scriptsdir)/$$p"; fi; \
	done
