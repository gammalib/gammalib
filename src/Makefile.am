## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Set sub directories
SUBDIRS = support \
          linalg \
          numerics \
          fits \
          xml \
          sky \
          opt \
          obs \
          model \
          app

# Define sources for this directory
sources =

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional instrument specific libraries
if WITH_INST_MWL
  INST_MWL = ../inst/mwl/libmwl.la
endif
if WITH_INST_CTA
  INST_CTA = ../inst/cta/libcta.la
endif
if WITH_INST_LAT
  INST_LAT = ../inst/lat/liblat.la
endif
if WITH_INST_COM
  INST_COM = ../inst/com/libcom.la
endif

# Scripts to be distributed
scripts = gammalib-setup \
          gammalib-init.sh \
          gammalib-init.csh

# Build libtool library
lib_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = $(sources)
libgamma_la_LIBADD = support/libsupport.la \
                     linalg/liblinalg.la \
                     numerics/libnumerics.la \
                     fits/libfits.la \
                     xml/libxml.la \
                     sky/libsky.la \
                     opt/libopt.la \
                     obs/libobs.la \
                     model/libmodel.la \
                     app/libapp.la \
                     $(INST_MWL) \
                     $(INST_CTA) \
                     $(INST_LAT) \
                     $(INST_COM)
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

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