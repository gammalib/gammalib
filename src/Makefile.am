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

# Optional instrument specific libraries
if WITH_INST_CTA
  INST_CTA = ../inst/cta/libcta.la
endif
if WITH_INST_LAT
  INST_LAT = ../inst/lat/liblat.la
endif

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
                     $(INST_CTA) \
                     $(INST_LAT)
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
