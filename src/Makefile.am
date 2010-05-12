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
          fits \
          sky \
          opt \
          obs \
          model \
          app \
          lat

# Define sources for this directory
sources =

# Build libtool library
lib_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = $(sources)
libgamma_la_LIBADD = support/libsupport.la \
                     linalg/liblinalg.la \
                     fits/libfits.la \
                     sky/libsky.la \
                     opt/libopt.la \
                     obs/libobs.la \
                     model/libmodel.la \
		     app/libapp.la \
                     lat/liblat.la
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
