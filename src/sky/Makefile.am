## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GSkyDir.cpp \
	  GSkymap.cpp \
	  GWcs.cpp \
	  GWcsHPX.cpp

# Build libtool library
noinst_LTLIBRARIES = libsky.la
libsky_la_SOURCES = $(sources)
libsky_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
