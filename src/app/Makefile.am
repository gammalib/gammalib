## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GApplication.cpp \
          GLog.cpp \
          GPars.cpp \
          GPar.cpp \
          GException_app.cpp

# Build libtool library
noinst_LTLIBRARIES = libapp.la
libapp_la_SOURCES = $(sources)
libapp_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
