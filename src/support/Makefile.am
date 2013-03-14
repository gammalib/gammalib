## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Define sources for this directory
sources = GException.cpp \
          GTools.cpp \
          GNodeArray.cpp \
          GCsv.cpp \
          GRan.cpp \
          GUrl.cpp \
          GUrlFile.cpp

# Build libtool library
noinst_LTLIBRARIES = libsupport.la
libsupport_la_SOURCES = $(sources)
libsupport_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
