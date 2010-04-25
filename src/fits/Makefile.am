## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GFitsHeaderCard.cpp \
	  GFitsHeader.cpp \
	  GFitsData.cpp \
	  GFitsImage.cpp \
	  GFitsImageFlt.cpp \
	  GFitsImageDbl.cpp \
	  GFitsTable.cpp \
	  GFitsAsciiTable.cpp \
	  GFitsBinTable.cpp \
	  GFitsTableCol.cpp \
	  GFitsTableLogCol.cpp \
	  GFitsTableStrCol.cpp \
	  GFitsTableShtCol.cpp \
	  GFitsTableLngCol.cpp \
	  GFitsTableFltCol.cpp \
	  GFitsTableDblCol.cpp \
	  GFitsHDU.cpp \
	  GFits.cpp

# Build libtool library
noinst_LTLIBRARIES = libfits.la
libfits_la_SOURCES = $(sources)
libfits_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
