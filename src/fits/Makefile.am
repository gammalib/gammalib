## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GFitsHeaderCard.cpp \
          GFitsHeader.cpp \
          GFitsImage.cpp \
          GFitsByteImage.cpp \
          GFitsImageFlt.cpp \
          GFitsImageDbl.cpp \
          GFitsTable.cpp \
          GFitsAsciiTable.cpp \
          GFitsBinTable.cpp \
          GFitsTableCol.cpp \
	  GFitsTableBitCol.cpp \
          GFitsTableBoolCol.cpp \
          GFitsTableStringCol.cpp \
          GFitsTableUShortCol.cpp \
          GFitsTableShortCol.cpp \
          GFitsTableULongCol.cpp \
          GFitsTableLongCol.cpp \
	  GFitsTableLongLongCol.cpp \
          GFitsTableFloatCol.cpp \
          GFitsTableDoubleCol.cpp \
          GFitsHDU.cpp \
          GFits.cpp \
          GException_fits.cpp

# Build libtool library
noinst_LTLIBRARIES = libfits.la
libfits_la_SOURCES = $(sources)
libfits_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
