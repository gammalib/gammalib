## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

files = GException.cpp \
	GTools.cpp \
	GSkyDir.cpp \
	GNodeArray.cpp \
	GImage.cpp \
	GHealpix.cpp \
	obs/GData.cpp \
	obs/GObservation.cpp \
	obs/GEvents.cpp \
	obs/GEventList.cpp \
	obs/GEventCube.cpp \
	obs/GEvent.cpp \
	obs/GEventAtom.cpp \
	obs/GEventBin.cpp \
	obs/GGti.cpp \
	obs/GResponse.cpp \
	lat/GLATObservation.cpp \
	lat/GLATEventList.cpp \
	lat/GLATEventCube.cpp \
	lat/GLATEventAtom.cpp \
	lat/GLATEventBin.cpp \
	lat/GLATResponse.cpp \
	lat/GLATResponseTable.cpp \
	lat/GLATAeff.cpp \
	lat/GLATPsf.cpp \
	lat/GLATEdisp.cpp \
	fits/GFitsHeaderCard.cpp \
	fits/GFitsHeader.cpp \
	fits/GFitsData.cpp \
	fits/GFitsImage.cpp \
	fits/GFitsDblImage.cpp \
	fits/GFitsTable.cpp \
	fits/GFitsAsciiTable.cpp \
	fits/GFitsBinTable.cpp \
	fits/GFitsTableCol.cpp \
	fits/GFitsTableLogCol.cpp \
	fits/GFitsTableStrCol.cpp \
	fits/GFitsTableShtCol.cpp \
	fits/GFitsTableLngCol.cpp \
	fits/GFitsTableFltCol.cpp \
	fits/GFitsTableDblCol.cpp \
	fits/GFitsHDU.cpp \
	fits/GFits.cpp \
	fit/GFitPar.cpp \
	linalg/GVector.cpp \
	linalg/GMatrixBase.cpp \
	linalg/GMatrix.cpp \
	linalg/GSymMatrix.cpp \
	linalg/GSparseMatrix.cpp \
	linalg/GSparseSymbolic.cpp \
	linalg/GSparseNumeric.cpp
	
# Build libtool library
lib_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = $(files)
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
