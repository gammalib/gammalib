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
	data/GData.cpp \
	rsp/GResponse.cpp \
	rsp/GLATResponse.cpp \
	rsp/GLATResponseTable.cpp \
	rsp/GLATAeff.cpp \
	rsp/GLATPsf.cpp \
	rsp/GLATEdisp.cpp \
	obs/GObservation.cpp \
	fits/GFitsHeaderCard.cpp \
	fits/GFitsHeader.cpp \
	fits/GFitsData.cpp \
	fits/GFitsImage.cpp \
	fits/GFitsDblImage.cpp \
	fits/GFitsTable.cpp \
	fits/GFitsAsciiTable.cpp \
	fits/GFitsBinTable.cpp \
	fits/GFitsTableCol.cpp \
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
