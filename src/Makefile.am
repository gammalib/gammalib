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
	wcs/GWcs.cpp \
	wcs/GWcsHPX.cpp \
	obs/GData.cpp \
	obs/GData_iterator.cpp \
	obs/GData_optimizer.cpp \
	obs/GObservation.cpp \
	obs/GEvents.cpp \
	obs/GEventList.cpp \
	obs/GEventCube.cpp \
	obs/GEvent.cpp \
	obs/GEventAtom.cpp \
	obs/GEventBin.cpp \
	obs/GGti.cpp \
	obs/GEbounds.cpp \
	obs/GResponse.cpp \
	model/GModels.cpp \
	model/GModel.cpp \
	model/GModelPar.cpp \
	model/GModelSpatial.cpp \
	model/GModelSpatialPtsrc.cpp \
	model/GModelSpectral.cpp \
	model/GModelSpectralPlaw.cpp \
	model/GModelTemporal.cpp \
	model/GModelTemporalConst.cpp \
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
	fits/GFitsImageFlt.cpp \
	fits/GFitsImageDbl.cpp \
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
	opt/GOptimizer.cpp \
	opt/GOptimizerLM.cpp \
	opt/GOptimizerPars.cpp \
	opt/GOptimizerFunction.cpp \
	linalg/GVector.cpp \
	linalg/GMatrixTools.cpp \
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
