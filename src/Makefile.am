## Process this file with automake to produce Makefile.in

files = GammaLib.hpp \
        GException.hpp GException.cpp \
        GVector.hpp GVector.cpp \
	GMatrixBase.hpp GMatrixBase.cpp \
	GMatrix.hpp GMatrix.cpp \
	GSymMatrix.hpp GSymMatrix.cpp \
	GSparseMatrix.hpp GSparseMatrix.cpp \
	GSparseSymbolic.hpp GSparseSymbolic.cpp \
	GSparseNumeric.hpp GSparseNumeric.cpp

# Define package included headers (to be installed)
pkginclude_HEADERS = GammaLib.hpp \
                     GException.hpp \
		     GVector.hpp \
		     GMatrix.hpp \
		     GSymMatrix.hpp \
		     GSparseMatrix.hpp \
		     GSparseSymbolic.hpp \
		     GSparseNumeric.hpp

# Build libtool library
lib_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = $(files)
libgamma_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
