## Process this file with automake to produce Makefile.in

files = GException.hpp GException.cpp \
        GVector.hpp GVector.cpp \
	GMatrix.hpp GMatrix.cpp \
	GSymMatrix.hpp GSymMatrix.cpp \
	GSparseMatrix.hpp GSparseMatrix.cpp \
	GSparseSymbolic.hpp GSparseSymbolic.cpp \
	GSparseNumeric.hpp GSparseNumeric.cpp

# Build dynamic library
#lib_LTLIBRARIES = libgamma.la
#libgamma_la_SOURCES = $(files)
#libgamma_la_LDFLAGS =

# Build static library
noinst_LIBRARIES = libgamma.a
libgamma_a_SOURCES = $(files)
#libgamma_a_LDFLAGS =
