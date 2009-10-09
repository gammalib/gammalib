## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

check_PROGRAMS = test_GVector \
		 test_GMatrix \
		 test_GSymMatrix \
		 test_GSparseMatrix \
		 test_GFits \
		 test_GData \
		 test_GResponse

TESTS = test_GVector \
	test_GMatrix \
	test_GSymMatrix \
	test_GSparseMatrix \
	test_GFits \
	test_GData \
	test_GResponse

test_GVector_SOURCES = test_GVector.cpp
test_GVector_LDADD = ../src/libgamma.la
test_GMatrix_SOURCES = test_GMatrix.cpp
test_GMatrix_LDADD = ../src/libgamma.la
test_GSymMatrix_SOURCES = test_GSymMatrix.cpp
test_GSymMatrix_LDADD = ../src/libgamma.la
test_GSparseMatrix_SOURCES = test_GSparseMatrix.cpp
test_GSparseMatrix_LDADD = ../src/libgamma.la
test_GFits_SOURCES = test_GFits.cpp
test_GFits_LDADD = ../src/libgamma.la
test_GData_SOURCES = test_GData.cpp
test_GData_LDADD = ../src/libgamma.la
test_GResponse_SOURCES = test_GResponse.cpp
test_GResponse_LDADD = ../src/libgamma.la
