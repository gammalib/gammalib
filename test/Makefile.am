## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Optional instrument specific tests
if WITH_INST_CTA
  INST_CTA = test_CTA
  INCL_CTA = $(top_srcdir)/inst/cta/include
  test_CTA_SOURCES = $(top_srcdir)/inst/cta/test/test_CTA.cpp
  test_CTA_LDADD = $(top_srcdir)/src/libgamma.la
endif
if WITH_INST_LAT
  INST_LAT = test_LAT
  INCL_LAT = $(top_srcdir)/inst/cta/include
  test_LAT_SOURCES = $(top_srcdir)/inst/cta/test/test_LAT.cpp
  test_LAT_LDADD = $(top_srcdir)/src/libgamma.la
endif

# Include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src/support \
	   -I$(INCL_CTA) \
           -I$(INCL_LAT)

# Test programs
check_PROGRAMS = test_GVector \
                 test_GMatrix \
                 test_GSymMatrix \
                 test_GSparseMatrix \
                 test_GFits \
                 test_GModel \
                 test_GSky \
		 $(INST_CTA) $(INST_LAT)

# Tests
TESTS = test_GVector \
        test_GMatrix \
        test_GSymMatrix \
        test_GSparseMatrix \
        test_GFits \
        test_GModel \
        test_GSky \
	$(INST_CTA) $(INST_LAT)

# Standard test sources and links
test_GVector_SOURCES = test_GVector.cpp
test_GVector_LDADD = $(top_srcdir)/src/libgamma.la
test_GMatrix_SOURCES = test_GMatrix.cpp
test_GMatrix_LDADD = $(top_srcdir)/src/libgamma.la
test_GSymMatrix_SOURCES = test_GSymMatrix.cpp
test_GSymMatrix_LDADD = $(top_srcdir)/src/libgamma.la
test_GSparseMatrix_SOURCES = test_GSparseMatrix.cpp
test_GSparseMatrix_LDADD = $(top_srcdir)/src/libgamma.la
test_GFits_SOURCES = test_GFits.cpp
test_GFits_LDADD = $(top_srcdir)/src/libgamma.la
test_GModel_SOURCES = test_GModel.cpp
test_GModel_LDADD = $(top_srcdir)/src/libgamma.la
test_GSky_SOURCES = test_GSky.cpp
test_GSky_LDADD = $(top_srcdir)/src/libgamma.la
