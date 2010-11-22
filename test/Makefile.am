## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Optional instrument specific tests
if WITH_INST_CTA
  INST_CTA = test_CTA
  INCL_CTA = -I$(top_srcdir)/inst/cta/include
  test_CTA_SOURCES = $(top_srcdir)/inst/cta/test/test_CTA.cpp
  test_CTA_LDFLAGS = @LDFLAGS@
  test_CTA_CPPFLAGS = @CPPFLAGS@
  test_CTA_LDADD = $(top_srcdir)/src/libgamma.la
endif
if WITH_INST_LAT
  INST_LAT = test_LAT
  INCL_LAT = -I$(top_srcdir)/inst/lat/include
  test_LAT_SOURCES = $(top_srcdir)/inst/lat/test/test_LAT.cpp
  test_LAT_LDFLAGS = @LDFLAGS@
  test_LAT_CPPFLAGS = @CPPFLAGS@
  test_LAT_LDADD = $(top_srcdir)/src/libgamma.la
endif

# Optional python bindings
# TODO: The build directory is hard wired for the moment. This has to be changed!!!
if HAVE_SWIG
if HAVE_PYTHON
if ENABLE_PYTHON_BINDING
  TEST_PYTHON_SCRIPT = test_python.py
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext/build/lib.linux-i686-2.5:$PYTHONPATH
endif # ENABLE_PYTHON_BINDING
endif # HAVE_PYTHON
endif # HAVE_SWIG

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           $(INCL_CTA) \
           $(INCL_LAT)

# Test programs (those will be compiled)
check_PROGRAMS = test_GVector \
                 test_GMatrix \
                 test_GSymMatrix \
                 test_GSparseMatrix \
                 test_GNumerics \
                 test_GFits \
                 test_GXml \
                 test_GApplication \
                 test_GModel \
                 test_GSky \
                 $(INST_CTA) $(INST_LAT)

# Set test environment (needed for linking with cfitsio)
# TODO: The cfitsio directory is harwdired. This has to be changed!!!
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=/usr/local/lib $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = test_GVector \
        test_GMatrix \
        test_GSymMatrix \
        test_GSparseMatrix \
        test_GNumerics \
        test_GFits \
        test_GXml \
        test_GApplication \
        test_GModel \
        test_GSky \
        $(INST_CTA) $(INST_LAT) \
        $(TEST_PYTHON_SCRIPT)

# Standard test sources and links
test_GVector_SOURCES = test_GVector.cpp
test_GVector_LDFLAGS = @LDFLAGS@
test_GVector_CPPFLAGS = @CPPFLAGS@
test_GVector_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GMatrix_SOURCES = test_GMatrix.cpp
test_GMatrix_LDFLAGS = @LDFLAGS@
test_GMatrix_CPPFLAGS = @CPPFLAGS@
test_GMatrix_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GSymMatrix_SOURCES = test_GSymMatrix.cpp
test_GSymMatrix_LDFLAGS = @LDFLAGS@
test_GSymMatrix_CPPFLAGS = @CPPFLAGS@
test_GSymMatrix_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GSparseMatrix_SOURCES = test_GSparseMatrix.cpp
test_GSparseMatrix_LDFLAGS = @LDFLAGS@
test_GSparseMatrix_CPPFLAGS = @CPPFLAGS@
test_GSparseMatrix_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GNumerics_SOURCES = test_GNumerics.cpp
test_GNumerics_LDFLAGS = @LDFLAGS@
test_GNumerics_CPPFLAGS = @CPPFLAGS@
test_GNumerics_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GFits_SOURCES = test_GFits.cpp
test_GFits_LDFLAGS = @LDFLAGS@
test_GFits_CPPFLAGS = @CPPFLAGS@
test_GFits_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GXml_SOURCES = test_GXml.cpp
test_GXml_LDFLAGS = @LDFLAGS@
test_GXml_CPPFLAGS = @CPPFLAGS@
test_GXml_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GApplication_SOURCES = test_GApplication.cpp
test_GApplication_LDFLAGS = @LDFLAGS@
test_GApplication_CPPFLAGS = @CPPFLAGS@
test_GApplication_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GModel_SOURCES = test_GModel.cpp
test_GModel_LDFLAGS = @LDFLAGS@
test_GModel_CPPFLAGS = @CPPFLAGS@
test_GModel_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GSky_SOURCES = test_GSky.cpp
test_GSky_LDFLAGS = @LDFLAGS@
test_GSky_CPPFLAGS = @CPPFLAGS@
test_GSky_LDADD = $(top_srcdir)/src/libgamma.la
