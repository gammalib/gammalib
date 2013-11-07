## Process this file with automake to produce Makefile.in

# we use srcdir/src for "private" includes
# we use srcdir/include for the static global include files
# we use builddir/src for the generated include files to get found

# Optional instrument specific tests
if WITH_INST_MWL
  INST_MWL = test_MWL
  INCL_MWL = -I$(top_srcdir)/inst/mwl/include
  TEST_MWL = :$(top_srcdir)/inst/mwl/test
  test_MWL_SOURCES = $(top_srcdir)/inst/mwl/test/test_MWL.cpp
  test_MWL_LDFLAGS = @LDFLAGS@
  test_MWL_CPPFLAGS = @CPPFLAGS@
  test_MWL_LDADD = $(top_srcdir)/src/libgamma.la
endif
if WITH_INST_CTA
  INST_CTA = test_CTA
  INCL_CTA = -I$(top_srcdir)/inst/cta/include
  TEST_CTA = :$(top_srcdir)/inst/cta/test
  test_CTA_SOURCES = $(top_srcdir)/inst/cta/test/test_CTA.cpp
  test_CTA_LDFLAGS = @LDFLAGS@
  test_CTA_CPPFLAGS = @CPPFLAGS@
  test_CTA_LDADD = $(top_srcdir)/src/libgamma.la
endif
if WITH_INST_LAT
  INST_LAT = test_LAT
  INCL_LAT = -I$(top_srcdir)/inst/lat/include
  TEST_LAT = :$(top_srcdir)/inst/lat/test
  test_LAT_SOURCES = $(top_srcdir)/inst/lat/test/test_LAT.cpp
  test_LAT_LDFLAGS = @LDFLAGS@
  test_LAT_CPPFLAGS = @CPPFLAGS@
  test_LAT_LDADD = $(top_srcdir)/src/libgamma.la
endif
if WITH_INST_COM
  INST_COM = test_COM
  INCL_COM = -I$(top_srcdir)/inst/com/include
  TEST_COM = :$(top_srcdir)/inst/com/test
  test_COM_SOURCES = $(top_srcdir)/inst/com/test/test_COM.cpp
  test_COM_LDFLAGS = @LDFLAGS@
  test_COM_CPPFLAGS = @CPPFLAGS@
  test_COM_LDADD = $(top_srcdir)/src/libgamma.la
endif

# Optional python bindings
if WITH_PYTHON
  TEST_PYTHON_SCRIPT = test_python.py
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext:$(top_builddir)/pyext/build/gammalib$(TEST_MWL)$(TEST_CTA)$(TEST_LAT):$(PYTHONPATH)
endif

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           $(INCL_MWL) \
           $(INCL_CTA) \
           $(INCL_LAT) \
           $(INCL_COM)

# Test programs (those will be compiled)
check_PROGRAMS = test_GSupport \
                 test_GVector \
                 test_GMatrix \
                 test_GMatrixSparse \
                 test_GMatrixSymmetric \
                 test_GNumerics \
                 test_GFits \
                 test_GXml \
                 test_GVO \
                 test_GXspec \
                 test_GApplication \
                 test_GModel \
                 test_GSky \
                 test_GOptimizer \
                 test_GObservation \
                 $(INST_MWL) $(INST_CTA) $(INST_LAT) $(INST_COM)

# Set test environment (needed for linking with cfitsio and readline)
TESTS_ENVIRONMENT = @RUNSHARED@=$(top_builddir)/src/.libs$(TEST_ENV_DIR):$(@RUNSHARED@) \
                    $(TEST_PYTHON_ENV)

# Tests (these are the tests to be executed)
TESTS = test_GSupport \
        test_GVector \
        test_GMatrix \
        test_GMatrixSparse \
        test_GMatrixSymmetric \
        test_GNumerics \
        test_GFits \
        test_GXml \
        test_GVO \
        test_GXspec \
        test_GApplication \
        test_GModel \
        test_GSky \
        test_GOptimizer \
        test_GObservation \
        $(INST_MWL) $(INST_CTA) $(INST_LAT) $(INST_COM) \
        $(TEST_PYTHON_SCRIPT)

# Standard test sources and links
test_GSupport_SOURCES = test_GSupport.cpp
test_GSupport_LDFLAGS = @LDFLAGS@
test_GSupport_CPPFLAGS = @CPPFLAGS@
test_GSupport_LDADD = $(top_srcdir)/src/libgamma.la
#
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
test_GMatrixSparse_SOURCES = test_GMatrixSparse.cpp
test_GMatrixSparse_LDFLAGS = @LDFLAGS@
test_GMatrixSparse_CPPFLAGS = @CPPFLAGS@
test_GMatrixSparse_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GMatrixSymmetric_SOURCES = test_GMatrixSymmetric.cpp
test_GMatrixSymmetric_LDFLAGS = @LDFLAGS@
test_GMatrixSymmetric_CPPFLAGS = @CPPFLAGS@
test_GMatrixSymmetric_LDADD = $(top_srcdir)/src/libgamma.la
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
test_GVO_SOURCES = test_GVO.cpp
test_GVO_LDFLAGS = @LDFLAGS@
test_GVO_CPPFLAGS = @CPPFLAGS@
test_GVO_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GXspec_SOURCES = test_GXspec.cpp
test_GXspec_LDFLAGS = @LDFLAGS@
test_GXspec_CPPFLAGS = @CPPFLAGS@
test_GXspec_LDADD = $(top_srcdir)/src/libgamma.la
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
#
test_GOptimizer_SOURCES = test_GOptimizer.cpp
test_GOptimizer_LDFLAGS = @LDFLAGS@
test_GOptimizer_CPPFLAGS = @CPPFLAGS@
test_GOptimizer_LDADD = $(top_srcdir)/src/libgamma.la
#
test_GObservation_SOURCES = test_GObservation.cpp
test_GObservation_LDFLAGS = @LDFLAGS@
test_GObservation_CPPFLAGS = @CPPFLAGS@
test_GObservation_LDADD = $(top_srcdir)/src/libgamma.la

# Add Valgrind rule
valgrind:
	@if type valgrind >/dev/null 2>&1; then \
		cd .. ; \
		for prog in $(check_PROGRAMS) ; do \
			echo "cd test; ./$$prog; cd .." > test.sh ; chmod u+x test.sh ; \
			$(SHELL) $(abs_top_builddir)/libtool --mode=execute valgrind --leak-check=full \
			--xml=yes --xml-file=tmp.xml --show-reachable=yes ./test.sh; \
			sed '/<\/valgrindoutput>/d' tmp.xml > test/valgrind_$$prog.xml; \
			echo '</valgrindoutput>' >> test/valgrind_$$prog.xml; \
		done ; \
		cd test ; \
	else \
		echo "valgrind  not found in PATH."; \
	fi
