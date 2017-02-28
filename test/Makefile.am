#############################################################################
# Makefile.am for unit test module                                          #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016 Juergen Knoedlseder                                    #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Optional instrument specific tests
if WITH_INST_MWL
  INST_MWL     = $(top_builddir)/inst/mwl/test/test_MWL
  TEST_MWL     = :$(top_srcdir)/inst/mwl/test
  TEST_MWL_ENV = TEST_MWL_DATA=$(top_srcdir)/inst/mwl/test/data
endif
if WITH_INST_CTA
  INST_CTA     = $(top_builddir)/inst/cta/test/test_CTA
  TEST_CTA     = :$(top_srcdir)/inst/cta/test
  TEST_CTA_ENV = TEST_CTA_DATA=$(top_srcdir)/inst/cta/test/data
endif
if WITH_INST_LAT
  INST_LAT     = $(top_builddir)/inst/lat/test/test_LAT
  TEST_LAT     = :$(top_srcdir)/inst/lat/test
  TEST_LAT_ENV = TEST_LAT_DATA=$(top_srcdir)/inst/lat/test/data
endif
if WITH_INST_COM
  INST_COM     = $(top_builddir)/inst/com/test/test_COM
  TEST_COM     = :$(top_srcdir)/inst/com/test
  TEST_COM_ENV = TEST_COM_DATA=$(top_srcdir)/inst/com/test/data
endif

# Optional python bindings
if WITH_PYTHON
  TEST_PYTHON_SCRIPT = test_python.sh test_examples.sh
  TEST_PYTHON_ENV    = PYTHONPATH=$(top_builddir)/pyext:$(top_builddir)/pyext/build/gammalib$(TEST_MWL)$(TEST_CTA)$(TEST_LAT)$(TEST_COM):$(PYTHONPATH)
endif

# Include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

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
                 test_GObservation

# Set test environment (needed for linking with cfitsio and readline)
TESTS_ENVIRONMENT = @RUNSHARED@=$(top_builddir)/src/.libs$(TEST_ENV_DIR):$(@RUNSHARED@) \
                    $(TEST_PYTHON_ENV) \
                    TEST_SRCDIR=$(abs_top_srcdir) \
                    TEST_BUILDDIR=$(abs_top_builddir) \
                    TEST_DATA=$(abs_top_srcdir)/test/data \
                    $(TEST_MWL_ENV) \
                    $(TEST_CTA_ENV) \
                    $(TEST_LAT_ENV) \
                    $(TEST_COM_ENV)

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

# Extra files to be distributed in the tar file
EXTRA_DIST = data reports testinst \
             README.md \
             test_GApplication.py \
             test_GFits.py \
             test_GLinalg.py \
             test_GModel.py \
             test_GNumerics.py \
             test_GObservation.py \
             test_GOptimizer.py \
             test_GSky.py \
             test_GSupport.py \
             test_GTest.py \
             test_GVO.py \
             test_GXml.py \
             test_GXspec.py \
             test_examples.py \
             test_python.py

# Standard test sources and links
test_GSupport_SOURCES = test_GSupport.cpp test_GSupport.hpp
test_GSupport_LDFLAGS = @LDFLAGS@
test_GSupport_CPPFLAGS = @CPPFLAGS@
test_GSupport_LDADD = $(top_builddir)/src/libgamma.la
#
test_GVector_SOURCES = test_GVector.cpp test_GVector.hpp
test_GVector_LDFLAGS = @LDFLAGS@
test_GVector_CPPFLAGS = @CPPFLAGS@
test_GVector_LDADD = $(top_builddir)/src/libgamma.la
#
test_GMatrix_SOURCES = test_GMatrix.cpp test_GMatrix.hpp
test_GMatrix_LDFLAGS = @LDFLAGS@
test_GMatrix_CPPFLAGS = @CPPFLAGS@
test_GMatrix_LDADD = $(top_builddir)/src/libgamma.la
#
test_GMatrixSparse_SOURCES = test_GMatrixSparse.cpp test_GMatrixSparse.hpp
test_GMatrixSparse_LDFLAGS = @LDFLAGS@
test_GMatrixSparse_CPPFLAGS = @CPPFLAGS@
test_GMatrixSparse_LDADD = $(top_builddir)/src/libgamma.la
#
test_GMatrixSymmetric_SOURCES = test_GMatrixSymmetric.cpp test_GMatrixSymmetric.hpp
test_GMatrixSymmetric_LDFLAGS = @LDFLAGS@
test_GMatrixSymmetric_CPPFLAGS = @CPPFLAGS@
test_GMatrixSymmetric_LDADD = $(top_builddir)/src/libgamma.la
#
test_GNumerics_SOURCES = test_GNumerics.cpp test_GNumerics.hpp
test_GNumerics_LDFLAGS = @LDFLAGS@
test_GNumerics_CPPFLAGS = @CPPFLAGS@
test_GNumerics_LDADD = $(top_builddir)/src/libgamma.la
#
test_GFits_SOURCES = test_GFits.cpp test_GFits.hpp
test_GFits_LDFLAGS = @LDFLAGS@
test_GFits_CPPFLAGS = @CPPFLAGS@
test_GFits_LDADD = $(top_builddir)/src/libgamma.la
#
test_GXml_SOURCES = test_GXml.cpp test_GXml.hpp
test_GXml_LDFLAGS = @LDFLAGS@
test_GXml_CPPFLAGS = @CPPFLAGS@
test_GXml_LDADD = $(top_builddir)/src/libgamma.la
#
test_GVO_SOURCES = test_GVO.cpp test_GVO.hpp
test_GVO_LDFLAGS = @LDFLAGS@
test_GVO_CPPFLAGS = @CPPFLAGS@
test_GVO_LDADD = $(top_builddir)/src/libgamma.la
#
test_GXspec_SOURCES = test_GXspec.cpp test_GXspec.hpp
test_GXspec_LDFLAGS = @LDFLAGS@
test_GXspec_CPPFLAGS = @CPPFLAGS@
test_GXspec_LDADD = $(top_builddir)/src/libgamma.la
#
test_GApplication_SOURCES = test_GApplication.cpp test_GApplication.hpp
test_GApplication_LDFLAGS = @LDFLAGS@
test_GApplication_CPPFLAGS = @CPPFLAGS@
test_GApplication_LDADD = $(top_builddir)/src/libgamma.la
#
test_GModel_SOURCES = test_GModel.cpp test_GModel.hpp
test_GModel_LDFLAGS = @LDFLAGS@
test_GModel_CPPFLAGS = @CPPFLAGS@
test_GModel_LDADD = $(top_builddir)/src/libgamma.la
#
test_GSky_SOURCES = test_GSky.cpp test_GSky.hpp
test_GSky_LDFLAGS = @LDFLAGS@
test_GSky_CPPFLAGS = @CPPFLAGS@
test_GSky_LDADD = $(top_builddir)/src/libgamma.la
#
test_GOptimizer_SOURCES = test_GOptimizer.cpp test_GOptimizer.hpp
test_GOptimizer_LDFLAGS = @LDFLAGS@
test_GOptimizer_CPPFLAGS = @CPPFLAGS@
test_GOptimizer_LDADD = $(top_builddir)/src/libgamma.la
#
test_GObservation_SOURCES = test_GObservation.cpp test_GObservation.hpp
test_GObservation_LDFLAGS = @LDFLAGS@
test_GObservation_CPPFLAGS = @CPPFLAGS@
test_GObservation_LDADD = $(top_builddir)/src/libgamma.la

# Create report directory
all-local:
	mkdir -p reports

# Files to be removed in build directory
clean-local:
	rm -rf $(top_builddir)/test/*.fits
	rm -rf $(top_builddir)/test/*.log
	rm -rf $(top_builddir)/test/*.html
	rm -rf $(top_builddir)/test/*.xml
	rm -rf $(top_builddir)/test/*.reg
	rm -rf $(top_builddir)/test/*.dat
	rm -rf $(top_builddir)/test/reports/*.xml

# Tests to be done on installed version
if WITH_PYTHON
installcheck-local:
	@export GAMMALIB=$(DESTDIR)$(prefix); \
	. $(DESTDIR)$(prefix)/bin/gammalib-init.sh; \
	python -c 'import gammalib; gammalib.test()'
endif

# Add Valgrind rule
valgrind:
	@if type valgrind >/dev/null 2>&1; then \
		cd .. ; \
		ln -s test/data data ; \
		mkdir reports ; \
		for prog in $(check_PROGRAMS) ; do \
			$(SHELL) $(abs_top_builddir)/libtool --mode=execute valgrind --leak-check=full \
			--xml=yes --xml-file=tmp.xml --show-reachable=yes test/$$prog; \
			sed '/<\/valgrindoutput>/d' tmp.xml > test/valgrind_$$prog.xml; \
			echo '</valgrindoutput>' >> test/valgrind_$$prog.xml; \
		done ; \
		mv reports/* test/reports/ ; \
		cd test ; \
	else \
		echo "valgrind  not found in PATH."; \
	fi
