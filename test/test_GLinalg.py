# ==========================================================================
# This module performs unit tests for the GammaLib linalg module.
#
# Copyright (C) 2012-2013 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
from gammalib import *
from math import *
import os


# ===================================== #
# Test class for GammaLib linalg module #
# ===================================== #
class Test(GPythonTestSuite):
    """
    Test class for GammaLib linalg module.
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        GPythonTestSuite.__init__(self)

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("linalg")

        # Append tests
        self.append(self.test_matrix,           "Test GMatrix")
        self.append(self.test_matrix_sparse,    "Test GMatrixSparse")
        self.append(self.test_matrix_symmetric, "Test GMatrixSymmetric")

        # Return
        return

    # Test GMatrix class
    def test_matrix(self):
        """
        Test GMatrix class.
        """
        # Set matrix size
        nrows = 3
        ncols = 5

        # Allocate matrix
        m = GMatrix(nrows, ncols)

        # Fill elements
        for i in range(nrows):
            for j in range(ncols):
                m[i, j] = i + j

        # Check elements
        for i in range(nrows):
            for j in range(ncols):
                ref = i + j
                self.test_value(m[i, j], ref, 0.0, "Test matrix element access")

        # Return
        return

    # Test GMatrixSparse class
    def test_matrix_sparse(self):
        """
        Test GMatrixSparse class.
        """
        # Set matrix size
        nrows = 5
        ncols = 3

        # Allocate matrix
        m = GMatrixSparse(nrows, ncols)

        # Fill elements
        for i in range(nrows):
            for j in range(ncols):
                m[i, j] = i + j

        # Check elements
        for i in range(nrows):
            for j in range(ncols):
                ref = i + j
                self.test_value(m[i, j], ref, 0.0, "Test matrix element access")

        # Return
        return

    # Test GMatrixSymmetric class
    def test_matrix_symmetric(self):
        """
        Test GMatrixSymmetric class.
        """
        # Set matrix size
        nrows = 3
        ncols = 3

        # Allocate matrix
        m = GMatrixSymmetric(nrows, ncols)

        # Fill elements
        for i in range(nrows):
            for j in range(ncols):
                m[i, j] = i + j

        # Check elements
        for i in range(nrows):
            for j in range(ncols):
                ref = i + j
                self.test_value(m[i, j], ref, 0.0, "Test matrix element access")

        # Return
        return

