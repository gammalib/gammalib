# ==========================================================================
# This module performs unit tests for the GammaLib linalg module.
#
# Copyright (C) 2012-2014 Juergen Knoedlseder
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
        self.append(self.test_vector,           "Test GVector")
        self.append(self.test_matrix,           "Test GMatrix")
        self.append(self.test_matrix_sparse,    "Test GMatrixSparse")
        self.append(self.test_matrix_symmetric, "Test GMatrixSymmetric")

        # Return
        return

    # Test GVector class
    def test_vector(self):
        """
        Test GVector class.
        """
        # Setup vectors
        vector    = GVector(3)
        vector[0] = 1.0
        vector[1] = 2.0
        vector[2] = 3.0
        vector_b  = vector.copy()

        # Unary addition operator
        vector += vector_b
        self.test_value(vector[0], 2.0);
        self.test_value(vector[1], 4.0);
        self.test_value(vector[2], 6.0);

        # Unary multiplication operator
        vector *= 2.0
        self.test_value(vector[0],  4.0);
        self.test_value(vector[1],  8.0);
        self.test_value(vector[2], 12.0);

        # Unary subtraction operator
        vector -= vector_b
        self.test_value(vector[0], 3.0);
        self.test_value(vector[1], 6.0);
        self.test_value(vector[2], 9.0);

        # Unary division operator
        vector /= 3.0
        self.test_value(vector[0], 1.0);
        self.test_value(vector[1], 2.0);
        self.test_value(vector[2], 3.0);

        # Unary scalar addition operator
        vector += 1.0
        self.test_value(vector[0], 2.0);
        self.test_value(vector[1], 3.0);
        self.test_value(vector[2], 4.0);

        # Unary scalar subtraction operator
        vector -= 1.0
        self.test_value(vector[0], 1.0);
        self.test_value(vector[1], 2.0);
        self.test_value(vector[2], 3.0);

        # Negation operator
        vector = -vector
        self.test_value(vector[0], -1.0);
        self.test_value(vector[1], -2.0);
        self.test_value(vector[2], -3.0);
        vector = -vector

        # Binary addition operator
        vector = vector + vector_b
        self.test_value(vector[0], 2.0);
        self.test_value(vector[1], 4.0);
        self.test_value(vector[2], 6.0);

        # Binary multiplication operator
        vector = vector * 2.0
        self.test_value(vector[0],  4.0);
        self.test_value(vector[1],  8.0);
        self.test_value(vector[2], 12.0);

        # Binary subtraction operator
        vector = vector - vector_b
        self.test_value(vector[0], 3.0);
        self.test_value(vector[1], 6.0);
        self.test_value(vector[2], 9.0);

        # Binary division operator
        vector = vector / 3.0
        self.test_value(vector[0], 1.0);
        self.test_value(vector[1], 2.0);
        self.test_value(vector[2], 3.0);

        # Binary scalar addition operator
        vector = vector + 1.0
        self.test_value(vector[0], 2.0);
        self.test_value(vector[1], 3.0);
        self.test_value(vector[2], 4.0);

        # Binary scalar subtraction operator
        vector = vector - 1.0
        self.test_value(vector[0], 1.0);
        self.test_value(vector[1], 2.0);
        self.test_value(vector[2], 3.0);

        # Scalar product
        product = vector * vector
        self.test_value(product, 14.0);

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

