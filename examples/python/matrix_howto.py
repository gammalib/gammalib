#! /usr/bin/env python
# ==========================================================================
# Scope
#
#   This script provides an example for matrix handling with GammaLib
#
# Usage
#   ./matrix_howto.py
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2013 Juergen Knoedlseder
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
import gammalib


# =========================== #
# Illustrates matrix handling #
# =========================== #
def handle_matrix():
    """
    Illustrates the handling of matrix. Although we use in this example the
    GMatrix class, the same operations can be performed on a sparse matrix
    GMatrixSparse or a symmetric matrix GMatrixSymmetric.
    """
    # Create matrix with 4 rows and 5 columns
    matrix = gammalib.GMatrix(4,5)

    # Get matrix attributes
    number_of_rows     = matrix.rows()     # Will be 4
    number_of_columns  = matrix.columns()  # Will be 5
    number_of_elements = matrix.size()     # Will be 20=4*5
    matrix_fill        = matrix.fill()     # Will be 0
    matrix_minimum     = matrix.min()      # Will be 0
    matrix_maximum     = matrix.max()      # Will be 0
    matrix_sum         = matrix.sum()      # Will be 0

    # Dump some information
    print("Attributes of an empty GMatrix:")
    print(gammalib.parformat("Matrix rows")+str(number_of_rows))
    print(gammalib.parformat("Matrix columns")+str(number_of_columns))
    print(gammalib.parformat("Matrix elements")+str(number_of_elements))
    print(gammalib.parformat("Matrix fill")+str(matrix_fill))
    print(gammalib.parformat("Smallest matrix element")+str(matrix_minimum))
    print(gammalib.parformat("Largest matrix element")+str(matrix_maximum))
    print(gammalib.parformat("Sum of all element")+str(matrix_sum))

    # Set all matrix elements to 2.0
    matrix.set(2.0)

    # Now set all matrix elements to specific values
    for row in range(matrix.rows()):
        for column in range(matrix.columns()):
            matrix[row,column] = row + column*matrix.rows() + 1.5
    
    # Get matrix attributes
    matrix_fill    = matrix.fill()             # Will be 1
    matrix_minimum = matrix.min()              # Will be 1.5
    matrix_maximum = matrix.max()              # Will be 20.5
    matrix_sum     = matrix.sum()              # Will be 220.0

    # Dump some information
    print("Attributes of a filled GMatrix:")
    print(gammalib.parformat("Matrix fill")+str(matrix_fill))
    print(gammalib.parformat("Smallest matrix element")+str(matrix_minimum))
    print(gammalib.parformat("Largest matrix element")+str(matrix_maximum))
    print(gammalib.parformat("Sum of all element")+str(matrix_sum))

    # Perform some matrix operations 
    other_matrix  = matrix.copy()              # Use copy() to create a deep copy
    sum_matrix    = matrix + other_matrix      # Add two matrices
    zero_matrix   = matrix - other_matrix      # Subtract two matrices
    other_matrix += matrix                     # Add-on matrix
    other_matrix -= matrix                     # Subtract off matrix
    other_matrix *= 2.0                        # Multiply elements by 2
    other_matrix /= 2.0                        # Divide elements by 2
    neg_matrix    = -matrix                    # Negate matrix
    abs_matrix    = neg_matrix.abs()           # Absolute value of elements

    # And now some more tricky operations
    transposed    = matrix.transpose()         # First transpose matrix
    mult_matrix   = matrix * transposed        # Now we can multiply it
    other_matrix *= transposed                 # Multiply on matrix
    #inverted      = matrix.invert()           # NOT YET IMPLEMENTED
    vector        = gammalib.GVector(5)        # Allocate vector
    multiplied    = matrix * vector            # Vector multiplication

    # Compare matrices
    if matrix == abs_matrix:
        print("Yes, got it!")
    if matrix != transposed:
        print("Got it again!")

    # Copy columns one-by-one from one matrix into another
    destination = gammalib.GMatrix(4,5)
    for column in range(matrix.columns()):
        vector = matrix.column(column)            # Extract column
        destination.column(column, vector)        # Set column

    # Add columns one-by-one from one matrix into another
    destination = matrix.copy()
    for column in range(matrix.columns()):
        vector = matrix.column(column)            # Extract column
        destination.add_to_column(column, vector) # Add column

    # Copy rows one-by-one from one matrix into another
    destination = gammalib.GMatrix(4,5)
    for row in range(matrix.rows()):
        vector = matrix.row(row)                  # Extract row
        destination.row(row, vector)              # Set row

    # Add rows one-by-one from one matrix into another
    destination = matrix.copy()
    for row in range(matrix.rows()):
        vector = matrix.row(row)                  # Extract row
        destination.add_to_row(row, vector)       # Add row
    
    # Return
    return


# ======================== #
# Illustrates matrix usage #
# ======================== #
def matrix_usage():
    """
    Illustrates the usage of the different matrix storage classes.
    """
    # Create matrices with 4 rows and 4 columns
    matrix_generic   = gammalib.GMatrix(4,4)
    matrix_sparse    = gammalib.GMatrixSparse(4,4)
    matrix_symmetric = gammalib.GMatrixSymmetric(4,4)

    # Now set all matrix elements to specific values
    for row in range(matrix_generic.rows()):
        for column in range(matrix_generic.columns()):
            matrix_generic[row,column]   = row * column
            matrix_sparse[row,column]    = row * column
            matrix_symmetric[row,column] = row * column

    # Convert from one storage class to another
    matrix = gammalib.GMatrix(matrix_sparse)
    matrix = gammalib.GMatrix(matrix_symmetric)
    matrix = gammalib.GMatrixSparse(matrix_generic)
    matrix = gammalib.GMatrixSparse(matrix_symmetric)
    matrix = gammalib.GMatrixSymmetric(matrix_generic)
    matrix = gammalib.GMatrixSymmetric(matrix_sparse)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    Example script for matrix handling.
    """
    # Dump header
    print("")
    print("*******************************")
    print("* Example for matrix handling *")
    print("*******************************")

    # Illustrate matrix handling
    handle_matrix()

    # Illustrate usage of different matrix classes
    matrix_usage()
    