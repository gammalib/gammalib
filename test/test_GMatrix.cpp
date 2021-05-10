/***************************************************************************
 *             test_GMatrix.cpp - Test generic matrix class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2021 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file test_GMatrix.cpp
 * @brief Implementation of unit tests for generic matrices
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "test_GMatrix.hpp"

/* __ Globals ____________________________________________________________ */
double g_matrix[] = {1.0, 2.0, 3.0, 4.0,
                     5.0, 6.0, 7.0, 8.0,
                     9.0, 10., 11., 12.};
double g_vector[] = {1.0, 2.0, 3.0, 4.0};
int    g_rows     = 3;
int    g_cols     = 4;

/* __ Test macros ________________________________________________________ */
#define TEST_FAILURE(WHAT,TEST, MSG, EXCEPTION) \
    test_try(WHAT); \
    try { \
        TEST; \
        test_try_failure(MSG); \
    } \
    catch (EXCEPTION &e) { \
        test_try_success(); \
    } \
    catch (std::exception &e) { \
        test_try_failure(e); \
    }


/***************************************************************************
 * @brief Set test matrix
 ***************************************************************************/
GMatrix TestGMatrix::set_matrix(void) const
{
    // Allocate matrix
    GMatrix matrix(g_rows,g_cols);

    // Set matrix values
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            matrix(row,col) = g_matrix[col+row*g_cols];
        }
    }

    // Return matrix
    return matrix;
}


/***************************************************************************
 * @brief Set test matrix with zero lines
 ***************************************************************************/
GMatrix TestGMatrix::set_matrix_zero(void) const
{
    // Allocate matrix
    GMatrix matrix(g_rows+1,g_cols+1);

    // Set matrix values
	int i = 0;
    for (int row = 0; row < g_rows+1; ++row) {
        if (row == 2) continue;
        int j = 0;
        for (int col = 0; col < g_cols+1; ++col) {
            if (col == 2) continue;
            matrix(row,col) = g_matrix[j+i*g_cols];
            j++;
        }
        i++;
    }

    // Return matrix
    return matrix;
}


/***************************************************************************
 * @brief Set test vector
 ***************************************************************************/
GVector TestGMatrix::set_vector(void) const
{
    // Allocate vector
    GVector vector(g_cols);

    // Set vector values
	for (int col = 0; col < g_cols; ++col) {
        vector[col] = g_vector[col];
    }

    // Return vector
    return vector;
}


/***************************************************************************
 * @brief Check if matrix corresponds to test matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] scale scale factor (default: 1.0).
 * @param[in] offset offset value (default: 0.0).
 *
 * Checks if a matrix corresponds to the test matrix. Optionally, the test
 * matrix can be scaled and an offset can be added.
 ***************************************************************************/
bool TestGMatrix::check_matrix(const GMatrix& matrix,
                               const double&  scale,
                               const double&  offset) const
{
    // Initialise check with true
    bool result = true;
    
    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            double value = g_matrix[col+row*g_cols] * scale + offset;
            if (std::abs(matrix(row,col)-value) > 1.0e-15) {
                result = false;
                break;
            }
        }
	}

    // Return result
    return result;
}


/***************************************************************************
 * @brief Check if transposed matrix corresponds to test matrix
 *
 * @param[in] matrix Trasposed matrix.
 * @param[in] scale scale factor (default: 1.0).
 * @param[in] offset offset value (default: 0.0).
 *
 * Checks if a matrix corresponds to the test matrix. Optionally, the test
 * matrix can be scaled and an offset can be added.
 ***************************************************************************/
bool TestGMatrix::check_matrix_trans(const GMatrix& matrix,
                                     const double&  scale,
                                     const double&  offset) const
{
    // Initialise check with true
    bool result = true;
    
    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            double value = g_matrix[col+row*g_cols] * scale + offset;
            if (std::abs(matrix(col,row)-value) > 1.0e-15) {
                result = false;
                break;
            }
        }
	}

    // Return result
    return result;
}


/***************************************************************************
 * @brief Check lower triangle general matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] ref Reference matrix.
 ***************************************************************************/
bool TestGMatrix::check_matrix_lt(const GMatrix& matrix,
                                  const GMatrix& ref) const
{
    // Initialise check with true
    bool result = true;
    
    // Compare all matrix elements
    for (int row = 0; row < ref.rows(); ++row) {
        for (int col = 0; col < ref.columns(); ++col) {
            double value = (col <= row) ? ref(row,col) : 0.0;
            if (matrix(row,col) != value) {
                result = false;
                break;
            }
        }
	}

    // Return result
    return result;
}


/***************************************************************************
 * @brief Check upper triangle general matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] ref Reference matrix.
 ***************************************************************************/
bool TestGMatrix::check_matrix_ut(const GMatrix& matrix,
                                  const GMatrix& ref) const
{
    // Initialise check with true
    bool result = true;
    
    // Compare all matrix elements
    for (int row = 0; row < ref.rows(); ++row) {
        for (int col = 0; col < ref.columns(); ++col) {
            double value = (col >= row) ? ref(row,col) : 0.0;
            if (matrix(row,col) != value) {
                result = false;
                break;
            }
        }
	}

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGMatrix::set(void)
{
    // Set test name
    name("GMatrix");

    // Append tests
    append(static_cast<pfunction>(&TestGMatrix::empty), "Test empty matrix");
    append(static_cast<pfunction>(&TestGMatrix::alloc_matrix), "Test matrix allocation");
    append(static_cast<pfunction>(&TestGMatrix::assign_values), "Test value assignment");
    append(static_cast<pfunction>(&TestGMatrix::copy_matrix), "Test matrix copying");
    append(static_cast<pfunction>(&TestGMatrix::matrix_operations), "Test matrix operations");
    append(static_cast<pfunction>(&TestGMatrix::matrix_arithmetics), "Test matrix arithmetics");
    append(static_cast<pfunction>(&TestGMatrix::matrix_functions), "Test matrix functions");
    append(static_cast<pfunction>(&TestGMatrix::matrix_compare), "Test matrix comparisons");
    //append(static_cast<pfunction>(&TestGMatrix::matrix_cholesky), "Test matrix Cholesky decomposition");
    append(static_cast<pfunction>(&TestGMatrix::matrix_print), "Test matrix printing");

    // Set members
    m_test   = set_matrix();
    v_test   = set_vector();
    m_bigger = GMatrix(g_rows+1, g_cols+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGMatrix* TestGMatrix::clone(void) const
{
    // Clone test suite
    return new TestGMatrix(*this);
}


/***********************************************************************//**
 * @brief Test empty matrix
 ***************************************************************************/
void TestGMatrix::empty(void)
{
    // Test void constructor
    GMatrix matrix1;
    test_assert(matrix1.is_empty(), "Test that void matrix is empty");
    test_value(matrix1.size(), 0, "Test that void matrix has zero size");
    test_value(matrix1.columns(), 0, "Test that void matrix has zero columns");
    test_value(matrix1.rows(), 0, "Test that void matrix has zero rows");
    test_value(matrix1.classname(), "GMatrix", "Test classname of void matrix");

    // Test allocation constructor with zero rows and columns
    GMatrix matrix2(0,0);
    test_assert(matrix2.is_empty(), "Test that zero-size matrix is empty");
    test_value(matrix2.size(), 0, "Test that zero-size matrix has zero size");
    test_value(matrix2.columns(), 0, "Test that zero-size matrix has zero columns");
    test_value(matrix2.rows(), 0, "Test that zero-size matrix has zero rows");
    test_value(matrix2.classname(), "GMatrix", "Test classname of zero-size matrix");

    // Test allocation constructor with zero rows
    GMatrix matrix3(0,10);
    test_assert(matrix3.is_empty(), "Test that zero-rows matrix is empty");
    test_value(matrix3.size(), 0, "Test that zero-rows matrix has zero size");
    test_value(matrix3.columns(), 0, "Test that zero-rows matrix has zero columns");
    test_value(matrix3.rows(), 0, "Test that zero-rows matrix has zero rows");
    test_value(matrix3.classname(), "GMatrix", "Test classname of zero-rows matrix");

    // Test allocation constructor with zero columns
    GMatrix matrix4(5,0);
    test_assert(matrix4.is_empty(), "Test that zero-columns matrix is empty");
    test_value(matrix4.size(), 0, "Test that zero-columns matrix has zero size");
    test_value(matrix4.columns(), 0, "Test that zero-columns matrix has zero columns");
    test_value(matrix4.rows(), 0, "Test that zero-columns matrix has zero rows");
    test_value(matrix4.classname(), "GMatrix", "Test classname of zero-columns matrix");

    // Test comparison operators
    test_assert(matrix1 == matrix2, "Test that void and zero-size matrix are equal");
    test_assert(!(matrix1 != matrix2), "Test that void and zero-size matrix are not unequal");
    test_assert(matrix1 == matrix3, "Test that void and zero-rows matrix are equal");
    test_assert(!(matrix1 != matrix3), "Test that void and zero-rows matrix are not unequal");
    test_assert(matrix1 == matrix4, "Test that void and zero-columns matrix are equal");
    test_assert(!(matrix1 != matrix4), "Test that void and columns-rows matrix are not unequal");

    // Test matrix attributes and printing
    test_value(matrix1.fill(), 0.0, "Test fill of empty matrix");
    test_value(matrix1.min(), 0.0, "Test minimum of empty matrix");
    test_value(matrix1.max(), 0.0, "Test maximum of empty matrix");
    test_value(matrix1.sum(), 0.0, "Test sum of empty matrix");
    test_value(matrix1.print(), "=== GMatrix ==="
                                "\n Number of rows ............: 0"
                                "\n Number of columns .........: 0"
                                "\n Number of elements ........: 0"
                                "\n Number of allocated cells .: 0",
               "Test printing of empty matrix");

    // Test vector multiplication operator
    GVector vector1;
    GVector vector2 = matrix1 * vector1;
    test_value(vector2.size(), 0, "Test that vector multiplication with empty matrix produces zero size vector");

    // Test matrix assignment operator
    GMatrix matrix5 = matrix1;
    test_assert(matrix5.is_empty(), "Test that assignment of empty matrix produces empty matrix");

    // Test matrix addition operators
    matrix5 = matrix1 + matrix1;
    test_assert(matrix5.is_empty(), "Test that matrix addition of empty matrix produces empty matrix");
    matrix5  = matrix1;
    matrix5 += matrix1;
    test_assert(matrix5.is_empty(), "Test that unary matrix addition of empty matrix produces empty matrix");
    matrix5 = matrix1 + 1.0;
    test_assert(matrix5.is_empty(), "Test that addition of scalar produces empty matrix");
    matrix5  = matrix1;
    matrix5 += 1.0;
    test_assert(matrix5.is_empty(), "Test that unary addition of scalar produces empty matrix");

    // Test matrix subtraction operators
    matrix5 = matrix1 - matrix1;
    test_assert(matrix5.is_empty(), "Test that matrix subtraction of empty matrix produces empty matrix");
    matrix5  = matrix1;
    matrix5 -= matrix1;
    test_assert(matrix5.is_empty(), "Test that unary matrix subtraction of empty matrix produces empty matrix");
    matrix5 = matrix1 - 1.0;
    test_assert(matrix5.is_empty(), "Test that subtraction of scalar produces empty matrix");
    matrix5  = matrix1;
    matrix5 -= 1.0;
    test_assert(matrix5.is_empty(), "Test that unary subtraction of scalar produces empty matrix");

    // Test matrix multiplication operators
    matrix5 = matrix1 * matrix1;
    test_assert(matrix5.is_empty(), "Test that matrix multiplication with empty matrix produces empty matrix");
    matrix5  = matrix1;
    matrix1 *= matrix1;
    test_assert(matrix5.is_empty(), "Test that unary matrix multiplication of empty matrix produces empty matrix");
    matrix5 = matrix1 * 3.0;
    test_assert(matrix5.is_empty(), "Test that multiplication by scalar of empty matrix produces empty matrix");
    matrix5  = matrix1;
    matrix5 *= 3.0;
    test_assert(matrix5.is_empty(), "Test that unary multiplication by scalar of empty matrix produces empty matrix");

    // Test matrix division operator
    matrix5  = matrix1;
    matrix1 /= 3.0;
    test_assert(matrix5.is_empty(), "Test that unary division by scalar of empty matrix produces empty matrix");

    // Test matrix negation operator
    matrix5 = -matrix1;
    test_assert(matrix5.is_empty(), "Test that matrix negation produces empty matrix");

    // Test matrix clearing
    matrix5.clear();
    test_assert(matrix5.is_empty(), "Test that cleared matrix produces empty matrix");

    // Test transpose method
    matrix5 = matrix1.transpose();
    test_assert(matrix5.is_empty(), "Test that transpose() method produces empty matrix");

    // Test invert method (feature not implemented)
    //matrix5 = matrix1.invert();
    //test_assert(matrix5.is_empty(), "Test that invert() method produces empty matrix");

    // Test solve method (feature not implemented)
    //GVector vector3 = matrix1.solve(GVector(0));
    //test_value(vector3.size(), 0, "Test that solve() method produces a zero-size vector");

    // Test abs method
    matrix5 = matrix1.abs();
    test_assert(matrix5.is_empty(), "Test that abs() method produces empty matrix");

    // Test extract_lower_triangle method
    matrix5 = matrix1.extract_lower_triangle();
    test_assert(matrix5.is_empty(), "Test that extract_lower_triangle() method produces empty matrix");

    // Test extract_upper_triangle method
    matrix5 = matrix1.extract_upper_triangle();
    test_assert(matrix5.is_empty(), "Test that extract_upper_triangle() method produces empty matrix");

    // Test matrix access
    TEST_FAILURE("Element access with at() of empty matrix",
                 double& element = matrix1.at(0,0),
                 "Exception expected for element access with at() of empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Element access with at() of empty matrix (const version)",
                 const double& element = matrix1.at(0,0),
                 "Exception expected for element access with at() of empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Row access of empty matrix",
                 GVector vector = matrix1.row(0),
                 "Exception expected for row access of empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Column access of empty matrix",
                 GVector vector = matrix1.column(0),
                 "Exception expected for column access of empty matrix.",
                 GException::out_of_range)

    // Test row and column insertion and addition
    TEST_FAILURE("Row insertion into empty matrix",
                 matrix1.row(0, GVector()),
                 "Exception expected for row insertion into empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Row addition to empty matrix",
                 matrix1.add_to_row(0, GVector()),
                 "Exception expected for row addition to empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Column insertion into empty matrix",
                 matrix1.column(0, GVector()),
                 "Exception expected for column insertion into empty matrix.",
                 GException::out_of_range)
    TEST_FAILURE("Column addition to empty matrix",
                 matrix1.add_to_column(0, GVector()),
                 "Exception expected for column addition to empty matrix.",
                 GException::out_of_range)

    // Test failure models
    TEST_FAILURE("Matrix inversion",
                 GMatrix matrix = matrix1.invert(),
                 "Exception expected for matrix inversion.",
                 GException::feature_not_implemented)
    TEST_FAILURE("Matrix solution with non-zero vector",
                 GVector vector = matrix1.solve(GVector(3)),
                 "Exception expected for matrix solution with non-zero vector.",
                 GException::feature_not_implemented)

    // Test Euler angle method
    matrix5.eulerx(10.0);
    test_assert(!matrix5.is_empty(), "Test that eulerx() method produces non-empty matrix");
    test_value(matrix5.size(), 9, "Test that matrix returned by eulerx() has 9 elements");
    test_value(matrix5.columns(), 3, "Test that matrix returned by eulerx() has 3 columns");
    test_value(matrix5.rows(), 3, "Test that matrix returned by eulerx() has 3 rows");
    matrix5.eulery(10.0);
    test_assert(!matrix5.is_empty(), "Test that eulery() method produces non-empty matrix");
    test_value(matrix5.size(), 9, "Test that matrix returned by eulerx() has 9 elements");
    test_value(matrix5.columns(), 3, "Test that matrix returned by eulerx() has 3 columns");
    test_value(matrix5.rows(), 3, "Test that matrix returned by eulerx() has 3 rows");
    matrix5.eulerz(10.0);
    test_assert(!matrix5.is_empty(), "Test that eulerz() method produces non-empty matrix");
    test_value(matrix5.size(), 9, "Test that matrix returned by eulerx() has 9 elements");
    test_value(matrix5.columns(), 3, "Test that matrix returned by eulerx() has 3 columns");
    test_value(matrix5.rows(), 3, "Test that matrix returned by eulerx() has 3 rows");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix allocation
 ***************************************************************************/
void TestGMatrix::alloc_matrix(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test value assignment
 ***************************************************************************/
void TestGMatrix::assign_values(void)
{
    // Setup 3x3 matrix
    GMatrix test(3,3);
    
    // Assignment individual values
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            double value = i*2.0 + k*2.0;
            test(i,k)   = value;
        }
    }

    // Check assignment of individual values
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            double value = i*2.0 + k*2.0;
            test_value(test(i,k), value, 1.0e-10, "Test matrix element assignment");
        }
    }

    // Check value assignment
    const double ref = 37.89;
    test = ref;
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            test_value(test(i,k), ref, 1.0e-10, "Test matrix element assignment");
        }
    }

    // Verify range checking
    #ifdef G_RANGE_CHECK
    test_try("Verify range checking");
    try {
        test.at(3,3) = 1.0;
        test_try_failure("Expected GException::out_of_range exception.");
    }
    catch (GException::out_of_range &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix copy
 ***************************************************************************/
void TestGMatrix::copy_matrix(void)
{
    // Copy matrix
	GMatrix test = m_test;
    
    // Test if original and compied matrices are correct
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test), "Test matrix copy operator",
                "Found:\n"+test.print()+"\nExpected:\n"+m_test.print());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix operations
 *
 * Tests matrix*vector and matrix*matrix multiplication operations.
 ***************************************************************************/
void TestGMatrix::matrix_operations(void)
{
    // Perform vector multiplication
	GVector test1 = m_test * v_test;

    // Check if the result vector is as expected
    GVector ref1 = test1;
    bool result = true;
    for (int row = 0; row < g_rows; ++row) {
        double value = 0.0;
        for (int col = 0; col < g_cols; ++col) {
            value += g_matrix[col+row*g_cols] * g_vector[col];
        }
        ref1[row] = value;
        if (test1[row] != value) {
            result = false;
            break;
        }
    }

    // Test if original matrix and result vector are correct
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(result, "Test matrix*vector multiplication",
                "Found:\n"+test1.print()+"\nExpected:\n"+ref1.print());

    // Test incompatible vector multiplication
    test_try("Test incompatible matrix*vector multiplication");
    try {
        GVector test2 = m_bigger * v_test;
        test_try_failure("Expected GException::matrix_vector_mismatch exception.");
    }
    catch (GException::matrix_vector_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test matrix multiplication
	GMatrix test3 = m_test * m_test.transpose();
    
    // Check if the result matrix is as expected
    GMatrix ref3 = test3;
    result = true;
    for (int row = 0; row < test3.rows(); ++row) {
        for (int col = 0; col < test3.columns(); ++col) {
            double value = 0.0;
            for (int i = 0; i < g_cols; ++i) {
                value += g_matrix[i+row*g_cols] * g_matrix[i+col*g_cols];
            }
            ref3(row,col) = value;
            if (test3(row,col) != value) {
                result = false;
                break;
            }
        }
    }

    // Test if original matrix and result matrix are correct
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(result, "Test matrix multiplication",
                "Found:\n"+test3.print()+"\nExpected:\n"+ref3.print());
    test_value(test3.rows(), g_rows, "Test number of rows of result matrix");
    test_value(test3.columns(), g_rows, "Test number of columns of result matrix");

    // Test incompatible matrix multiplication
    test_try("Test incompatible matrix multiplication");
    try {
        GMatrix test4 = m_bigger * m_test;
        test_try_failure("Expected GException::matrix_mismatch exception.");
    }
    catch (GException::matrix_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test another incompatible matrix multiplication
    test_try("Test incompatible matrix multiplication");
    try {
        GMatrix test5 = m_bigger * m_test;
        test_try_failure("Expected GException::matrix_mismatch exception.");
    }
    catch (GException::matrix_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix arithmetics
 *
 * Tests matrix arithmetics.
 ***************************************************************************/
void TestGMatrix::matrix_arithmetics(void)
{
	// -GMatrix
	GMatrix test = -m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, -1.0, 0.0), "Test -GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix += GMatrix
	test  = m_test;
	test += m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 2.0, 0.0), "Test GMatrix += GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix -= GMatrix
	test  = m_test;
	test -= m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 0.0, 0.0), "Test GMatrix -= GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix *= 3.0
	test  = m_test;
	test *= 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test GMatrix *= 3.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix /= 3.0
	test  = m_test;
	test /= 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0/3.0, 0.0), "Test GMatrix /= 3.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix + GMatrix
	test = m_test + m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 2.0, 0.0), "Test GMatrix + GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix - GMatrix
	test = m_test - m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 0.0, 0.0), "Test GMatrix - GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix * 3.0
	test = m_test * 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test GMatrix * 3.0",
                "Unexpected result matrix:\n"+test.print());

	// 3.0 * GMatrix
	test = 3.0 * m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test 3.0 * GMatrix",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix / 3.0
	test = m_test / 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0/3.0, 0.0), "Test GMatrix / 3.0",
                "Unexpected result matrix:\n"+test.print());
   
	// GMatrix + 3.0
	test = m_test + 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0, 3.0), "Test GMatrix + 3.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix += 3.0
    test = m_test;
    test += 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0, 3.0), "Test GMatrix + 3.0",
                "Unexpected result matrix:\n"+test.print());
    
	// GMatrix - 5.0
	test = m_test - 5.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0, -5.0), "Test GMatrix -= 5.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrix -= 3.0
    test = m_test;
    test -= 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0, -3.0), "Test GMatrix -= 3.0",
                "Unexpected result matrix:\n"+test.print());

    // Test invalid matrix addition
    test_try("Test invalid matrix addition");
    try {
        test  = m_test;
        test += m_bigger;
        test_try_failure("Expected GException::matrix_mismatch exception.");
    }
    catch (GException::matrix_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix functions
 *
 * Tests matrix functions.
 ***************************************************************************/
void TestGMatrix::matrix_functions(void)
{
    // Minimum
	double min = m_test.min();

    // Check mimimum
    double value = g_matrix[0];
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            if (g_matrix[col+row*g_cols] < value) {
                value = g_matrix[col+row*g_cols];
            }
        }
    }
    test_value(min, value, 0.0, "Test minimum function");

    // Maximum
	double max = m_test.max();

    // Check maximum
    value = g_matrix[0];
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            if (g_matrix[col+row*g_cols] > value) {
                value = g_matrix[col+row*g_cols];
            }
        }
    }
    test_value(max, value, 0.0, "Test maximum function");

	// Sum
	double sum = m_test.sum();

    // Check sum
    value = 0.0;
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
            value += g_matrix[col+row*g_cols];
        }
    }
    test_value(sum, value, 1.0e-20, "Test sum function");

    // Transpose function
	GMatrix test1 = m_test.transpose();
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix_trans(test1, 1.0, 0.0),
                "Test transpose(GMatrix) function",
                "Unexpected transposed matrix:\n"+test1.print());

    // Convert to general matrix
    GMatrix test2 = GMatrix(m_test);
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test2, 1.0, 0.0), 
                "Test GMatrix(GMatrix) constructor",
                "Unexpected GMatrix:\n"+test2.print());

    // Try extracting of lower triangle
    test_try("Try extracting of lower triangle");
    try {
        test2 = m_test.extract_lower_triangle();
        test_try_failure("Expected GException::matrix_not_square exception.");
    }
    catch (GException::matrix_not_square &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Try extracting of upper triangle
    test_try("Try extracting of lower triangle");
    try {
        test2 = m_test.extract_upper_triangle();
        test_try_failure("Expected GException::matrix_not_square exception.");
    }
    catch (GException::matrix_not_square &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Extract lower triangle for square matrix
    GMatrix test3 = m_test * m_test.transpose();
    GMatrix test4 = test3.extract_lower_triangle();
    test_assert(check_matrix_lt(test4, test3), 
                "Test GMatrix.extract_lower_triangle() method",
                "Unexpected lower triangle:\n"+test4.print());

    // Extract upper triangle
    test4 = test3.extract_upper_triangle();
    test_assert(check_matrix_ut(test4, test3), 
                "Test GMatrix.extract_upper_triangle() method",
                "Unexpected lower triangle:\n"+test4.print());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix comparisons
 ***************************************************************************/
void TestGMatrix::matrix_compare(void)
{
    // Allocate an empty matrix
    GMatrix empty(g_rows,g_cols);

    // Test operators
    test_assert((m_test == m_test), "Test == operator");
    test_assert(!(m_test == empty), "Test == operator");
    test_assert(!(m_test == m_bigger), "Test == operator");
    test_assert(!(m_test != m_test), "Test != operator");
    test_assert((m_test != empty), "Test != operator");
    test_assert((m_test != m_bigger), "Test != operator");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test Cholesky decomposition
 ***************************************************************************/
/*
void TestGMatrix::matrix_cholesky(void)
{
    // Test Cholesky decomposition
	GMatrix cd           = cholesky_decompose(m_test);
	GMatrix    cd_lower     = cd.extract_lower_triangle();
	GMatrix    cd_upper     = transpose(cd_lower);
	GMatrix    cd_product   = cd_lower * cd_upper;
	GMatrix    cd_residuals = matrix(m_test) - cd_product;
	double res = (abs(cd_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test cholesky_decompose() method");

    // Test compressed Cholesky decomposition
    GMatrix test_zero         = set_matrix_zero();
	GMatrix cd_zero           = cholesky_decompose(test_zero);
	GMatrix    cd_zero_lower     = cd_zero.extract_lower_triangle();
	GMatrix    cd_zero_upper     = transpose(cd_zero_lower);
	GMatrix    cd_zero_product   = cd_zero_lower * cd_zero_upper;
	GMatrix    cd_zero_residuals = matrix(test_zero) - cd_zero_product;
	res = (abs(cd_zero_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_decompose() method");

	// Test Cholesky inplace decomposition
	GMatrix test = m_test;
    test.cholesky_decompose();
	GMatrix cd_lower2 = test.extract_lower_triangle();
    test_assert((cd_lower2 == cd_lower), "Test inplace cholesky_decompose() method");

    // Test Cholesky solver (first test)
	GVector e0(g_rows);
	GVector a0(g_rows);
	e0[0] = 1.0;
	e0[1] = 0.0;
	e0[2] = 0.0;
	a0[0] = g_matrix[0];
	a0[1] = g_matrix[3];
	a0[2] = g_matrix[6];
	GVector s0 = cd.cholesky_solver(a0) - e0;
	res = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method");

    // Test Cholesky solver (second test)
	e0[0] = 0.0;
	e0[1] = 1.0;
	e0[2] = 0.0;
	a0[0] = g_matrix[1];
	a0[1] = g_matrix[4];
	a0[2] = g_matrix[7];
	s0 = cd.cholesky_solver(a0) - e0;
	res = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method");

    // Test Cholesky solver (third test)
	e0[0] = 0.0;
	e0[1] = 0.0;
	e0[2] = 1.0;
	a0[0] = g_matrix[2];
	a0[1] = g_matrix[5];
	a0[2] = g_matrix[8];
	s0 = cd.cholesky_solver(a0) - e0;
	res = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method");

    // Test compressed Cholesky solver (first test)
	e0 = GVector(g_rows+1);
	a0 = GVector(g_rows+1);
	e0[0] = 1.0;
	e0[1] = 0.0;
	e0[2] = 0.0;
	e0[3] = 0.0;
	a0[0] = g_matrix[0];
	a0[1] = g_matrix[3];
	a0[2] = 0.0;
	a0[3] = g_matrix[6];
	s0    = cd_zero.cholesky_solver(a0) - e0;
	res   = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method");

    // Test compressed Cholesky solver (second test)
	e0[0] = 0.0;
	e0[1] = 1.0;
	e0[2] = 0.0;
	e0[3] = 0.0;
	a0[0] = g_matrix[1];
	a0[1] = g_matrix[4];
	a0[2] = 0.0;
	a0[3] = g_matrix[7];
	s0    = cd_zero.cholesky_solver(a0) - e0;
	res   = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method");

    // Test compressed Cholesky solver (third test)
	e0[0] = 0.0;
	e0[1] = 0.0;
	e0[2] = 0.0;
	e0[3] = 1.0;
	a0[0] = g_matrix[2];
	a0[1] = g_matrix[5];
	a0[2] = 0.0;
	a0[3] = g_matrix[8];
	s0    = cd_zero.cholesky_solver(a0) - e0;
	res   = max(abs(s0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method");

	// Test Cholesky inverter
	GMatrix unit(g_rows,g_cols);
	unit(0,0) = unit(1,1) = unit(2,2) = 1.0;
	GMatrix test_inv = m_test;
	test_inv.cholesky_invert();
    GMatrix ci_product   = m_test * test_inv;
    GMatrix ci_residuals = ci_product - unit;
	res = (abs(ci_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test cholesky_invert method");

	// Test Cholesky inverter for compressed matrix
	unit = GMatrix(4,4);
	unit(0,0) = unit(1,1) = unit(3,3) = 1.0;
	GMatrix test_zero_inv = test_zero;
	test_zero_inv.cholesky_invert();
    GMatrix ciz_product   = test_zero * test_zero_inv;
    GMatrix ciz_residuals = ciz_product - unit;
	res = (abs(ciz_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_invert method");

    // Return
    return;
}
*/

/***************************************************************************
 * @brief Test matrix printing
 ***************************************************************************/
void TestGMatrix::matrix_print(void)
{
    // Set reference
    std::string reference;
    reference.append("=== GMatrix ===\n");
    reference.append(" Number of rows ............: 3\n");
    reference.append(" Number of columns .........: 4\n");
    reference.append(" Number of elements ........: 12\n");
    reference.append(" Number of allocated cells .: 12\n");
    reference.append(" 1, 2, 3, 4\n");
    reference.append(" 5, 6, 7, 8\n");
    reference.append(" 9, 10, 11, 12");

    // Test the print method
    std::string output = m_test.print();
    test_assert((output == reference), "Test matrix printing",
                "Unexpected print() output:\n"+output);

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("GMatrix class testing");

    // Create a test suite
    TestGMatrix test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    bool success = testsuites.run();

    // Save test report
    testsuites.save("reports/GMatrix.xml");

    // Return success status
    return (success ? 0 : 1);
}
