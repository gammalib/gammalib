/***************************************************************************
 *            test_GMatrixSparse.cpp - Test sparse matrix class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
 * @file test_GMatrixSparse.cpp
 * @brief Implementation of unit tests for sparse matrices
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "test_GMatrixSparse.hpp"

/* __ Globals ____________________________________________________________ */
double g_matrix[] = {1.0, 7.0, 3.0, 2.0, 4.0, 8.0, 5.0, 6.0, 9.0};
int    g_row[]    = {  0,   0,   1,   2,   2,   2,   3,   3,   3};
int    g_col[]    = {  0,   4,   1,   0,   2,   4,   2,   3,   4};
double g_vector[] = {1.0, 2.0, 3.0, 4.0, 5.0};
int    g_elements = 9;
int    g_rows     = 4;
int    g_cols     = 5;


/***************************************************************************
 * @brief Set test matrix
 ***************************************************************************/
GMatrixSparse TestGMatrixSparse::set_matrix(void) const
{
    // Allocate matrix
    GMatrixSparse matrix(g_rows, g_cols);

    // Set matrix values
    for (int i = 0; i < g_elements; ++i) {
        matrix(g_row[i],g_col[i]) = g_matrix[i];
    }

    // Return matrix
    return matrix;
}


/***************************************************************************
 * @brief Set test matrix with zero lines
 ***************************************************************************/
/*
GMatrixSparse TestGMatrixSparse::set_matrix_zero(void) const
{
    // Allocate matrix
    GMatrixSparse matrix(g_rows+1,g_cols+1);

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
*/

/***************************************************************************
 * @brief Set test vector
 ***************************************************************************/
GVector TestGMatrixSparse::set_vector(void) const
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
bool TestGMatrixSparse::check_matrix(const GMatrixSparse& matrix,
                                     const double&        scale,
                                     const double&        offset) const
{
    // Initialise check with true
    bool result = true;
    
    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {

            // Set reference value
            double ref_value = offset;
            for (int i = 0; i < g_elements; ++i) {
                if (g_row[i] == row && g_col[i] == col) {
                    ref_value = g_matrix[i] * scale + offset;
                    break;
                }
            }

            // Check matrix
            if (std::abs(matrix(row,col)-ref_value) > 1.0e-15) {
                result = false;
                break;
            }
            
        } // endfor: looped over columns
	} // endfor: looped over rows

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
bool TestGMatrixSparse::check_matrix_trans(const GMatrixSparse& matrix,
                                           const double&        scale,
                                           const double&        offset) const
{
    // Initialise check with true
    bool result = true;

    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {

            // Set reference value
            double ref_value = offset;
            for (int i = 0; i < g_elements; ++i) {
                if (g_row[i] == row && g_col[i] == col) {
                    ref_value = g_matrix[i] * scale + offset;
                    break;
                }
            }

            // Check matrix
            if (std::abs(matrix(col,row)-ref_value) > 1.0e-15) {
                result = false;
                break;
            }
            
        } // endfor: looped over columns
	} // endfor: looped over rows

    // Return result
    return result;
}


/***************************************************************************
 * @brief Check lower triangle general matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] ref Reference matrix.
 ***************************************************************************/
bool TestGMatrixSparse::check_matrix_lt(const GMatrixSparse& matrix,
                                        const GMatrixSparse& ref) const
{
    // Initialise check with true
    bool result = true;

    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {

            // Set reference value
            double ref_value = 0.0;
            if (col <= row) {
                for (int i = 0; i < g_elements; ++i) {
                    if (g_row[i] == row && g_col[i] == col) {
                        ref_value = g_matrix[i];
                        break;
                    }
                }
            }

            // Check matrix
            if (std::abs(matrix(row,col)-ref_value) > 1.0e-15) {
                result = false;
                break;
            }
            
        } // endfor: looped over columns
	} // endfor: looped over rows

    // Return result
    return result;
}


/***************************************************************************
 * @brief Check upper triangle general matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] ref Reference matrix.
 ***************************************************************************/
bool TestGMatrixSparse::check_matrix_ut(const GMatrixSparse& matrix,
                                        const GMatrixSparse& ref) const
{
    // Initialise check with true
    bool result = true;

    // Compare all matrix elements
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {

            // Set reference value
            double ref_value = 0.0;
            if (col >= row) {
                for (int i = 0; i < g_elements; ++i) {
                    if (g_row[i] == row && g_col[i] == col) {
                        ref_value = g_matrix[i];
                        break;
                    }
                }
            }

            // Check matrix
            if (std::abs(matrix(row,col)-ref_value) > 1.0e-15) {
                result = false;
                break;
            }
            
        } // endfor: looped over columns
	} // endfor: looped over rows

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGMatrixSparse::set(void)
{
    // Set test name
    name("GMatrixSparse");

    // Add tests
    append(static_cast<pfunction>(&TestGMatrixSparse::alloc_matrix), "Test matrix allocation");
    append(static_cast<pfunction>(&TestGMatrixSparse::assign_values), "Test value assignment");
    append(static_cast<pfunction>(&TestGMatrixSparse::copy_matrix), "Test matrix copying");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_operations), "Test matrix operations");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_arithmetics), "Test matrix arithmetics");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_functions), "Test matrix functions");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_compare), "Test matrix comparisons");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_cholesky), "Test matrix Cholesky decomposition");
    append(static_cast<pfunction>(&TestGMatrixSparse::matrix_print), "Test matrix printing");

    // Set members
    m_test   = set_matrix();
    v_test   = set_vector();
    m_bigger = GMatrixSparse(g_rows+1, g_cols+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix allocation
 ***************************************************************************/
void TestGMatrixSparse::alloc_matrix(void)
{
    // Allocate zero matrix. The allocation should fail.
    test_try("Allocate zero matrix");
    try {
        GMatrixSparse test(0,0);
        test_try_failure("Expected GException::empty exception.");
    }
    catch (GException::empty &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Setup a symmetric sparse matrix
    int size = 30;
    GMatrixSparse symmetric(size,size);
    for (int i = 0; i < size; i+=2) {
        for (int j = 0; j < size; j+=2) {
            symmetric(i,j) = 1.0+i+j;
        }
    }

    // Convert to GMatrix
    test_try("Test symmetric GMatrix conversion");
    try {
        GMatrix cnv_matrix        = GMatrix(symmetric);
        GMatrixSparse back_matrix = GMatrixSparse(cnv_matrix);
        test_assert((symmetric == back_matrix),
                    "Test symmetric GMatrixSparse - GMatrix conversion",
                    "Found:\n"+back_matrix.print()+"\nExpected:\n"+symmetric.print());
        test_try_success();
    }
    catch (GException::empty &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GMatrix conversion
    test_try("Test GMatrix conversion");
    try {
        GMatrix       cnv_matrix  = GMatrix(m_test);
        GMatrixSparse back_matrix = GMatrixSparse(cnv_matrix);
        test_assert((m_test == back_matrix),
                    "Test GMatrixSparse - GMatrix conversion",
                    "Found:\n"+back_matrix.print()+"\nExpected:\n"+m_test.print());
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GMatrixSparse <-> GMatrixSymmetric conversion
    test_try("Test GMatrixSymmetric conversion");
    try {
        GMatrixSymmetric cnv_sym  = GMatrixSymmetric(symmetric);
        GMatrixSparse    back_sym = GMatrixSparse(cnv_sym);
        test_assert((symmetric == back_sym),
                    "Test GMatrixSparse - GMatrixSymmetric conversion",
                    "Found:\n"+back_sym.print()+"\nExpected:\n"+symmetric.print());
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test invalid GMatrixSparse <-> GMatrixSymmetric conversion
    test_try("Test invalid GMatrixSymmetric conversion");
    try {
        GMatrixSymmetric bad_sym = GMatrixSymmetric(m_test);
        test_try_failure("Expected GException::matrix_not_symmetric exception.");
    }
    catch (GException::matrix_not_symmetric &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test value assignment
 ***************************************************************************/
void TestGMatrixSparse::assign_values(void)
{
    // Setup 3x3 matrix
    GMatrixSparse test(3,3);
    
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

    // Setup 10x10 matrix and keep for reference
	GMatrixSparse sparse(10,10);
	for (int i = 3; i < 5; ++i) {
        sparse(i,i) = 5.0;
    }
	GMatrixSparse initial   = sparse;
	GMatrixSparse reference = sparse;

    // Insert column into 10 x 10 matrix using large matrix stack and the
    // add_col(GVector) method
	sparse.stack_init(100,50);
	GVector column(10);
    for (int j = 0; j < 10; ++j) {
        int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
        if (col > 9) col -= 10;          // This avoids overflow
        int i_min = (j < 2) ?  0 : j-2;
        int i_max = (j > 8) ? 10 : j+2;
        column = 0.0;
        for (int i = i_min; i < i_max; ++i) {
            column[i]         = (i+1)*1;
            reference(i,col) += column[i];
        }
        sparse.add_to_column(col, column);
	}
	sparse.stack_destroy();
    test_assert((sparse == reference),
                "Test stack fill with large stack using add_col(GVector) method",
                "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

    // Insert column into 10 x 10 matrix using small matrix stack and the
    // add_col(GVector) method
    sparse = initial;
	sparse.stack_init(100,3);
    for (int j = 0; j < 10; ++j) {
        int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
        if (col > 9) col -= 10;          // This avoids overflow
        int i_min = (j < 2) ?  0 : j-2;
        int i_max = (j > 8) ? 10 : j+2;
        column = 0.0;
        for (int i = i_min; i < i_max; ++i) {
            column[i]         = (i+1)*1;
        }
        sparse.add_to_column(col, column);
	}
	sparse.stack_destroy();
    test_assert((sparse == reference),
                "Test stack fill with small stack using add_col(GVector) method",
                "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

    // Insert column into 10 x 10 matrix using tiny matrix stack and the
    // add_col(GVector) method
    sparse = initial;
	sparse.stack_init(8,3);
    for (int j = 0; j < 10; ++j) {
        int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
        if (col > 9) col -= 10;          // This avoids overflow
        int i_min = (j < 2) ?  0 : j-2;
        int i_max = (j > 8) ? 10 : j+2;
        column = 0.0;
        for (int i = i_min; i < i_max; ++i) {
            column[i]         = (i+1)*1;
        }
        sparse.add_to_column(col, column);
	}
	sparse.stack_destroy();
    test_assert((sparse == reference),
                "Test stack fill with tiny stack using add_col(GVector) method",
                "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

    // Insert column into 10 x 10 matrix using no matrix stack and the
    // add_col(GVector) method
    sparse = initial;
    for (int j = 0; j < 10; ++j) {
        int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
        if (col > 9) col -= 10;          // This avoids overflow
        int i_min = (j < 2) ?  0 : j-2;
        int i_max = (j > 8) ? 10 : j+2;
        column = 0.0;
        for (int i = i_min; i < i_max; ++i) {
            column[i]         = (i+1)*1;
        }
        sparse.add_to_column(col, column);
	}
	sparse.stack_destroy();
    test_assert((sparse == reference),
                "Test fill using add_col(GVector) method",
                "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

	// Set-up workspace for compressed column adding
	double* wrk_data = new double[10];
	int*    wrk_row  = new int[10];

    // Compressed tests
    test_try("Verify compressed column add_col() method");
    try {

        // Insert column into 10 x 10 matrix using large matrix stack and the
        // compressed column add_col() method
        sparse    = initial;
        reference = initial;
        sparse.stack_init(100,50);
        for (int j = 0; j < 10; ++j) {
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            int inx   = 0;
            int i_min = (j < 3) ?  0 : j-3;
            int i_max = (j > 8) ? 10 : j+2;
            for (int i = i_min; i < i_max; ++i) {
                wrk_data[inx]     = (i+1)*3.7;
                wrk_row[inx]      = i;
                reference(i,col) += wrk_data[inx];
                inx++;
            }
            sparse.add_to_column(col, wrk_data, wrk_row, inx);
        }
        sparse.stack_destroy();
        test_assert((sparse == reference),
                    "Test stack fill with large stack using compressed add_col() method",
                    "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

        // Insert column into 10 x 10 matrix using small matrix stack and the
        // compressed column add_col() method
        sparse = initial;
        sparse.stack_init(100,2);
        for (int j = 0; j < 10; ++j) {
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            int inx   = 0;
            int i_min = (j < 3) ?  0 : j-3;
            int i_max = (j > 8) ? 10 : j+2;
            for (int i = i_min; i < i_max; ++i) {
                wrk_data[inx] = (i+1)*3.7;
                wrk_row[inx]  = i;
                inx++;
            }
            sparse.add_to_column(col, wrk_data, wrk_row, inx);
        }
        sparse.stack_destroy();
        test_assert((sparse == reference),
                    "Test stack fill with small stack using compressed add_col() method",
                    "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

        // Insert column into 10 x 10 matrix using tiny matrix stack and the
        // compressed column add_col() method
        sparse = initial;
        sparse.stack_init(3,2);
        for (int j = 0; j < 10; ++j) {
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            int inx   = 0;
            int i_min = (j < 3) ?  0 : j-3;
            int i_max = (j > 8) ? 10 : j+2;
            for (int i = i_min; i < i_max; ++i) {
                wrk_data[inx] = (i+1)*3.7;
                wrk_row[inx]  = i;
                inx++;
            }
            sparse.add_to_column(col, wrk_data, wrk_row, inx);
        }
        sparse.stack_destroy();
        test_assert((sparse == reference),
                    "Test stack fill with tiny stack using compressed add_col() method",
                    "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

        // Insert column into 10 x 10 matrix using no matrix stack and the
        // compressed column add_col() method
        sparse = initial;
        for (int j = 0; j < 10; ++j) {
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            int inx   = 0;
            int i_min = (j < 3) ?  0 : j-3;
            int i_max = (j > 8) ? 10 : j+2;
            for (int i = i_min; i < i_max; ++i) {
                wrk_data[inx] = (i+1)*3.7;
                wrk_row[inx]  = i;
                inx++;
            }
            sparse.add_to_column(col, wrk_data, wrk_row, inx);
        }
        sparse.stack_destroy();
        test_assert((sparse == reference),
                    "Test fill using compressed add_col() method",
                    "Found:\n"+sparse.print()+"\nExpected:\n"+reference.print());

        // Signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

	// Free workspace
	delete [] wrk_data;
	delete [] wrk_row;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix copy
 ***************************************************************************/
void TestGMatrixSparse::copy_matrix(void)
{
    // Copy matrix
	GMatrixSparse test = m_test;
    
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
void TestGMatrixSparse::matrix_operations(void)
{
    // Perform vector multiplication
	GVector test1 = m_test * v_test;

    // Check result
    GVector ref1(g_rows);
    for (int i = 0; i < g_elements; ++i) {
        ref1[g_row[i]] += g_matrix[i] * v_test[g_col[i]];
    }
    bool result = true;
    for (int i = 0; i < g_rows; ++i) {
        if (ref1[i] != test1[i]) {
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
	GMatrixSparse test3 = m_test * transpose(m_test);

    // Check if the result matrix is as expected
    GMatrixSparse ref3(g_rows, g_rows);
    for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_rows; ++col) {
            double value = 0.0;
            for (int i = 0; i < g_cols; ++i) {
                double ref_value_1 = 0.0;
                double ref_value_2 = 0.0;
                for (int k = 0; k < g_elements; ++k) {
                    if (g_row[k] == row && g_col[k] == i) {
                        ref_value_1 = g_matrix[k];
                        break;
                    }
                }
                for (int k = 0; k < g_elements; ++k) {
                    if (g_row[k] == col && g_col[k] == i) {
                        ref_value_2 = g_matrix[k];
                        break;
                    }
                }
                value += ref_value_1 * ref_value_2;
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
        GMatrixSparse test4 = m_bigger * m_test;
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
        GMatrixSparse test5 = m_bigger * m_test;
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
void TestGMatrixSparse::matrix_arithmetics(void)
{
	// -GMatrixSparse
	GMatrixSparse test = -m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, -1.0, 0.0), "Test -GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse += GMatrixSparse
	test  = m_test;
	test += m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 2.0, 0.0), "Test GMatrixSparse += GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse -= GMatrixSparse
	test  = m_test;
	test -= m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 0.0, 0.0), "Test GMatrixSparse -= GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse *= 3.0
	test  = m_test;
	test *= 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test GMatrixSparse *= 3.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse /= 3.0
	test  = m_test;
	test /= 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0/3.0, 0.0), "Test GMatrixSparse /= 3.0",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse + GMatrixSparse
	test = m_test + m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 2.0, 0.0), "Test GMatrixSparse + GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse - GMatrixSparse
	test = m_test - m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 0.0, 0.0), "Test GMatrixSparse - GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse * 3.0
	test = m_test * 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test GMatrixSparse * 3.0",
                "Unexpected result matrix:\n"+test.print());

	// 3.0 * GMatrixSparse
	test = 3.0 * m_test;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 3.0, 0.0), "Test 3.0 * GMatrixSparse",
                "Unexpected result matrix:\n"+test.print());

	// GMatrixSparse / 3.0
	test = m_test / 3.0;
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test, 1.0/3.0, 0.0), "Test GMatrixSparse / 3.0",
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
void TestGMatrixSparse::matrix_functions(void)
{
    // Minimum
	double min = m_test.min();

    // Check mimimum
    double value = 0.0;
    test_value(min, value, 0.0, "Test minimum function");

    // Maximum
	double max = m_test.max();

    // Check maximum
    value = g_matrix[0];
    for (int i = 1; i < g_elements; ++i) {
        if (g_matrix[i] > value) {
            value = g_matrix[i];
        }
    }
    test_value(max, value, 0.0, "Test maximum function");

	// Sum
	double sum = m_test.sum();

    // Check sum
    value = 0.0;
    for (int i = 0; i < g_elements; ++i) {
        value += g_matrix[i];
    }
    test_value(sum, value, 1.0e-20, "Test sum function");

    // Transpose function
	GMatrixSparse test1 = transpose(m_test);
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix_trans(test1, 1.0, 0.0),
                "Test transpose(GMatrixSparse) function",
                "Unexpected transposed matrix:\n"+test1.print());

    // Transpose method
	test1 = m_test;
	test1.transpose();
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix_trans(test1, 1.0, 0.0), 
                "Test GMatrixSparse.transpose() method",
                "Unexpected transposed matrix:\n"+test1.print());

    // Convert to general matrix
    GMatrixSparse test2 = GMatrixSparse(m_test);
    test_assert(check_matrix(m_test), "Test source matrix");
    test_assert(check_matrix(test2, 1.0, 0.0), 
                "Test GMatrixSparse(GMatrixSparse) constructor",
                "Unexpected GMatrixSparse:\n"+test2.print());
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test matrix comparisons
 ***************************************************************************/
void TestGMatrixSparse::matrix_compare(void)
{
    // Allocate an empty matrix
    GMatrixSparse empty(g_rows, g_cols);

    // Allocate a full matrix
    GMatrixSparse full(g_rows, g_cols);
    for (int i = 0; i < g_rows; ++i) {
        for (int j = 0; j < g_rows; ++j) {
            full(i,j) = 1.0 + i*1000.0 + j*10.0;
        }
    }

    // Test == operator
    test_assert((m_test == m_test), "Test == operator",
                "Identity test faileded on:\n"+m_test.print());
    test_assert(!(m_test == empty), "Test == operator",
                "Matrix\n"+m_test.print()+"\nshould not be equal to\n"+empty.print());
    test_assert(!(m_test == full), "Test == operator",
                "Matrix\n"+m_test.print()+"\nshould not be equal to\n"+full.print());
    test_assert(!(m_test == m_bigger), "Test == operator",
                "Matrix\n"+m_test.print()+"\nshould not be equal to\n"+m_bigger.print());

    // Test != operator
    test_assert(!(m_test != m_test), "Test != operator",
                "Unequality test faileded on:\n"+m_test.print());
    test_assert((m_test != empty), "Test != operator",
                "Matrix\n"+m_test.print()+"\nshould differ from\n"+empty.print());
    test_assert((m_test != full), "Test != operator",
                "Matrix\n"+m_test.print()+"\nshould differ from\n"+full.print());
    test_assert((m_test != m_bigger), "Test != operator",
                "Matrix\n"+m_test.print()+"\nshould differ from\n"+m_bigger.print());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test Cholesky decomposition
 ***************************************************************************/
void TestGMatrixSparse::matrix_cholesky(void)
{
    // Setup matrix for Cholesky decomposition
    GMatrixSparse chol_test(5,5);
    chol_test(0,0) = 1.0;
    chol_test(0,1) = 0.2;
    chol_test(0,2) = 0.2;
    chol_test(0,3) = 0.2;
    chol_test(0,4) = 0.2;
    chol_test(1,0) = 0.2;
    chol_test(2,0) = 0.2;
    chol_test(3,0) = 0.2;
    chol_test(4,0) = 0.2;
    chol_test(1,1) = 1.0;
    chol_test(2,2) = 1.0;
    chol_test(3,3) = 1.0;
    chol_test(4,4) = 1.0;

    // Try to solve now (should not work)
    test_try("Try Cholesky solver without factorisation");
    try {
        GVector vector(5);
        vector = chol_test.cholesky_solver(vector);
        test_try_failure("Expected GException::matrix_not_factorised exception.");
    }
    catch (GException::matrix_not_factorised &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Perform Cholesky decomposition
    GMatrixSparse cd = cholesky_decompose(chol_test);

    // Perform inplace Cholesky decomposition
    cd = chol_test;
    cd.cholesky_decompose();

    // Test Cholesky solver (first test)
    GVector e0(5);
    GVector a0(5);
    e0[0] = 1.0;
    a0[0] = 1.0;
    a0[1] = 0.2;
    a0[2] = 0.2;
    a0[3] = 0.2;
    a0[4] = 0.2;
    GVector s0 = cd.cholesky_solver(a0);
    double res = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method - 1");

    // Test Cholesky solver (second test)
	e0[0] = 0.0;
	e0[1] = 1.0;
	a0[0] = 0.2;
	a0[1] = 1.0;
	a0[2] = 0.0;
	a0[3] = 0.0;
	a0[4] = 0.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method - 2");

    // Test Cholesky solver (third test)
	e0[1] = 0.0;
	e0[2] = 1.0;
	a0[1] = 0.0;
	a0[2] = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method - 3");

    // Test Cholesky solver (forth test)
	e0[2] = 0.0;
	e0[3] = 1.0;
	a0[2] = 0.0;
	a0[3] = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method - 4");

    // Test Cholesky solver (fifth test)
	e0[3] = 0.0;
	e0[4] = 1.0;
	a0[3] = 0.0;
	a0[4] = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test cholesky_solver() method - 5");

	// Setup matrix for Cholesky decomposition with zero row/col
	GMatrixSparse chol_test_zero(6,6);
	chol_test_zero(0,0) = 1.0;
	chol_test_zero(0,1) = 0.2;
	chol_test_zero(0,2) = 0.2;
	chol_test_zero(0,4) = 0.2;
	chol_test_zero(0,5) = 0.2;
	chol_test_zero(1,0) = 0.2;
	chol_test_zero(2,0) = 0.2;
	chol_test_zero(4,0) = 0.2;
	chol_test_zero(5,0) = 0.2;
	chol_test_zero(1,1) = 1.0;
	chol_test_zero(2,2) = 1.0;
	chol_test_zero(4,4) = 1.0;
	chol_test_zero(5,5) = 1.0;

    // Test compressed Cholesky decomposition
	GMatrixSparse cd_zero = chol_test_zero;
	cd_zero.cholesky_decompose();

    // Test compressed Cholesky solver (first test)
	e0 = GVector(6);
	a0 = GVector(6);
	e0[0] = 1.0;
	a0[0] = 1.0;
	a0[1] = 0.2;
	a0[2] = 0.2;
	a0[4] = 0.2;
	a0[5] = 0.2;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method - 1");

    // Test compressed Cholesky solver (second test)
	e0[0] = 0.0;
	e0[1] = 1.0;
	a0[0] = 0.2;
	a0[1] = 1.0;
	a0[2] = 0.0;
	a0[4] = 0.0;
	a0[5] = 0.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method - 2");

    // Test compressed Cholesky solver (third test)
	e0[1] = 0.0;
	e0[2] = 1.0;
	a0[1] = 0.0;
	a0[2] = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method - 3");

    // Test compressed Cholesky solver (forth test)
	e0[2] = 0.0;
	e0[4] = 1.0;
	a0[2] = 0.0;
	a0[4] = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method - 4");

    // Test compressed Cholesky solver (fifth test)
	e0[4] = 0.0;
	e0[5] = 1.0;
	a0[4] = 0.0;
	a0[5] = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, "Test compressed cholesky_solver() method - 5");

	// Setup matrix for Cholesky decomposition with zero row/col (unsymmetric case)
	GMatrixSparse chol_test_zero2(6,5);
	chol_test_zero2(0,0) = 1.0;
	chol_test_zero2(0,1) = 0.2;
	chol_test_zero2(0,2) = 0.2;
	chol_test_zero2(0,3) = 0.2;
	chol_test_zero2(0,4) = 0.2;
	chol_test_zero2(1,0) = 0.2;
	chol_test_zero2(2,0) = 0.2;
	chol_test_zero2(4,0) = 0.2;
	chol_test_zero2(5,0) = 0.2;
	chol_test_zero2(1,1) = 1.0;
	chol_test_zero2(2,2) = 1.0;
	chol_test_zero2(4,3) = 1.0;
	chol_test_zero2(5,4) = 1.0;

    // Test compressed Cholesky decomposition (unsymmetric case)
	GMatrixSparse cd_zero2 = chol_test_zero2;
	cd_zero2.cholesky_decompose();

    // Test compressed Cholesky solver (unsymmetric case)
	e0 = GVector(5);
	a0 = GVector(6);
	e0[0] = 1.0;
	a0[0] = 1.0;
	a0[1] = 0.2;
	a0[2] = 0.2;
	a0[4] = 0.2;
	a0[5] = 0.2;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, 
               "Test unsymmetric compressed cholesky_solver() method - 1");

    // Test compressed Cholesky solver (unsymmetric case)
	e0[0] = 0.0;
	e0[1] = 1.0;
	a0[0] = 0.2;
	a0[1] = 1.0;
	a0[2] = 0.0;
	a0[4] = 0.0;
	a0[5] = 0.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, 
               "Test unsymmetric compressed cholesky_solver() method - 2");

    // Test compressed Cholesky solver (unsymmetric case)
	e0[1] = 0.0;
	e0[2] = 1.0;
	a0[1] = 0.0;
	a0[2] = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, 
               "Test unsymmetric compressed cholesky_solver() method - 3");

    // Test compressed Cholesky solver (unsymmetric case)
	e0[2] = 0.0;
	e0[3] = 1.0;
	a0[2] = 0.0;
	a0[4] = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, 
               "Test unsymmetric compressed cholesky_solver() method - 4");

    // Test compressed Cholesky solver (unsymmetric case)
	e0[3] = 0.0;
	e0[4] = 1.0;
	a0[4] = 0.0;
	a0[5] = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(abs(s0-e0));
    test_value(res, 0.0, 1.0e-15, 
               "Test unsymmetric compressed cholesky_solver() method - 5");

	// Test Cholesky inverter (inplace)
	GMatrixSparse unit(5,5);
	unit(0,0) = 1.0;
	unit(1,1) = 1.0;
	unit(2,2) = 1.0;
	unit(3,3) = 1.0;
	unit(4,4) = 1.0;
	GMatrixSparse chol_test_inv = chol_test;
	chol_test_inv.cholesky_invert();
    GMatrixSparse ci_product   = chol_test * chol_test_inv;
    GMatrixSparse ci_residuals = ci_product - unit;
	res = (abs(ci_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test inplace Cholesky inverter");
            
	// Test Cholesky inverter
	chol_test_inv = cholesky_invert(chol_test);
    ci_product    = chol_test * chol_test_inv;
    ci_residuals  = ci_product - unit;
	res           = (abs(ci_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test Cholesky inverter");

    // Test Cholesky inverter for compressed matrix
    unit = GMatrixSparse(6,6);
    unit(0,0) = 1.0;
    unit(1,1) = 1.0;
    unit(2,2) = 1.0;
    unit(4,4) = 1.0;
    unit(5,5) = 1.0;
    GMatrixSparse chol_test_zero_inv = chol_test_zero;
    chol_test_zero_inv.cholesky_invert();
    GMatrixSparse ciz_product   = chol_test_zero * chol_test_zero_inv;
    GMatrixSparse ciz_residuals = ciz_product - unit;
    res = (abs(ciz_residuals)).max();
    test_value(res, 0.0, 1.0e-15, "Test compressed matrix Cholesky inverter");

    // Return
    return;
}


/***************************************************************************
 * @brief Test matrix printing
 ***************************************************************************/
void TestGMatrixSparse::matrix_print(void)
{
    // Set reference
    std::string reference;
    reference.append("=== GMatrixSparse ===\n");
    reference.append(" Number of rows ............: 4\n");
    reference.append(" Number of columns .........: 5\n");
    reference.append(" Number of nonzero elements : 9\n");
    reference.append(" Pending element ...........: (3,4)=9\n");
    reference.append(" Number of allocated cells .: 512\n");
    reference.append(" Memory block size .........: 512\n");
    reference.append(" Sparse matrix fill ........: 0.45\n");
    reference.append(" Pending element ...........: 9\n");
    reference.append(" Fill stack size ...........: 0 (none)\n");
    reference.append(" 1, 0, 0, 0, 7\n");
    reference.append(" 0, 3, 0, 0, 0\n");
    reference.append(" 2, 0, 4, 0, 8\n");
    reference.append(" 0, 0, 5, 6, 9");
 
    // Test the print method
    std::string output = m_test.print();
    test_assert((output == reference), "Test matrix printing",
                "Unexpected print() output:\n"+output+"\nExpected:\n"+reference);

    // Allocate big matrix
    GMatrixSparse big(30, 40);
    for (int i = 0; i < 30; ++i) {
        for (int j = 0; j < 40; ++j) {
            big(i,j) = i+j;
        }
    }

    // Set reference for big matrix
    std::string ref_big;
    ref_big.append("=== GMatrixSparse ===\n");
    ref_big.append(" Number of rows ............: 30\n");
    ref_big.append(" Number of columns .........: 40\n");
    ref_big.append(" Number of nonzero elements : 1199\n");
    ref_big.append(" Pending element ...........: (29,39)=68\n");
    ref_big.append(" Number of allocated cells .: 1536\n");
    ref_big.append(" Memory block size .........: 512\n");
    ref_big.append(" Sparse matrix fill ........: 0.999167\n");
    ref_big.append(" Pending element ...........: 68\n");
    ref_big.append(" Fill stack size ...........: 0 (none)\n");
    ref_big.append(" 0, 1, 2, 3, 4, ... 35, 36, 37, 38, 39\n");
    ref_big.append(" 1, 2, 3, 4, 5, ... 36, 37, 38, 39, 40\n");
    ref_big.append(" 2, 3, 4, 5, 6, ... 37, 38, 39, 40, 41\n");
    ref_big.append(" 3, 4, 5, 6, 7, ... 38, 39, 40, 41, 42\n");
    ref_big.append(" 4, 5, 6, 7, 8, ... 39, 40, 41, 42, 43\n");
    ref_big.append(" ... ... ... ... ... ... ... ... ... ... ... \n");
    ref_big.append(" 25, 26, 27, 28, 29, ... 60, 61, 62, 63, 64\n");
    ref_big.append(" 26, 27, 28, 29, 30, ... 61, 62, 63, 64, 65\n");
    ref_big.append(" 27, 28, 29, 30, 31, ... 62, 63, 64, 65, 66\n");
    ref_big.append(" 28, 29, 30, 31, 32, ... 63, 64, 65, 66, 67\n");
    ref_big.append(" 29, 30, 31, 32, 33, ... 64, 65, 66, 67, 68");
 
    // Test the print method
    output = big.print();
    test_assert((output == ref_big), "Test big matrix printing",
                "Unexpected print() output:\n"+output+"\nExpected:\n"+ref_big);

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("GMatrixSparse class testing");

    // Initially assume that we pass all tests
    bool success = true;

    // Create a test suite
    TestGMatrixSparse test;

    // Append test suite to the container
    testsuites.append(test);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GMatrixSparse.xml");

    // Return success status
    return (success ? 0 : 1);
}
