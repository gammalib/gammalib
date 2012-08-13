/***************************************************************************
 *                 test_GMatrix.cpp  -  test GMatrix class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include "GammaLib.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Globals ____________________________________________________________ */
double   g_matrix[] = {1.0, 2.0, 3.0, 4.0,
                       5.0, 6.0, 7.0, 8.0,
                       9.0, 10., 11., 12.};
double   g_vector[] = {1.0, 2.0, 3.0, 4.0};
int g_rows    = 3;
int g_cols    = 4;


/***************************************************************************
 *                              Set test matrix                            *
 ***************************************************************************/
GMatrix set_matrix(GTestSuite& testsuite)
{
    testsuite.test_try("Set test matrix");
    try {
        GMatrix matrix(g_rows,g_cols);
        for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col)
            matrix(row,col) = g_matrix[col+row*g_cols];
        }

        testsuite.test_try_success();
        return matrix;
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return GMatrix(g_rows,g_cols);
}


/***************************************************************************
 *                              Set test vector                            *
 ***************************************************************************/
GVector set_vector(GTestSuite& testsuite)
{
    testsuite.test_try("Set test vector");
    try {
        GVector vector(g_cols);
            for (int col = 0; col < g_cols; ++col)
                vector[col] = g_vector[col];

        testsuite.test_try_success();
        return vector;
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    return GVector(g_cols);
}


/***************************************************************************
 *                             Check test matrix                           *
 ***************************************************************************/
int check_matrix(const GMatrix& m, const double scale, const double add, GTestSuite& testsuite)
{
    int result = 1;
    testsuite.test_try("Check test matrix");
    try {
        for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
                double value = g_matrix[col+row*g_cols] * scale + add;
                if (abs(m(row,col)-value) > 1.0e-15) {
                    result = 0;
                    break;
                    }
            }
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                        Check transposed test matrix                     *
 ***************************************************************************/
int check_transpose_matrix(const GMatrix& m, const double scale, const double add, GTestSuite& testsuite)
{
    int result = 1;
    testsuite.test_try("Check transposed test matrix");
    try {
        for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
                double value = g_matrix[col+row*g_cols] * scale + add;
                if (m(col,row) != value) {
                    result = 0;
                    break;
                    }
            }
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                     Check full test matrix (lower triangle)             *
 ***************************************************************************/
int check_matrix_lt(const GMatrix& m, const double scale, const double add, GTestSuite& testsuite)
{
    int result = 1;
    testsuite.test_try("Check full test matrix (lower triangle)");
    try {
        for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
                double value = (col <= row) ? g_matrix[col+row*g_cols] * scale + add : 0.0;
                if (m(row,col) != value) {
                    result = 0;
                    break;
                    }
            }
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                     Check full test matrix (upper triangle)             *
 ***************************************************************************/
int check_matrix_ut(const GMatrix& m, const double scale, const double add, GTestSuite& testsuite)
{
    int result = 1;
    testsuite.test_try("Check full test matrix (upper triangle)");
    try {
        for (int row = 0; row < g_rows; ++row) {
        for (int col = 0; col < g_cols; ++col) {
                double value = (col >= row) ? g_matrix[col+row*g_cols] * scale + add : 0.0;
                if (m(row,col) != value) {
                    result = 0;
                    break;
                    }
            }
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                             Check matrix*vector                         *
 ***************************************************************************/
int check_matrix_vector(const GVector& v, GTestSuite& testsuite)
{
    int result = 1;
    testsuite.test_try(" Check matrix*vector");
    try {
        for (int row = 0; row < v.size(); ++row) {
        double value = 0.0;
            for (int col = 0; col < g_cols; ++col)
                value += g_matrix[col+row*g_cols] * g_vector[col];
            if (v[row] != value) {
                result = 0;
                break;
            }
        }
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                             Check matrix*matrix                         *
 ***************************************************************************/
int check_matrix_matrix(const GMatrix& m, GTestSuite& testsuite)
{
  if (m.rows() != g_rows || m.cols() != g_rows)
    return 0;
  int result = 1;
  testsuite.test_try("Check matrix*matrix");
  try {
    for (int row = 0; row < m.rows(); ++row) {
	  for (int col = 0; col < m.cols(); ++col) {
        double value = 0.0;
		for (int i = 0; i < g_cols; ++i)
	      value += g_matrix[i+row*g_cols] * g_matrix[i+col*g_cols];
	    if (m(row,col) != value) {
	      result = 0;
	      break;
	    }
      }
	}
        testsuite.test_try_success();
  }
  catch (exception &e) {
      testsuite.test_try_failure(e);
  }
  return result;
}


/***************************************************************************
 *                              Check matrix min                           *
 ***************************************************************************/
int check_matrix_min(const double min, GTestSuite& testsuite)
{
  double value = g_matrix[0];
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col) {
	  if (g_matrix[col+row*g_cols] < value)
		value = g_matrix[col+row*g_cols];
	}
  }
  testsuite.test_assert(min==value, "Check matrix min");
  return (min == value);
}


/***************************************************************************
 *                              Check matrix max                           *
 ***************************************************************************/
int check_matrix_max(const double max, GTestSuite& testsuite)
{
  double value = g_matrix[0];
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col) {
	  if (g_matrix[col+row*g_cols] > value)
		value = g_matrix[col+row*g_cols];
	}
  }
  testsuite.test_assert(max==value,"Check matrix max");
  return (max == value);
}


/***************************************************************************
 *                              Check matrix sum                           *
 ***************************************************************************/
int check_matrix_sum(const double sum, GTestSuite& testsuite)
{
  double value = 0.0;
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col)
      value += g_matrix[col+row*g_cols];
  }
  testsuite.test_assert(sum==value,"Check matrix sum");
  return (sum == value);
}


/***************************************************************************
 *                                Test: Output                             *
 ***************************************************************************/
void test_output(const GMatrix& m_test, GTestSuite& testsuite)
{
  testsuite.test_try("Output test matrix");
  try {
      cout << "Output test matrix"<<endl<< m_test << endl;
    testsuite.test_try_success();
  }
  catch (exception &e) {
      testsuite.test_try_failure(e);
  }
}


/***************************************************************************
 *                     Test: Conversion between matrix types               *
 ***************************************************************************/
void test_conversion(GTestSuite& testsuite)
{
    testsuite.test_try("Matrix conversions");
    try {
        //
        // Setup a symmetric matrix
        int     num = 10;
        GMatrix symmetric(num,num);
        for (int i = 0; i < num; ++i) {
        for (int j = 0; j < num; ++j)
            symmetric(i,j) = (i+j+1)/2.0;
        }

        //
        // Convert symmetric matrix into GSymMatrix object
        //    GSymMatrix converted = sym_matrix(symmetric);
        GSymMatrix converted = GSymMatrix(symmetric);

        //
        // Convert GSymMatrix back to full matrix
        GMatrix back_convert = matrix(converted);

        //
        // Compare back converted matrix to original one. They should be identical
        testsuite.test_assert(symmetric == back_convert,"Compare back converted matrix to original one");
        //
        // Determine the fill of the matrix. It should be 1.0
        double fill = back_convert.fill();
        testsuite.test_assert(abs(fill-1.0) <= 1.0e-15,"Test the fill of the matrix","Bad fill "+str(fill)+" determined (expected 1.0).");
        //
        // Extract lower triangle and check values
        GMatrix lower = symmetric.extract_lower_triangle();
        int ok = 1;
        for (int j = 1; j < num; ++j) {
        for (int i = 0; i < j; ++i) {
            if (lower(i,j) != 0.0)
            ok = 0;
        }
        }
        for (int j = 0; j < num; ++j) {
        for (int i = j; i < num; ++i) {
            if (lower(i,j) != symmetric(i,j))
            ok = 0;
        }
        }
        testsuite.test_assert(ok,"Extract lower triangle");

        //
        // Extract upper triangle and check values
        GMatrix upper = symmetric.extract_upper_triangle();
        ok = 1;
        for (int j = 0; j < num; ++j) {
        for (int i = j+1; i < num; ++i) {
            if (upper(i,j) != 0.0)
            ok = 0;
        }
        }
        for (int j = 0; j < num; ++j) {
        for (int i = 0; i <= j; ++i) {
            if (upper(i,j) != symmetric(i,j))
            ok = 0;
        }
        }
    
        testsuite.test_assert(ok,"Extract upper triangle","Corrupt extract_upper_triangle.");

        //
        // Now make the matrix unsymmetric
        symmetric(0,num-1) = 1000.0;
        //
        // Try converting now into GSymMatrix object (this should fail)
        testsuite.test_try("Try converting now into GSymMatrix object");
        try {
            converted = sym_matrix(symmetric);
            testsuite.test_try_failure();
        }
        catch (GException::matrix_not_symmetric &e) {
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
        //
        // Now zero some elements to emulate a sparse matrix
        symmetric(0,num-1) = 0.0;
        for (int i = 5; i < 8; ++i) {
        for (int j = i-2; j < i+2; ++j)
            symmetric(i,j) = 0.0;
        }

        //
        // Convert symmetric matrix into GSparseMatrix object
        GSparseMatrix sparse = sparse_matrix(symmetric);
        //
        // Convert GSparseMatrix back to full matrix
        back_convert = matrix(sparse);
        //
        // Compare back converted matrix to original one. They should be identical
        testsuite.test_assert(symmetric == back_convert,
                              "Compare back converted matrix to original one",
                              "Unable to convert matrixes (sparse)");
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
 
    // Return
    return;
}


/***************************************************************************
 *                        Test: extraction and insertion                   *
 ***************************************************************************/
void test_extract(GTestSuite& testsuite)
{
    testsuite.test_try("Test GMatrix: Vector extraction, insertion and addition:");
    try {
    //
	// Set-up test matrix
    int rows = 10;
    int cols = 20;
    GMatrix test(rows, cols);
	//
	// Add and extract column vectors
        testsuite.test_try("Add and extract column vectors");
        try
        {
            for (int col = 0; col < cols; ++col) {
                GVector column(rows);
                for (int row = 0; row < rows; ++row)
                    column[row] = (col+1)*100.0 + (row+1)*1.0;
                test.add_col(column, col);
                GVector check = test.extract_col(col);
                if (check != column) {
                    throw testsuite.exception_failure("Unable to add and extract columns.\nAdded column ...: "+column.print()+"\nExtracted column: "+check.print());
                }
            }
                testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }
	
        //
	// Insert and extract column vectors
        testsuite.test_try("Insert and extract column vectors");
        try
        {
            for (int col = 0; col < cols; ++col) {
            GVector column(rows);
            for (int row = 0; row < rows; ++row)
                column[row] = (col+1)*100.0 + (row+1)*1.0;
            test.insert_col(column, col);
            GVector check = test.extract_col(col);
            if (check != column) {
                throw testsuite.exception_failure("Unable to insert and extract columns.\nIserted column ...: "+column.print()+"\nExtracted column: "+check.print());
            }
            }
            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        //
	// Extract rows
        testsuite.test_try("Extract rows");
        try
        {
            for (int row = 0; row < rows; ++row) {
            GVector v_row(cols);
            for (int col = 0; col < cols; ++col)
                v_row[col] = (col+1)*100.0 + (row+1)*1.0;
            GVector check = test.extract_row(row);
            if (check != v_row) {
                throw testsuite.exception_failure("Unable to extract rows.\nIserted row ...: "+v_row.print()+"\nExtracted row: "+check.print());
            }
            }

            testsuite.test_try_success();
        }
        catch(std::exception& e)
        {
            testsuite.test_try_failure(e);
        }

        testsuite.test_try_success();
    }
    catch (exception &e)
    {
            testsuite.test_try_failure(e);
    }

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Create a test suites
    GTestSuites testsuites("GMatrix class testing");
    
    // Create a test suite pointer
    GTestSuite testsuite("GMatrix");
    testsuites.append(testsuite);
    
    // Set test matrix and vector
    GMatrix m_test = set_matrix(testsuite);
    GVector v_test = set_vector(testsuite);
    
    // Set bigger matrix (only used for collision, don't care about content)
    GMatrix bigger(g_rows+1,g_cols+1);
    
    // Prepare result matrix
    GMatrix result = m_test;
    
    // Execute the tests
    test_output(m_test,testsuite);
    test_conversion(testsuite);
    test_extract(testsuite);
    
    // Test 1: Allocate zero matrix
    testsuite.test_try("Allocate zero matrix");
    try {
        GMatrix test1(0,0);
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 2: Allocate too large matrix
    /*
    try {
        cout << "GMatrix - Test 2: Allocate too large matrix: ";
        GMatrix test2(100000,100000);
    }
    catch (bad_alloc &e) {
        cout << "ok." << endl;
    }
    catch (exception &e) {
            cout << endl << "TEST ERROR: Did not signal out of memory condition." << endl;
        cout << e.what() << endl;
            throw;
    }
    */
    
    // Test 3: Assign values

    GMatrix test3(3,3);
    test3(1,1) = 1.0;
    testsuite.test_assert(test3(0,0) == 0.0 && test3(1,0) == 0.0 && test3(2,0) == 0.0 &&
                          test3(0,1) == 0.0 && test3(1,1) == 1.0 && test3(2,1) == 0.0 &&
                          test3(0,2) == 0.0 && test3(1,2) == 0.0 && test3(2,2) == 0.0,"Assign value");

        #if defined(G_RANGE_CHECK)
        testsuite.test_try("Range check");
        try {
            GMatrix test3(3,3);
            test3(3,3) = 1.0;
            testsuite.test_try_failure();
        }
        catch (GException::out_of_range &e) {
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
        #endif

    // Test 4: Matrix copy constructor
    testsuite.test_try("Define matrix using copy constructor");
    try {
            GMatrix test4 = m_test;
            if (!check_matrix(test4, 1.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt copy constructor");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 5: Matrix assignment
    testsuite.test_try("Matrix assignment");
    try {
            //
            // GMatrix = GMatrix
            result = m_test;
            if (!check_matrix(result, 1.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt assignment");
            }
            //
            // GMatrix = GMatrix (bigger matrix)
            result = bigger;
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 6: Transposition
    testsuite.test_try("Transposition");
    try {
            //
            // transpose(GMatrix)
            result = transpose(m_test);
            if (!check_transpose_matrix(result, 1.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt transpose(GMatrix) function.");
            }
            //
            // GMatrix.transpose()
            result = m_test;
            result.transpose();
            if (!check_transpose_matrix(result, 1.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix.transpose() function.");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 7: Matrix*Vector multiplication
    testsuite.test_try("Matrix*Vector multiplication");
    try {
            GVector v_test7 = m_test*v_test;
            if (!check_matrix_vector(v_test7,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
                throw testsuite.exception_failure("Corrupt Matrix*Vector multiplication");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("Bigger*Vector");
    try {
            GVector v_test7 = bigger*v_test;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_vector_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 8: Matrix*Matrix multiplication
    testsuite.test_try("Matrix*Matrix multiplication");
    try {
            GMatrix m_test8 = m_test * transpose(m_test);
            if (!check_matrix_matrix(m_test8,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) 
                throw testsuite.exception_failure("Corrupt Matrix*Matrix multiplication");
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
            throw;
    }
    
    testsuite.test_try("Bigger*Matrix multiplication");
    try {
            GMatrix m_test8 = m_test*bigger;
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("Bigger*Matrix multiplication 2");
    try {
            GMatrix m_test8 = bigger*m_test;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 9: Assignment and arithmetics
    testsuite.test_try("Assignment and arithmetics");
    try {
            //
            // -GMatrix
            result = -m_test;
            if (!check_matrix(result, -1.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
                throw testsuite.exception_failure("Corrupt -GMatrix operator.");
            }
            //
            // GMatrix += GMatrix
            result  = m_test;
            result += m_test;
            if (!check_matrix(result, 2.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
                throw testsuite.exception_failure("Corrupt GMatrix += GMatrix operator.");
            }
            //
            // GMatrix -= GMatrix
            result  = m_test;
            result -= m_test;
            if (!check_matrix(result, 0.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
                throw testsuite.exception_failure("Corrupt GMatrix -= GMatrix operator.");
            }
            //
            // GMatrix *= 3.0
            result  = m_test;
            result *= 3.0;
            if (!check_matrix(result, 3.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix *= double operator.");
        }

            //
            // GMatrix /= 3.0
            result  = m_test;
            result /= 3.0;
            if (!check_matrix(result, 1.0/3.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix /= double operator.");
        }
            //
            // GMatrix + GMatrix
            result = m_test + m_test;
            if (!check_matrix(result, 2.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix + GMatrix operator.");
        }

            //
            // GMatrix - GMatrix
        result = m_test - m_test;
        if (!check_matrix(result, 0.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)){
            throw testsuite.exception_failure("Corrupt GMatrix - GMatrix operator.");
        }
            //
            // GMatrix - GMatrix
            result = m_test - m_test;
            if (!check_matrix(result, 0.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix - GMatrix operator.");
        }
            //
            // GMatrix * 3.0
            result = m_test * 3.0;
            if (!check_matrix(result, 3.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix * double operator.");
        }
            //
            // 3.0 * GMatrix
            result = 3.0 * m_test;
            if (!check_matrix(result, 3.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt double * GMatrix operator.");
        }
            //
            // GMatrix / 3.0
            result = m_test / 3.0;
            if (!check_matrix(result, 1.0/3.0, 0.0,testsuite) || !check_matrix(m_test, 1.0, 0.0,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix / double operator.");
        }
    testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("Test += with bigger");
    try {
            result  = m_test;
            result += bigger;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 10: Matrix functions
    testsuite.test_try("Matrix functions");
    try {
            //
            double min = m_test.min();
        if (!check_matrix_min(min,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix.min() function.");
        }

            //
            double max = m_test.max();
            if (!check_matrix_max(max,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix.max() function.");
        }
            //
            double sum = m_test.sum();
            if (!check_matrix_sum(sum,testsuite)) {
            throw testsuite.exception_failure("Corrupt GMatrix.sum() function.");
        }
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test 11: Matrix comparison
    testsuite.test_try("Matrix comparison");
    try {
            //
            if ((m_test == m_test) != 1) {
                throw testsuite.exception_failure("Corrupt GMatrix == GMatrix operator.");
            }
            //
            GMatrix m_test10(g_rows,g_cols);
            if ((m_test == m_test10) != 0) {
                throw testsuite.exception_failure("Corrupt GMatrix == GMatrix operator.");
            }

            //
            if ((m_test == bigger) != 0) {
                throw testsuite.exception_failure("Corrupt GMatrix == GMatrix operator.");
            }

            //
            if ((m_test != m_test) != 0) {
                throw testsuite.exception_failure("Corrupt GMatrix != GMatrix operator.");
            }

            //
            if ((m_test != m_test10) != 1) {
                throw testsuite.exception_failure("Corrupt GMatrix != GMatrix operator.");
            }

            //
            if ((m_test != bigger) != 1) {
                throw testsuite.exception_failure("Corrupt GMatrix != GMatrix operator.");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.end_test();

    //save xml report
    testsuites.save("reports/GMatrix.xml");

    // Return
    return 0;
}
