/***************************************************************************
 *            test_GSparseMatrix.cpp  -  test GSparseMatrix class          *
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
/**
 * @file test_GSparseMatrix.cpp
 * @brief Testing of GSparseMatrix class
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cstdlib>
#include <cmath>
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include "GammaLib.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Definitions ________________________________________________________ */
//#define DUMP_RESIDUALS                                    // Dump residuals
#define DUMP_TIMING                                          // Dump timing

/* __ Globals ____________________________________________________________ */
double g_matrix[] = {1.0, 7.0, 3.0, 2.0, 4.0, 8.0, 5.0, 6.0, 9.0};
int    g_row[]    = {  0,   0,   1,   2,   2,   2,   3,   3,   3};
int    g_col[]    = {  0,   4,   1,   0,   2,   4,   2,   3,   4};
double g_vector[] = {1.0, 2.0, 3.0, 4.0, 5.0};
int    g_elements = 9;
int    g_rows     = 4;
int    g_cols     = 5;


/***************************************************************************
 *                              Set test matrix                            *
 ***************************************************************************/
GSparseMatrix set_matrix(GTestSuite& testsuite)
{
    testsuite.test_try("Set test matrix");
    try {
        GSparseMatrix matrix(g_rows,g_cols);
        for (int i = 0; i < g_elements; ++i)
            matrix(g_row[i],g_col[i]) = g_matrix[i];

        testsuite.test_try_success();
        return matrix;
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    return GSparseMatrix(g_rows,g_cols);
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
int check_matrix(const GSparseMatrix& m, const double scale, const double add)
{
  int result = 1;
  try {
    for (int row = 0; row < g_rows; ++row) {
      for (int col = 0; col < g_cols; ++col) {
	    int i;
		double   ref_value = add;
		for (i = 0; i < g_elements; ++i) {
		  if (g_row[i] == row && g_col[i] == col) {
		    ref_value = g_matrix[i] * scale + add;
			break;
		  }
		}
	    if (abs(m(row,col)-ref_value) > 1.0e-15) {
		  result = 0;
		  break;
		}
	  }
	}
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check test matrix." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                             Check matrix*vector                         *
 ***************************************************************************/
int check_matrix_vector(const GVector& v)
{
  int result = 1;
  try {
    for (int row = 0; row < v.size(); ++row) {
      double value = 0.0;
	  for (int col = 0; col < g_cols; ++col) {
	    int i;
		double   ref_value = 0.0;
		for (i = 0; i < g_elements; ++i) {
		  if (g_row[i] == row && g_col[i] == col) {
		    ref_value = g_matrix[i];
			break;
		  }
		}
	    value += ref_value * g_vector[col];
	  }
	  if (v[row] != value) {
	    result = 0;
	    break;
	  }
    }
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check matrix*vector product." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                             Check matrix*matrix                         *
 ***************************************************************************/
int check_matrix_matrix(const GSparseMatrix& m)
{
  if (m.rows() != g_rows || m.cols() != g_rows)
    return 0;
  int result = 1;
  try {
    for (int row = 0; row < m.rows(); ++row) {
	  for (int col = 0; col < m.cols(); ++col) {
        double value = 0.0;
		for (int i = 0; i < g_cols; ++i) {
	      int k;
		  double   ref_value_1 = 0.0;
		  double   ref_value_2 = 0.0;
		  for (k = 0; k < g_elements; ++k) {
		    if (g_row[k] == row && g_col[k] == i) {
		      ref_value_1 = g_matrix[k];
			  break;
		    }
		  }
		  for (k = 0; k < g_elements; ++k) {
		    if (g_row[k] == col && g_col[k] == i) {
		      ref_value_2 = g_matrix[k];
			  break;
		    }
		  }
		  value += ref_value_1 * ref_value_2;
		}
//cout << row << "," << col << ": " << abs(m(row,col)-value) << endl;
		if (m(row,col) != value) {
		  result = 0;
		  break;
		}
	  }
	}
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check matrix*matrix product." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                        Check transposed test matrix                     *
 ***************************************************************************/
int check_transpose_matrix(const GSparseMatrix& m, const double scale, const double add)
{
  int result = 1;
  try {
    for (int row = 0; row < g_rows; ++row) {
      for (int col = 0; col < g_cols; ++col) {
	    int i;
		double   ref_value = 0.0;
		for (i = 0; i < g_elements; ++i) {
		  if (g_row[i] == row && g_col[i] == col) {
		    ref_value = g_matrix[i] * scale + add;
			break;
		  }
		}
	    if (abs(m(col,row)-ref_value) > 1.0e-15) {
		  result = 0;
		  break;
		}
	  }
	}
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check test matrix." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                              Check matrix min                           *
 ***************************************************************************/
int check_matrix_min(const double min)
{
  double value = g_matrix[0];
  for (int i = 0; i < g_elements; ++i) {
	if (g_matrix[i] < value)
	  value = g_matrix[i];
  }
//  return (min == value);
  return (min == 0.0);
}


/***************************************************************************
 *                              Check matrix max                           *
 ***************************************************************************/
int check_matrix_max(const double max)
{
  double value = g_matrix[0];
  for (int i = 0; i < g_elements; ++i) {
	if (g_matrix[i] > value)
	  value = g_matrix[i];
  }
  return (max == value);
}


/***************************************************************************
 *                              Check matrix sum                           *
 ***************************************************************************/
int check_matrix_sum(const double sum)
{
  double value = 0.0;
  for (int i = 0; i < g_elements; ++i)
	  value += g_matrix[i];
  return (sum == value);
}


/***************************************************************************
 *                                Test: Output                             *
 ***************************************************************************/
void test_output(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    testsuite.test_try("Test output");
    try {
        cout << m_test << endl;
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                          Test: Allocate zero matrix                     *
 ***************************************************************************/
void test_allocate(GTestSuite& testsuite)
{
    testsuite.test_try("Allocate zero matrix");
    try {
        GSparseMatrix test(0,0);
        testsuite.test_try_failure();
    }
    catch (GException::empty &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                              Test: Assign values                        *
 ***************************************************************************/
void test_assign_values(const GSparseMatrix& m_test, GTestSuite& testsuite)
{

    //
    GSparseMatrix test(3,3);
    test(1,1) = 1.0; 
    testsuite.test_assert((test(0,0) == 0.0 && test(1,0) == 0.0 && test(2,0) == 0.0 &&
                           test(0,1) == 0.0 && test(1,1) == 1.0 && test(2,1) == 0.0 &&
                           test(0,2) == 0.0 && test(1,2) == 0.0 && test(2,2) == 0.0),
                           "Test values");
    // 
    GSparseMatrix result = m_test;
    result(0,0) = 0.0;
    result(0,0) = 1.0;
    testsuite.test_assert(check_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                          "Test assignement 1");
    //
    result      = m_test;
    result(0,0) = 0.0;
    result(0,4) = 0.0;
    result(0,0) = 1.0;
    result(0,4) = 0.0;
    result(0,4) = 7.0;
    testsuite.test_assert(check_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                          "Test assignement 2");

    //
    result      = m_test;
    result(0,1) = 2.0;
    result(0,2) = 2.0;
    result(0,3) = 2.0;
    result(0,1) = 0.0;
    result(0,2) = 0.0;
    result(0,3) = 0.0;
    testsuite.test_assert(check_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                          "Test assignement 3");

    //
    result.clear();
    testsuite.test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                          "Test clear");

    #if defined(G_RANGE_CHECK)
    testsuite.test_try("Test range check");
    try {
        GSparseMatrix test(3,3);
        test(3,3) = 1.0;
        testsuite.test_try_failure("Should failed");
    }
    catch (GException::out_of_range &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    #endif

    //
    // Insert column into 10 x 10 matrix using matrix stack and the
    // compressed vector format
    testsuite.test_try("Insert column into 10 x 10 matrix stack");
    try {
        //
            // Set-up sparse matrix for adding  
            GSparseMatrix sparse(10,10);
            for (int i = 3; i < 5; ++i)
            sparse(i,i) = 5.0;
            GSparseMatrix reference = sparse;
            //
            // Initialise vector for column adding
            GVector column(10);
            //
            // Initialise matrix stack  
            sparse.stack_init(100,50);
            //
            // Fill matrix using add_col
        for (int j = 0; j < 10; ++j) {
            // 
            // Set column index
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            //
            // Set-up vector
            int i_min = (j < 2) ?  0 : j-2;
            int i_max = (j > 8) ? 10 : j+2;
            column    = 0.0;
            for (int i = i_min; i < i_max; ++i) {
                    column[i]         = (i+1)*1;
                    reference(i,col) += column[i];
            }
            //
            // Add vector
            sparse.add_col(column, col);
        }

        //
        // Destroy stack
        sparse.stack_destroy();

        if (sparse != reference) {
            testsuite.exception_failure("Corrupt add_col function (vector version).");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    //
    // Insert column into 10 x 10 matrix using matrix stack and the
    // compressed vector format
    testsuite.test_try("Insert column into 10 x 10 matrix stack 2");
    try {
        //
            // Set-up sparse matrix for adding  
            GSparseMatrix sparse(10,10);
            sparse.clear();
            for (int i = 3; i < 5; ++i)
            sparse(i,i) = 5.0;
            GSparseMatrix reference = sparse;
            //
            // Set-up workspace for column adding
            double* wrk_data = new double[10];
            int*    wrk_row  = new int[10];
            //
            // Initialise matrix stack  
            sparse.stack_init(100,50);
            //
            // Fill matrix using add_col
        for (int j = 0; j < 10; ++j) {
            // 
            // Set column index
            int col = int(0.8 * j + 0.5);    // This allows that some columns are twice
            if (col > 9) col -= 10;          // This avoids overflow
            //
            // Set-up vector
            int inx   = 0;
            int i_min = (j < 2) ?  0 : j-2;
            int i_max = (j > 8) ? 10 : j+2;
            for (int i = i_min; i < i_max; ++i) {
                    wrk_data[inx]     = (i+1)*1;
                    wrk_row[inx]      = i;
                    reference(i,col) += wrk_data[inx];
                    inx++;
            }
            //
            // Add vector
            sparse.add_col(wrk_data, wrk_row, inx, col);
        }

        //
        // Destroy stack
        sparse.stack_destroy();
        if (sparse != reference) {
            testsuite.exception_failure("Corrupt add_col function (compressed array version).");
        }
        //
        // Free workspace
        delete [] wrk_data;
        delete [] wrk_row;

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

}


/***************************************************************************
 *                        Test: Matrix copy constructor                    *
 ***************************************************************************/
void test_copy_constructor(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    testsuite.test_try("Test copy constructor");
    try {
            //
            GSparseMatrix test = m_test;
        if (!check_matrix(test, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt copy constructor.");
            }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                           Test: Matrix assignment                       *
 ***************************************************************************/
void test_assign_matrix(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    testsuite.test_try("GSparseMatrix = GSparseMatrix");
    GSparseMatrix result(1,1);
    try {
        //
        // GSparseMatrix = GSparseMatrix
        result = m_test;
        if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt assignment.");
        }
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("GSparseMatrix = GSparseMatrix (bigger matrix)");
    try {
        //
        // GSparseMatrix = GSparseMatrix (bigger matrix)
        GSparseMatrix bigger(1000,1000);
        for (int i = 0; i < 1000; ++i)
        bigger(i,i) = (i+1)*1.1;
        result = bigger;

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                           Test: Matrix transpose                        *
 ***************************************************************************/
void test_transpose(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    GSparseMatrix result(1,1);

    testsuite.test_try("transpose(GSparseMatrix)");
    try {
        //
        // transpose(GSparseMatrix)
        result = transpose(m_test);
        if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt transpose(GSparseMatrix) function.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("GSparseMatrix.transpose()");
    try {
        //
        // GSparseMatrix.transpose()
        result = m_test;
        result.transpose();
         if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
             throw testsuite.exception_failure("Corrupt transpose() function.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                     Test: Matrix*Vector multiplication                  *
 ***************************************************************************/
void test_matrix_vector(const GSparseMatrix& m_test, const GVector& v_test, GTestSuite& testsuite)
{
    testsuite.test_try("matrix*vector");
    try {
        //
        GVector test = m_test*v_test;
        if (!check_matrix_vector(test) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt Matrix*Vector multiplication.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("bigger*vector");
    try {
        GSparseMatrix bigger(1000,1000);
        for (int i = 0; i < 1000; ++i)
        bigger(i,i) = (i+1)*1.1;
        GVector test = bigger*v_test;

        testsuite.test_try_failure();
    }
    catch (GException::matrix_vector_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

}


/***************************************************************************
 *                     Test: Matrix*Vector multiplication                  *
 ***************************************************************************/
void test_matrix_matrix(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    testsuite.test_try("matrix*vector");
    try {

        GSparseMatrix test = m_test * transpose(m_test);
        if (!check_matrix_matrix(test) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt Matrix*Matrix multiplication.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    GSparseMatrix bigger(1000,1000);
    for (int i = 0; i < 1000; ++i)
    bigger(i,i) = (i+1)*1.1;
    testsuite.test_try("matrix*bigger");
    try {
        GSparseMatrix test = m_test * bigger;
        testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_try("bigger*matrix");
    try {
        GSparseMatrix test7 = bigger * m_test;
        testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

}


/***************************************************************************
 *                           Test: Arithmetics                             *
 ***************************************************************************/
void test_arithmetics(const GSparseMatrix& m_test, GTestSuite& testsuite)
{

    //
    // -GSparseMatrix
    GSparseMatrix result(1,1);
    result = -m_test;
    testsuite.test_assert(check_matrix(result, -1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "-GSparseMatrix",
                            "Corrupt -GSparseMatrix operator.");

    //
    // abs(GSparseMatrix)
    result = abs(m_test);
    testsuite.test_assert(check_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "abs(GSparseMatrix)",
                            "Corrupt abs(GSparseMatrix) operator.");
    //
    // GSparseMatrix += GSparseMatrix
    result  = m_test;
    result += m_test;
    testsuite.test_assert(check_matrix(result, 2.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix += GSparseMatrix",
                            "Corrupt GSparseMatrix += GSparseMatrix operator.");
    //
    // GSparseMatrix -= GSparseMatrix
    result  = m_test;
    result -= m_test;
    testsuite.test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix -= GSparseMatrix",
                            "Corrupt GSparseMatrix -= GSparseMatrix operator.");
    //
    // GSparseMatrix *= 3.0
    result  = m_test;
    result *= 3.0;
    testsuite.test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix *= 3.0",
                            "Corrupt GSparseMatrix *= double operator.");
    //
    // GSparseMatrix /= 3.0
    result  = m_test;
    result /= 3.0;
    testsuite.test_assert(check_matrix(result, 1.0/3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix /= 3.0",
                            "Corrupt GSparseMatrix /= double operator.");
    //
    // GSparseMatrix + GSparseMatrix
    result = m_test + m_test;
    testsuite.test_assert(check_matrix(result, 2.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix + GSparseMatrix",
                            "Corrupt GSparseMatrix + GSparseMatrix operator.");
    //
    // GSparseMatrix - GSparseMatrix
    result = m_test - m_test;
    testsuite.test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix - GSparseMatrix",
                            "Corrupt GSparseMatrix - GSparseMatrix operator.");
    //
    // GSparseMatrix - GSparseMatrix
    result = m_test - m_test;
    testsuite.test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix - GSparseMatrix 2",
                            "Corrupt GSparseMatrix - GSparseMatrix operator.");
    //
    // GSparseMatrix * 3.0
    result = m_test * 3.0;
    testsuite.test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix * 3.0",
                            "Corrupt GSparseMatrix * double operator.");
    //
    // GSparseMatrix * -7.12
    result = m_test * -7.12;
    testsuite.test_assert(check_matrix(result, -7.12, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix * -7.12",
                            "Corrupt GSparseMatrix * double operator.");
    //
    // 3.0 * GSparseMatrix
    result = 3.0 * m_test;
    testsuite.test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "3.0 * GSparseMatrix",
                            "Corrupt double * GSparseMatrix operator.");
    //
    // GSparseMatrix / 3.0
    result = m_test / 3.0;
    testsuite.test_assert(check_matrix(result, 1.0/3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix / 3.0",
                            "Corrupt GSparseMatrix / double operator.");

    //
    // GSparseMatrix / 0.0
    double zero = 0.0;
    result = m_test / zero;
    testsuite.test_assert(check_matrix(result, 1.0/zero, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix / 0.0",
                            "Corrupt GSparseMatrix / 0.0");

    //
    // GSparseMatrix.add_col()
    result = m_test;
    GVector v_add(g_rows);
    for (int i = 0; i < g_rows; ++i)
    v_add[i] = 7.0;
    for (int col = 0; col < g_cols; ++col)
    result.add_col(v_add, col);
    testsuite.test_assert(check_matrix(result, 1.0, 7.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix.add_col()",
                            "Corrupt GSparseMatrix.add_col(+v).");
    //
    // GSparseMatrix.add_col()
    for (int col = 0; col < g_cols; ++col)
    result.add_col(-v_add, col);
    testsuite.test_assert(check_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),
                            "GSparseMatrix.add_col() 2",
                            "CCorrupt GSparseMatrix.add_col(-v).");
    //
    // GSparseMatrix.add_col() using a compressed array
    GSparseMatrix compare = m_test;
    GVector v_sparse(g_rows);
    for (int i = 0; i < g_rows; i += 2)  // Build vector for comparison
    v_sparse[i] = 7.0;
    for (int col = 0; col < g_cols; ++col)
        compare.add_col(v_sparse, col);
    result = m_test;
    double* values = new double[g_rows];
    int*    rows   = new int[g_rows];
    int     number = 0;
    for (int i = 0; i < g_rows; i += 2, ++number) {
        values[number] = 7.0;
        rows[number]   = i;
    }
    for (int col = 0; col < g_cols; ++col)
        result.add_col(values, rows, number, col);
    delete [] values;
    delete [] rows;

    testsuite.test_assert(result == compare,
                            "GSparseMatrix.add_col() using a compressed array",
                            "Corrupt GSparseMatrix.add_col(compressed array).");


    testsuite.test_try("matrix+bigger");
    try {
        GSparseMatrix bigger(1000,1000);
        GSparseMatrix result(1,1);
        for (int i = 0; i < 1000; ++i)
        bigger(i,i) = (i+1)*1.1;
        result  = m_test;
        result += bigger;

        testsuite.test_try_failure(); 
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e); 
    }

}


/***************************************************************************
 *                           Test: Matrix functions                        *
 ***************************************************************************/
void test_functions(const GSparseMatrix& m_test, GTestSuite& testsuite)
{
    //Test min
    testsuite.test_try("Test min()");
    try {
        //
        double min = m_test.min();
        if (!check_matrix_min(min)) {
            throw testsuite.exception_failure("Corrupt GSparseMatrix.min() function.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e); 
    }

    //Test max
    testsuite.test_try("Test max()");
    try {
        //
        double max = m_test.max();
        if (!check_matrix_max(max)) {
            throw testsuite.exception_failure("TEST ERROR: Corrupt GSparseMatrix.max() function.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e); 
    }

    //Test sum
    testsuite.test_try("Test sum()");
    try {
        //
        double sum = m_test.sum();
        if (!check_matrix_sum(sum)) {
            throw testsuite.exception_failure("TEST ERROR: Corrupt GSparseMatrix.sum() function.");
        }

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e); 
    }
}


/***************************************************************************
 *                           Test: Matrix comparison                       *
 ***************************************************************************/
void test_comparison(const GSparseMatrix& m_test, GTestSuite& testsuite)
{

    GSparseMatrix bigger(1000,1000);
    for (int i = 0; i < 1000; ++i)
        bigger(i,i) = (i+1)*1.1;
    //
    testsuite.test_assert((m_test == m_test) == 1,
                            "GSparseMatrix == GSparseMatrix",
                            "Corrupt GSparseMatrix == GSparseMatrix operator.");
    //
    GSparseMatrix test(g_rows,g_cols);
    testsuite.test_assert((m_test == test) == 0,
                            "GSparseMatrix == GSparseMatrix",
                            "Corrupt GSparseMatrix == GSparseMatrix operator.");
    //
    testsuite.test_assert((m_test == bigger) == 0,
                            "GSparseMatrix == GSparseMatrix",
                            "Corrupt GSparseMatrix == GSparseMatrix operator.");
    //
    testsuite.test_assert((m_test != m_test) == 0,
                            "GSparseMatrix != GSparseMatrix",
                            "Corrupt GSparseMatrix != GSparseMatrix operator.");
    //
    testsuite.test_assert((m_test != test) == 1,
                            "GSparseMatrix != GSparseMatrix",
                            "Corrupt GSparseMatrix != GSparseMatrix operator.");
    //
    testsuite.test_assert((m_test != bigger) == 1,
                            "GSparseMatrix != GSparseMatrix",
                            "Corrupt GSparseMatrix != GSparseMatrix operator.");
}


/***************************************************************************
 *                     Test: Conversion between matrix types               *
 ***************************************************************************/
void test_conversion(GTestSuite& testsuite)
{

    // Test 1
    testsuite.test_try("Conversion");
    try {
        //
        // Setup a symmetric sparse matrix
        int           num = 10;
        GSparseMatrix sparse(num,num);
        for (int i = 0; i < num; ++i)
            sparse(i,i) = (i+1);
        for (int j = 1; j < num; ++j) {
            sparse(0,j) = (j+1);
            sparse(j,0) = (j+1);
        }
    
        //
        // Convert GSparseMatrix matrix into GSymMatrix object
        // GSymMatrix converted = sym_matrix(sparse);
        GSymMatrix converted(1,1);
        testsuite.test_try("Convert GSparseMatrix matrix into GSymMatrix object");
        try {
            converted = GSymMatrix(sparse);
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
    
        //
        // Convert GSymMatrix back to GSparseMatrix matrix
        GSparseMatrix back_convert(1,1);
        testsuite.test_try("Convert GSymMatrix back to GSparseMatrix matrix");
        try {
            back_convert = sparse_matrix(converted);
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
        
        //
        // Compare back converted matrix to original one. They should be identical
        testsuite.test_assert(sparse == back_convert,"Compare back converted matrix to original one","Unable to convert matrixes (symmetric).");
    
        //
        // Determine the fill of the matrix. It should be 0.28
        double fill = back_convert.fill();
        testsuite.test_assert(abs(fill-0.28) <= 1.0e-15,"Determine the fill of the matrix. It should be 0.28","Bad fill "+str(fill)+" determined (expected 1.0)");
    
        //
        // Extract lower triangle and check values
        #if 0
        GMatrix lower = sparse.extract_lower_triangle();
        int ok = 1;
        for (int j = 1; j < num; ++j) {
        for (int i = 0; i < j; ++i) {
            if (lower(i,j) != 0.0)
            ok = 0;
        }
        }
        for (int j = 0; j < num; ++j) {
        for (int i = j; i < num; ++i) {
            if (lower(i,j) != sparse(i,j))
            ok = 0;
        }
        }
        if (!ok) {
        cout << endl << "TEST ERROR: Corrupt extract_lower_triangle." << endl;
        cout << "Original matrix " << sparse << endl;
        cout << "Lower triangle matrix " << lower << endl;
        throw;
        }
        cout << ".";
        //
        // Extract upper triangle and check values
        GMatrix upper = sparse.extract_upper_triangle();
        ok = 1;
        for (int j = 0; j < num; ++j) {
        for (int i = j+1; i < num; ++i) {
            if (upper(i,j) != 0.0)
            ok = 0;
        }
        }
        for (int j = 0; j < num; ++j) {
        for (int i = 0; i <= j; ++i) {
            if (upper(i,j) != sparse(i,j))
            ok = 0;
        }
        }
        if (!ok) {
        cout << endl << "TEST ERROR: Corrupt extract_upper_triangle." << endl;
        cout << "Original matrix " << sparse << endl;
        cout << "Upper triangle matrix " << upper << endl;
        throw;
        }
        cout << ".";
        #endif
        //
        // Now make the matrix unsymmetric
        sparse(0,num-1) = 1000.0;
        //
        // Try converting now into GSymMatrix object (this should fail)
        testsuite.test_try("Try converting now into GSymMatrix object");
        try {
            converted = sym_matrix(sparse);
            testsuite.test_try_failure();
        }
        catch (GException::matrix_not_symmetric &e) {
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
    
        //
        // Convert matrix into GMatrix object
        GMatrix full(1,1);
        testsuite.test_try("Convert matrix into GMatrix object");
        try {
            full = matrix(sparse);
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }

        //
        // Convert full matrix back to GSparseMatrix
        testsuite.test_try("Convert full matrix back to GSparseMatrix");
        try {
            back_convert = sparse_matrix(full);
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }

        //
        // Compare back converted matrix to original one. They should be identical
        testsuite.test_assert(sparse == back_convert,
                            "Compare back conUnable to convert matrixes (sparse).verted matrix to original one. They should be identical",
                            "Unable to convert matrixes (sparse).");

        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
        // Return
    return;
}


/***************************************************************************
 *                        Test: Cholesky factorisation                     *
 ***************************************************************************/
void test_cholesky(GTestSuite& testsuite)
{
    // Test 1
    testsuite.test_try( "Test 1");
    try {
        //
        // Setup matrix for Cholesky decomposition
        GSparseMatrix m_chol_test(5,5);
        m_chol_test(0,0) = 1.0;
        m_chol_test(0,1) = 0.2;
        m_chol_test(0,2) = 0.2;
        m_chol_test(0,3) = 0.2;
        m_chol_test(0,4) = 0.2;
        m_chol_test(1,0) = 0.2;
        m_chol_test(2,0) = 0.2;
        m_chol_test(3,0) = 0.2;
        m_chol_test(4,0) = 0.2;
        m_chol_test(1,1) = 1.0;
        m_chol_test(2,2) = 1.0;
        m_chol_test(3,3) = 1.0;
        m_chol_test(4,4) = 1.0;
        //
        // Try to solve now (should not work)
        testsuite.test_try("Try to solve now");
        try {
            GVector v_test12(5);
            v_test12 = m_chol_test.cholesky_solver(v_test12);
            testsuite.test_try_failure();
        }
        catch (GException::matrix_not_factorised) {
            testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }

        //
        // Test Cholesky decomposition
        GSparseMatrix cd(1,1);
        testsuite.test_try("Test Cholesky decomposition");
        try {
            cd = cholesky_decompose(m_chol_test);
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        //
        // Test inplace Cholesky decomposition
        testsuite.test_try("Test inplace Cholesky decomposition");
        try {
        cd = m_chol_test;
        cd.cholesky_decompose();
        testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GVector e0(5);
        GVector a0(5);
        GVector s0;
        double res;
        //
        // Test Cholesky solver
        testsuite.test_try("Test Cholesky solver");
        try {
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
            if (res < 1.0e-15) {
                #if defined(DUMP_RESIDUALS)
                std::cout << " Res(S0)=" << res << " ";
                #else
                std::cout << ".";
                #endif
            }
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        //
        testsuite.test_try("Test Cholesky solver 2");
        try {
            e0[0] = 0.0;
            e0[1] = 1.0;
            a0[0] = 0.2;
            a0[1] = 1.0;
            a0[2] = 0.0;
            a0[3] = 0.0;
            a0[4] = 0.0;
            s0    = cd.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S1)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test Cholesky solver 3");
        try {
            e0[1] = 0.0;
            e0[2] = 1.0;
            a0[1] = 0.0;
            a0[2] = 1.0;
            s0    = cd.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S2)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test Cholesky solver 4");
        try {
            e0[2] = 0.0;
            e0[3] = 1.0;
            a0[2] = 0.0;
            a0[3] = 1.0;
            s0    = cd.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S3)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test Cholesky solver 4");
        try {
            e0[3] = 0.0;
            e0[4] = 1.0;
            a0[3] = 0.0;
            a0[4] = 1.0;
            s0    = cd.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S4)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GSparseMatrix m_chol_test_zero(6,6);
	//
	// Setup matrix for Cholesky decomposition with zero row/col
        testsuite.test_try("Setup matrix for Cholesky decomposition with zero row/col");
        try {
            m_chol_test_zero(0,0) = 1.0;
            m_chol_test_zero(0,1) = 0.2;
            m_chol_test_zero(0,2) = 0.2;
            m_chol_test_zero(0,4) = 0.2;
            m_chol_test_zero(0,5) = 0.2;
            m_chol_test_zero(1,0) = 0.2;
            m_chol_test_zero(2,0) = 0.2;
            m_chol_test_zero(4,0) = 0.2;
            m_chol_test_zero(5,0) = 0.2;
            m_chol_test_zero(1,1) = 1.0;
            m_chol_test_zero(2,2) = 1.0;
            m_chol_test_zero(4,4) = 1.0;
            m_chol_test_zero(5,5) = 1.0;

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GSparseMatrix cd_zero = m_chol_test_zero;
	//
        // Test compressed Cholesky decomposition
        testsuite.test_try("Test compressed Cholesky decomposition");
        try {
            cd_zero.cholesky_decompose();
    
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        // Test compressed Cholesky solver
        testsuite.test_try("Test compressed Cholesky solver 5");
        try {
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
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S0Z)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver 6");
        try {
            e0[0] = 0.0;
            e0[1] = 1.0;
            a0[0] = 0.2;
            a0[1] = 1.0;
            a0[2] = 0.0;
            a0[4] = 0.0;
            a0[5] = 0.0;
            s0    = cd_zero.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S1Z)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver 7");
        try {
            e0[1] = 0.0;
            e0[2] = 1.0;
            a0[1] = 0.0;
            a0[2] = 1.0;
            s0    = cd_zero.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
                #if defined(DUMP_RESIDUALS)
                cout << "Res(S2Z)=" << res << " ";
                #else
                cout << ".";
                #endif
            else {
                    throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver 8");
        try {
            e0[2] = 0.0;
            e0[4] = 1.0;
            a0[2] = 0.0;
            a0[4] = 1.0;
            s0    = cd_zero.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S3Z)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }
            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver 8");
        try {
            e0[4] = 0.0;
            e0[5] = 1.0;
            a0[4] = 0.0;
            a0[5] = 1.0;
            s0    = cd_zero.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S4Z)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GSparseMatrix m_chol_test_zero2(6,5);
	//
	// Setup matrix for Cholesky decomposition with zero row/col (unsymmetric case)
        testsuite.test_try("Setup matrix for Cholesky decomposition with zero row/col (unsymmetric case)");
        try {
            m_chol_test_zero2(0,0) = 1.0;
            m_chol_test_zero2(0,1) = 0.2;
            m_chol_test_zero2(0,2) = 0.2;
            m_chol_test_zero2(0,3) = 0.2;
            m_chol_test_zero2(0,4) = 0.2;
            m_chol_test_zero2(1,0) = 0.2;
            m_chol_test_zero2(2,0) = 0.2;
            m_chol_test_zero2(4,0) = 0.2;
            m_chol_test_zero2(5,0) = 0.2;
            m_chol_test_zero2(1,1) = 1.0;
            m_chol_test_zero2(2,2) = 1.0;
            m_chol_test_zero2(4,3) = 1.0;
            m_chol_test_zero2(5,4) = 1.0;

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GSparseMatrix cd_zero2 = m_chol_test_zero2;
        //
        // Test compressed Cholesky decomposition (unsymmetric case)
        testsuite.test_try("Test compressed Cholesky decomposition (unsymmetric case)");
        try {
            cd_zero2.cholesky_decompose();

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        // Test compressed Cholesky solver (unsymmetric case)
        testsuite.test_try("Test compressed Cholesky solver (unsymmetric case)");
        try {
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
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S0ZU)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        testsuite.test_try("Test compressed Cholesky solver (unsymmetric case) 2");
        try {
            e0[0] = 0.0;
            e0[1] = 1.0;
            a0[0] = 0.2;
            a0[1] = 1.0;
            a0[2] = 0.0;
            a0[4] = 0.0;
            a0[5] = 0.0;
            s0    = cd_zero2.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S1ZU)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver (unsymmetric case) 3");
        try {
            e0[1] = 0.0;
            e0[2] = 1.0;
            a0[1] = 0.0;
            a0[2] = 1.0;
            s0    = cd_zero2.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S2ZU)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver (unsymmetric case) 4");
        try {
            e0[2] = 0.0;
            e0[3] = 1.0;
            a0[2] = 0.0;
            a0[4] = 1.0;
            s0    = cd_zero2.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S3ZU)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        testsuite.test_try("Test compressed Cholesky solver (unsymmetric case) 5");
        try {
            e0[3] = 0.0;
            e0[4] = 1.0;
            a0[4] = 0.0;
            a0[5] = 1.0;
            s0    = cd_zero2.cholesky_solver(a0);
            res   = max(abs(s0-e0));
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(S4ZU)=" << res << " "; 
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function.\nResidual vector (all elements should be zero):"+s0.print());
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        GSparseMatrix m_chol_test_inv = m_chol_test;
        GSparseMatrix ci_product(1,1);
        GSparseMatrix ci_residuals(1,1);
        GSparseMatrix unit(5,5);
	//
	// Test Cholesky inverter (inplace)
        testsuite.test_try("Test Cholesky inverter (inplace)");
        try {
            unit(0,0) = 1.0;
            unit(1,1) = 1.0;
            unit(2,2) = 1.0;
            unit(3,3) = 1.0;
            unit(4,4) = 1.0;
            m_chol_test_inv.cholesky_invert();
            ci_product   = m_chol_test * m_chol_test_inv;
            ci_residuals = ci_product - unit;
            res = (abs(ci_residuals)).max();
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(CI)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("TEST ERROR: Corrupt GSparseMatrix.cholesky_invert() function.");
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
	// Test Cholesky inverter
        testsuite.test_try("Test Cholesky inverter");
        try {
            m_chol_test_inv = cholesky_invert(m_chol_test);
            ci_product      = m_chol_test * m_chol_test_inv;
            ci_residuals    = ci_product - unit;
            res             = (abs(ci_residuals)).max();
            if (res < 1.0e-15)
            #if defined(DUMP_RESIDUALS)
            cout << "Res(CI2)=" << res << " ";
            #else
            cout << ".";
            #endif
            else {
                throw testsuite.exception_failure("TEST ERROR: Corrupt cholesky_invert(GSparseMatrix) function.");
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

	//
        // Test Cholesky inverter for compressed matrix
        testsuite.test_try("Test Cholesky inverter for compressed matrix");
        try {
            unit = GSparseMatrix(6,6);
            unit(0,0) = 1.0;
            unit(1,1) = 1.0;
            unit(2,2) = 1.0;
            unit(4,4) = 1.0;
            unit(5,5) = 1.0;
            GSparseMatrix m_chol_test_zero_inv = m_chol_test_zero;
            m_chol_test_zero_inv.cholesky_invert();
            GSparseMatrix ciz_product   = m_chol_test_zero * m_chol_test_zero_inv;
            GSparseMatrix ciz_residuals = ciz_product - unit;
            res = (abs(ciz_residuals)).max();
            if (res < 1.0e-15) {
                #if defined(DUMP_RESIDUALS)
                std::cout << "Res(CIZ)=" << res << " ";
                #else
                std::cout << ".";
                #endif
            }
            else {
                throw testsuite.exception_failure("TEST ERROR: Corrupt compressed GSparseMatrix.cholesky_invert()");
            }

            testsuite.test_try_success();
        }
        catch(exception& e)
        {
            testsuite.test_try_failure(e); 
        }

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e); 
    }

    // Return
    return;
}


/***************************************************************************
 *                            Test: Heavy calculus                         *
 ***************************************************************************/
void test_heavy(GTestSuite& testsuite,int block = 512)
{
  // Allocate some variables
  clock_t t_start;
  double  t_elapse;
  int     number;
  int     num_cols;
  int     i_min, i_max;

  GSparseMatrix large(1,1);
    //
    // Fill 10000 x 10000 matrix with 10000 values
    testsuite.test_try(" Fill 10000 x 10000 matrix with 10000 values");
    try {
            t_start = clock();
            number  = 10000;
            large=GSparseMatrix (number,number);
            large.set_mem_block(block);
            for (int i = 0; i < number; ++i)
            large(i,i) = (i+1)*3.14;
            t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
            #if defined(DUMP_TIMING)
            cout << endl << " - Fill of 10000 values needed " << t_elapse << 
                            " sec (reference ~ 0.2 sec)";
            #endif

          testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e); 
    }

    //
    // Modify 10000 x 10000 matrix with 10000 (exisiting) values
    testsuite.test_try("Modify 10000 x 10000 matrix with 10000 (exisiting) values");
    try {
        t_start = clock();
        for (int i = 0; i < number; ++i)
        large(i,i) = (i+1)*1.0;
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        #if defined(DUMP_TIMING)
        cout << endl << " - Modification of 10000 values needed " << t_elapse << 
                        " sec (reference ~ 0 sec)";
        #endif

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e); 
    }

    GSparseMatrix diag(1,1);
    GVector column;
    //
    // Insert columns into 10000 x 10000 matrix
    testsuite.test_try("Insert columns into 10000 x 10000 matrix");
    try {
        number   = 10000;
        num_cols = 10000;
        diag = GSparseMatrix(number,number);
        for (int j = 10; j < number-10; ++j)
        diag(j,j) = 3.14;
        large = diag;
        large.set_mem_block(block);
        column= GVector(number);
        t_start = clock();
        for (int j = 0; j < num_cols; ++j) {
            i_min = (j < 2) ? 0 : j-2;
            i_max = (j > num_cols-2) ? num_cols : j+2;
            column = 0.0;
            for (int i = i_min; i < i_max; ++i)
                column[i] = (i+1)*0.01;
            large.add_col(column, j);
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        #if defined(DUMP_TIMING)
        cout << endl << " - " << num_cols << " columns adding needed " << t_elapse << 
                        " sec (reference ~ 4.6 sec)";
        #endif

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e); 
    }

    GSparseMatrix stack(1,1);
    //
    // Insert columns into 10000 x 10000 matrix using matrix stack
    testsuite.test_try("Insert columns into 10000 x 10000 matrix using matrix stack");
    try {
        stack=GSparseMatrix(number,number);
        stack = diag;
        stack.set_mem_block(block);
        for (int j = 200; j < 400; ++j)
        stack(j,j) = 3.14;
        t_start = clock();
        column = GVector(number);
        stack.stack_init(10000, 1000);
        for (int j = 0; j < num_cols; ++j) {
            i_min = (j < 2)          ? 0 : j-2;
            i_max = (j > num_cols-2) ? num_cols : j+2;
            column = 0.0;
            for (int i = i_min; i < i_max; ++i)
                column[i] = (i+1)*0.01;
            stack.add_col(column, j);
        }
        stack.stack_destroy();
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        #if defined(DUMP_TIMING)
        cout << endl << " - " << num_cols << " columns stack-adding needed " << t_elapse << 
                        " sec (reference ~ 1 sec)";
        #endif

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e); 
    }

    //
    // Check that both are identical
    testsuite.test_assert(large==stack,"Check that both are identical","Stack-based insertion corrupted.");
    std::cout<<"large : "<<large<<" et stack : "<<stack<<std::endl;
    GSparseMatrix sum(1,1);
    //
    // Subsequent addition and subtractions
    testsuite.test_try("Subsequent addition and subtractions");
    try {
        number  = 300;
        t_start = clock();
        sum=GSparseMatrix(number,number);
        GSparseMatrix add(number,number);
        GSparseMatrix sub(number,number);
        sum.set_mem_block(block);
        add.set_mem_block(block);
        sub.set_mem_block(block);
        for (int i = 0; i < number; ++i) {
        add(i,number-i-1) = (i*i+1)*0.01;
        sum += add;
        sub(i,i) = (i+1)*1.0;
        sum -= sub;
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        #if defined(DUMP_TIMING)
        cout << endl << " - Matrix adding and subtraction 300 times " << t_elapse << 
                        " sec (reference ~ 1.4 sec)";
        #endif

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e);
    }

    //
    // Subsequent multiplications
    testsuite.test_try("Subsequent multiplications");
    try {
        GSparseMatrix fac(number,number);
        fac.set_mem_block(block);
        t_start = clock();
        for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < number; ++j)
            fac(j,j) = (i+1)*1.0;
        sum *= fac;
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
        #if defined(DUMP_TIMING)
        cout << endl << " - Matrix multiplication 1 times " << t_elapse << 
                        " sec (reference ~ 3.1 sec)";
        #endif

        testsuite.test_try_success();
    }
    catch(exception& e)
    {
        testsuite.test_try_failure(e);
    }

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
        //Create a test suites
    GTestSuites testsuites("GSparseMatrix class testing");

        // Create a test suite
    GTestSuite * testsuite = NULL;

    // Set test matrix and vector
    testsuite = new GTestSuite("Set test matrix and vector");
    testsuites.append(*testsuite);
    GSparseMatrix m_test = set_matrix(*testsuite);
    GVector       v_test = set_vector(*testsuite);
    testsuite->end_test();

    // Execute the tests

    //Test output
    testsuite = new GTestSuite("Output test matrix");
    testsuites.append(*testsuite);
    test_output(m_test,*testsuite);
    testsuite->end_test();

    //Test allocate
    testsuite = new GTestSuite("Allocate zero matrix");
    testsuites.append(*testsuite);
    test_allocate(*testsuite);
    testsuite->end_test();

    //Test assign values
    testsuite = new GTestSuite("Assign matrix values");
    testsuites.append(*testsuite);
    test_assign_values(m_test,*testsuite);
    testsuite->end_test();

    //Test copy constructor
    testsuite = new GTestSuite("Define matrix using copy constructor");
    testsuites.append(*testsuite);
    test_copy_constructor(m_test,*testsuite);
    testsuite->end_test();

    //Test assign matrix
    testsuite = new GTestSuite("Matrix assignment");
    testsuites.append(*testsuite);
    test_assign_matrix(m_test,*testsuite);
    testsuite->end_test();

    //Test transpose
    testsuite = new GTestSuite("Matrix transpose");
    testsuites.append(*testsuite);
    test_transpose(m_test,*testsuite);
    testsuite->end_test();

    //Test matrix*vector
    testsuite = new GTestSuite("Matrix*Vector multiplication");
    testsuites.append(*testsuite);
    test_matrix_vector(m_test, v_test, *testsuite);
    testsuite->end_test();

    //Test matrix*matrix
    testsuite = new GTestSuite("Matrix*Matrix multiplication");
    testsuites.append(*testsuite);
    test_matrix_matrix(m_test,*testsuite);
    testsuite->end_test();

    //Test arithmetics
    testsuite = new GTestSuite("Arithmetics");
    testsuites.append(*testsuite);
    test_arithmetics(m_test,*testsuite);
    testsuite->end_test();

    //Test functions
    testsuite = new GTestSuite("Matrix functions");
    testsuites.append(*testsuite);
    test_functions(m_test,*testsuite);
    testsuite->end_test();

    //Test comparison
    testsuite = new GTestSuite("Comparison");
    testsuites.append(*testsuite);
    test_comparison(m_test,*testsuite);
    testsuite->end_test();

    //Test conversion
    testsuite = new GTestSuite("Conversion");
    testsuites.append(*testsuite);
    test_conversion(*testsuite);
    testsuite->end_test();

    //Test cholesky
    testsuite = new GTestSuite("Cholesky decomposition, solver and inverter");
    testsuites.append(*testsuite);
    test_cholesky(*testsuite);
    testsuite->end_test();

    //Test heavy
    testsuite = new GTestSuite("Heavy calculus");
    testsuites.append(*testsuite);
    test_heavy(*testsuite);
    testsuite->end_test();

    //save xml report
    testsuites.save("reports/GSparseMatrix.xml");
    // Return
    return 0;
}
