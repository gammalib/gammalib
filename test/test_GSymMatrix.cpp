/***************************************************************************
 *              test_GSymMatrix.cpp  -  test GSymMatrix class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
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
 * @file test_GSymMatrix.cpp
 * @brief Testing of GSymMatrix class
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include <cmath>
#include "GammaLib.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Globals ____________________________________________________________ */
double g_matrix[] = {4.0, 1.0, 2.0, 1.0, 5.0, 3.0, 2.0, 3.0, 6.0};
double g_vector[] = {1.0, 2.0, 3.0};
int    g_rows    = 3;
int    g_cols    = 3;


/***************************************************************************
 *                              Set test matrix                            *
 ***************************************************************************/
GSymMatrix set_matrix(GTestSuite& testsuite)
{
    testsuite.test_try("Set test matrix");
    try {
        GSymMatrix matrix(g_rows,g_cols);
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
    return GSymMatrix(g_rows,g_cols);
}


/***************************************************************************
 *                       Set test matrix with zero lines                   *
 ***************************************************************************/
GSymMatrix set_matrix_zero(void)
{
  try {
    GSymMatrix matrix(g_rows+1,g_cols+1);
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
    return matrix;
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to set test matrix." << endl;
	throw;
  }
  return GSymMatrix(g_rows+1,g_cols+1);
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
int check_matrix(const GSymMatrix& m, const double scale, const double add)
{
  int result = 1;
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
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check test matrix." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                          Check full test matrix                         *
 ***************************************************************************/
int check_matrix(const GMatrix& m, const double scale, const double add)
{
  int result = 1;
  try {
    for (int row = 0; row < g_rows; ++row) {
      for (int col = 0; col < g_cols; ++col) {
	    double value = g_matrix[col+row*g_cols] * scale + add;
	    if (m(row,col) != value) {
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
 *                     Check full test matrix (lower triangle)             *
 ***************************************************************************/
int check_matrix_lt(const GMatrix& m, const double scale, const double add)
{
  int result = 1;
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
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check test matrix." << endl;
	throw;
  }
  return result;
}


/***************************************************************************
 *                     Check full test matrix (upper triangle)             *
 ***************************************************************************/
int check_matrix_ut(const GMatrix& m, const double scale, const double add)
{
  int result = 1;
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
	  for (int col = 0; col < g_cols; ++col)
	    value += g_matrix[col+row*g_cols] * g_vector[col];
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
int check_matrix_matrix(const GSymMatrix& m)
{
  if (m.rows() != g_rows || m.cols() != g_rows)
    return 0;
  int result = 1;
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
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to check matrix*matrix product." << endl;
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
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col) {
	  if (g_matrix[col+row*g_cols] < value)
		value = g_matrix[col+row*g_cols];
	}
  }
  return (min == value);
}


/***************************************************************************
 *                              Check matrix max                           *
 ***************************************************************************/
int check_matrix_max(const double max)
{
  double value = g_matrix[0];
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col) {
	  if (g_matrix[col+row*g_cols] > value)
		value = g_matrix[col+row*g_cols];
	}
  }
  return (max == value);
}


/***************************************************************************
 *                              Check matrix sum                           *
 ***************************************************************************/
int check_matrix_sum(const double sum)
{
  double value = 0.0;
  for (int row = 0; row < g_rows; ++row) {
    for (int col = 0; col < g_cols; ++col)
      value += g_matrix[col+row*g_cols];
  }
  return (sum == value);
}


/***************************************************************************
 *                                Test: Output                             *
 ***************************************************************************/
void test_output(const GSymMatrix& m_test,GTestSuite& testsuite)
{
    testsuite.test_try("Output test matrix");
    try {
        cout << m_test << endl;
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{ 
    //Create a test suites
    GTestSuites testsuites("GSymMatrix class testing");
    
    // Create a test suite
    GTestSuite testsuite("GSymMatrix");
    
    testsuites.append(testsuite);
    
    // Set test matrix and vector
    GSymMatrix m_test = set_matrix(testsuite);
    GVector    v_test = set_vector(testsuite);
    
    // Set bigger matrix (only used for collision, don't care about content)
    GSymMatrix bigger(g_rows+1,g_cols+1);
    
    // Prepare result matrix
    GSymMatrix result = m_test;
    
    // Execute the tests
    test_output(m_test,testsuite);
    
    // Test 1: Allocate zero matrix
    testsuite.test_try("Allocate zero matrix");
    try {
        GSymMatrix test1(0,0);
    
        testsuite.test_try_failure(); 
    }
    catch (GException::empty &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 2: Allocate too large matrix
    /*
    try {
        cout << "GSymMatrix - Test 2: Allocate too large matrix: ";
        GSymMatrix test2(100000000,100000000);
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
        GSymMatrix test3(3,3);	
        test3(1,1) = 1.0;
        testsuite.test_assert((test3(0,0) == 0.0 && test3(1,0) == 0.0 && test3(2,0) == 0.0 &&
                            test3(0,1) == 0.0 && test3(1,1) == 1.0 && test3(2,1) == 0.0 &&
                            test3(0,2) == 0.0 && test3(1,2) == 0.0 && test3(2,2) == 0.0),
                            "Test 3: Assign matrix values"," Unable to set matrix value.");
    
        #if defined(G_RANGE_CHECK)
        try {
            testsuite.test_try("Range check");
            GSymMatrix test3(3,3);
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
        testsuite.test_try("Test 4: Matrix copy constructor");
        try {
                //
                GSymMatrix test4 = m_test;
            if (!check_matrix(test4, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt copy constructor");
                }
    
                testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
    
        // Test 5: Matrix assignment
        testsuite.test_try("Test 5: Matrix assignment");
        try {
                //
                // GSymMatrix = GSymMatrix
                result = m_test;
            if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0))
            {
                throw testsuite.exception_failure("Corrupt assignment");
            }
                //
                // GSymMatrix = GSymMatrix (bigger matrix)
                result = bigger;
    
                testsuite.test_try_success();
        }
        catch (exception &e) {
            testsuite.test_try_failure(e);
        }
    
    // Test 6: Matrix*Vector multiplication
        testsuite.test_try("Test 6: Matrix*Vector multiplication");
    try {
        //
        GVector v_test6 = m_test*v_test;
        if (!check_matrix_vector(v_test6) || !check_matrix(m_test, 1.0, 0.0)) {
            throw testsuite.exception_failure("Corrupt Matrix*Vector multiplication");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    testsuite.test_try("Bigger*Vector");
    try {
            GVector v_test6 = bigger*v_test;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_vector_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 7: Matrix*Matrix multiplication
    testsuite.test_try("Test 7: Matrix*Matrix multiplication");
    try {
        //
        GSymMatrix m_test7 = m_test*m_test;
        if (!check_matrix_matrix(m_test7) || !check_matrix(m_test, 1.0, 0.0))
        {
            throw testsuite.exception_failure("Corrupt Matrix*Matrix multiplication");
        }
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
        throw;
    }
    
    testsuite.test_try("Bigger*Matrix multiplication");
    try {
            GSymMatrix m_test7 = m_test*bigger;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    testsuite.test_try("Bigger*Matrix multiplication 2");
    try {
            GSymMatrix m_test7 = bigger*m_test;
            testsuite.test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 8: Assignment and arithmetics
    testsuite.test_try("Test 8: Assignment and arithmetics");
    try {
            //
            // -GSymMatrix
            result = -m_test;
            if (!check_matrix(result, -1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt -GSymMatrix operator.");
            }
            //
            // GSymMatrix += GSymMatrix
            result  = m_test;
            result += m_test;
            if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix += GSymMatrix operator.");
            }
            //
            // GSymMatrix -= GSymMatrix
            result  = m_test;
            result -= m_test;
            if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix -= GSymMatrix operator.");
            }
            //
            // GSymMatrix *= 3.0
            result  = m_test;
            result *= 3.0;
            if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix *= double operator.");
            }
            //
            // GSymMatrix /= 3.0
            result  = m_test;
            result /= 3.0;
            if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix /= double operator.");
            }
    
            //
            // GSymMatrix + GSymMatrix
            result = m_test + m_test;
            if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix + GSymMatrix operator.");
            }
            //
            // GSymMatrix - GSymMatrix
            result = m_test - m_test;
            if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix - GSymMatrix operator.");
            }
    
            //
            // GSymMatrix - GSymMatrix
            result = m_test - m_test;
            if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix - GSymMatrix operator.");
            }
    
            //
            // GSymMatrix * 3.0
            result = m_test * 3.0;
            if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0))
            {
                throw testsuite.exception_failure("Corrupt GSymMatrix * double operator.");
            }
            //
            // 3.0 * GSymMatrix
            result = 3.0 * m_test;
            if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt double * GSymMatrix operator.");
            }
            //
            // GSymMatrix / 3.0
            result = m_test / 3.0;
            if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix / double operator.");
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
    
    // Test 9: Matrix functions
    testsuite.test_try("Test 9: Matrix functions");
    try {
            //
            double min = m_test.min();
            if (!check_matrix_min(min)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix.min() function.");
            }
            //
            double max = m_test.max();
            if (!check_matrix_max(max)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix.max() function.");
            }
            //
            double sum = m_test.sum();
            if (!check_matrix_sum(sum)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix.sum() function.");
            }
    
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 10: Matrix comparison
    testsuite.test_try("Test 10: Matrix comparison");
    try {
            //
            if ((m_test == m_test) != 1) {
                throw testsuite.exception_failure("Corrupt GSymMatrix == GSymMatrix operator.");
            }
            //
            GSymMatrix m_test10(g_rows,g_cols);
            if ((m_test == m_test10) != 0) {
                throw testsuite.exception_failure("Corrupt GSymMatrix == GSymMatrix operator.");
            }
            //
            if ((m_test == bigger) != 0) {
                throw testsuite.exception_failure("Corrupt GSymMatrix == GSymMatrix operator.");
            }
            //
            if ((m_test != m_test) != 0) {
                throw testsuite.exception_failure("Corrupt GSymMatrix != GSymMatrix operator.");
            }
        //
            if ((m_test != m_test10) != 1) {
                throw testsuite.exception_failure("Corrupt GSymMatrix != GSymMatrix operator.");
            }
            //
            if ((m_test != bigger) != 1) {
                throw testsuite.exception_failure("Corrupt GSymMatrix != GSymMatrix operator.");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 11: Transformations  
    testsuite.test_try("Test 11: Transformations");
    try {
        //
            // transpose(GSymMatrix)
            result = transpose(m_test);
            if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt transpose(GSymMatrix) function.");
            }
            //
            // GSymMatrix.transpose()
            result = m_test;
            result.transpose();
    
            if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt GSymMatrix.transpose() function.");
            }
            //
            // convert_to_full()
            GMatrix m_test11 = matrix(m_test);
            if (!check_matrix(m_test11, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt convert_to_full() operator.");
            }
    
            //
            // extract_lower_triangle()
            m_test11 = m_test.extract_lower_triangle();
            if (!check_matrix_lt(m_test11, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt extract_lower_triangle() function.");
            }
    
            //
            // extract_upper_triangle()
            m_test11 = m_test.extract_upper_triangle();
            if (!check_matrix_ut(m_test11, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
                throw testsuite.exception_failure("Corrupt extract_upper_triangle() function.");
            }
    
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }
    
    // Test 12: Perform Cholesky decomposition
    testsuite.test_try("Test 12: Perform Cholesky decomposition");
    try {
            //
            // Test Cholesky decomposition
            GSymMatrix cd           = cholesky_decompose(m_test);
            GMatrix    cd_lower     = cd.extract_lower_triangle();
            GMatrix    cd_upper     = transpose(cd_lower);
            GMatrix    cd_product   = cd_lower * cd_upper;
            GMatrix    cd_residuals = matrix(m_test) - cd_product;
            //
            double res = (abs(cd_residuals)).max();
            if (res < 1.0e-15){
            //cout << "Res(CD)=" << res << " ";
            }
            else{
                throw testsuite.exception_failure("Corrupt cholesky_decompose(GSymMatrix) function.");
            }
    
            //
            // Test compressed Cholesky decomposition
            GSymMatrix m_test_zero       = set_matrix_zero();
            GSymMatrix cd_zero           = cholesky_decompose(m_test_zero);
            GMatrix    cd_zero_lower     = cd_zero.extract_lower_triangle();
            GMatrix    cd_zero_upper     = transpose(cd_zero_lower);
            GMatrix    cd_zero_product   = cd_zero_lower * cd_zero_upper;
            GMatrix    cd_zero_residuals = matrix(m_test_zero) - cd_zero_product;
            //
            res = (abs(cd_zero_residuals)).max();
            if (res < 1.0e-15){
            //cout << "Res(CDZ)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt cholesky_decompose(GSymMatrix) function. 2");
            }
    
            //
            // Test Cholesky inplace decomposition
            GSymMatrix m_test12 = m_test;
            m_test12.cholesky_decompose();
            GMatrix cd_lower2 = m_test12.extract_lower_triangle();
            if (cd_lower2 != cd_lower) {
                throw testsuite.exception_failure("Corrupt cholesky_decompose() function. 3");
            }
    
            //
            // Test Cholesky solver
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
            if (res < 1.0e-15){
            //cout << "Res(S0)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function. 4");
            }
            //
            e0[0] = 0.0;
            e0[1] = 1.0;
            e0[2] = 0.0;
            a0[0] = g_matrix[1];
            a0[1] = g_matrix[4];
            a0[2] = g_matrix[7];
            s0 = cd.cholesky_solver(a0) - e0;
            res = max(abs(s0));
            if (res < 1.0e-15){
            //cout << "Res(S1)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function. 5");
            }
            //
            e0[0] = 0.0;
            e0[1] = 0.0;
            e0[2] = 1.0;
            a0[0] = g_matrix[2];
            a0[1] = g_matrix[5];
            a0[2] = g_matrix[8];
            s0 = cd.cholesky_solver(a0) - e0;
            res = max(abs(s0));
            if (res < 1.0e-15){
            //cout << "Res(S2)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt cholesky_solver(GVector) function. 6");
            }
    
            //
            // Test compressed Cholesky solver
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
            if (res < 1.0e-15){
            //cout << "Res(S0Z)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt compressed cholesky_solver(GVector) function. 7");
            }
    
            //
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
            if (res < 1.0e-15){
            //cout << "Res(S1Z)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt compressed cholesky_solver(GVector) function. 8");
            }
            //
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
            if (res < 1.0e-15){
            // cout << "Res(S2Z)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt compressed cholesky_solver(GVector) function. 9");
            }
        //
            // Test Cholesky inverter
            GSymMatrix unit(g_rows,g_cols);
            unit(0,0) = unit(1,1) = unit(2,2) = 1.0;
            GSymMatrix m_test12_inv = m_test;
            m_test12_inv.cholesky_invert();
        GSymMatrix ci_product   = m_test * m_test12_inv;
        GSymMatrix ci_residuals = ci_product - unit;
            //
            res = (abs(ci_residuals)).max();
            if (res < 1.0e-15){
            //cout << "Res(CI)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt GSymMatrix.cholesky_invert(GVector) function. 10");
            }
            //
            // Test Cholesky inverter for compressed matrix
            unit = GSymMatrix(4,4);
            unit(0,0) = unit(1,1) = unit(3,3) = 1.0;
            GSymMatrix m_test12_zero_inv = m_test_zero;
            m_test12_zero_inv.cholesky_invert();
            GSymMatrix ciz_product   = m_test_zero * m_test12_zero_inv;
            GSymMatrix ciz_residuals = ciz_product - unit;
            //
            res = (abs(ciz_residuals)).max();
            if (res < 1.0e-15){
            //cout << "Res(CIZ)=" << res << " ";
            }
            else {
                throw testsuite.exception_failure("Corrupt compressed GSymMatrix.cholesky_invert() function. 11");
            }
            testsuite.test_try_success();
    }
    catch (exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.end_test();

        //save xml report
    testsuites.save("reports/GSymMatrix.xml");
    // Return
    return 0;
}
