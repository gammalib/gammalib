/***************************************************************************
 *                 test_GMatrix.cpp  -  test GMatrix class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
GMatrix set_matrix(void)
{
  try {
    GMatrix matrix(g_rows,g_cols);
    for (int row = 0; row < g_rows; ++row) {
      for (int col = 0; col < g_cols; ++col)
        matrix(row,col) = g_matrix[col+row*g_cols];
    }
    return matrix;
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to set test matrix." << endl;
    throw;
  }
}


/***************************************************************************
 *                              Set test vector                            *
 ***************************************************************************/
GVector set_vector(void)
{
  try {
    GVector vector(g_cols);
	for (int col = 0; col < g_cols; ++col)
	  vector[col] = g_vector[col];
    return vector;
  }
  catch (exception &e) {
    cout << "TEST ERROR: Unable to set test vector." << endl;
	throw;
  }
}


/***************************************************************************
 *                             Check test matrix                           *
 ***************************************************************************/
int check_matrix(const GMatrix& m, const double scale, const double add)
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
 *                        Check transposed test matrix                     *
 ***************************************************************************/
int check_transpose_matrix(const GMatrix& m, const double scale, const double add)
{
  int result = 1;
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
int check_matrix_matrix(const GMatrix& m)
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
void test_output(const GMatrix& m_test)
{
  cout << "Test GMatrix: Output test matrix: " << endl;
  try {
    cout << m_test << endl;
  }
  catch (exception &e) {
	cout << endl << "TEST ERROR: Unable to output matrix." << endl;
    cout << e.what() << endl;
	throw;
  }
}


/***************************************************************************
 *                     Test: Conversion between matrix types               *
 ***************************************************************************/
void test_conversion(void)
{
  cout << "Test GMatrix: Matrix conversions: ";
  try {
    //
    // Setup a symmetric matrix
    int     num = 10;
    GMatrix symmetric(num,num);
    for (int i = 0; i < num; ++i) {
      for (int j = 0; j < num; ++j)
        symmetric(i,j) = (i+j+1)/2.0;
    }
    cout << ".";
    //
    // Convert symmetric matrix into GSymMatrix object
//    GSymMatrix converted = sym_matrix(symmetric);
    GSymMatrix converted = GSymMatrix(symmetric);
    cout << ".";
    //
    // Convert GSymMatrix back to full matrix
    GMatrix back_convert = matrix(converted);
    cout << ".";
    //
    // Compare back converted matrix to original one. They should be identical
    if (symmetric != back_convert) {
      cout << endl << "TEST ERROR: Unable to convert matrixes (symmetric)." << endl;
      cout << "Original matrix " << symmetric << endl;
      cout << "GSymMatrix matrix " << converted << endl;
      cout << "Back converted matrix " << back_convert << endl;
      throw;
    }
    //
    // Determine the fill of the matrix. It should be 1.0
    double fill = back_convert.fill();
    if (abs(fill-1.0) > 1.0e-15) {
      cout << endl << "TEST ERROR: Bad fill " << fill << " determined (expected 1.0)." <<
           endl;
      throw;
    }
    cout << ".";
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
    if (!ok) {
      cout << endl << "TEST ERROR: Corrupt extract_lower_triangle." << endl;
      cout << "Original matrix " << symmetric << endl;
      cout << "Lower triangle matrix " << lower << endl;
      throw;
    }
    cout << ".";
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
    if (!ok) {
      cout << endl << "TEST ERROR: Corrupt extract_upper_triangle." << endl;
      cout << "Original matrix " << symmetric << endl;
      cout << "Upper triangle matrix " << upper << endl;
      throw;
    }
    cout << ".";
    //
    // Now make the matrix unsymmetric
    symmetric(0,num-1) = 1000.0;
    //
    // Try converting now into GSymMatrix object (this should fail)
    try {
      converted = sym_matrix(symmetric);
      cout << endl << "TEST ERROR: Conversion to symmetric should have failed." << endl;
      throw;
    }
    catch (GException::matrix_not_symmetric &e) {
      cout << ".";
    }
    catch (exception &e) {
      cout << e.what() << endl;
      throw;
    }
    //
    // Now zero some elements to emulate a sparse matrix
    symmetric(0,num-1) = 0.0;
    for (int i = 5; i < 8; ++i) {
      for (int j = i-2; j < i+2; ++j)
        symmetric(i,j) = 0.0;
    }
    cout << ".";
    //
    // Convert symmetric matrix into GSparseMatrix object
    GSparseMatrix sparse = sparse_matrix(symmetric);
    cout << ".";
    //
    // Convert GSparseMatrix back to full matrix
    back_convert = matrix(sparse);
    cout << ".";
    //
    // Compare back converted matrix to original one. They should be identical
    if (symmetric != back_convert) {
      cout << endl << "TEST ERROR: Unable to convert matrixes (sparse)." << endl;
      cout << "Original matrix " << symmetric << endl;
      cout << "GSparseMatrix matrix " << sparse << endl;
      cout << "Back converted matrix " << back_convert << endl;
      throw;
    }
    cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
    throw;
  }
  cout << " ok." << endl;

  // Return
  return;
}


/***************************************************************************
 *                        Test: extraction and insertion                   *
 ***************************************************************************/
void test_extract(void)
{
  cout << "Test GMatrix: Vector extraction, insertion and addition: ";
  try {
    //
	// Set-up test matrix
    int rows = 10;
    int cols = 20;
    GMatrix test(rows, cols);
	//
	// Add and extract column vectors
	for (int col = 0; col < cols; ++col) {
	  GVector column(rows);
	  for (int row = 0; row < rows; ++row)
	    column[row] = (col+1)*100.0 + (row+1)*1.0;
	  test.add_col(column, col);
	  GVector check = test.extract_col(col);
	  if (check != column) {
	    cout << endl << "TEST ERROR: Unable to add and extract columns." << endl;
	    cout << "Added column ...: " << column << endl;
	    cout << "Extracted column: " << check << endl;
        throw;
	  }
	}
	cout << ".";
	//
	// Insert and extract column vectors
	for (int col = 0; col < cols; ++col) {
	  GVector column(rows);
	  for (int row = 0; row < rows; ++row)
	    column[row] = (col+1)*100.0 + (row+1)*1.0;
	  test.insert_col(column, col);
	  GVector check = test.extract_col(col);
	  if (check != column) {
	    cout << endl << "TEST ERROR: Unable to insert and extract columns." << endl;
	    cout << "Inserted column : " << column << endl;
	    cout << "Extracted column: " << check << endl;
        throw;
	  }
	}
	cout << ".";
	//
	// Extract rows
    for (int row = 0; row < rows; ++row) {
	  GVector v_row(cols);
	  for (int col = 0; col < cols; ++col)
	    v_row[col] = (col+1)*100.0 + (row+1)*1.0;
	  GVector check = test.extract_row(row);
	  if (check != v_row) {
	    cout << endl << "TEST ERROR: Unable to extract rows." << endl;
	    cout << "Inserted row : " << v_row << endl;
	    cout << "Extracted row: " << check << endl;
        throw;
	  }
	}
	cout << ".";
  }
  catch (exception &e) {
	cout << endl << "TEST ERROR: Corrupt vector extraction, insertion, or addition." << endl;
    cout << e.what() << endl;
	throw;
  }
  cout << ". ok." << endl;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
  // Dump header
  cout << endl;
  cout << "*************************" << endl;
  cout << "* GMatrix class testing *" << endl;
  cout << "*************************" << endl;

  // Set test matrix and vector
  GMatrix m_test = set_matrix();
  GVector v_test = set_vector();

  // Set bigger matrix (only used for collision, don't care about content)
  GMatrix bigger(g_rows+1,g_cols+1);

  // Prepare result matrix
  GMatrix result = m_test;

  // Execute the tests
  test_output(m_test);
  test_conversion();
  test_extract();

  // Test 1: Allocate zero matrix
  try {
    cout << "GMatrix - Test 1: Allocate zero matrix: ";
    GMatrix test1(0,0);
  }
  catch (GException::empty &e) {
    cout << "ok." << endl;
  }
  catch (exception &e) {
	cout << endl << "TEST ERROR: Did not signal zero matrix allocation." << endl;
    cout << e.what() << endl;
	throw;
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
  try {
    cout << "GMatrix - Test 3: Assign matrix values: ";
    //
    GMatrix test3(3,3);
	test3(1,1) = 1.0;
    if (test3(0,0) != 0.0 || test3(1,0) != 0.0 || test3(2,0) != 0.0 ||
        test3(0,1) != 0.0 || test3(1,1) != 1.0 || test3(2,1) != 0.0 ||
        test3(0,2) != 0.0 || test3(1,2) != 0.0 || test3(2,2) != 0.0) {
      cout << endl << "TEST ERROR: Unable to set matrix value." << endl;
	  cout << test3 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
    #if defined(G_RANGE_CHECK)
    try {
        GMatrix test3(3,3);
        test3(3,3) = 1.0;
    }
    catch (GException::out_of_range &e) {
    }
    catch (exception &e) {
        std::cout << e.what() << std::endl;
        throw;
    }
    #endif
    cout << "ok." << endl;

  // Test 4: Matrix copy constructor
  try {
    cout << "GMatrix - Test 4: Define matrix using copy constructor: ";
	//
	GMatrix test4 = m_test;
    if (!check_matrix(test4, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt copy constructor." << endl;
	  cout << test4 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 5: Matrix assignment
  try {
    cout << "GMatrix - Test 5: Matrix assignment: ";
	//
	// GMatrix = GMatrix
	result = m_test;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt assignment." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix = GMatrix (bigger matrix)
	result = bigger;
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 6: Transposition
  try {
    cout << "GMatrix - Test 6: Transposition: ";
    //
	// transpose(GMatrix)
	result = transpose(m_test);
    if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt transpose(GMatrix) function." << endl;
	  cout << result << endl;
	  throw;
	}
    //
	// GMatrix.transpose()
	result = m_test;
	result.transpose();
    if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix.transpose() function." << endl;
	  cout << result << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 7: Matrix*Vector multiplication
  try {
    cout << "GMatrix - Test 7: Matrix*Vector multiplication: ";
    //
	GVector v_test7 = m_test*v_test;
	if (!check_matrix_vector(v_test7) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Vector multiplication." << endl;
	  cout << v_test7 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GVector v_test7 = bigger*v_test;
  }
  catch (GException::matrix_vector_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 8: Matrix*Matrix multiplication
  try {
    cout << "GMatrix - Test 8: Matrix*Matrix multiplication: ";
    //
	GMatrix m_test8 = m_test * transpose(m_test);
	if (!check_matrix_matrix(m_test8) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Matrix multiplication." << endl;
	  cout << m_test8 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GMatrix m_test8 = m_test*bigger;
  }
  catch (GException::matrix_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GMatrix m_test8 = bigger*m_test;
  }
  catch (GException::matrix_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 9: Assignment and arithmetics
  try {
    cout << "GMatrix - Test 9: Assignment and arithmetics: ";
	//
	// -GMatrix
	result = -m_test;
    if (!check_matrix(result, -1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt -GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix += GMatrix
	result  = m_test;
	result += m_test;
    if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix += GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix -= GMatrix
	result  = m_test;
	result -= m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix -= GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix *= 3.0
	result  = m_test;
	result *= 3.0;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix *= double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix /= 3.0
	result  = m_test;
	result /= 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix /= double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix + GMatrix
	result = m_test + m_test;
    if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix + GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix - GMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix - GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix - GMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix - GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix * 3.0
	result = m_test * 3.0;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix * double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// 3.0 * GMatrix
	result = 3.0 * m_test;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt double * GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix / 3.0
	result = m_test / 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix / double operator." << endl;
	  cout << result << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	result  = m_test;
	result += bigger;
  }
  catch (GException::matrix_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 10: Matrix functions
  try {
    cout << "GMatrix - Test 10: Matrix functions: ";
	//
	double min = m_test.min();
    if (!check_matrix_min(min)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix.min() function." << endl;
	  cout << min << endl;
	  throw;
	}
	//
	double max = m_test.max();
    if (!check_matrix_max(max)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix.max() function." << endl;
	  cout << max << endl;
	  throw;
	}
	//
	double sum = m_test.sum();
    if (!check_matrix_sum(sum)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix.sum() function." << endl;
	  cout << sum << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  cout << "ok." << endl;

  // Test 11: Matrix comparison
  try {
    cout << "GMatrix - Test 11: Comparison: ";
	//
	if ((m_test == m_test) != 1) {
      cout << endl << "TEST ERROR: Corrupt GMatrix == GMatrix operator." << endl;
	  throw;
	}
	//
	GMatrix m_test10(g_rows,g_cols);
	if ((m_test == m_test10) != 0) {
      cout << endl << "TEST ERROR: Corrupt GMatrix == GMatrix operator." << endl;
	  throw;
	}
	//
	if ((m_test == bigger) != 0) {
      cout << endl << "TEST ERROR: Corrupt GMatrix == GMatrix operator." << endl;
	  throw;
	}
	//
	if ((m_test != m_test) != 0) {
      cout << endl << "TEST ERROR: Corrupt GMatrix != GMatrix operator." << endl;
	  throw;
	}
    //
	if ((m_test != m_test10) != 1) {
      cout << endl << "TEST ERROR: Corrupt GMatrix != GMatrix operator." << endl;
	  throw;
	}
	//
	if ((m_test != bigger) != 1) {
      cout << endl << "TEST ERROR: Corrupt GMatrix != GMatrix operator." << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Return
  return 0;
}
