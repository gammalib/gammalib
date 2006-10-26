/***************************************************************************
 *                 test_GMatrix.cpp  -  test GMatrix class                 *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "test_GMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Globals ____________________________________________________________ */
double   g_matrix[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double   g_vector[] = {1.0, 2.0, 3.0};
unsigned g_rows    = 3;
unsigned g_cols    = 3;


/***************************************************************************
 *                              Set test matrix                            *
 ***************************************************************************/
GMatrix set_matrix(void)
{
  try {
    GMatrix matrix(g_rows,g_cols);
    for (unsigned row = 0; row < g_rows; ++row) {
      for (unsigned col = 0; col < g_cols; ++col)
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
	for (unsigned col = 0; col < g_cols; ++col)
	  vector(col) = g_vector[col];
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
    for (unsigned row = 0; row < g_rows; ++row) {
      for (unsigned col = 0; col < g_cols; ++col) {
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
    for (unsigned row = 0; row < g_rows; ++row) {
      for (unsigned col = 0; col < g_cols; ++col) {
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
    for (unsigned row = 0; row < g_rows; ++row) {
      for (unsigned col = 0; col < g_cols; ++col) {
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
    for (unsigned row = 0; row < v.size(); ++row) {
      double value = 0.0;
	  for (unsigned col = 0; col < g_cols; ++col)
	    value += g_matrix[col+row*g_cols] * g_vector[col];
	  if (v(row) != value) {
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
  if (m.rows() != g_rows || m.cols() != g_cols)
    return 0;
  int result = 1;
  try {
    for (unsigned row = 0; row < g_rows; ++row) {
	  for (unsigned col = 0; col < g_cols; ++col) {
        double value = 0.0;
		for (unsigned i = 0; i < g_cols; ++i)
	      value += g_matrix[i+row*g_cols] * g_matrix[row+i*g_cols];
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
  for (unsigned row = 0; row < g_rows; ++row) {
    for (unsigned col = 0; col < g_cols; ++col) {
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
  for (unsigned row = 0; row < g_rows; ++row) {
    for (unsigned col = 0; col < g_cols; ++col) {
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
  for (unsigned row = 0; row < g_rows; ++row) {
    for (unsigned col = 0; col < g_cols; ++col)
      value += g_matrix[col+row*g_cols];
  }
  return (sum == value);
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

  // Test 1: Allocate zero matrix
  try {
    cout << "GMatrix - Test 1: Allocate zero matrix: ";
    GMatrix test1(0,0);
  }
  catch (empty &e) {
    cout << "ok." << endl;
  }
  catch (exception &e) {
	cout << endl << "TEST ERROR: Did not signal zero matrix allocation." << endl;
    cout << e.what() << endl;
	throw;
  }

  // Test 2: Allocate too large matrix
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
  #if G_RANGE_CHECK
  try {
    GMatrix test3(3,3);
    test3(3,3) = 1.0;
  }
  catch (GMatrix::out_of_range &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  #endif  
  cout << "ok." << endl;

  // Test 4: Matrix copy constructor
  try {
    cout << "GMatrix - Test 4: Define matrix using copy constructor: ";
	//
	GMatrix test4 = m_test;
    if (!check_matrix(test4, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
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
    if (!check_matrix(result, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt copy constructor." << endl;
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

  // Test 6: Matrix*Vector multiplication
  try {
    cout << "GMatrix - Test 6: Matrix*Vector multiplication: ";
    //
	GVector v_test6 = m_test*v_test;
	if (!check_matrix_vector(v_test6) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Vector multiplication." << endl;
	  cout << v_test6 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GVector v_test6 = bigger*v_test;
  }
  catch (GMatrix::vec_mat_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 7: Matrix*Matrix multiplication
  try {
    cout << "GMatrix - Test 7: Matrix*Matrix multiplication: ";
    //
	GMatrix m_test7 = m_test*m_test;
	if (!check_matrix_matrix(m_test7) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Matrix multiplication." << endl;
	  cout << m_test7 << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GMatrix m_test7 = m_test*bigger;
  }
  catch (GMatrix::dim_mult_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GMatrix m_test7 = bigger*m_test;
  }
  catch (GMatrix::dim_mult_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 8: Assignment and arithmetics
  try {
    cout << "GMatrix - Test 8: Assignment and arithmetics: ";
	//
	// -GMatrix
	result = -m_test;
    if (!check_matrix(result, -1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt -GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix = 3.0
	result = m_test;
	result = 3.0;
    if (!check_matrix(result, 0.0, 3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix = double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix += GMatrix
	result  = m_test;
	result += m_test;
    if (!check_matrix(result, 2.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix += GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix += 3.0
	result  = m_test;
	result += 3.0;
    if (!check_matrix(result, 1.0, 3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix += double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix -= GMatrix
	result  = m_test;
	result -= m_test;
    if (!check_matrix(result, 0.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix -= GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix -= 3.0
	result  = m_test;
	result -= 3.0;
    if (!check_matrix(result, 1.0, -3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix -= double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix *= 3.0
	result  = m_test;
	result *= 3.0;
    if (!check_matrix(result, 3.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix *= double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix /= 3.0
	result  = m_test;
	result /= 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix /= double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix + GMatrix
	result = m_test + m_test;
    if (!check_matrix(result, 2.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix + GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix - GMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix - GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix - GMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix - GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix + 3.0
	result = m_test + 3.0;
    if (!check_matrix(result, 1.0, 3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix + double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// 3.0 + GMatrix
	result = 3.0 + m_test;
    if (!check_matrix(result, 1.0, 3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt double + GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix - 3.0
	result = m_test - 3.0;
    if (!check_matrix(result, 1.0, -3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix - double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// 3.0 - GMatrix
	result = 3.0 - m_test;
    if (!check_matrix(result, -1.0, 3.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt double - GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix * 3.0
	result = m_test * 3.0;
    if (!check_matrix(result, 3.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix * double operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// 3.0 * GMatrix
	result = 3.0 * m_test;
    if (!check_matrix(result, 3.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt double * GMatrix operator." << endl;
	  cout << result << endl;
	  throw;
	}
	//
	// GMatrix / 3.0
	result = m_test / 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
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
  catch (GMatrix::dim_add_mismatch) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

  // Test 9: Matrix functions
  try {
    cout << "GMatrix - Test 9: Matrix functions: ";
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
	  cout << min << endl;
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

  // Test 10: Matrix comparison
  try {
    cout << "GMatrix - Test 10: Comparison: ";
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

  // Test 11: Transformations
  try {
    cout << "GMatrix - Test 11: Transformations: ";
    //
	// transpose(GMatrix)
	result = transpose(m_test);
    if (!check_matrix(result, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt transpose(GMatrix) function." << endl;
	  cout << result << endl;
	  throw;
	}
    //
	// GMatrix.transpose()
	result = m_test;
	result.transpose();
    if (!check_matrix(result, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GMatrix.transpose() function." << endl;
	  cout << result << endl;
	  throw;
	}
/*
	//
	// lower_triangle(GMatrix)
	m_test11 = lower_triangle(m_test);
    if (!check_matrix_lt(m_test11, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt lower_triangle(GMatrix) function." << endl;
	  cout << m_test11 << endl;
	  throw;
	}
	//
	// upper_triangle(GMatrix)
	m_test11 = upper_triangle(m_test);
    if (!check_matrix_ut(m_test11, 1.0, 0.0) && !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt upper_triangle(GMatrix) function." << endl;
	  cout << m_test11 << endl;
	  throw;
	}
*/
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;

}