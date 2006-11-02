/***************************************************************************
 *            test_GSparseMatrix.cpp  -  test GSparseMatrix class          *
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
#include "test_GSparseMatrix.hpp"

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
GSparseMatrix set_matrix(void)
{
  try {
    GSparseMatrix matrix(g_rows,g_cols);
    for (int i = 0; i < g_elements; ++i)
	  matrix(g_row[i],g_col[i]) = g_matrix[i];
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
	    if (fabs(m(row,col)-ref_value) > 1.0e-15) {
//cout << row << "," << col << ": " << fabs(m(row,col)-ref_value) << endl;
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
//cout << row << "," << col << ": " << fabs(m(row,col)-value) << endl;
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
	    if (fabs(m(col,row)-ref_value) > 1.0e-15) {
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
void test_output(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Output test matrix: " << endl;
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
 *                          Test: Allocate zero matrix                     *
 ***************************************************************************/
void test_allocate(void)
{
  cout << "Test GSparseMatrix: Allocate zero matrix: ";
  try {
    GSparseMatrix test(0,0);
  }
  catch (empty &e) {
    cout << ". ok." << endl;
  }
  catch (exception &e) {
	cout << endl << "TEST ERROR: Did not signal zero matrix allocation." << endl;
    cout << e.what() << endl;
	throw;
  }
  cout << ". ok." << endl;
}


/***************************************************************************
 *                              Test: Assign values                        *
 ***************************************************************************/
void test_assign_values(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Assign matrix values: ";
  try {
    //
    GSparseMatrix test(3,3);
	test(1,1) = 1.0; 
    if (test(0,0) != 0.0 || test(1,0) != 0.0 || test(2,0) != 0.0 ||
        test(0,1) != 0.0 || test(1,1) != 1.0 || test(2,1) != 0.0 ||
        test(0,2) != 0.0 || test(1,2) != 0.0 || test(2,2) != 0.0) {
      cout << endl << "TEST ERROR: Unable to set matrix values." << endl;
	  cout << test << endl;
	  throw;
	}
    cout << ".";
	//
	GSparseMatrix result = m_test;
	result(0,0) = 0.0;
	result(0,0) = 1.0;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Bad M(row,col) assignment." << endl;
	  cout << "result:" << endl << result << endl;
	  cout << "m_test:" << endl << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	result      = m_test;
	result(0,0) = 0.0;
	result(0,4) = 0.0;
	result(0,0) = 1.0;
	result(0,4) = 0.0;
	result(0,4) = 7.0;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Bad M(row,col) assignment." << endl;
	  cout << "result:" << endl << result << endl;
	  cout << "m_test:" << endl << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	result      = m_test;
	result(0,1) = 2.0;
	result(0,2) = 2.0;
	result(0,3) = 2.0;
	result(0,1) = 0.0;
	result(0,2) = 0.0;
	result(0,3) = 0.0;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Bad M(row,col) assignment." << endl;
	  cout << "result:" << endl << result << endl;
	  cout << "m_test:" << endl << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	result.clear();
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.clear()." << endl;
	  cout << "result:" << endl << result << endl;
	  cout << "m_test:" << endl << m_test << endl;
	  throw;
	}
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  #if defined(G_RANGE_CHECK)
  try {
    GSparseMatrix test(3,3);
    test(3,3) = 1.0;
  }
  catch (GMatrix::out_of_range &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << ".";
  #endif
  //
  // Insert column into 10 x 10 matrix using matrix stack
  try {
	GSparseMatrix stack(10,10);
	stack.clear();
	for (int i = 3; i < 5; ++i)
	  stack(i,i) = 5.0;
	GVector column(10);
	stack.stack_init(100,50);
	for (int j = 0; j < 10; ++j) {
	  int i_min = (j < 2) ?  0 : j-2;
	  int i_max = (j > 8) ? 10 : j+2;
	  column = 0.0;
	  for (int i = i_min; i < i_max; ++i)
	    column(i) = (i+1)*1;
	  stack.add_col(column, j);
    }
	stack.stack_destroy();
//	cout << stack << endl;
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                        Test: Matrix copy constructor                    *
 ***************************************************************************/
void test_copy_constructor(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Define matrix using copy constructor: ";
  try {
	//
	GSparseMatrix test = m_test;
    if (!check_matrix(test, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt copy constructor." << endl;
	  cout << test << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << ". ok." << endl;
}


/***************************************************************************
 *                           Test: Matrix assignment                       *
 ***************************************************************************/
void test_assign_matrix(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Matrix assignment: ";
  try {
	//
	// GSparseMatrix = GSparseMatrix
	GSparseMatrix result(1,1);
	result = m_test;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt assignment." << endl;
	  cout << result << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix = GSparseMatrix (bigger matrix)
	GSparseMatrix bigger(1000,1000);
	for (int i = 0; i < 1000; ++i)
	  bigger(i,i) = (i+1)*1.1;
	result = bigger;
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                           Test: Matrix transpose                        *
 ***************************************************************************/
void test_transpose(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Matrix transpose: ";
  try {
    //
	// transpose(GSparseMatrix)
	GSparseMatrix result(1,1);
	result = transpose(m_test);
    if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt transpose(GSparseMatrix) function." << endl;
	  cout << result << endl;
	  throw;
	}
	cout << ".";
    //
	// GSparseMatrix.transpose()
	result = m_test;
	result.transpose();
    if (!check_transpose_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.transpose() function." << endl;
	  cout << result << endl;
	  throw;
	}
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                     Test: Matrix*Vector multiplication                  *
 ***************************************************************************/
void test_matrix_vector(const GSparseMatrix& m_test, const GVector& v_test)
{
  cout << "Test GSparseMatrix: Matrix*Vector multiplication: ";
  try {
    //
	GVector test = m_test*v_test;
	if (!check_matrix_vector(test) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Vector multiplication." << endl;
	  cout << test << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << ".";
  try {
	GSparseMatrix bigger(1000,1000);
	for (int i = 0; i < 1000; ++i)
	  bigger(i,i) = (i+1)*1.1;
	GVector test = bigger*v_test;
  }
  catch (GMatrix::matrix_vector_mismatch &e) {
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << ". ok." << endl;
}


/***************************************************************************
 *                     Test: Matrix*Vector multiplication                  *
 ***************************************************************************/
void test_matrix_matrix(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Matrix*Matrix multiplication: ";
  try {
    //
	GSparseMatrix test = m_test * transpose(m_test);
	if (!check_matrix_matrix(test) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt Matrix*Matrix multiplication." << endl;
	  cout << "test:"   << endl << test << endl;
	  cout << "m_test:" << endl << m_test << endl;
	  throw;
	}
    cout << ".";
    //
	GSparseMatrix bigger(1000,1000);
	for (int i = 0; i < 1000; ++i)
	  bigger(i,i) = (i+1)*1.1;
    try {
	  GSparseMatrix test = m_test * bigger;
    }
    catch (GMatrix::matrix_mismatch &e) {
      cout << ".";
    }
    catch (exception &e) {
      cout << e.what() << endl;
	  throw;
    }
	//
    try {
	  GSparseMatrix test7 = bigger * m_test;
    }
    catch (GMatrix::matrix_mismatch &e) {
      cout << ".";
    }
    catch (exception &e) {
      cout << e.what() << endl;
	  throw;
    }
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                           Test: Arithmetics                             *
 ***************************************************************************/
void test_arithmetics(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Arithmetics: ";
  try {
	//
	// -GSparseMatrix
	GSparseMatrix result(1,1);
	result = -m_test;
    if (!check_matrix(result, -1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt -GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// fabs(GSparseMatrix)
	result = fabs(m_test);
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt fabs(GSparseMatrix) function." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix += GSparseMatrix
	result  = m_test;
	result += m_test;
    if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix += GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix -= GSparseMatrix
	result  = m_test;
	result -= m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix -= GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix *= 3.0
	result  = m_test;
	result *= 3.0;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix *= double operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix /= 3.0
	result  = m_test;
	result /= 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix /= double operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix + GSparseMatrix
	result = m_test + m_test;
    if (!check_matrix(result, 2.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix + GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix - GSparseMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix - GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix - GSparseMatrix
	result = m_test - m_test;
    if (!check_matrix(result, 0.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix - GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix * 3.0
	result = m_test * 3.0;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix * double operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix * -7.12
	result = m_test * -7.12;
    if (!check_matrix(result, -7.12, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix * double operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// 3.0 * GSparseMatrix
	result = 3.0 * m_test;
    if (!check_matrix(result, 3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt double * GSparseMatrix operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix / 3.0
	result = m_test / 3.0;
    if (!check_matrix(result, 1.0/3.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix / double operator." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix / 0.0
	double zero = 0.0;
	result = m_test / zero;
    if (!check_matrix(result, 1.0/zero, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix / 0.0." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	cout << ".";
	}
	//
	// GSparseMatrix.add_col()
	result = m_test;
 	GVector v_add(g_rows);
	for (int i = 0; i < g_rows; ++i)
	  v_add(i) = 7.0;
	for (int col = 0; col < g_cols; ++col)
	  result.add_col(v_add, col);
    if (!check_matrix(result, 1.0, 7.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.add_col(+v)." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix.add_col()
	for (int col = 0; col < g_cols; ++col)
	  result.add_col(-v_add, col);
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.add_col(-v)." << endl;
	  cout << "result " << result << endl;
	  cout << "m_test " << m_test << endl;
	  throw;
	}
	cout << ".";
	//
	// GSparseMatrix.add_col() using a compressed array
	GSparseMatrix compare = m_test;
 	GVector v_sparse(g_rows);
	for (int i = 0; i < g_rows; i += 2)  // Build vector for comparison
	  v_sparse(i) = 7.0;
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
    if (result != compare) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.add_col(compressed array)." << endl;
	  cout << "m_test " << m_test << endl;
	  cout << "result " << result << endl;
	  cout << "compare " << compare << endl;
	  throw;
	}
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  try {
	GSparseMatrix bigger(1000,1000);
	GSparseMatrix result(1,1);
	for (int i = 0; i < 1000; ++i)
	  bigger(i,i) = (i+1)*1.1;
	result  = m_test;
	result += bigger;
  }
  catch (GMatrix::matrix_mismatch) {
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                           Test: Matrix functions                        *
 ***************************************************************************/
void test_functions(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Matrix functions: ";
  try {
	//
	double min = m_test.min();
    if (!check_matrix_min(min)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.min() function." << endl;
	  cout << "min:" << min << endl;
	  throw;
	}
	cout << ".";
	//
	double max = m_test.max();
    if (!check_matrix_max(max)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.max() function." << endl;
	  cout << "max:" << max << endl;
	  throw;
	}
	cout << ".";
	//
	double sum = m_test.sum();
    if (!check_matrix_sum(sum)) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.sum() function." << endl;
	  cout << "sum:" << sum << endl;
	  throw;
	}
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                           Test: Matrix comparison                       *
 ***************************************************************************/
void test_comparison(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Comparison: ";
  try {
	GSparseMatrix bigger(1000,1000);
	for (int i = 0; i < 1000; ++i)
	  bigger(i,i) = (i+1)*1.1;
	//
	if ((m_test == m_test) != 1) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix == GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
	//
	GSparseMatrix test(g_rows,g_cols);
	if ((m_test == test) != 0) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix == GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
	//
	if ((m_test == bigger) != 0) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix == GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
	//
	if ((m_test != m_test) != 0) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix != GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
    //
	if ((m_test != test) != 1) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix != GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
	//
	if ((m_test != bigger) != 1) {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix != GSparseMatrix operator." << endl;
	  throw;
	}
	cout << ".";
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << " ok." << endl;
}


/***************************************************************************
 *                           Test: Matrix conversions                      *
 ***************************************************************************/
void test_conversion(const GSparseMatrix& m_test)
{
  cout << "Test GSparseMatrix: Conversions: ";
  cout << " ok." << endl;
}


/***************************************************************************
 *                        Test: Cholesky factorisation                     *
 ***************************************************************************/
void test_cholesky(void)
{
  cout << "Test GSparseMatrix: Cholesky decomposition, solver and inverter: ";
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
	try {
	  GVector v_test12(5);
	  v_test12 = m_chol_test.cholesky_solver(v_test12);
    }
    catch (GMatrix::matrix_not_factorised) {
	  cout << ".";
    }
    catch (exception &e) {
      cout << e.what() << endl;
	  throw;
    }
	//
    // Test Cholesky decomposition
	GSparseMatrix cd = cholesky_decompose(m_chol_test);
	cout << ".";
	//
    // Test inplace Cholesky decomposition
	cd = m_chol_test;
	cd.cholesky_decompose();
	cout << ".";
	//
    // Test Cholesky solver
	GVector e0(5);
	GVector a0(5);
	e0(0) = 1.0;
	a0(0) = 1.0;
	a0(1) = 0.2;
	a0(2) = 0.2;
	a0(3) = 0.2;
	a0(4) = 0.2;
	GVector s0 = cd.cholesky_solver(a0);
	double res = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << " Res(S0)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(0) = 0.0;
	e0(1) = 1.0;
	a0(0) = 0.2;
	a0(1) = 1.0;
	a0(2) = 0.0;
	a0(3) = 0.0;
	a0(4) = 0.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S1)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(1) = 0.0;
	e0(2) = 1.0;
	a0(1) = 0.0;
	a0(2) = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S2)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(2) = 0.0;
	e0(3) = 1.0;
	a0(2) = 0.0;
	a0(3) = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S3)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(3) = 0.0;
	e0(4) = 1.0;
	a0(3) = 0.0;
	a0(4) = 1.0;
	s0    = cd.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S4)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	// Setup matrix for Cholesky decomposition with zero row/col
	GSparseMatrix m_chol_test_zero(6,6);
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
	//
    // Test compressed Cholesky decomposition
	GSparseMatrix cd_zero = m_chol_test_zero;
	cd_zero.cholesky_decompose();
	//
    // Test compressed Cholesky solver
	e0 = GVector(6);
	a0 = GVector(6);
	e0(0) = 1.0;
	a0(0) = 1.0;
	a0(1) = 0.2;
	a0(2) = 0.2;
	a0(4) = 0.2;
	a0(5) = 0.2;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S0Z)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(0) = 0.0;
	e0(1) = 1.0;
	a0(0) = 0.2;
	a0(1) = 1.0;
	a0(2) = 0.0;
	a0(4) = 0.0;
	a0(5) = 0.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S1Z)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(1) = 0.0;
	e0(2) = 1.0;
	a0(1) = 0.0;
	a0(2) = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S2Z)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(2) = 0.0;
	e0(4) = 1.0;
	a0(2) = 0.0;
	a0(4) = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S3Z)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(4) = 0.0;
	e0(5) = 1.0;
	a0(4) = 0.0;
	a0(5) = 1.0;
	s0    = cd_zero.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S4Z)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}	
	//
	// Setup matrix for Cholesky decomposition with zero row/col (unsymmetric case)
	GSparseMatrix m_chol_test_zero2(6,5);
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
	//
    // Test compressed Cholesky decomposition (unsymmetric case)
	GSparseMatrix cd_zero2 = m_chol_test_zero2;
	cd_zero2.cholesky_decompose();
	//
    // Test compressed Cholesky solver (unsymmetric case)
	e0 = GVector(5);
	a0 = GVector(6);
	e0(0) = 1.0;
	a0(0) = 1.0;
	a0(1) = 0.2;
	a0(2) = 0.2;
	a0(4) = 0.2;
	a0(5) = 0.2;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S0ZU)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(0) = 0.0;
	e0(1) = 1.0;
	a0(0) = 0.2;
	a0(1) = 1.0;
	a0(2) = 0.0;
	a0(4) = 0.0;
	a0(5) = 0.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S1ZU)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(1) = 0.0;
	e0(2) = 1.0;
	a0(1) = 0.0;
	a0(2) = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S2ZU)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(2) = 0.0;
	e0(3) = 1.0;
	a0(2) = 0.0;
	a0(4) = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S3ZU)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	e0(3) = 0.0;
	e0(4) = 1.0;
	a0(4) = 0.0;
	a0(5) = 1.0;
	s0    = cd_zero2.cholesky_solver(a0);
	res   = max(fabs(s0-e0));
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(S4ZU)=" << res << " "; 
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_solver(GVector) function." << endl;
	  cout << "Residual vector (all elements should be zero):" << endl;
	  cout << s0 << endl;
	  throw;
	}
	//
	// Test Cholesky inverter (inplace)
	GSparseMatrix unit(5,5);
	unit(0,0) = 1.0;
	unit(1,1) = 1.0;
	unit(2,2) = 1.0;
	unit(3,3) = 1.0;
	unit(4,4) = 1.0;
	GSparseMatrix m_chol_test_inv = m_chol_test;
	m_chol_test_inv.cholesky_invert();
    GSparseMatrix ci_product   = m_chol_test * m_chol_test_inv;
    GSparseMatrix ci_residuals = ci_product - unit;
	res = (fabs(ci_residuals)).max();
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(CI)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt GSparseMatrix.cholesky_invert() function." << endl;
	  cout << "Residual matrix (all elements should be zero):" << endl;
	  cout << ci_residuals << endl;
	  throw;
	}
	//
	// Test Cholesky inverter
	m_chol_test_inv = cholesky_invert(m_chol_test);
    ci_product      = m_chol_test * m_chol_test_inv;
    ci_residuals    = ci_product - unit;
	res             = (fabs(ci_residuals)).max();
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(CI2)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt cholesky_invert(GSparseMatrix) function." << endl;
	  cout << "Residual matrix (all elements should be zero):" << endl;
	  cout << ci_residuals << endl;
	  throw;
	}
	//
	// Test Cholesky inverter for compressed matrix
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
	res = (fabs(ciz_residuals)).max();
	if (res < 1.0e-15)
	  #if defined(DUMP_RESIDUALS)
	  cout << "Res(CIZ)=" << res << " ";
	  #else
	  cout << ".";
	  #endif
	else {
      cout << endl << "TEST ERROR: Corrupt compressed GSparseMatrix.cholesky_invert() function." << endl;
	  cout << "Residual matrix (all elements should be zero):" << endl;
	  cout << ciz_residuals << endl;
	  throw;
	}
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  cout << "ok." << endl;
}


/***************************************************************************
 *                            Test: Heavy calculus                         *
 ***************************************************************************/
void test_heavy(int block = 512)
{
  // Allocate some variables
  clock_t t_start;
  double  t_elapse;
  int     number;
  int     num_cols;
  int     i_min, i_max;
  
  // Perform test
  cout << "Test GSparseMatrix: Heavy calculus: ";
  try {
    //
	// Fill 10000 x 10000 matrix with 10000 values
    t_start = clock();
	number  = 10000;
	GSparseMatrix large(number,number);
	large.set_mem_block(block);
	for (int i = 0; i < number; ++i)
	  large(i,i) = (i+1)*3.14;
	t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    #if defined(DUMP_TIMING)
	cout << endl << " - Fill of 10000 values needed " << t_elapse << 
	                " sec (reference ~ 0.2 sec)";
	#endif   
    //
	// Modify 10000 x 10000 matrix with 10000 (exisiting) values
    t_start = clock();
	for (int i = 0; i < number; ++i)
	  large(i,i) = (i+1)*1.0;
	t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    #if defined(DUMP_TIMING)
	cout << endl << " - Modification of 10000 values needed " << t_elapse << 
	                " sec (reference ~ 0 sec)";
	#endif   
    //
	// Insert columns into 10000 x 10000 matrix
	number   = 10000;
	num_cols = 10000;
	GSparseMatrix diag(number,number);
	for (int j = 10; j < number-10; ++j)
	  diag(j,j) = 3.14;
	large = diag;
	large.set_mem_block(block);
	GVector column(number);
    t_start = clock();
	for (int j = 0; j < num_cols; ++j) {
	  i_min = (j < 2) ? 0 : j-2;
	  i_max = (j > num_cols-2) ? num_cols : j+2;
	  column = 0.0;
	  for (int i = i_min; i < i_max; ++i)
	    column(i) = (i+1)*0.01;
	  large.add_col(column, j);
    }
	t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    #if defined(DUMP_TIMING)
	cout << endl << " - " << num_cols << " columns adding needed " << t_elapse << 
	                " sec (reference ~ 4.8 sec)";
	#endif   
    //
	// Insert columns into 10000 x 10000 matrix using matrix stack
	GSparseMatrix stack(number,number);
	stack = diag;
	stack.set_mem_block(block);
	for (int j = 200; j < 400; ++j)
	  stack(j,j) = 3.14;
    t_start = clock();
	column = GVector(number);
	stack.stack_init(10000,10000);
	for (int j = 0; j < num_cols; ++j) {
	  i_min = (j < 2)          ? 0 : j-2;
	  i_max = (j > num_cols-2) ? num_cols : j+2;
	  column = 0.0;
	  for (int i = i_min; i < i_max; ++i)
	    column(i) = (i+1)*0.01;
	  stack.add_col(column, j);
    }
	stack.stack_destroy();
	t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);
    #if defined(DUMP_TIMING)
	cout << endl << " - " << num_cols << " columns stack-adding needed " << t_elapse << 
	                " sec (reference ~ 1.3 sec)";
	#endif
	//
	// Check that both are identical
	if (large != stack) {
      cout << endl << "TEST ERROR: Stack-based insertion corrupted." << endl;
	  throw;
	}
    //
	// Subsequent addition and subtractions
	number  = 300;
    t_start = clock();
	GSparseMatrix sum(number,number);
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
	                " sec (reference ~ 1.5 sec)";
	#endif   
    //
	// Subsequent multiplications
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
  }
  catch (exception &e) {
    cout << e.what() << endl;
	throw;
  }
  #if defined(DUMP_TIMING)
  cout << endl;
  #endif
  cout << "ok." << endl;
}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
  // Dump header
  cout << endl;
  cout << "*******************************" << endl;
  cout << "* GSparseMatrix class testing *" << endl;
  cout << "*******************************" << endl;
  
  // Set test matrix and vector
  GSparseMatrix m_test = set_matrix();
  GVector       v_test = set_vector();

  // Execute the tests
  test_output(m_test);
  test_allocate();
  test_assign_values(m_test);
  test_copy_constructor(m_test);
  test_assign_matrix(m_test);
  test_transpose(m_test);
  test_matrix_vector(m_test, v_test);
  test_matrix_matrix(m_test);
  test_arithmetics(m_test);
  test_functions(m_test);
  test_comparison(m_test);
  test_conversion(m_test);
  test_cholesky();
  test_heavy();

  // Return
  return 0;
}
