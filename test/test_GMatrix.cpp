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
#include <iostream>                           // std::cout, cerr
#include <stdexcept>                          // std::exception
#include "GammaLib.hpp"
#include "test_GMatrix.hpp"
/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

void TestGMatrix::set(void){
    // Test name
    name("GMatrix");

    //Initialize m_test
    init_matrix();

    //Initialize m_vector
    init_vector();

    //Set parameters
    m_g_rows = 3;
    m_g_cols = 4;

    //add tests
    //TODO: slipt tests

    add_test(static_cast<pfunction>(&TestGMatrix::set_matrix),"Set matrix");
    add_test(static_cast<pfunction>(&TestGMatrix::set_vector),"Set Vector");
    add_test(static_cast<pfunction>(&TestGMatrix::set_bigger),"Set matrix bigger");
    add_test(static_cast<pfunction>(&TestGMatrix::test_output),"Output test matrix");
    add_test(static_cast<pfunction>(&TestGMatrix::test_conversion),"Matrix conversions");
    add_test(static_cast<pfunction>(&TestGMatrix::test_extract),"Vector extraction, insertion and addition");

    add_test(static_cast<pfunction>(&TestGMatrix::test1),"Test 1: Allocate zero matrix");
    //add_test(static_cast<pfunction>(&TestGMatrix::test2),"Test 2:  Allocate too large matrix");
    add_test(static_cast<pfunction>(&TestGMatrix::test3),"Test3: Assign values");
    add_test(static_cast<pfunction>(&TestGMatrix::test4),"Test 4: Matrix copy constructor");
    add_test(static_cast<pfunction>(&TestGMatrix::test5),"Test 5: Matrix assignment");
    add_test(static_cast<pfunction>(&TestGMatrix::test6),"Test 6: Transposition");
    add_test(static_cast<pfunction>(&TestGMatrix::test7),"Test 7: Matrix*Vector multiplication");
    add_test(static_cast<pfunction>(&TestGMatrix::test8),"Test 8: Matrix*Matrix multiplication");
    add_test(static_cast<pfunction>(&TestGMatrix::test9),"Test 9: Assignment and arithmetics");
    add_test(static_cast<pfunction>(&TestGMatrix::test10),"Test 10: Matrix functions");
    add_test(static_cast<pfunction>(&TestGMatrix::test11),"Test 11: Matrix comparison");

    return;
}

/***************************************************************************
*                              Initializer m_test                          *
***************************************************************************/
void TestGMatrix::init_matrix(void)
{
    double g_matrix[] = {1.0, 2.0, 3.0, 4.0,
                         5.0, 6.0, 7.0, 8.0,
                         9.0, 10., 11., 12.};

    for(int i=0;i<12;i++){
        m_g_matrix[i]=g_matrix[i];
    }

    return;
}

/***************************************************************************
*                              Initializer m_vector                        *
***************************************************************************/
void TestGMatrix::init_vector(void)
{
    double   g_vector[] = {1.0, 2.0, 3.0, 4.0};

    for(int i=0;i<4;i++){
        m_g_vector[i]=g_vector[i];
    }

    return;
}

/***************************************************************************
 *                              Set test matrix                            *
 ***************************************************************************/
void TestGMatrix::set_matrix(void)
{
    GMatrix matrix(m_g_rows,m_g_cols);
    for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col)
        matrix(row,col) = m_g_matrix[col+row*m_g_cols];
    }

    m_test = matrix;

    return;
}


/***************************************************************************
 *                              Set test vector                            *
 ***************************************************************************/
void TestGMatrix::set_vector(void)
{
    GVector vector(m_g_cols);
    for (int col = 0; col < m_g_cols; ++col){
            vector[col] = m_g_vector[col];
    }

    v_test = vector;

    return;
}

void TestGMatrix::set_bigger(void)
{
    bigger= GMatrix(m_g_rows+1,m_g_cols+1); 
}
/***************************************************************************
 *                             Check test matrix                           *
 ***************************************************************************/
int TestGMatrix::check_matrix(const GMatrix& m, const double scale, const double add)
{
    int result = 1;

    for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col) {
            double value = m_g_matrix[col+row*m_g_cols] * scale + add;
            if (fabs(m(row,col)-value) > 1.0e-15) {
                result = 0;
                break;
                }
        }
        }

    return result;
}


/***************************************************************************
 *                        Check transposed test matrix                     *
 ***************************************************************************/
int TestGMatrix::check_transpose_matrix(const GMatrix& m, const double scale, const double add)
{
    int result = 1;

    for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col) {
            double value = m_g_matrix[col+row*m_g_cols] * scale + add;
            if (m(col,row) != value) {
                result = 0;
                break;
                }
        }
        }

    return result;
}


/***************************************************************************
 *                     Check full test matrix (lower triangle)             *
 ***************************************************************************/
int TestGMatrix::check_matrix_lt(const GMatrix& m, const double scale, const double add)
{
    int result = 1;
    test_try("Check full test matrix (lower triangle)");
    try {
        for (int row = 0; row < m_g_rows; ++row) {
        for (int col = 0; col < m_g_cols; ++col) {
                double value = (col <= row) ? m_g_matrix[col+row*m_g_cols] * scale + add : 0.0;
                if (m(row,col) != value) {
                    result = 0;
                    break;
                    }
            }
            }
            test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                     Check full test matrix (upper triangle)             *
 ***************************************************************************/
int TestGMatrix::check_matrix_ut(const GMatrix& m, const double scale, const double add)
{
    int result = 1;
    test_try("Check full test matrix (upper triangle)");
    try {
        for (int row = 0; row < m_g_rows; ++row) {
        for (int col = 0; col < m_g_cols; ++col) {
                double value = (col >= row) ? m_g_matrix[col+row*m_g_cols] * scale + add : 0.0;
                if (m(row,col) != value) {
                    result = 0;
                    break;
                    }
            }
            }
            test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                             Check matrix*vector                         *
 ***************************************************************************/
int TestGMatrix::check_matrix_vector(const GVector& v)
{
    int result = 1;
    test_try(" Check matrix*vector");
    try {
        for (int row = 0; row < v.size(); ++row) {
        double value = 0.0;
            for (int col = 0; col < m_g_cols; ++col)
                value += m_g_matrix[col+row*m_g_cols] * m_g_vector[col];
            if (v[row] != value) {
                result = 0;
                break;
            }
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    return result;
}


/***************************************************************************
 *                             Check matrix*matrix                         *
 ***************************************************************************/
int TestGMatrix::check_matrix_matrix(const GMatrix& m)
{
  if (m.rows() != m_g_rows || m.cols() != m_g_rows)
    return 0;
  int result = 1;
  test_try("Check matrix*matrix");
  try {
    for (int row = 0; row < m.rows(); ++row) {
	  for (int col = 0; col < m.cols(); ++col) {
        double value = 0.0;
		for (int i = 0; i < m_g_cols; ++i)
	      value += m_g_matrix[i+row*m_g_cols] * m_g_matrix[i+col*m_g_cols];
	    if (m(row,col) != value) {
	      result = 0;
	      break;
	    }
      }
	}
        test_try_success();
  }
  catch (std::exception &e) {
      test_try_failure(e);
  }
  return result;
}


/***************************************************************************
 *                              Check matrix min                           *
 ***************************************************************************/
int TestGMatrix::check_matrix_min(const double min)
{
  double value = m_g_matrix[0];
  for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col) {
	  if (m_g_matrix[col+row*m_g_cols] < value)
		value = m_g_matrix[col+row*m_g_cols];
	}
  }
  test_assert(min==value, "Check matrix min");
  return (min == value);
}


/***************************************************************************
 *                              Check matrix max                           *
 ***************************************************************************/
int TestGMatrix::check_matrix_max(const double max)
{
  double value = m_g_matrix[0];
  for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col) {
	  if (m_g_matrix[col+row*m_g_cols] > value)
		value = m_g_matrix[col+row*m_g_cols];
	}
  }
  test_assert(max==value,"Check matrix max");
  return (max == value);
}


/***************************************************************************
 *                              Check matrix sum                           *
 ***************************************************************************/
int TestGMatrix::check_matrix_sum(const double sum)
{
  double value = 0.0;
  for (int row = 0; row < m_g_rows; ++row) {
    for (int col = 0; col < m_g_cols; ++col)
      value += m_g_matrix[col+row*m_g_cols];
  }
  test_assert(sum==value,"Check matrix sum");
  return (sum == value);
}


/***************************************************************************
 *                                Test: Output                             *
 ***************************************************************************/
void TestGMatrix::test_output(void)
{
      std::cout << "Output test matrix"<<std::endl<< m_test << std::endl;
}


/***************************************************************************
 *                     Test: Conversion between matrix types               *
 ***************************************************************************/
void TestGMatrix::test_conversion(void)
{
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
    test_assert(symmetric == back_convert,"Compare back converted matrix to original one");
    //
    // Determine the fill of the matrix. It should be 1.0
    double fill = back_convert.fill();
    test_assert(fabs(fill-1.0) <= 1.0e-15,"Test the fill of the matrix","Bad fill "+str(fill)+" determined (expected 1.0).");
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
    test_assert(ok,"Extract lower triangle");

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

    test_assert(ok,"Extract upper triangle","Corrupt extract_upper_triangle.");

    //
    // Now make the matrix unsymmetric
    symmetric(0,num-1) = 1000.0;
    //
    // Try converting now into GSymMatrix object (this should fail)
    test_try("Try converting now into GSymMatrix object");
    try {
        converted = sym_matrix(symmetric);
        test_try_failure();
    }
    catch (GException::matrix_not_symmetric &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
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
    test_assert(symmetric == back_convert,
                            "Compare back converted matrix to original one",
                            "Unable to convert matrixes (sparse)");

    // Return
    return;
}


/***************************************************************************
 *                        Test: extraction and insertion                   *
 ***************************************************************************/
void TestGMatrix::test_extract(void)
{

    //
    // Set-up test matrix
    int rows = 10;
    int cols = 20;
    GMatrix test(rows, cols);

    //
    // Add and extract column vectors
    test_try("Add and extract column vectors");
    try
    {
        for (int col = 0; col < cols; ++col) {
            GVector column(rows);
            for (int row = 0; row < rows; ++row)
                column[row] = (col+1)*100.0 + (row+1)*1.0;
            test.add_col(column, col);
            GVector check = test.extract_col(col);
            if (check != column) {
                throw exception_failure("Unable to add and extract columns.\nAdded column ...: "+column.print()+"\nExtracted column: "+check.print());
            }
        }
            test_try_success();
    }
    catch(std::exception& e)
    {
        test_try_failure(e);
    }

    //
    // Insert and extract column vectors
    test_try("Insert and extract column vectors");
    try
    {
        for (int col = 0; col < cols; ++col) {
        GVector column(rows);
        for (int row = 0; row < rows; ++row)
            column[row] = (col+1)*100.0 + (row+1)*1.0;
        test.insert_col(column, col);
        GVector check = test.extract_col(col);
        if (check != column) {
            throw exception_failure("Unable to insert and extract columns.\nIserted column ...: "+column.print()+"\nExtracted column: "+check.print());
        }
        }
        test_try_success();
    }
    catch(std::exception& e)
    {
        test_try_failure(e);
    }

    //
    // Extract rows
    test_try("Extract rows");
    try
    {
        for (int row = 0; row < rows; ++row) {
        GVector v_row(cols);
        for (int col = 0; col < cols; ++col)
            v_row[col] = (col+1)*100.0 + (row+1)*1.0;
        GVector check = test.extract_row(row);
        if (check != v_row) {
            throw exception_failure("Unable to extract rows.\nIserted row ...: "+v_row.print()+"\nExtracted row: "+check.print());
        }
        }

        test_try_success();
    }
    catch(std::exception& e)
    {
        test_try_failure(e);
    }

}

// Test 1:  Allocate zero matrix
void TestGMatrix::test1(void)
{
    GMatrix test1(0,0);

    return;
}

// Test 2:  Allocate too large matrix
void TestGMatrix::test2(void)
{
    //GMatrix test2(100000,100000);

    return;
}

// Test 3: Assign values
void TestGMatrix::test3(void)
{
    GMatrix test3(3,3);
    test3(1,1) = 1.0;
    test_assert(test3(0,0) == 0.0 && test3(1,0) == 0.0 && test3(2,0) == 0.0 &&
            test3(0,1) == 0.0 && test3(1,1) == 1.0 && test3(2,1) == 0.0 &&
            test3(0,2) == 0.0 && test3(1,2) == 0.0 && test3(2,2) == 0.0,"Assign value");

        #if defined(G_RANGE_CHECK)
        test_try("Range check");
        try {
            GMatrix test3(3,3);
            test3(3,3) = 1.0;
            test_try_failure();
        }
        catch (GException::out_of_range &e) {
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }
        #endif

    return;
}

// Test 4: Matrix copy constructor
void TestGMatrix::test4(void)
{
    GMatrix test4 = m_test;
    if (!check_matrix(test4, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
        throw exception_failure("Corrupt copy constructor");
    }

    return;
}

// Test 5: Matrix assignment
void TestGMatrix::test5(void)
{
    //
    // GMatrix = GMatrix
    result = m_test;
    if (!check_matrix(result, 1.0, 0.0) || !check_matrix(m_test, 1.0, 0.0)) {
        throw exception_failure("Corrupt assignment");
    }
    //
        // GMatrix = GMatrix (bigger matrix)
    result = bigger;

    return;
}

// Test 6: Transposition
void TestGMatrix::test6(void)
{
    //
    // transpose(GMatrix)
    result = transpose(m_test);
    test_assert(check_transpose_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"transpose(GMatrix)","Corrupt transpose(GMatrix) function.");
    //
    // GMatrix.transpose()
    result = m_test;
    result.transpose();
    test_assert(check_transpose_matrix(result, 1.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix.transpose()","Corrupt GMatrix.transpose() function.");
}

// Test 7: Matrix*Vector multiplication
void TestGMatrix::test7(void)
{
    GVector v_test7 = m_test*v_test;
    test_assert(check_matrix_vector(v_test7) && check_matrix(m_test, 1.0, 0.0),"","Corrupt Matrix*Vector multiplication");

    test_try("Bigger*Vector");
    try {
        GVector v_test7 = bigger*v_test;
        test_try_failure();
    }
    catch (GException::matrix_vector_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    return;
}

// Test 8: Matrix*Matrix multiplication
void TestGMatrix::test8(void)
{
    GMatrix m_test8 = m_test * transpose(m_test);
    test_assert(check_matrix_matrix(m_test8) && check_matrix(m_test, 1.0, 0.0),"","Corrupt Matrix*Matrix multiplication");

    test_try("Bigger*Matrix multiplication");
    try {
        GMatrix m_test8 = m_test*bigger;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    test_try("Bigger*Matrix multiplication 2");
    try {
        GMatrix m_test8 = bigger*m_test;
        test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    return;
}

// Test 9: Assignment and arithmetics
void TestGMatrix::test9(void)
{
    //
    // -GMatrix
    result = -m_test;
    test_assert(check_matrix(result, -1.0, 0.0)&&check_matrix(m_test, 1.0, 0.0),"-GMatrix","Corrupt -GMatrix operator.");

    //
    // GMatrix += GMatrix
    result  = m_test;
    result += m_test;
    test_assert(check_matrix(result, 2.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix += GMatrix","Corrupt GMatrix += GMatrix operator.");

    //
    // GMatrix -= GMatrix
    result  = m_test;
    result -= m_test;
    test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix -= GMatrix","Corrupt GMatrix -= GMatrix operator.");

    //
    // GMatrix *= 3.0
    result  = m_test;
    result *= 3.0;
    test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix *= 3.0","Corrupt GMatrix *= double operator.");

    //
    // GMatrix /= 3.0
    result  = m_test;
    result /= 3.0;
    test_assert(check_matrix(result, 1.0/3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix /= 3.0","Corrupt GMatrix /= double operator.");

    //
    // GMatrix + GMatrix
    result = m_test + m_test;
    test_assert(check_matrix(result, 2.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix + GMatrix","Corrupt GMatrix + GMatrix operator.");

    //
    // GMatrix - GMatrix
    result = m_test - m_test;
    test_assert(check_matrix(result, 0.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix - GMatrix","Corrupt GMatrix - GMatrix operator.");

    // GMatrix * 3.0
    result = m_test * 3.0;
    test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix * 3.0","Corrupt GMatrix * double operator.");

    //
    // 3.0 * GMatrix
    result = 3.0 * m_test;
    test_assert(check_matrix(result, 3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"3.0 * GMatrix","Corrupt double * GMatrix operator.");

    //
    // GMatrix / 3.0
    result = m_test / 3.0;
    test_assert(check_matrix(result, 1.0/3.0, 0.0) && check_matrix(m_test, 1.0, 0.0),"GMatrix / 3.0","Corrupt GMatrix / double operator.");


    test_try("Test += with bigger");
    try {
        result  = m_test;
        result += bigger;
        test_try_failure();
    }
    catch (GException::matrix_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    return;
}

// Test 10: Matrix functions
void TestGMatrix::test10(void){
    //
    double min = m_test.min();
    test_assert(check_matrix_min(min),"Test min","Corrupt GMatrix.min() function.");

    //
    double max = m_test.max();
    test_assert(check_matrix_max(max),"Test max","Corrupt GMatrix.max() function.");

    //
    double sum = m_test.sum();
    test_assert(check_matrix_sum(sum),"Test sum","Corrupt GMatrix.sum() function.");

    return;
}

// Test 11: Matrix comparison
void TestGMatrix::test11(void){

    test_assert(((m_test == m_test) == 1),"GMatrix == GMatrix","Corrupt GMatrix == GMatrix operator.");
    //
    GMatrix m_test10(m_g_rows,m_g_cols);
    test_assert(((m_test == m_test10) == 0),"GMatrix == GMatrix 2","Corrupt GMatrix == GMatrix operator.");
    //
    test_assert(((m_test == bigger) == 0),"GMatrix == GMatrix 3","Corrupt GMatrix == GMatrix operator.");
    //
    test_assert(((m_test != m_test) == 0),"GMatrix != GMatrix","Corrupt GMatrix != GMatrix operator.");
    //
    test_assert(((m_test != m_test10) == 1),"GMatrix != GMatrix 2","Corrupt GMatrix != GMatrix operator.");
    //
    test_assert(((m_test != bigger) == 1),"GMatrix != GMatrix 3","Corrupt GMatrix != GMatrix operator.");

    return;
}

int main(void)
{
    GTestSuites testsuite("GMatrix");

    bool was_successful=true;

    //Create a test suite
    TestGMatrix test;

    //Append to the container
    testsuite.append(test);

    //Run
    was_successful=testsuite.run();

    //save xml report
    testsuite.save("reports/GMatrix.xml");

    // Return
    return was_successful ? 0:1;
}
