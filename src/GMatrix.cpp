/***************************************************************************
 *                       GMatrix.cpp  -  matrix class                      *
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
#include "GMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/*==========================================================================
 =                                                                         =
 =                      GMatrix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *   GMatrix private constructor (used for derived class initialization)   *
 ***************************************************************************/
GMatrix::GMatrix(void)
{
  // Initialise private members for clean destruction
  m_rows       = 0;
  m_cols       = 0;
  m_elements   = 0;
  m_num_rowsel = 0;
  m_num_colsel = 0;
  m_data       = NULL;
  m_colstart   = NULL;
  m_rowsel     = NULL;
  m_colsel     = NULL;
	
  // Return
  return;
}


/***************************************************************************
 *                           GMatrix constructor                           *
 ***************************************************************************/
GMatrix::GMatrix(int rows, int cols)
{
  // Initialise private members for clean destruction
  m_rows       = 0;
  m_cols       = 0;
  m_elements   = 0;
  m_num_rowsel = 0;
  m_num_colsel = 0;
  m_data       = NULL;
  m_colstart   = NULL;
  m_rowsel     = NULL;
  m_colsel     = NULL;

  // Determine number of physical elements in matrix
  int elements = rows*cols;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw empty("GMatrix constructor");
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new int[cols+1];
  if (m_data == NULL || m_colstart == NULL)
	throw mem_alloc("GMatrix constructor", elements);
	
  // Store matrix size (logical and physical)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;
  
  // Set-up column start indices
  m_colstart[0] = 0;
  for (int col = 1; col <= m_cols; ++col)
    m_colstart[col] = m_colstart[col-1] + m_rows;
	
  // Initialise matrix elements to 0.0
  for (int i = 0; i < m_elements; ++i)
    m_data[i] = 0.0;
	
  // Return
  return;
}


/***************************************************************************
 *                           GMatrix copy constructor                      *
 * ----------------------------------------------------------------------- *
 * The copy constructor is sufficiently general to provide the base        *
 * constructor for all derived classes, including sparse matrices.         *
 ***************************************************************************/
GMatrix::GMatrix(const GMatrix& m)
{
  // Initialise private members for clean destruction
  m_rows       = 0;
  m_cols       = 0;
  m_elements   = 0;
  m_num_rowsel = 0;
  m_num_colsel = 0;
  m_data       = NULL;
  m_colstart   = NULL;
  m_rowsel     = NULL;
  m_colsel     = NULL;
  
  // Allocate memory for column start array and copy content
  m_colstart = new int[m.m_cols+1];
  if (m_colstart == NULL)
    throw mem_alloc("GMatrix copy constructor", m.m_cols+1);
  for (int i = 0; i <= m.m_cols; ++i)
    m_colstart[i] = m.m_colstart[i];
  
  // Allocate memory for elements and copy them (only if there are elements)
  if (m.m_elements > 0) {
    m_data = new double[m.m_elements];
    if (m_data == NULL)
	  throw mem_alloc("GMatrix copy constructor", m.m_elements);
    for (int i = 0; i < m.m_elements; ++i)
      m_data[i] = m.m_data[i];
  }
  
  // If there is a row selection then copy it
  if (m.m_rowsel != NULL && m.m_num_rowsel > 0) {
    m_rowsel = new int[m.m_num_rowsel];
    if (m_rowsel == NULL)
	  throw mem_alloc("GMatrix copy constructor", m.m_num_rowsel);
    for (int i = 0; i < m.m_num_rowsel; ++i)
      m_rowsel[i] = m.m_rowsel[i];
	m_num_rowsel = m.m_num_rowsel;
  }

  // If there is a column selection then copy it
  if (m.m_colsel != NULL && m.m_num_colsel > 0) {
    m_colsel = new int[m.m_num_colsel];
    if (m_colsel == NULL)
	  throw mem_alloc("GMatrix copy constructor", m.m_num_colsel);
    for (int i = 0; i < m.m_num_colsel; ++i)
      m_colsel[i] = m.m_colsel[i];
	m_num_colsel = m.m_num_colsel;
  }

  // Copy matrix attributes
  m_rows     = m.m_rows;
  m_cols     = m.m_cols;
  m_elements = m.m_elements;
  
  // Return
  return;
}


/***************************************************************************
 *                           GMatrix destructor                            *
 ***************************************************************************/
GMatrix::~GMatrix()
{
  // De-allocate only if memory has indeed been allocated
  if (m_colstart != NULL) delete [] m_colstart;
  if (m_rowsel   != NULL) delete [] m_rowsel;
  if (m_colsel   != NULL) delete [] m_colsel;
  if (m_data     != NULL) delete [] m_data;

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                             GMatrix operators                           =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        GMatrix assignment operator                      *
 ***************************************************************************/
GMatrix& GMatrix::operator= (const GMatrix& m)
{
  // Execute only if object is not identical
  if (this != &m) { 
  
    // De-allocate memory if it has indeed been allocated
    if (m_colstart != NULL) delete [] m_colstart;
    if (m_rowsel   != NULL) delete [] m_rowsel;
    if (m_colsel   != NULL) delete [] m_colsel;
    if (m_data     != NULL) delete [] m_data;

    // Initialise private members for clean destruction
    m_rows       = 0;
    m_cols       = 0;
    m_elements   = 0;
    m_num_rowsel = 0;
    m_num_colsel = 0;
    m_data       = NULL;
    m_colstart   = NULL;
    m_rowsel     = NULL;
    m_colsel     = NULL;
  
    // Allocate memory for column start array and copy content
    m_colstart = new int[m.m_cols+1];
    if (m_colstart == NULL)
      throw mem_alloc("GMatrix::operator= (const GMatrix&)", m.m_cols+1);
    for (int i = 0; i <= m.m_cols; ++i)
      m_colstart[i] = m.m_colstart[i];
  
    // Allocate memory for elements and copy them (only if there are elements)
    if (m.m_elements > 0) {
      m_data = new double[m.m_elements];
      if (m_data == NULL)
	    throw mem_alloc("GMatrix::operator= (const GMatrix&)", m.m_elements);
      for (int i = 0; i < m.m_elements; ++i)
        m_data[i] = m.m_data[i];
    }

    // If there is a row selection then copy it
    if (m.m_rowsel != NULL && m.m_num_rowsel > 0) {
      m_rowsel = new int[m.m_num_rowsel];
      if (m_rowsel == NULL)
  	    throw mem_alloc("GMatrix::operator= (const GMatrix&)", m.m_num_rowsel);
      for (int i = 0; i < m.m_num_rowsel; ++i)
        m_rowsel[i] = m.m_rowsel[i];
	  m_num_rowsel = m.m_num_rowsel;
    }

    // If there is a column selection then copy it
    if (m.m_colsel != NULL && m.m_num_colsel > 0) {
      m_colsel = new int[m.m_num_colsel];
      if (m_colsel == NULL)
	    throw mem_alloc("GMatrix::operator= (const GMatrix&)", m.m_num_colsel);
      for (int i = 0; i < m.m_num_colsel; ++i)
        m_colsel[i] = m.m_colsel[i];
	  m_num_colsel = m.m_num_colsel;
    }

    // Copy matrix attributes
    m_rows     = m.m_rows;
    m_cols     = m.m_cols;
    m_elements = m.m_elements;
  } // endif: object was not identical

  // Return this object
  return *this;
}


/***************************************************************************
 *                           Vector multiplication                         *
 ***************************************************************************/
GVector GMatrix::operator* (const GVector& v) const
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_cols != v.m_num)
    throw matrix_vector_mismatch("GMatrix::operator* (const GVector&)", 
	                             v.m_num, m_rows, m_cols);

  // Perform vector multiplication
  GVector result(m_rows);
  for (int row = 0; row < m_rows; ++row) {
    double sum = 0.0;
	for (int col = 0; col < m_cols; ++col)
	  sum += (*this)(row,col) * v(col);
	result(row) = sum;
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                   GMatrix comparison (equality operator)                *
 ***************************************************************************/
int GMatrix::operator== (const GMatrix &m) const
{
  // Initalise the result to 'equal matrices'
  int result = 1;
  
  // Perform comparison (only if matrix dimensions and number of physical
  // elements are identical; otherwise set result to 'false')
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (int i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 0;
		break;
	  }
	}
  }
  else
    result = 0;

  // Return result
  return result;
}


/***************************************************************************
 *                 GMatrix comparison (non-equality operator)              *
 ***************************************************************************/
int GMatrix::operator!= (const GMatrix &m) const
{
  // Initalise the result to 'equal matrices'
  int result = 0;

  // Perform comparison (only if matrix dimensions and number of physical
  // elements are identical; otherwise set result to 'true')
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (int i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 1;
		break;
	  }
	}
  }
  else
    result = 1;
	
  // Return result
  return result;
}


/***************************************************************************
 *                  GMatrix negation (unary minus operator)                *
 ***************************************************************************/
GMatrix GMatrix::operator- ( ) const
{
  // Copy argument into result
  GMatrix result = *this;
  
  // Negate result
  for (int i = 0; i < m_elements; ++i)
    result.m_data[i] = -result.m_data[i];
	
  // Return result	
  return result;
}


/***************************************************************************
 *                  GMatrix addition (unary addition operator)             *
 ***************************************************************************/
GMatrix& GMatrix::operator+= (const GMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw matrix_mismatch("GMatrix::operator+= (const GMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Add matrices
  for (int i = 0; i < m_elements; ++i)
    m_data[i] += m.m_data[i];

  // Return result
  return *this;
}


/***************************************************************************
 *             GMatrix subtraction (unary subtraction operator)            *
 ***************************************************************************/
GMatrix& GMatrix::operator-= (const GMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw matrix_mismatch("GMatrix::operator-= (const GMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Subtract matrices
  for (int i = 0; i < m_elements; ++i)
    m_data[i] -= m.m_data[i];

  // Return result
  return *this;
}


/***************************************************************************
 *          Matrix multiplication (unary multiplication operator)          *
 * ----------------------------------------------------------------------- *
 * In case of rectangular matrices the result matrix does not change and   *
 * the operations can be performed 'quasi' inplace (except of a vector).   *
 * For the general case the result matrix changes the size so for          *
 * simplicity a new matrix is allocated to hold the results.               *
 ***************************************************************************/
GMatrix& GMatrix::operator*= (const GMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows)
    throw matrix_mismatch("GMatrix::operator*= (const GMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Case A: Matrices are rectangular, so perform 'inplace' multiplication
  if (m_rows == m_cols) {
    for (int row = 0; row < m_rows; ++row) {
      GVector v_row = extract_row(row);
      for (int col = 0; col < m_cols; ++col) {
	    double sum = 0.0;
	    for (int i = 0; i < m_cols; ++i)
	      sum += v_row(i) * m(i,col);
        (*this)(row,col) = sum;
      }
    }
  }
  
  // Case B: Matrices are not rectangular, so we cannot work inplace
  else {
  
    // Allocate result matrix
	GMatrix result(m_rows, m.m_cols);
	
	// Loop over all elements of result matrix
    for (int row = 0; row < m_rows; ++row) {
      for (int col = 0; col < m.m_cols; ++col) {
	    double sum = 0.0;
	    for (int i = 0; i < m_cols; ++i)
	      sum += (*this)(row,i) * m(i,col);
		result(row,col) = sum;
	  }
	}
	
	// Assign result
	*this = result;
	
  } // endelse: matrices were not rectangular

  // Return result
  return *this;
}


/***************************************************************************
 *  GMatrix scalar multiplication (unary scalar multiplication operator)   *
 ***************************************************************************/
GMatrix& GMatrix::operator*= (const double& d)
{
  // Case A: If multiplicator is 0.0 then set entire matrix to 0.0
  if (d == 0.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = 0.0;
  }
  
  // Case B: If multiplicator is not +/- 1.0 then perform multiplication
  else if (fabs(d) != 1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] *= d;
  }
  
  // Case C: If multiplication is -1.0 then negate
  else if (d == -1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = -m_data[i];
  }
  
  // Return result
  return *this;
}


/***************************************************************************
 *         Matrix scalar division (unary scalar division operator)         *
 ***************************************************************************/
GMatrix& GMatrix::operator/= (const double& d)
{
  // Case A: If divider is not +/- 1.0 then perform division
  if (fabs(d) != 1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] /= d;
  }
  
  // Case B: If divider is -1.0 then negate
  else if (d == -1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = -m_data[i];
  }
  
  // Return result
  return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GMatrix member functions                       =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                       Set all matrix elements to zero                   *
 ***************************************************************************/
void GMatrix::clear()
{
  // Loop over all matrix elements
  for (int i = 0; i < m_elements; ++i)
    m_data[i] = 0.0;
  
  // Return
  return;
}


/***************************************************************************
 *                    Get minimum physical matrix element                  *
 ***************************************************************************/
double GMatrix::min() const
{
  // Initialise minimum with first element
  double result = m_data[0];
  
  // Search all elements for the smallest one
  for (int i = 1; i < m_elements; ++i) {
    if (m_data[i] < result)
	  result = m_data[i];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                    Get maximum physical matrix element                  *
 ***************************************************************************/
double GMatrix::max() const
{
  // Initialise minimum with first element
  double result = m_data[0];
  
  // Search all elements for the largest one
  for (int i = 1; i < m_elements; ++i) {
    if (m_data[i] > result)
	  result = m_data[i];
  }

  // Return result
  return result;
}


/***************************************************************************
 *                              Get matrix sum                             *
 ***************************************************************************/
double GMatrix::sum() const
{
  // Initialise matrix sum
  double result = 0.0;
  
  // Add all elements  
  for (int i = 0; i < m_elements; ++i)
    result += m_data[i];
  
  // Return result
  return result;
}


/***************************************************************************
 *                         Inplace transpose of matrix                     *
 * ----------------------------------------------------------------------- *
 * The transpose operation exchanges the number of rows against the number *
 * of columns. The element transformation is done inplace.                 * 
 ***************************************************************************/
void GMatrix::transpose()
{
  // Case A: matrix is rectangular then simply swap the elements
  if (m_rows == m_cols) {
    double  swap;
    double* ptr_dst;
    double* ptr_src;
    for (int row = 0; row < m_rows; ++row) {
      for (int col = row; col < m_cols; ++col) {
	    ptr_dst  = m_data + m_cols*row + col;
	    ptr_src  = m_data + m_rows*col + row;
	    swap     = *ptr_dst;
	    *ptr_dst = *ptr_src;
	    *ptr_src = swap;
	  }
    }
  }
  
  // Case B: Dummy for non-rectangular transpose
  else {
    GMatrix result(m_cols, m_rows);
    for (int row = 0; row < m_rows; ++row) {
      for (int col = 0; col < m_cols; ++col)
	    result(col, row) = (*this)(row, col);
	}
	*this = result;
  }
  
}


/***************************************************************************
 *                     Extract row from matrix into vector                 *
 ***************************************************************************/
GVector GMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw out_of_range("GMatrix::extract_row(int)", 
	                   row, 0, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_cols);
  
  // Extract row into vector
  for (int col = 0; col < m_cols; ++col)
    result(col) = (*this)(row,col);
	
  // Return vector
  return result;

}


/***************************************************************************
 *                   Extract column from matrix into vector                *
 ***************************************************************************/
GVector GMatrix::extract_col(int col) const
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw out_of_range("GMatrix::extract_col(int)", 
	                   0, col, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_rows);
  
  // Extract column into vector
  for (int row = 0; row < m_rows; ++row)
    result(row) = (*this)(row,col);
	
  // Return vector
  return result;

}


/*==========================================================================
 =                                                                         =
 =                     GMatrix private member functions                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                              select_non_zero                            *
 * ----------------------------------------------------------------------- *
 * Determines the non-zero rows and columns in matrix and set up index     *
 * arrays that points to these rows/columns. These arraya are used for     *
 * compressed matrix factorisations.                                       *
 ***************************************************************************/
void GMatrix::select_non_zero(void)
{
  // Free existing selection arrays
  if (m_rowsel != NULL) delete [] m_rowsel;
  if (m_colsel != NULL) delete [] m_colsel;
  
  // Allocate selection arrays
  m_rowsel = new int[m_rows];
  if (m_rowsel == NULL)
    throw mem_alloc("GMatrix::select_non_zero()", m_rows);
  m_colsel = new int[m_cols];
  if (m_colsel == NULL)
    throw mem_alloc("GMatrix::select_non_zero()", m_cols);

  // Initialise non-zero row and column counters
  m_num_rowsel = 0;
  m_num_colsel = 0;
  
  // Declare loop variables
  int row;
  int col;
  
  // Find all non-zero rows
  for (row = 0; row < m_rows; ++row) {
    for (col = 0; col < m_cols; ++col) {
	  if ((*this)(row,col) != 0.0)
	    break;
	}
	if (col < m_cols)      // found a non-zero element in row
      m_rowsel[m_num_rowsel++] = row;
  }

  // Find all non-zero columns
  for (col = 0; col < m_cols; ++col) {
    for (row = 0; row < m_rows; ++row) {
	  if ((*this)(row,col) != 0.0)
	    break;
	}
	if (row < m_rows)      // found a non-zero element in column
      m_colsel[m_num_colsel++] = col;
  }
  
  // Return
  return;
}

/*==========================================================================
 =                                                                         =
 =                             GMatrix friends                             =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GMatrix& m)
{
  // Loop over all matrix elements
  for (int row = 0; row < m.m_rows; ++row) {
    for (int col = 0; col < m.m_cols; ++col) {
      os << m(row,col);
	  if (col != m.m_cols-1)
	    os << ", ";
	}
	if (row != m.m_rows-1)
	  os << endl;
  }
  
  // Return output stream
  return os;
}


/***************************************************************************
 *                          Matrix absolute values                         *
 ***************************************************************************/
GMatrix fabs(const GMatrix& m)
{
  // Define result matrix
  GMatrix result = m;
  
  // Convert all elements to absolute values  
  for (int i = 0; i < result.m_elements; ++i)
    result.m_data[i] = fabs(result.m_data[i]);
  
  // Return result
  return result;
}


/*==========================================================================
 =                                                                         =
 =                        GMatrix exception classes                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                      Matrix row or column out of range                  *
 ***************************************************************************/
GMatrix::out_of_range::out_of_range(string origin, int row, int col, 
                                    int rows, int cols)
{
  m_origin = origin;
  ostringstream s_row;
  ostringstream s_col;
  ostringstream s_rows;
  ostringstream s_cols;
  s_row  << row;
  s_col  << col;
  s_rows << rows-1;
  s_cols << cols-1;
  m_message = "Matrix element (" + s_row.str() + "," + s_col.str() +
              ") out of range ([0," + s_rows.str() + "], [0," +
			  s_cols.str() + "])";
}


/***************************************************************************
 *                     Mismatch between matrix and vector                  *
 ***************************************************************************/
GMatrix::matrix_vector_mismatch::matrix_vector_mismatch(string origin, 
                                                int num, int rows, int cols)
{
  m_origin = origin;
  ostringstream s_num;
  ostringstream s_rows;
  ostringstream s_cols;
  s_num  << num;
  s_rows << rows;
  s_cols << cols;
  m_message = "Vector dimension [" + s_num.str() + 
              "] is incompatible with matrix size [" + 
			  s_rows.str() + "," + s_cols.str() + "]";
}


/***************************************************************************
 *                       Matrix mismatch in operation                      *
 ***************************************************************************/
GMatrix::matrix_mismatch::matrix_mismatch(string origin, int rows1, 
                                            int cols1, int rows2, int cols2)
{
  m_origin = origin;
  ostringstream s_rows1;
  ostringstream s_rows2;
  ostringstream s_cols1;
  ostringstream s_cols2;
  s_rows1 << rows1;
  s_cols1 << cols1;
  s_rows2 << rows2;
  s_cols2 << cols2;
  m_message = "Matrix mismatch: M1(" + s_rows1.str() + "," + s_cols1.str() +
              ") incompatible with M2(" + s_rows2.str() + "," + 
			  s_cols2.str() + ")";
}


/***************************************************************************
 *                           Matrix not rectangular                        *
 ***************************************************************************/
GMatrix::matrix_not_rect::matrix_not_rect(string origin, int rows1, 
                                            int cols1, int rows2, int cols2)
{
  m_origin = origin;
  ostringstream s_rows1;
  ostringstream s_rows2;
  ostringstream s_cols1;
  ostringstream s_cols2;
  s_rows1 << rows1;
  s_cols1 << cols1;
  s_rows2 << rows2;
  s_cols2 << cols2;
  m_message = "Operation only valid for rectangular matrices: M[" + 
              s_rows1.str() + "," + s_cols1.str() +
              "] *= M[" + s_rows2.str() + "," + s_cols2.str() + 
			  "] not allowed";
}


/***************************************************************************
 *                     Matrix is not positive definite                     *
 ***************************************************************************/
GMatrix::matrix_not_pos_definite::matrix_not_pos_definite(string origin, 
                                                        int row, double sum)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_sum;
  s_rows << row;
  s_sum << scientific << sum;
  m_message = "Matrix is not positive definite (sum " + s_sum.str() + 
              " occured in row/column " + s_rows.str() + ")";
}


/***************************************************************************
 *                         Matrix is not symmetric                         *
 ***************************************************************************/
GMatrix::matrix_not_symmetric::matrix_not_symmetric(string origin, int cols, 
                                                                   int rows)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_cols;
  s_rows << rows;
  s_cols << cols;
  m_message = "Matrix is not symmetric [" + s_rows.str() + "," +
			  s_cols.str() + "]";
}


/***************************************************************************
 *                       Matrix has not been factorised                    *
 ***************************************************************************/
GMatrix::matrix_not_factorised::matrix_not_factorised(string origin, string type)
{
  m_origin  = origin;
  m_message = "Matrix has not been factorised using " + type;
}


/***************************************************************************
 *                       All matrix elements are zero                      *
 ***************************************************************************/
GMatrix::matrix_zero::matrix_zero(string origin)
{
  m_origin = origin;
  m_message = "All matrix elements are zero";
}
