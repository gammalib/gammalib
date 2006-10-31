/***************************************************************************
 *                  GSparseMatrix.cpp  -  sparse matrix class              *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * This class implements the compressed sparse column format. The          *
 * following arrays are allocated:                                         *
 * m_data     : holds all 'elements' non-zero values, in column order (Ax) *
 * m_colstart : holds the index of the first element of each column   (Ap) *
 * m_rowinx   : holds the row indices for all elements                (Ai) *
 * Column 'col' covers therefore [m_colstart[j], ..., m_colstart[j+1]-1]   *
 * Note that 'm_colstart' has m_cols+1 elements.                           *
/***************************************************************************

/* __ Includes ___________________________________________________________ */
#include "GSparseMatrix.hpp"
#include "GSparseSymbolic.hpp"
#include "GSparseNumeric.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Definitions ________________________________________________________ */
//#define G_DEBUG_SPARSE_PENDING                    // Analyse pending values
//#define G_DEBUG_SPARSE_INSERTION                 // Analyse value insertion
//#define G_DEBUG_SPARSE_COMPRESSION   // Analyse zero row/column compression


/*==========================================================================
 =                                                                         =
 =                   GSparseMatrix constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        GSparseMatrix constructor                        *
 ***************************************************************************/
GSparseMatrix::GSparseMatrix(int rows, int cols, int elements) : GMatrix()
{
  // Initialise private members for clean destruction
  m_rows     = 0;
  m_cols     = 0;
  m_elements = 0;
  m_data     = NULL;
  m_colstart = NULL;
  m_rowinx   = NULL;
  m_zero     = 0.0;
  m_fill_val = 0.0;
  m_fill_row = 0;
  m_fill_col = 0;
  m_symbolic = NULL;
  m_numeric  = NULL;
	
  // Allocate column start array. This is the only array that we can
  // allocate at this time. The other arrays can only be allocated during
  // filling of the matrix
  m_colstart = new int[cols+1];
  if (m_colstart == NULL)
	throw mem_alloc("GSparseMatrix constructor", cols+1);

  // Store (logical) matrix size
  m_rows = rows;
  m_cols = cols;

  // Initialise column start indices to 0
  for (int col = 0; col <= m_cols; ++col)
    m_colstart[col] = 0;
	
  // Optionally allocate elements
  if (elements > 0)
    alloc_elements(0, elements);
  
  // Return
  return;
}


/***************************************************************************
 *                       GSparseMatrix copy constructor                    *
 * ----------------------------------------------------------------------- *
 * First invoques GMatrix copy constructor (which sets the logical matrix  *
 * dimensions and copies the elements and the column start array). It      *
 * remains then to this constructor to copy the row indices.               *
 * NOTE: The copy constructor does NOT fill a pending element.             *
 ***************************************************************************/
GSparseMatrix::GSparseMatrix(const GSparseMatrix& m) : GMatrix(m)
{ 
  // Initialise private members for clean destruction
  m_rowinx   = NULL;
  m_symbolic = NULL;
  m_numeric  = NULL;

  // Copy data members
  m_zero     = m.m_zero;
  m_fill_val = m.m_fill_val;
  m_fill_row = m.m_fill_row;
  m_fill_col = m.m_fill_col;
  
  // Allocate memory for row indices and copy them (only if there are elements)
  if (m_elements > 0) {
    m_rowinx = new int[m_elements];
    if (m_rowinx == NULL)
	  throw mem_alloc("GSparseMatrix copy constructor", m_elements);
    for (int i = 0; i < m_elements; ++i)
      m_rowinx[i] = m.m_rowinx[i];
  }

  // If there is a symbolic decomposition then copy it
  if (m.m_symbolic != NULL) {
    GSparseSymbolic* symbolic = new GSparseSymbolic();
    if (symbolic == NULL)
	  throw mem_alloc("GSparseMatrix copy constructor", 1);
    *symbolic  = *((GSparseSymbolic*)m.m_symbolic);
    m_symbolic = (void*)symbolic;
  }

  // If there is a numeric decomposition then copy it
  if (m.m_numeric != NULL) {
    GSparseNumeric* numeric = new GSparseNumeric();
    if (numeric == NULL)
	  throw mem_alloc("GSparseMatrix copy constructor", 1);
    *numeric  = *((GSparseNumeric*)m.m_numeric);
    m_numeric = (void*)numeric;
  }
  
  // Return
  return;
}


/***************************************************************************
 *                        GSparseMatrix destructor                         *
 * ----------------------------------------------------------------------- *
 * First invoques this destructor to de-allocate row indices, then invoque *
 * ~GMatrix destructor to de-allocate the other arrays.                    *
 ***************************************************************************/
GSparseMatrix::~GSparseMatrix()
{
  // De-allocate only if memory has indeed been allocated by derived class
  if (m_rowinx   != NULL) delete [] m_rowinx;
  if (m_symbolic != NULL) delete (GSparseSymbolic*)m_symbolic;
  if (m_numeric  != NULL) delete (GSparseNumeric*)m_numeric;

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                         GSparseMatrix operators                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                     GSparseMatrix assignment operator                   *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator= (const GSparseMatrix& m)
{ 
  // Execute only if object is not identical
  if (this != &m) {

    // Invoque base class assignment operator
    this->GMatrix::operator=(m);

    // De-allocate memory if it has indeed been allocated
    if (m_rowinx   != NULL) delete [] m_rowinx;
    if (m_symbolic != NULL) delete (GSparseSymbolic*)m_symbolic;
    if (m_numeric  != NULL) delete (GSparseNumeric*)m_numeric;

    // Initialise private members for clean destruction
    m_rowinx   = NULL;
	m_symbolic = NULL;
	m_numeric  = NULL;

    // Copy data members
    m_zero     = m.m_zero;
    m_fill_val = m.m_fill_val;
    m_fill_row = m.m_fill_row;
    m_fill_col = m.m_fill_col;

    // Allocate memory for row indices and copy them (only if there are elements)
    if (m_elements > 0) {
      m_rowinx = new int[m_elements];
      if (m_rowinx == NULL)
	    throw mem_alloc("GSparseMatrix::operator= (const GSparseMatrix&)", 
		                m_elements);
      for (int i = 0; i < m_elements; ++i)
        m_rowinx[i] = m.m_rowinx[i];
    }

    // If there is a symbolic decomposition then copy it
    if (m.m_symbolic != NULL) {
      GSparseSymbolic* symbolic = new GSparseSymbolic();
      if (symbolic == NULL)
	    throw mem_alloc("GSparseMatrix copy constructor", 1);
      *symbolic  = *((GSparseSymbolic*)m.m_symbolic);
      m_symbolic = (void*)symbolic;
    }

    // If there is a numeric decomposition then copy it
    if (m.m_numeric != NULL) {
      GSparseNumeric* numeric = new GSparseNumeric();
      if (numeric == NULL)
	    throw mem_alloc("GSparseMatrix copy constructor", 1);
      *numeric  = *((GSparseNumeric*)m.m_numeric);
      m_numeric = (void*)numeric;
    }

  } // endif: object was not identical
  
  // Return
  return *this;
}


/***************************************************************************
 *                           Vector multiplication                         *
 * ----------------------------------------------------------------------- *
 * Performs sparse matrix*vector multiplication including an eventually    *
 * pending fill.                                                           *
 ***************************************************************************/
GVector GSparseMatrix::operator* (const GVector& v) const
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_cols != v.m_num)
    throw matrix_vector_mismatch("GSparseMatrix::operator* (const GVector&)", 
	                             v.m_num, m_rows, m_cols);

  // Initialise result vector
  GVector result(m_rows);

  // Multiply only if there are elements in matrix
  if (m_elements > 0) {

    // Perform vector multiplication
    for (int col = 0; col < m_cols; ++col) {
	  int  i_start    = m_colstart[col];
	  int  i_stop     = m_colstart[col+1];
      double*   ptr_data   = m_data   + i_start;
	  int* ptr_rowinx = m_rowinx + i_start;
      for (int i = i_start; i < i_stop; ++i)
	    result(*ptr_rowinx++) += *ptr_data++ * v(col);
    }
  
    // If fill is pending then add-in also this element into the product
    if (m_fill_val != 0.0) {
      result(m_fill_row) += m_fill_val * v(m_fill_col);
      #if defined(G_DEBUG_SPARSE_PENDING)
      cout << "GSparseMatrix::operator* (const GVector&) const: pending value " << 
			  m_fill_val << " for location (" <<  m_fill_row << "," << m_fill_col <<
			  ") has been used." << endl;
	  #endif
    }
	
  } // endif: there were elements in matrix
  
  // Return result
  return result;
}


/***************************************************************************
 *               GSparseMatrix comparison (equality operator)              *
 ***************************************************************************/
int GSparseMatrix::operator== (const GSparseMatrix &m) const
{
  // Initalise the result to 'equal matrices'
  int result = 1;
  
  // Perform comparison (only if matrix dimensions are identical)
  if (m_rows == m.m_rows && m_cols == m.m_cols) {
  
    // Loop over all matrix elements
    for (int row = 0; row < m_rows; ++row) {
	  for (int col = 0; col < m_cols; ++col) {
	    if ((*this)(row,col) != m(row,col)) {
	      result = 0;
		  break;
	    }
	  }
	  if (!result)
	    break;
    }	  
  } 
  else
	result = 0;

  // Return result
  return result;
}


/***************************************************************************
 *             GSparseMatrix comparison (non-equality operator)            *
 ***************************************************************************/
int GSparseMatrix::operator!= (const GSparseMatrix &m) const
{
  // Initalise the result to 'equal matrices'
  int result = 0;

  // Perform comparison (only if matrix dimensions are identical)
  if (m_rows == m.m_rows && m_cols == m.m_cols) {
  
    // Loop over all matrix elements
    for (int row = 0; row < m_rows; ++row) {
	  for (int col = 0; col < m_cols; ++col) {
	    if ((*this)(row,col) != m(row,col)) {
	      result = 1;
		  break;
	    }
	  }
	  if (!result)
	    break;
    }
  } 
  else
	result = 1;

  // Return result
  return result;
}


/***************************************************************************
 *                  GSparseMatrix negation (unary minus operator)          *
 ***************************************************************************/
GSparseMatrix GSparseMatrix::operator- ( ) const
{
  // Copy argument into result
  GSparseMatrix result = *this;

  // Fill pending matrix element
  result.fill_pending();
  
  // Negate result
  for (int i = 0; i < result.m_elements; ++i)
    result.m_data[i] = -result.m_data[i];
	
  // Return result	
  return result;
}


/***************************************************************************
 *             GSparseMatrix addition (unary addition operator)            *
 * ----------------------------------------------------------------------- *
 * NOT YET IN PURE NATIVE CODE !!!!!!                                      *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator+= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw matrix_mismatch("GSparseMatrix::operator+= (const GSparseMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Perform inplace matrix addition
  for (int col = 0; col < m_cols; ++col) {
    GVector v_result  = extract_col(col);
    GVector v_operand = m.extract_col(col);
	v_result += v_operand;
	insert_col(v_result, col);
  }

  // Return result
  return *this;
}


/***************************************************************************
 *            GSparseMatrix subtraction (unary addition operator)          *
 * ----------------------------------------------------------------------- *
 * NOT YET IN PURE NATIVE CODE !!!!!!                                      *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator-= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw matrix_mismatch("GSparseMatrix::operator-= (const GSparseMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Perform inplace matrix subtraction
  for (int col = 0; col < m_cols; ++col) {
    GVector v_result  = extract_col(col);
    GVector v_operand = m.extract_col(col);
	v_result -= v_operand;
	insert_col(v_result, col);
  }

  // Return result
  return *this;
}


/***************************************************************************
 *          Matrix multiplication (unary multiplication operator)          *
 * ----------------------------------------------------------------------- *
 * NOT YET IN PURE NATIVE CODE !!!!!!                                      *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator*= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows || m_rows != m.m_cols)
    throw matrix_mismatch("GSparseMatrix::operator*= (const GSparseMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Allocate result matrix
  GSparseMatrix result(m_rows, m.m_cols);

  // Multiply only if there are elements in both matrices
  if (m_elements > 0 && m.m_elements > 0) {

	// Loop over all elements of result matrix
    for (int row = 0; row < m_rows; ++row) {
      for (int col = 0; col < m.m_cols; ++col) {
	    double sum = 0.0;
	    for (int i = 0; i < m_cols; ++i)
	      sum += (*this)(row,i) * m(i,col);
		result(row,col) = sum;
	  }
	}

  }

  // Assign result
  *this = result;

  // Return result
  return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GSparseMatrix member functions                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        Get minimum matrix element                       *
 ***************************************************************************/
double GSparseMatrix::min() const
{
  // Initialise minimum with fill value
  double result = m_fill_val;
  
  // Search all elements for the smallest one
  for (int i = 0; i < m_elements; ++i) {
    if (m_data[i] < result)
	  result = m_data[i];
  }
  
  // If minimum is > 0.0 and there are zero elements then set the minimum
  // to 0.0
  if ((result > 0.0) && (m_elements < (m_rows*m_cols)))
    result = 0.0;
  
  // Return result
  return result;
}


/***************************************************************************
 *                          Get maximum matrix element                     *
 ***************************************************************************/
double GSparseMatrix::max() const
{
  // Initialise maximum with fill value
  double result = m_fill_val;

  // Search all elements for the largest one
  for (int i = 0; i < m_elements; ++i) {
    if (m_data[i] > result)
	  result = m_data[i];
  }

  // If maximum is < 0.0 and there are zero elements then set the maximum
  // to 0.0
  if ((result < 0.0) && (m_elements < (m_rows*m_cols)))
    result = 0.0;

  // Return result
  return result;
}


/***************************************************************************
 *                       Compute sum of non-zero elements                  *
 ***************************************************************************/
double GSparseMatrix::sum() const
{
  // Initialise matrix sum with fill value
  double result = m_fill_val;
  
  // Add all elements  
  for (int i = 0; i < m_elements; ++i)
    result += m_data[i];
  
  // Return result
  return result;
}


/***************************************************************************
 *                         Inplace transpose of matrix                     *
 * ----------------------------------------------------------------------- *
 * Call cs function                                                        *
 ***************************************************************************/
void GSparseMatrix::transpose()
{
  *this = cs_transpose(this, 1);
}


/***************************************************************************
 *                     Extract row from matrix into vector                 *
 ***************************************************************************/
GVector GSparseMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw out_of_range("GSparseMatrix::extract_row(int)", 
	                   row, 0, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_cols);
  
  // Loop over all columns to extract data
  for (int col = 0; col < m_cols; ++col) {
  
    // Get the start and stop of the elements
    int i_start = m_colstart[col];
    int i_stop  = m_colstart[col+1];
	
	// Search requested row in elements
	int i;
    for (i = i_start; i < i_stop; ++i) {
	  if (m_rowinx[i] == row)
	    break;
	}
	
	// Copy element if we found one
	if (i < i_stop)
	  result(col) = m_data[i];
	  
  } // endfor: looped over all columns
  
  // If there is a pending element then put it in the vector
  if (m_fill_val != 0.0 && m_fill_row == row)
    result(m_fill_col) = m_fill_val;
  
  // Return vector
  return result;

}


/***************************************************************************
 *                   Extract column from matrix into vector                *
 ***************************************************************************/
GVector GSparseMatrix::extract_col(int col) const
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw out_of_range("GSparseMatrix::extract_col(int)", 
	                   0, col, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_rows);

  // Get the start and stop of the elements
  int i_start = m_colstart[col];
  int i_stop  = m_colstart[col+1];
  
  // Extract elements into vector
  for (int i = i_start; i < i_stop; ++i)
    result(m_rowinx[i]) = m_data[i];

  // If there is a pending element then put it in the vector
  if (m_fill_val != 0.0 && m_fill_col == col)
    result(m_fill_row) = m_fill_val;
	
  // Return vector
  return result;

}


/***************************************************************************
 *                      Insert vector column into matrix                   *
 ***************************************************************************/
void GSparseMatrix::insert_col(const GVector& v, int col)
{
  // Debug header
  #if defined(G_DEBUG_SPARSE_INSERTION)
  cout << "GSparseMatrix::insert_col([" << v << "], " << col << "):" << endl;
  cout << " In Data : ";
  for (int i = 0; i < m_elements; ++i)
	cout << m_data[i] << " ";
  cout << endl << " In Row .: ";
  for (int i = 0; i < m_elements; ++i)
	cout << m_rowinx[i] << " ";
  cout << endl << " In Col .: ";
  for (int i = 0; i < m_cols+1; ++i)
	cout << m_colstart[i] << " ";
  cout << endl;
  #endif
  
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw out_of_range("GSparseMatrix::insert_col(const GVector&, int)", 
	                   0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw matrix_vector_mismatch("GSparseMatrix::insert_col(const GVector&, int)", 
	                             v.m_num, m_rows, m_cols);

  // If there is a pending element for this column then delete it since
  // the vector overwrites this element
  if (m_fill_val != 0.0 && m_fill_col == col) {
    #if defined(G_DEBUG_SPARSE_PENDING)
    cout << "GSparseMatrix::insert_col(const GVector, int): pending value " <<
	        m_fill_val << " for location (" << m_fill_row << "," << m_fill_col << 
		    ") became obsolete" << endl;
	#endif
    m_fill_val = 0.0;
	m_fill_row = 0;
	m_fill_col = 0;
  }

  // Determine the number of non-zero elements in the vector
  int n_vector = 0;
  for (int row = 0; row < m_rows; ++row) {
    if (v(row) != 0.0)
	  n_vector++;
  }
  
  // Get the start and stop indices of the actual column and compute
  // the number of exisiting elements in the column
  int i_start = m_colstart[col];
  int i_stop  = m_colstart[col+1];
  int n_exist = i_stop - i_start;
	
  // Compute the size difference for the new matrix. It is positive if
  // the number of non-zero entries in the vector is larger than the
  // number of non-zero entries in the matrix (in this case we have to
  // increase the matrix size).
  int n_diff = n_vector - n_exist;
	
  // If we need space then allocate it, if we have to much space the free it
  if (n_diff > 0) {
	alloc_elements(i_start, n_diff);
	#if defined(G_DEBUG_SPARSE_INSERTION)
    cout << " Insert .: " << n_diff << " elements at index " << i_start << endl;
    #endif
  }
  else if (n_diff < 0) {
	free_elements(i_start, -n_diff);
	#if defined(G_DEBUG_SPARSE_INSERTION)
    cout << " Remove .: " << -n_diff << " elements at index " << i_start << endl;
    #endif
  }
	  
  // Insert the vector elements in the matrix
  if (n_vector > 0) {
    for (int row = 0, i = i_start; row < m_rows; ++row) {
	  if (v(row) != 0.0) {
	    m_data[i]   = v(row);
	    m_rowinx[i] = row;
	    i++;
	  }
    }
  }
	
  // Update column start indices
  for (int i = col+1; i <= m_cols; ++i)
	m_colstart[i] += n_diff;

  // Debugging: show sparse matrix after insertion
  #if defined(G_DEBUG_SPARSE_INSERTION)
  cout << " Out Data: ";
  for (int i = 0; i < m_elements; ++i)
	cout << m_data[i] << " ";
  cout << endl << " Out Row : ";
  for (int i = 0; i < m_elements; ++i)
	cout << m_rowinx[i] << " ";
  cout << endl << " Out Col : ";
  for (int i = 0; i < m_cols+1; ++i)
	cout << m_colstart[i] << " ";
  cout << endl;
  #endif
  
  // Return
  return;
}


/***************************************************************************
 *                         Add vector to matrix column                     *
 ***************************************************************************/
void GSparseMatrix::add_col(const GVector& v, int col)
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw out_of_range("GSparseMatrix::add_col(const GVector&, int)", 
	                   0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are incompatible
  if (m_rows != v.m_num)
    throw matrix_vector_mismatch("GSparseMatrix::add_col(const GVector&, int)", 
	                             v.m_num, m_rows, m_cols);

  // Extract vector for column, add elements, and re-insert vector
  GVector column = extract_col(col);
  column += v;
  insert_col(column, col);
  
  // Return
  return;
}


/***************************************************************************
 *                       Convert matrix into full matrix                   *
 ***************************************************************************/
GMatrix GSparseMatrix::convert_to_full() const
{
  // Define result matrix
  GMatrix result(m_rows, m_cols);
  
  // Extract all elements
  for (int row = 0; row < m_rows; ++row) {
    for (int col = 0; col < m_cols; ++col)
	  result(row,col) = (*this)(row,col);
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                   GSparseMatrix Cholesky decomposition                  *
 * ----------------------------------------------------------------------- *
 * Cholesky decomposition of a sparse matrix. The decomposition is stored  *
 * within the GSparseMatrix without destroying the original matrix.        *
***************************************************************************/
void GSparseMatrix::cholesky_decompose(int compress)
{
  // Save original matrix size
  int matrix_rows = m_rows;
  int matrix_cols = m_cols;

  // Delete any existing symbolic analysis object
  if (m_symbolic != NULL) delete (GSparseSymbolic*)m_symbolic;

  // Allocate symbolic analysis object
  GSparseSymbolic* symbolic = new GSparseSymbolic();
  if (symbolic == NULL)
	throw mem_alloc("GSparseMatrix::cholesky_decompose(int)", 1);

  // Declare numeric analysis object. We don't allocate one since we'll
  // throw it away at the end of the function (the L matrix will be copied
  // in this object)
  GSparseNumeric numeric;
  
  // Fill pending element into matrix
  fill_pending();

  // Remove rows and columns containing only zeros if matrix compression
  // has been selected
  if (compress)
	remove_zero_row_col();
  
  // Ordering an symbolic analysis of matrix. This sets up an array 'pinv'
  // which contains the fill-in reducing permutations
  symbolic->cholesky_symbolic_analysis(1, *this);
  
  // Perform numeric Cholesky decomposition
  numeric.cholesky_numeric_analysis(*this, *symbolic);

  // Insert zero rows and columns if they have been removed previously.
  // Since no selection information is present in numeric.m_L we have to
  // copy it from this in order to allow for insertion.
  if (compress) {

    // If there is a row selection then copy it
    if (m_rowsel != NULL && m_num_rowsel > 0) {
      numeric.m_L->m_rowsel = new int[m_num_rowsel];
      if (numeric.m_L->m_rowsel == NULL)
	    throw mem_alloc("GSparseMatrix::cholesky_decompose(int)", m_num_rowsel);
      for (int i = 0; i < m_num_rowsel; ++i)
        numeric.m_L->m_rowsel[i] = m_rowsel[i];
	  numeric.m_L->m_num_rowsel = m_num_rowsel;
    }

    // If there is a column selection then copy it
    if (m_colsel != NULL && m_num_colsel > 0) {
      numeric.m_L->m_colsel = new int[m_num_colsel];
      if (numeric.m_L->m_colsel == NULL)
	    throw mem_alloc("GSparseMatrix::cholesky_decompose(int)", m_num_colsel);
      for (int i = 0; i < m_num_colsel; ++i)
        numeric.m_L->m_colsel[i] = m_colsel[i];
	  numeric.m_L->m_num_colsel = m_num_colsel;
    }

    // Now we can insert zero rows and columns
    numeric.m_L->insert_zero_row_col(matrix_rows, matrix_cols);
  }

  // Store Cholesky decomposition in matrix
  *this = *(numeric.m_L);

  // Store symbolic pointer in sparse matrix object
  m_symbolic = (void*)symbolic;

  // Return
  return;
}


/***************************************************************************
 *                      GSparseMatrix Cholesky solver                      *
 * ----------------------------------------------------------------------- *
 * Solves the linear equation A*x=b using a Cholesky decomposition of A.   *
 * This function is to be applied on a GSparseMatrix matrix for which a    *
 * Choleksy factorization has been produced using 'cholesky_decompose'.    *
 ***************************************************************************/
GVector GSparseMatrix::cholesky_solver(const GVector& v, int compress)
{
  // Declare loop variables
  int row;
  int col;
  int c_row;
  int c_col;

  // Dump header
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << "GSparseMatrix::cholesky_solver" << endl;
  cout << " Input vector .....: " << v << endl;
  #endif

  // Raise an exception if the matrix and vector dimensions are incompatible
  if (m_rows != v.m_num)
	throw matrix_vector_mismatch(
		  "GSparseMatrix::cholesky_solver(const GVector&, int)", 
		  v.m_num, m_rows, m_cols);
  
  // Raise an exception if there is no symbolic pointer
  if (!m_symbolic)
    throw matrix_not_factorised(
	      "GSparseMatrix::cholesky_solver(const GVector&, int)",
		  "Cholesky decomposition");

  // Get symbolic pointer
  GSparseSymbolic* symbolic = (GSparseSymbolic*)m_symbolic;
  
  // Raise an exception if there is no permutation
  if (!symbolic->m_pinv)
    throw matrix_not_factorised(
	      "GSparseMatrix::cholesky_solver(const GVector&, int)",
		  "Cholesky decomposition");

  // Flag row and column compression
  int row_compressed = (m_rowsel != NULL && m_num_rowsel < m_rows);
  int col_compressed = (m_colsel != NULL && m_num_colsel < m_cols);

  // Decide if we need a compression algorithm or not
  int no_zero = !(compress && (row_compressed || col_compressed));
  
  // Allocate vector for permutation and result vector
  GVector result(m_cols);

  // Setup pointers to L matrix and x vector
  int*    Lp = m_colstart;
  int*    Li = m_rowinx; 
  double* Lx = m_data;

  // Case A: no zero-row/col compression needed
  if (no_zero) {

    // Declare working vector
	GVector x(m_rows);

    // Perform inverse vector permutation
    for (int i = 0; i < v.m_num; ++i)
      x.m_data[symbolic->m_pinv[i]] = v.m_data[i];  

    // Inplace solve L\x=x
    for (col = 0; col < m_cols; ++col) {                // loop over columns
	  x.m_data[col] /= Lx[Lp[col]];                     // divide by diag. 
	  for (int p = Lp[col]+1; p < Lp[col+1]; p++)       // loop over elements
        x.m_data[Li[p]] -= Lx[p] * x.m_data[col];
    }
  
    // Inplace solve L'\x=x
    for (col = m_cols-1; col >= 0; --col) {             // loop over columns
	  for (int p = Lp[col]+1; p < Lp[col+1]; p++)       // loop over elements  
	      x.m_data[col] -= Lx[p] * x.m_data[Li[p]];
	  x.m_data[col] /= Lx[Lp[col]];
    }
	
	// Perform vector permutation
	for (int i = 0; i < m_cols; ++i)
      result.m_data[i] = x.m_data[symbolic->m_pinv[i]];
	  
  } // endif: Case A
  
  // Case B: zero-row/column compression requested
  else {
  
    // Allocate row and column mapping arrays
    int* row_map = new int[m_rows];
    if (row_map == NULL)
      throw mem_alloc("GSparseMatrix::cholesky_solver(const GVector&, int)", m_rows);
    int* col_map = new int[m_cols];
    if (col_map == NULL)
      throw mem_alloc("GSparseMatrix::cholesky_solver(const GVector&, int)", m_cols);

    // Setup row mapping array that maps original matrix rows into compressed
    // matrix rows. An entry of -1 indicates that the row should be dropped.
	// If no selection exists then setup an identity map.
	if (row_compressed) { 
      for (row = 0; row < m_rows; ++row)
        row_map[row] = -1;
      for (c_row = 0; c_row < m_num_rowsel; ++c_row)
        row_map[m_rowsel[c_row]] = c_row;
	}
	else {
      for (row = 0; row < m_rows; ++row)
        row_map[row] = row;
	}
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Row mapping ......:";
	for (row = 0; row < m_rows; ++row)
	  cout << " " << row_map[row];
	cout << endl;
    #endif

    // Setup column mapping array that maps original matrix column into compressed
    // matrix columns. An entry of -1 indicates that the column should be dropped
	// If no selection exists then setup an identity map.
	if (col_compressed) {
      for (col = 0; col < m_cols; ++col)
        col_map[col] = -1;
      for (c_col = 0; c_col < m_num_colsel; ++c_col)
        col_map[m_colsel[c_col]] = c_col;
	}
	else {
      for (col = 0; col < m_cols; ++col)
        col_map[col] = col;
	}
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Column mapping ...:";
	for (col = 0; col < m_cols; ++col)
	  cout << " " << col_map[col];
	cout << endl;
    #endif

    // Declare working vector
	GVector x(row_compressed ? m_num_rowsel : m_rows);
	
	// Compress input vector v -> c_v if required
    if (m_rowsel != NULL && m_num_rowsel < m_rows) {
      for (c_row = 0; c_row < m_num_rowsel; ++c_row)
	    x.m_data[c_row] = v.m_data[m_rowsel[c_row]];
	}
	else
	  x = v;
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Compressed vector : " << x << endl;
    #endif

    // Perform inverse permutation
    x = iperm(x, symbolic->m_pinv);
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Permutated vector : " << x << endl;
    #endif
   
    // Inplace solve L\x=x. The column and row maps are just use to see which 
	// columns or rows should be skipped in the calculations.
    for (col = 0; col < m_cols; ++col) {                // loop over columns
	  if ((c_col = col_map[col]) >= 0) {                // use only non-zero cols
	    x.m_data[c_col] /= Lx[Lp[col]];                 // divide by diag.
	    for (int p = Lp[col]+1; p < Lp[col+1]; p++) {   // loop over elements
		  if ((c_row = row_map[Li[p]]) >= 0)            // use only non-zero rows
            x.m_data[c_row] -= Lx[p] * x.m_data[c_col];
		}
	  }
    }
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Solve Lx=x .......: " << x << endl;
    #endif
  
    // Inplace solve L'\x=x. The column and row maps are just use to see which 
	// columns or rows should be skipped in the calculations.
    for (col = m_cols-1; col >= 0; --col) {             // loop over columns
	  if ((c_col = col_map[col]) >= 0) {                // use only non-zero cols
	    for (int p = Lp[col]+1; p < Lp[col+1]; p++) {   // loop over elements
		  if ((c_row = row_map[Li[p]]) >= 0)            // use only non-zero rows
	        x.m_data[c_col] -= Lx[p] * x.m_data[c_row];
		}
	    x.m_data[c_col] /= Lx[Lp[col]];
	  }
    }
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Solve L'x=x ......: " << x << endl;
    #endif

    // Perform vector permutation
    x = perm(x, symbolic->m_pinv);
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Permutated vector : " << x << endl;
    #endif
  
    // If column compression has been performed the expand the result vector
    // accordingly
    if (m_colsel != NULL && m_num_colsel < m_cols) {
	  for (c_col = 0; c_col < m_num_colsel; ++c_col)
	    result.m_data[m_colsel[c_col]] = x.m_data[c_col];
	}
	else
	  result = x;
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    cout << " Restored vector ..: " << result << endl;
    #endif

  // Free mapping arrays
  delete [] row_map;
  delete [] col_map;
	
  } // endelse: Case B

  // Return result vector
  return result;
}


/***************************************************************************
 *                       GSparseMatrix Cholesky invert                     *
 * ----------------------------------------------------------------------- *
 * Inplace matrix inversion using the Cholesky decomposition.              *
 ***************************************************************************/
void GSparseMatrix::cholesky_invert(int compress)
{
  // Generate Cholesky decomposition of matrix
  cholesky_decompose(compress);

  // Allocate result matrix and unit vector
  GSparseMatrix result(m_rows, m_cols);
  GVector       unit(m_rows);

  // Column-wise solving of the problem
  for (int col = 0; col < m_cols; ++col) {
  
    // Set unit vector
	unit.m_data[col] = 1.0;
	
	// Solve for column
	GVector x = cholesky_solver(unit, compress);
	 
	// Insert column in matrix
    result.insert_col(x, col);

    // Clear unit vector for next round
	unit.m_data[col] = 0.0;
    
  }
  
  *this = result;
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                     GSparseMatrix private functions                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                   Determines element index for (row,col)                *
 * ----------------------------------------------------------------------- *
 * Returns the index in the compressed array for (row,col). The following  *
 * special results exist:                                                  *
 *         -1: Requested index does not exist in the matrix                *
 * m_elements: Requested index is the pending element                      *
 ***************************************************************************/
int GSparseMatrix::get_index(int row, int col) const
{
  // Initialise element to 'not found'
  int index = -1;
  
  // If we have a pending element then check if this element is requested
  if (m_fill_val != 0.0) {
    if (row == m_fill_row && col == m_fill_col)
	  return m_elements;
  }
  
  // If requested element is not the pending element then check if it exists
  // in the matrix. Only if it is found its index is returned. Otherwise
  // the default index is -1, signalling that the element is absent.
  if (m_elements > 0) {
    int* ptr_colstart = m_colstart + col;
    int  i_start      = *ptr_colstart++;
    int  i_stop       = *ptr_colstart;
	int* ptr_rowinx   = m_rowinx + i_start;
    for (int i = i_start; i < i_stop; ++i) {
	  int row_test = *ptr_rowinx++;
	  if (row_test == row) {
	    index = i;
		break;
	  }
	}
  }
  
  // Return index
  return index;
}


/***************************************************************************
 *                       Fills pending matrix element                      *
 * ----------------------------------------------------------------------- *
 * If 'm_fill_val' is non-zero a pending matrix element exists that should *
 * be filled into (row,col)=(m_fill_row,m_fill_col). This routine performs *
 * the filling of the matrix with this element and resets 'm_fill_val' to  *
 * zero. This routine allows for element-by-element filling of a sparse    *
 * matrix. This is, of course, very time consuming and should in general   *
 * be avoided. However, it allows to design a sparse matrix class that     *
 * hides the matrix sparsity completely to the user.                       *
 ***************************************************************************/
void GSparseMatrix::fill_pending(void)
{
  // If we have a pending element then fill it into the matrix
  if (m_fill_val != 0.0) {

    // Debugging
    #if defined(G_DEBUG_SPARSE_PENDING) || defined(G_DEBUG_SPARSE_INSERTION)
    cout << "GSparseMatrix::fill_pending(): pending value " << 
			m_fill_val << " will be filled in location (" << 
			m_fill_row << "," << m_fill_col <<  ")" << endl;
	#endif

    // If there are so far no elements in the matrix then append element ...
	int inx_insert;
    if (m_elements == 0)
	  inx_insert = 0;
	  
	// ... otherwise search for index to insert
	else {
      int* ptr_colstart = m_colstart + m_fill_col;
      int  i_start      = *ptr_colstart++;
      int  i_stop       = *ptr_colstart;
	  int* ptr_rowinx   = m_rowinx + i_start;
      for (inx_insert = i_start; inx_insert < i_stop; ++inx_insert) {
	    int row_test = *ptr_rowinx++;
	    if (row_test > m_fill_row) {
		  break;
	    }
	  }
	}
	
	// Allocate memory for new element
	alloc_elements(inx_insert, 1);

	// Insert element
	m_data[inx_insert]   = m_fill_val;
	m_rowinx[inx_insert] = m_fill_row;
	
	// Update column start indices
	for (int col = m_fill_col+1; col <= m_cols; ++col)
	  m_colstart[col] += 1;
		
    // Reset fill value
	m_fill_val = 0.0;
	m_fill_row = 0;
	m_fill_col = 0;

    // Debugging: show sparse matrix after filling
    #if defined(G_DEBUG_SPARSE_PENDING) || defined(G_DEBUG_SPARSE_INSERTION)
    cout << " Data: ";
    for (int i = 0; i < m_elements; ++i)
      cout << m_data[i] << " ";
    cout << endl << " Row.: ";
    for (int i = 0; i < m_elements; ++i)
      cout << m_rowinx[i] << " ";
    cout << endl << " Col.: ";
    for (int i = 0; i < m_cols+1; ++i)
      cout << m_colstart[i] << " ";
    cout << endl;
	#endif

  } // endif: a pending matrix element was found
  
}


/***************************************************************************
 *                  Allocate memory for new matrix elements                *
 * ----------------------------------------------------------------------- *
 * Inserts a memory allocation for 'num' elements at the index 'start'.    *
 * The new elements are filled with 0.0 and the corresponding row indices  *
 * are set to 0.                                                           *
 * NOTE: This routine does not take care of updating 'm_colstart'. This    *
 * has to be done by the client.                                           *
 ***************************************************************************/
void GSparseMatrix::alloc_elements(int start, int num)
{
  // Continue only if we need memory
  if (num > 0) {
  
    // If start is after the end the append memory
	if (start > m_elements)
	  start = m_elements;
  
	// Allocate memory for new elements
	double* new_data   = new double[m_elements+num];
	int*    new_rowinx = new int[m_elements+num];
	if (new_data == NULL || new_rowinx == NULL)
	  throw mem_alloc("GSparseMatrix::alloc_elements(int, int)", 
	                  m_elements+num);

    // Copy all elements before index to insert
	for (int i = 0; i < start; ++i) {
	  new_data[i]   = m_data[i];
	  new_rowinx[i] = m_rowinx[i];
	}

	// Clear new elements (zero content for row 0)
	for (int i = start; i < start+num; ++i) {
	  new_data[i]   = 0.0;
	  new_rowinx[i] = 0;
	}
	
	// Copy all elements after index to insert
	for (int i = start; i < m_elements; ++i) {
	  new_data[i+num]   = m_data[i];
	  new_rowinx[i+num] = m_rowinx[i];
	}
	
    // Delete old memory
    if (m_data   != NULL) delete [] m_data;
    if (m_rowinx != NULL) delete [] m_rowinx;
	
	// Update pointers to new memory and update element counter
	m_data      = new_data;
	m_rowinx    = new_rowinx;
	m_elements += num;
	
  } // endif: needed new memory

}


/***************************************************************************
 *                  Free memory for obsolete matrix elements               *
 * ----------------------------------------------------------------------- *
 * Free memory for 'num' elements starting from index 'start'.             *
 * NOTE: This routine does not take care of updating 'm_colstart'. This    *
 * has to be done by the client.                                           *
 ***************************************************************************/
void GSparseMatrix::free_elements(int start, int num)
{
  // Continue only if we need to free memory and if start is within the
  // range
  if (num > 0 && start < m_elements) {

    // Determine the new size of the matrix. If there are no elements then
	// simply delete all matrix elements ...
	int new_size = m_elements - num;
	if (new_size < 1) {
      if (m_data   != NULL) delete [] m_data;
      if (m_rowinx != NULL) delete [] m_rowinx;
	  m_data     = NULL;
	  m_rowinx   = NULL;
	  m_elements = 0;
	}
	
	// ... otherwise shrink the array
	else {
	
	  // Allocate new memory
	  double* new_data   = new double[new_size];
	  int*    new_rowinx = new int[new_size];
	  if (new_data == NULL || new_rowinx == NULL)
	    throw mem_alloc("GSparseMatrix::free_elements(int, int)", 
	                    new_size);

	  // Copy all elements before the starting index
	  for (int i = 0; i < start; ++i) {
	    new_data[i]   = m_data[i];
	    new_rowinx[i] = m_rowinx[i];
	  }

	  // Copy all elements after the starting index
	  for (int i = start; i < new_size; ++i) {
	    new_data[i]   = m_data[i+num];
	    new_rowinx[i] = m_rowinx[i+num];
	  }
	
      // Delete old memory
      if (m_data   != NULL) delete [] m_data;
      if (m_rowinx != NULL) delete [] m_rowinx;
	
	  // Update pointers to new memory and update element counter
	  m_data     = new_data;
	  m_rowinx   = new_rowinx;
	  m_elements = new_size;

	} // endif: array shrinkage needed	
  } // endif: needed new memory

}


/***************************************************************************
 *                           remove_zero_row_col                           *
 * ----------------------------------------------------------------------- *
 * Remove all rows and columns that contain only zeros from matrix. This   *
 * function is needed for compressed matrix factorisation. The resulting   *
 * matrix has reduced size (# of rows and columns) and dimensions (# of    *
 * elements). Note that the physical memory is not reduced by the routine. *
 ***************************************************************************/
void GSparseMatrix::remove_zero_row_col(void)
{
  // Dump header
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << "GSparseMatrix::remove_zero_row_col" << endl;
  #endif

  // Fill pending value in matrix
  fill_pending();

  // Select non-zero rows and columns of matrix
  select_non_zero();
  
  // Stop if there is no compression
  if (m_rows == m_num_rowsel && m_cols == m_num_colsel)
    return;

  // Raise exception if all matrix elements are zero
  if (m_num_rowsel < 1 || m_num_colsel < 1)
    throw matrix_zero("GSparseMatrix::remove_zero_row_col()");

  // Allocate row mapping array
  int* row_map = new int[m_rows];
  if (row_map == NULL)
    throw mem_alloc("GSparseMatrix::remove_zero_row_col()", m_rows);
	
  // Setup row mapping array that maps original matrix rows into compressed
  // matrix rows. An entry of -1 indicates that the row should be dropped
  for (int row = 0; row < m_rows; ++row)
    row_map[row] = -1;
  for (int c_row = 0; c_row < m_num_rowsel; ++c_row)
    row_map[m_rowsel[c_row]] = c_row;
  
  // Initialise pointers to compressed array
  double* d_data     = m_data;
  int*    d_rowinx   = m_rowinx;
	
  // Initialise column start of first column to zero
  m_colstart[0] = 0;

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " before compression:";
  for (int col = 0; col <= m_cols; ++col)
    cout << " " << m_colstart[col];
  cout << endl;
  #endif

  // Loop over all columns of compressed matrix
  for (int c_col = 0; c_col < m_num_colsel; ++c_col) {
	
    // Get index range for elements in original matrix 
	int i_start = m_colstart[m_colsel[c_col]];
	int i_stop  = m_colstart[m_colsel[c_col]+1];
	  
	// Initialise number of element counter for the compressed column
	int num = 0;
	  
	// Loop over all elements in original column
	int c_row;
	for (int i = i_start; i < i_stop; ++i) {
	  if ((c_row = row_map[m_rowinx[i]]) >= 0) {
		*d_data++   = m_data[i];
		*d_rowinx++ = c_row;
		num++;
	  }
	}
	  
	// Update column stop
	m_colstart[c_col+1] = m_colstart[c_col] + num;
	  
  } // endfor: looped over all columns of compressed matrix

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " after compression :";
  for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
    cout << " " << m_colstart[c_col];
  cout << endl;
  #endif

  // Free row mapping array
  delete [] row_map;
	
  // Update matrix attributes
  m_rows     = m_num_rowsel;
  m_cols     = m_num_colsel;
  m_elements = m_colstart[m_num_colsel];
	 
  // Return
  return;
}


/***************************************************************************
 *                           insert_zero_row_col                           *
 * ----------------------------------------------------------------------- *
 * Insert zero rows and columns into matrix. Since for a sparse matrix     *
 * this does not require any allocation of additional memory, the data are *
 * not moved by this function, but the pointers are re-arranged.           *
 ***************************************************************************/
void GSparseMatrix::insert_zero_row_col(int rows, int cols)
{
  // Dump header
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << "GSparseMatrix::insert_zero_row_col(" << rows << "," << cols << 
          ")" << endl;
  #endif

  // Fill pending value in matrix
  fill_pending();

  // Stop if there is no insertion
  if (m_rows == rows && m_cols == cols)
    return;
  
  // If row selection exists then restore row indices
  if (m_rowsel != NULL) {
    for (int i = 0; i < m_elements; ++i)
	  m_rowinx[i] = m_rowsel[m_rowinx[i]];
  }

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " before restoration:";
  for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
    cout << " " << m_colstart[c_col];
  cout << endl;
  #endif
  
  // If column selection exists then restore column counters
  if (m_colsel != NULL) {
    int col_stop = cols - 1;
    for (int c_col = m_num_colsel-1; c_col > 0; --c_col) {
	  int col_start = m_colsel[c_col-1] + 1;
	  int col_value = m_colstart[c_col];
	  for (int col = col_start; col <= col_stop; ++col)
	    m_colstart[col] = col_value;
	  col_stop = col_start - 1;
	}
	m_colstart[cols] = m_elements;
  }

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " after restoration :";
  for (int col = 0; col <= cols; ++col)
    cout << " " << m_colstart[col];
  cout << endl;
  #endif  

  // Update matrix attributes
  m_rows = rows;
  m_cols = cols;
	
  // Return
  return;
}


/***************************************************************************
 *                              cs_transpose                               *
 * ----------------------------------------------------------------------- *
 * Transpose matrix. The flag 'values' allows to avoid copying the actual  *
 * data values. This allows to perform a logical matrix transposition, as  *
 * needed by the symbolic matrix analysis class.                           *
 ***************************************************************************/
GSparseMatrix GSparseMatrix::cs_transpose(GSparseMatrix* m, int values)
{
  // Fill pending value in input matrix
  m->fill_pending();

  // Declare and allocate result matrix 
  GSparseMatrix result(m->m_cols, m->m_rows, m->m_elements);

  // Allocate and initialise workspace
  int  wrk_size = m->m_rows;
  int* wrk_int  = new int[wrk_size];
  if (wrk_int == NULL)
	throw mem_alloc("GSparseMatrix::cs_transpose(GSparseMatrix*, int)", 
	                wrk_size);
  for (int i = 0; i < wrk_size; i++) 
    wrk_int[i] = 0;

  // Setup the number of non-zero elements in each row
  for (int p = 0; p < m->m_elements; p++) 
    wrk_int[m->m_rowinx[p]]++;
  
  // Set row pointers. To use a GSparseSymbolic function we have to
  // allocate and object (but this does not take memory)
  cs_cumsum(result.m_colstart, wrk_int, m->m_rows);

  // Case A: Normal transponse, including assignment of values
  if (values) {
    for (int col = 0; col < m->m_cols; col++) {
	  for (int p = m->m_colstart[col] ; p < m->m_colstart[col+1] ; p++) {
	    int i         = wrk_int[m->m_rowinx[p]]++;
	    result.m_rowinx[i] = col;
		result.m_data[i]   = m->m_data[p] ;
	  }
    }
  }
  
  // Case B: Logical transponse, no assignment of values is performed
  else {
    for (int col = 0; col < m->m_cols; col++) {
	  for (int p = m->m_colstart[col] ; p < m->m_colstart[col+1] ; p++)
	    result.m_rowinx[wrk_int[m->m_rowinx[p]]++] = col;
    }
  }
  
  // Return transponse matrix
  return result;
	
}


/*==========================================================================
 =                                                                         =
 =                           GSparseMatrix friends                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GSparseMatrix& m)
{
  // Put header in stream
  os << "=== GSparseMatrix ===" << endl;
  if (m.m_rowsel != NULL)
    os << " Number of rows ............: " << m.m_rows << " (compressed " <<
	      m.m_num_rowsel << ")" << endl;
  else
    os << " Number of rows ............: " << m.m_rows << endl;
  if (m.m_colsel != NULL)
    os << " Number of columns .........: " << m.m_cols << " (compressed " <<
	      m.m_num_colsel << ")" << endl;
  else
    os << " Number of columns .........: " << m.m_cols << endl;
  os << " Number of non-zero elements: " << m.m_elements << endl;
  os << " Sparse matrix fill ........: " << m.fill() << endl;

  // Loop over all matrix elements and put them into the stream
  for (int row = 0; row < m.m_rows; ++row) {
    os << " ";
    for (int col = 0; col < m.m_cols; ++col) {
      os << m(row,col);
	  if (col != m.m_cols-1)
	    os << ", ";
	}
	if (row != m.m_rows-1)
	  os << endl;
  }
  
  // If there is a row compression the show scheme
  if (m.m_rowsel != NULL) {
    os << endl << " Row selection ..:";
    for (int row = 0; row < m.m_num_rowsel; ++row)
      os << " " << m.m_rowsel[row];
  }

  // If there is a column compression the show scheme
  if (m.m_colsel != NULL) {
    os << endl << " Column selection:";
    for (int col = 0; col < m.m_num_colsel; ++col)
      os << " " << m.m_colsel[col];
  }
  
  // If there is a symbolic decomposition then put it also in the stream
  if (m.m_symbolic != NULL) {
    os << endl << *((GSparseSymbolic*)m.m_symbolic);
  }

  // If there is a numeric decomposition then put it also in the stream
  if (m.m_numeric != NULL) {
    os << endl << *((GSparseNumeric*)m.m_numeric);
  }

  // Return output stream
  return os;
}


/***************************************************************************
 *                          Matrix absolute values                         *
 ***************************************************************************/
GSparseMatrix fabs(const GSparseMatrix& m)
{
  // Define result matrix
  GSparseMatrix result = m;
  
  // Fill pending element
  result.fill_pending();
  
  // Convert all elements to absolute values  
  for (int i = 0; i < result.m_elements; ++i)
    result.m_data[i] = fabs(result.m_data[i]);
  
  // Return result
  return result;
}


/*==========================================================================
 =                                                                         =
 =                      GSparseMatrix exception classes                    =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GSparseMatrix                 =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                  cs_cumsum                              *
 * ----------------------------------------------------------------------- *
 * Evaluate p[0..n] = cumulative sum of c[0..n-1]                          *
 * ----------------------------------------------------------------------- *
 * Input:   c[0..n-1]      integer array of n elements                     *
 * Output:  nz2            cumulative sum of c[0..n-1]                     *
 *          p[0..n]        integer array of n+1 elements                   *
 *          c[0..n-1]      holds p[0..n-1] on output                       *
 ***************************************************************************/
double cs_cumsum(int* p, int* c, int n)
{
  // Signal error if one of the input pointers is NULL
  if (!p || !c) return (-1);
	
  // Initialise sums (integer and double)
  int    nz  = 0;
  double nz2 = 0.0;
  
  // Evaluate p[0..n] = cumulative sum of c[0..n-1]
  for (int i = 0; i < n; ++i) {
	p[i]  = nz ;
	nz   += c[i];
	nz2  += c[i];		    // also in double to avoid int overflow
	c[i]  = p[i];		    // also copy p[0..n-1] back into c[0..n-1]
  }
  p[n] = nz ;
  
  // Return cumulative sum of c[0..n-1]
  return nz2;
}

