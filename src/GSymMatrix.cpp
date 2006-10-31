/***************************************************************************
 *                  GSymMatrix.cpp  -  symmetric matrix class              *
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
#include "GSymMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/*==========================================================================
 =                                                                         =
 =                   GSymMatrix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         GSymMatrix constructor                          *
 ***************************************************************************/
GSymMatrix::GSymMatrix(int rows, int cols) : GMatrix()
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

  // Throw exception if number of rows and columns is not identical
  if (rows != cols)
    throw matrix_not_symmetric("GSymMatrix constructor", rows, cols);

  // Determine number of physical elements in matrix
  int elements = rows*(rows+1)/2;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw empty("GSymMatrix constructor");
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new int[cols+1];
  m_inx      = new int[cols];
  if (m_data == NULL || m_colstart == NULL || m_inx == NULL)
	throw mem_alloc("GSymMatrix constructor", elements);

  // Store matrix size (logical and physical)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;

  // Set-up column start indices
  m_colstart[0]   = 0;
  int offset = rows;
  for (int col = 1; col <= m_cols; ++col)
    m_colstart[col] = m_colstart[col-1] + offset--;

  // Initialise matrix elements to 0.0
  for (int i = 0; i < m_elements; ++i)
    m_data[i] = 0.0;

  // Return
  return;
}


/***************************************************************************
 *                        GSymMatrix copy constructor                      *
 * ----------------------------------------------------------------------- *
 * First invoques GMatrix copy constructor, then this function.            *
 ***************************************************************************/
GSymMatrix::GSymMatrix(const GSymMatrix& m) : GMatrix(m)
{ 
  // Copy members
  m_num_inx = m.m_num_inx;
  m_inx     = new int[m_cols];
  if (m_inx == NULL)
	throw mem_alloc("GSymMatrix copy constructor", m_cols);
  for (int i = 0; i < m_cols; ++i)
	m_inx[i] = m.m_inx[i];

  // Return
  return;
}


/***************************************************************************
 *                         GSymMatrix destructor                           *
 * ----------------------------------------------------------------------- *
 * First invoques this function, then ~GMatrix destructor.                 *
 ***************************************************************************/
GSymMatrix::~GSymMatrix()
{
  // De-allocate only if memory has indeed been allocated by derived class
  if (m_inx != NULL) delete [] m_inx;

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                           GSymMatrix operators                          =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                      GSymMatrix assignment operator                     *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator= (const GSymMatrix& m)
{ 
  // Execute only if object is not identical
  if (this != &m) { 

    // Invoque base class assignment operator
    this->GMatrix::operator=(m);

    // De-allocate memory if it has indeed been allocated
    if (m_inx != NULL) delete [] m_inx;

    // Copy data members
    m_num_inx = m.m_num_inx;

    // Allocate memory for indices and copy them
    m_inx = new int[m_cols];
    if (m_inx == NULL)
	  throw mem_alloc("GSymMatrix assignment operator", m_cols);
    for (int i = 0; i < m_cols; ++i)
	  m_inx[i] = m.m_inx[i];
  }

  // Return this object
  return *this;
}


/***************************************************************************
 *                           Vector multiplication                         *
 ***************************************************************************/
GVector GSymMatrix::operator* (const GVector& v) const
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_cols != v.m_num)
    throw matrix_vector_mismatch("GSymMatrix::operator* (const GVector&)", 
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
 *                GSymMatrix negation (unary minus operator)               *
 ***************************************************************************/
GSymMatrix GSymMatrix::operator- ( ) const
{
  // Copy argument into result
  GSymMatrix result = *this;
  
  // Negate result
  for (int i = 0; i < m_elements; ++i)
    result.m_data[i] = -result.m_data[i];
	
  // Return result	
  return result;
}


/***************************************************************************
 *          Matrix multiplication (unary multiplication operator)          *
 * ----------------------------------------------------------------------- *
 * Symmetric matrices must have the same dimensions to be multiplied.      *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator*= (const GSymMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows || m_rows != m.m_cols)
    throw matrix_mismatch("GSymMatrix::operator*= (const GSymMatrix&)", 
						  m_rows, m_cols, m.m_rows, m.m_cols);

  // Allocate result matrix
  GSymMatrix result(m_rows, m.m_cols);
	
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

  // Return result
  return *this;
}


/*==========================================================================
 =                                                                         =
 =                        GSymMatrix member functions                      =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Sum matrix elements                          *
 ***************************************************************************/
double GSymMatrix::sum() const
{
  // Initialise matrix sums (diagonal and off-diagonal)
  double diag     = 0.0;
  double off_diag = 0.0;

  // Calulate sum over diagonal elements
  for (int row = 0; row < m_rows; ++row)
    diag += m_data[m_colstart[row]];

  // Calulate sum over off-diagonal elements
  for (int row = 0; row < m_rows; ++row) {
    for (int col = row+1; col < m_cols; ++col)
      off_diag += m_data[m_colstart[row]+(col-row)];
  }

  // Calculate total
  double result = diag + 2.0 * off_diag;
  
  // Return result
  return result;
}


/***************************************************************************
 *                     Extract row from matrix into vector                 *
 ***************************************************************************/
GVector GSymMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw out_of_range("GSymMatrix::extract_row(int)", 
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
GVector GSymMatrix::extract_col(int col) const
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw out_of_range("GSymMatrix::extract_col(int)", 
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


/***************************************************************************
 *                 Convert symmetric matrix into full matrix               *
 ***************************************************************************/
GMatrix GSymMatrix::convert_to_full() const
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
 *        Extract lower triangle of symmetric matrix into full matrix      *
 ***************************************************************************/
GMatrix GSymMatrix::extract_lower_triangle() const
{
  // Define result matrix
  GMatrix result(m_rows, m_cols);

  // Extract all elements
  for (int row = 0; row < m_rows; ++row) {
    for (int col = 0; col <= row; ++col)
	  result(row,col) = m_data[m_colstart[col]+(row-col)];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *        Extract upper triangle of symmetric matrix into full matrix      *
 ***************************************************************************/
GMatrix GSymMatrix::extract_upper_triangle() const
{
  // Define result matrix
  GMatrix result(m_rows, m_cols);

  // Extract all elements
  for (int row = 0; row < m_rows; ++row) {
    for (int col = row; col < m_cols; ++col)
	  result(row,col) = m_data[m_colstart[row]+(col-row)];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                     GSymMatrix Cholesky decomposition                   *
 * ----------------------------------------------------------------------- *
 * Inplace Cholesky decomposition inspired by Numerical Recipes algorithm. *
 * The decomposition, which is a matrix occupying only the lower triange,  *
 * is stored in the elements of the symmetric matrix. To visualise the     *
 * matrix one has to use 'lower_triangle()' to extract the relevant part.  *
 * Case A operates on a full matrix, Case B operates on a (logically)      *
 * compressed matrix where zero rows/columns have been removed.            *
 ***************************************************************************/
void GSymMatrix::cholesky_decompose(int compress)
{
  // Set-up incides of non zero rows if matrix compression is requested
  if (compress)
    set_inx();

  // Check if zero-row/col compression is needed  
  int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);
  
  // Case A: no zero-row/col compression needed
  if (no_zeros) {
  
    // Loop over upper triangle (col >= row)
    for (int row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
      for (int col = row; col < m_cols; ++col, ++ptr) {
	    double sum = *ptr;                                // sum = M(row,col)
	    for (int k = 0; k < row; ++k) {
	      int offset = m_colstart[k] - k;            // is always positive
	      sum -= m_data[offset+row] * m_data[offset+col]; // sum -= M(row,k)*M(col,k)
	    }
	    double diag;
	    if (row == col) {
	      if (sum <= 0.0)
		    throw matrix_not_pos_definite("GSymMatrix::cholesky_decompose(int)", 
			                              row, sum);
          *ptr = sqrt(sum);                               // M(row,row) = sqrt(sum)
		  diag = 1.0/(*ptr);
	    }
	    else
          *ptr = sum*diag;                                // M(row,col) = sum/M(row,row)
	  }
    }
  } // endif: there were no zero rows/cols in matrix
  
  // Case B: zero-row/col compression needed
  else if (m_num_inx > 0) {

    // Allocate loop variables and pointers
	int  row;
	int  col;
	int  k;
	int* row_ptr;
	int* col_ptr;
	int* k_ptr;

    // Loop over upper triangle (col >= row)
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
	  double* ptr_0 = m_data + m_colstart[*row_ptr] - *row_ptr;
      for (col = row, col_ptr = m_inx + row; col < m_num_inx; ++col, ++col_ptr) {
		double* ptr = ptr_0 + *col_ptr;
	    double  sum = *ptr;                                         // sum = M(row,col)
	    for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr) {
	      int offset = m_colstart[*k_ptr] - *k_ptr;            // is always positive
	      sum -= m_data[offset+*row_ptr] * m_data[offset+*col_ptr]; // sum -= M(row,k)*M(col,k)
	    }
	    double diag;
	    if (*row_ptr == *col_ptr) {
	      if (sum <= 0.0)
		    throw matrix_not_pos_definite("GSymMatrix::cholesky_decompose(int)", 
			                              *row_ptr, sum);
          *ptr = sqrt(sum);                                         // M(row,row) = sqrt(sum)
		  diag = 1.0/(*ptr);
	    }
	    else
          *ptr = sum*diag;                                          // M(row,col) = sum/M(row,row)
	  }
    }
  } // endelse: zero-row/col compression needed
  
  // Case C: all matrix elements are zero
  else
    throw matrix_zero("GSymMatrix::cholesky_decompose(int)");
	
}


/***************************************************************************
 *                       GSymMatrix Cholesky solver                        *
 * ----------------------------------------------------------------------- *
 * Solves the linear equation A*x=b using a Cholesky decomposition of A.   *
 * This function is to be applied on a decomposition GSymMatrix matrix     *
 * that is produced by 'cholesky_decompose'. Case A operates on a full     *
 * matrix, Case B on a zero rows/columns (logically) compressed matrix.    *
 ***************************************************************************/
GVector GSymMatrix::cholesky_solver(const GVector& v, int compress)
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw matrix_vector_mismatch(
	      "GSymMatrix::cholesky_solver(const GVector&, int)", 
	      v.m_num, m_rows, m_cols);

  // Allocate result vector
  GVector x(m_rows);

  // Check if zero-row/col compression is needed  
  int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

  // Case A: no zero-row/col compression needed
  if (no_zeros) {
    
    // Solve L*y=b, storing y in x (row>k)
    for (int row = 0; row < m_rows; ++row) {
      double sum = v(row);
	  for (int k = 0; k < row; ++k)
	    sum -= m_data[m_colstart[k]+(row-k)] * x(k); // sum -= M(row,k) * x(k)
	  x(row) = sum/m_data[m_colstart[row]];          // x(row) = sum/M(row,row) 
    }
	  
    // Solve trans(L)*x=y (k>row)
    for (int row = m_rows-1; row >= 0; --row) {
      double  sum = x(row);
	  double* ptr = m_data + m_colstart[row] + 1;
	  for (int k = row+1; k < m_rows; ++k)
	    sum -= *ptr++ * x(k);                        // sum -= M(k,row) * x(k)
	  x(row) = sum/m_data[m_colstart[row]];          // x(row) = sum/M(row,row) 
    }
  } // endif: no zero-row/col compression needed
  
  // Case B: zero-row/col compression needed
  else if (m_num_inx > 0) {

    // Allocate loop variables and pointers
	int       row;
	int  k;
	int* row_ptr;
	int* k_ptr;

    // Solve L*y=b, storing y in x (row>k)
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
      double  sum = v(*row_ptr);
	  double* ptr = m_data + *row_ptr;
	  for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr)
	    sum -= *(ptr + m_colstart[*k_ptr] - *k_ptr) * x(*k_ptr); // sum -= M(row,k) * x(k)
	  x(*row_ptr) = sum/m_data[m_colstart[*row_ptr]];            // x(row) = sum/M(row,row) 
    }
	  
    // Solve trans(L)*x=y (k>row)
    for (row = m_num_inx-1, row_ptr = m_inx+m_num_inx-1; row >= 0; --row, --row_ptr) {
      double  sum      = x(*row_ptr);
	  double* ptr_diag = m_data + m_colstart[*row_ptr];
	  double* ptr      = ptr_diag - *row_ptr;
	  for (k = row+1, k_ptr = m_inx+row+1; k < m_num_inx; ++k, ++k_ptr)
	    sum -= *(ptr + *k_ptr) * x(*k_ptr);                     // sum -= M(k,row) * x(k)
	  x(*row_ptr) = sum/(*ptr_diag);                            // x(row) = sum/M(row,row) 
    }
  } // endelse: zero-row/col compression needed
  
  // Case C: all matrix elements are zero
  else
    throw matrix_zero("GSymMatrix::cholesky_solver(const GVector&, int)");

  // Return result vector
  return x;
}


/***************************************************************************
 *                        GSymMatrix Cholesky invert                       *
 * ----------------------------------------------------------------------- *
 * Inplace matrix inversion using the Cholesky decomposition. Case A       *
 * operates on a full matrix while Case B operates on a (logically)        *
 * compressed matrix where all zero rows/columns are skipped.              *
 ***************************************************************************/
void GSymMatrix::cholesky_invert(int compress)
{
  // Generate Cholesky decomposition of matrix
  this->cholesky_decompose(compress);

  // Check if zero-row/col compression is needed  
  int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

  // Case A: no zero-row/col compression needed
  if (no_zeros) {
  
    // Generate inverse of Cholesky decomposition (col>row)
    for (int row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
	  *ptr        = 1.0/(*ptr);                       // M(row,row) = 1/M(row,row)
	  for (int col = row+1; col < m_cols; ++col) {
	    double   sum = 0.0;
	    double* ptr1 = m_data + col - row;
	    double* ptr2 = ptr;
	    for (int k = row; k < col; ++k)
	      sum -= *(ptr1-- + m_colstart[k]) * *ptr2++; // sum -= M(col,k)*M(k,row)
	    *(ptr+col-row) = sum/m_data[m_colstart[col]]; // M(col,row) = sum/M(col,col)
	  }
    }
  
    // Matrix multiplication (col>=row)
    for (int row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
	  for (int col = row; col < m_cols; ++col) {
	    double   sum = 0.0;
	    double* ptr1 = ptr + col - row;
	    double* ptr2 = m_data + m_colstart[col];
	    for (int k = col; k < m_cols; ++k)
	      sum += *ptr1++ * *ptr2++;                   // sum += M(row,k)*M(k,col)
	    *(ptr+col-row) = sum;                         // M(row,col) = sum
	  }
    }
  } // endif: no zero-row/col compression needed
  
  // Case B: zero-row/col compression needed
  else if (m_num_inx > 0) {

    // Allocate loop variables and pointers
	int  row;
	int  col;
	int  k;
	int* row_ptr;
	int* col_ptr;
	int* k_ptr;

    // Generate inverse of Cholesky decomposition (col>row)
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
	  double* ptr_diag = m_data + m_colstart[*row_ptr];
	  double* ptr_2    = ptr_diag - *row_ptr;
	  *ptr_diag = 1.0/(*ptr_diag);                              // M(row,row) = 1/M(row,row)
	  for (col = row+1, col_ptr = m_inx+row+1; col < m_num_inx; ++col, ++col_ptr) {
	    double  sum   = 0.0;
	    double* ptr_1 = m_data + *col_ptr;
	    for (k = row, k_ptr = m_inx+row; k < col; ++k, ++k_ptr)
          sum -= *(ptr_1 + m_colstart[*k_ptr] - *k_ptr) * 
		         *(ptr_2 + *k_ptr);                             // sum -= M(col,k)*M(k,row)
	    *(ptr_2 + *col_ptr) = sum/m_data[m_colstart[*col_ptr]]; // M(col,row) = sum/M(col,col)
	  }
    }
  
    // Matrix multiplication (col>=row)
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
	  double* ptr_diag = m_data + m_colstart[*row_ptr];
	  double* ptr_1    = ptr_diag - *row_ptr;
	  for (col = row, col_ptr = m_inx+row; col < m_num_inx; ++col, ++col_ptr) {
	    double  sum   = 0.0;
	    double* ptr_2 = m_data + m_colstart[*col_ptr] - *col_ptr;
	    for (k = col, k_ptr = m_inx+col; k < m_num_inx; ++k, ++k_ptr)
		  sum += *(ptr_1 + *k_ptr) * *(ptr_2 + *k_ptr);         // sum += M(row,k)*M(k,col)
	    *(ptr_1 + *col_ptr) = sum;                              // M(row,col) = sum
	  }
    }
  } // endelse: zero-row/col compression needed
  
  // Case C: all matrix elements are zero
  else
    throw matrix_zero("GSymMatrix::cholesky_invert(int)");

}


/*==========================================================================
 =                                                                         =
 =                      GSymMatrix private functions                       =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            GSymMatrix set_inx                           *
 * ----------------------------------------------------------------------- *
 * Determines the non-zero rows/columns in matrix and set up index array   *
 * that points to these rows/columns. This array is used for compressed    *
 * matrix calculations (Case B in the above routines).                     *
 ***************************************************************************/
void GSymMatrix::set_inx(void)
{
  // Allocate loop variables and pointers
  int row;
  int col;
	
  // Find all nonzero rows/cols
  m_num_inx = 0;
  for (row = 0; row < m_rows; ++row) {
    // Check first for zero diagonal. If we have one then check the rest of
	// the row
	if (m_data[m_colstart[row]] == 0.0) {
      for (col = 0; col < row; ++col) {
        if (m_data[m_colstart[col]+(row-col)] != 0.0)
		  break;
	  }
	  if (col < row)      // found a non-zero element
	    m_inx[m_num_inx++] = row;
	  else {
        for (col = row+1; col < m_cols; ++col) {
          if (m_data[m_colstart[row]+(col-row)] != 0.0)
		    break;
	    }
	    if (col < m_cols) // found a non-zero element
	      m_inx[m_num_inx++] = row;
	  }
	}
	else
      m_inx[m_num_inx++] = row;
  }

}


/*==========================================================================
 =                                                                         =
 =                            GSymMatrix friends                           =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                          Matrix absolute values                         *
 ***************************************************************************/
GSymMatrix fabs(const GSymMatrix& m)
{
  // Define result matrix
  GSymMatrix result(m.m_rows,m.m_cols);
  
  // Convert all elements to absolute values  
  for (int i = 0; i < m.m_elements; ++i)
    result.m_data[i] = fabs(m.m_data[i]);
  
  // Return result
  return result;
}


/*==========================================================================
 =                                                                         =
 =                      GSymMatrix exception classes                       =
 =                                                                         =
 ==========================================================================*/
