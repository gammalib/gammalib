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


/***************************************************************************
 *                          GSymMatrix constructor                         *
 ***************************************************************************/
GSymMatrix::GSymMatrix(unsigned rows, unsigned cols) : GMatrix()
{
  // Initialise private members for clean destruction
  m_rows        = 0;
  m_cols        = 0;
  m_elements    = 0;
  m_num_inx     = 0;
  m_data        = NULL;
  m_colstart    = NULL;
  m_inx         = NULL;

  // Throw exception if number of rows and columns is not identical
  if (rows != cols)
    throw not_sym("GSymMatrix constructor", rows, cols);

  // Determine number of physical elements in matrix
  unsigned elements = rows*(rows+1)/2;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw empty("GSymMatrix constructor");
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new unsigned[cols+1];
  m_inx      = new unsigned[cols];
  if (m_data == NULL || m_colstart == NULL || m_inx == NULL)
	throw mem_alloc("GSymMatrix constructor", elements);

  // Store matrix size (logical and physical)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;

  // Set-up column start indices
  m_colstart[0]   = 0;
  unsigned offset = rows;
  for (unsigned col = 1; col <= m_cols; ++col)
    m_colstart[col] = m_colstart[col-1] + offset--;

  // Initialise matrix elements to 0.0
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] = 0.0;
}


/***************************************************************************
 *                        GSymMatrix copy constructor                      *
 * ----------------------------------------------------------------------- *
 * First invoques GMatrix copy constructor, then this function.            *
 ***************************************************************************/
GSymMatrix::GSymMatrix(const GSymMatrix& m) : GMatrix(m)
{ 
  m_num_inx = m.m_num_inx;
  m_inx     = new unsigned[m_cols];
  if (m_inx == NULL)
	throw mem_alloc("GSymMatrix copy constructor", m_cols);
  for (unsigned i = 0; i < m_cols; ++i)
	m_inx[i] = m.m_inx[i];
}


/***************************************************************************
 *                         GSymMatrix destructor                           *
 * ----------------------------------------------------------------------- *
 * First invoques this function, then ~GMatrix destructor.                 *
 ***************************************************************************/
GSymMatrix::~GSymMatrix()
{
  // De-allocate only if memory has indeed been allocated by derived class
  if (m_inx != NULL) delete[] m_inx;
}


/***************************************************************************
 *                      GSymMatrix assignment operator                     *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator= (const GSymMatrix& m)
{ 
  if (this != &m) { 
    this->GMatrix::operator=(m);
    m_num_inx = m.m_num_inx;
    m_inx     = new unsigned[m_cols];
    if (m_inx == NULL)
	  throw mem_alloc("GSymMatrix assignment operator", m_cols);
    for (unsigned i = 0; i < m_cols; ++i)
	  m_inx[i] = m.m_inx[i];
  }
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
    for (unsigned row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
      for (unsigned col = row; col < m_cols; ++col, ++ptr) {
	    double sum = *ptr;                                // sum = M(row,col)
	    for (unsigned k = 0; k < row; ++k) {
	      unsigned offset = m_colstart[k] - k;            // is always positive
	      sum -= m_data[offset+row] * m_data[offset+col]; // sum -= M(row,k)*M(col,k)
	    }
	    double diag;
	    if (row == col) {
	      if (sum <= 0.0)
		    throw not_pos_definite("GSymMatrix cholesky_decompose", row, sum);
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
	unsigned  row;
	unsigned  col;
	unsigned  k;
	unsigned* row_ptr;
	unsigned* col_ptr;
	unsigned* k_ptr;

    // Loop over upper triangle (col >= row)
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
	  double* ptr_0 = m_data + m_colstart[*row_ptr] - *row_ptr;
      for (col = row, col_ptr = m_inx + row; col < m_num_inx; ++col, ++col_ptr) {
		double* ptr = ptr_0 + *col_ptr;
	    double  sum = *ptr;                                         // sum = M(row,col)
	    for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr) {
	      unsigned offset = m_colstart[*k_ptr] - *k_ptr;            // is always positive
	      sum -= m_data[offset+*row_ptr] * m_data[offset+*col_ptr]; // sum -= M(row,k)*M(col,k)
	    }
	    double diag;
	    if (*row_ptr == *col_ptr) {
	      if (sum <= 0.0)
		    throw not_pos_definite("GSymMatrix cholesky_decompose", *row_ptr, sum);
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
    throw all_zero("GSymMatrix cholesky_decompose");
	
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
  // Allocate result vector
  GVector x(m_rows);

  // Check if zero-row/col compression is needed  
  int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

  // Case A: no zero-row/col compression needed
  if (no_zeros) {
    
    // Solve L*y=b, storing y in x (row>k)
    for (unsigned row = 0; row < m_rows; ++row) {
      double sum = v(row);
	  for (unsigned k = 0; k < row; ++k)
	    sum -= m_data[m_colstart[k]+(row-k)] * x(k); // sum -= M(row,k) * x(k)
	  x(row) = sum/m_data[m_colstart[row]];          // x(row) = sum/M(row,row) 
    }
	  
    // Solve trans(L)*x=y (k>row)
    for (int row = m_rows-1; row >= 0; --row) {
      double  sum = x(row);
	  double* ptr = m_data + m_colstart[row] + 1;
	  for (unsigned k = row+1; k < m_rows; ++k)
	    sum -= *ptr++ * x(k);                        // sum -= M(k,row) * x(k)
	  x(row) = sum/m_data[m_colstart[row]];          // x(row) = sum/M(row,row) 
    }
  } // endif: no zero-row/col compression needed
  
  // Case B: zero-row/col compression needed
  else if (m_num_inx > 0) {

    // Allocate loop variables and pointers
	int       row;
	unsigned  k;
	unsigned* row_ptr;
	unsigned* k_ptr;

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
    throw all_zero("GSymMatrix cholesky_solver");

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
    for (unsigned row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
	  *ptr        = 1.0/(*ptr);                       // M(row,row) = 1/M(row,row)
	  for (unsigned col = row+1; col < m_cols; ++col) {
	    double   sum = 0.0;
	    double* ptr1 = m_data + col - row;
	    double* ptr2 = ptr;
	    for (unsigned k = row; k < col; ++k)
	      sum -= *(ptr1-- + m_colstart[k]) * *ptr2++; // sum -= M(col,k)*M(k,row)
	    *(ptr+col-row) = sum/m_data[m_colstart[col]]; // M(col,row) = sum/M(col,col)
	  }
    }
  
    // Matrix multiplication (col>=row)
    for (unsigned row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
	  for (unsigned col = row; col < m_cols; ++col) {
	    double   sum = 0.0;
	    double* ptr1 = ptr + col - row;
	    double* ptr2 = m_data + m_colstart[col];
	    for (unsigned k = col; k < m_cols; ++k)
	      sum += *ptr1++ * *ptr2++;                   // sum += M(row,k)*M(k,col)
	    *(ptr+col-row) = sum;                         // M(row,col) = sum
	  }
    }
  } // endif: no zero-row/col compression needed
  
  // Case B: zero-row/col compression needed
  else if (m_num_inx > 0) {

    // Allocate loop variables and pointers
	unsigned  row;
	unsigned  col;
	unsigned  k;
	unsigned* row_ptr;
	unsigned* col_ptr;
	unsigned* k_ptr;

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
    throw all_zero("GSymMatrix cholesky_invert");

}


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
  unsigned row;
  unsigned col;
	
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


/***************************************************************************
 *                          Class exception handlers                       *
 ***************************************************************************/
// Matrix not symmetric
GSymMatrix::not_sym::not_sym(string origin, unsigned cols, unsigned rows)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_cols;
  s_rows << rows;
  s_cols << cols;
  m_message = "Matrix is not symmetric [" + s_rows.str() + "," +
			  s_cols.str() + "]";
}

// Matrix not positive definite
GSymMatrix::not_pos_definite::not_pos_definite(string origin, unsigned row,
                                               double sum)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_sum;
  s_rows << row;
  s_sum << scientific << sum;
  m_message = "Matrix is not positive definite (sum " + s_sum.str() + 
              " occured in row " + s_rows.str() + ")";
}

// All matrix elements are zero
GSymMatrix::all_zero::all_zero(string origin)
{
  m_origin = origin;
  m_message = "All matrix elements are zero";
}
