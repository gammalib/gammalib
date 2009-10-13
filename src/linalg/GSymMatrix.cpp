/***************************************************************************
 *                  GSymMatrix.cpp  -  symmetric matrix class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2009 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CAST_MATRIX       "GSymMatrix::GSymMatrix(const GMatrix&)"
#define G_CAST_SPARSEMATRIX "GSymMatrix::GSymMatrix(const GSparseMatrix&)"
#define G_OP_ADD            "GSymMatrix::operator+= (const GSymMatrix&)"
#define G_OP_SUB            "GSymMatrix::operator-= (const GSymMatrix&)"
#define G_OP_MUL_VEC        "GSymMatrix::operator* (const GVector&) const"
#define G_OP_MAT_MUL        "GSymMatrix::operator*= (const GSymMatrix&)"
#define G_ADD_COL           "GSymMatrix::add_col(const GVector&, int)"
#define G_CHOL_DECOMP       "GSymMatrix::cholesky_decompose(int)"
#define G_CHOL_SOLVE        "GSymMatrix::cholesky_solver(const GVector&, int)"
#define G_CHOL_INVERT       "GSymMatrix::cholesky_invert(int)"
#define G_EXTRACT_ROW       "GSymMatrix::extract_row(int)"
#define G_EXTRACT_COL       "GSymMatrix::extract_col(int)"
#define G_INSERT_COL        "GSymMatrix::insert_col(const GVector&, int)"
#define G_CONSTRUCTOR       "GSymMatrix::constructor(int, int)"
#define G_COPY_MEMBERS      "GSymMatrix::copy_members(const GSymMatrix&)"


/*==========================================================================
 =                                                                         =
 =                   GSymMatrix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GMatrix initialisation method and finally  *
 * construct the object.                                                   *
 ***************************************************************************/
GSymMatrix::GSymMatrix(int rows, int cols) : GMatrixBase()
{
  // Initialise class members for clean destruction
  init_members();

  // Construct full matrix
  constructor(rows, cols);

  // Return
  return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
 * First calls the void base class copy constructor that copys all base    *
 * class members, then call the GSymMatrix member copy method.             *
 ***************************************************************************/
GSymMatrix::GSymMatrix(const GSymMatrix& m) : GMatrixBase(m)
{ 
  // Initialise private members for clean destruction
  init_members();

  // Copy members
  copy_members(m);

  // Return
  return;
}


/***************************************************************************
 *                                Destructor                               *
 * ----------------------------------------------------------------------- *
 * First destroys class members, then destroy base class members.          *
 ***************************************************************************/
GSymMatrix::~GSymMatrix()
{
  // Free members
  free_members();

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                           GSymMatrix operators                          =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                           Assignment operator                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator= (const GSymMatrix& m)
{ 
  // Execute only if object is not identical
  if (this != &m) { 

    // Copy base class members
	this->GMatrixBase::operator=(m);

    // Free members
    free_members();
  
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(m);
	
  } // endif: object was not identical

  // Return this object
  return *this;
}


/***************************************************************************
 *                           Vector multiplication                         *
 * ----------------------------------------------------------------------- *
 * Multiply matrix by a vector. The result is a vector.                    *
 ***************************************************************************/
GVector GSymMatrix::operator* (const GVector& v) const
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_cols != v.m_num)
    throw GException::matrix_vector_mismatch(G_OP_MUL_VEC, v.m_num, m_rows, m_cols);

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
 *                        Unary matrix addition operator                   *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator+= (const GSymMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw GException::matrix_mismatch(G_OP_ADD, m_rows, m_cols, m.m_rows, m.m_cols);

  // Add matrices
  addition(m);

  // Return result
  return *this;
}


/***************************************************************************
 *           GSymMatrix subtraction (unary subtraction operator)           *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator-= (const GSymMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw GException::matrix_mismatch(G_OP_SUB, m_rows, m_cols, m.m_rows, m.m_cols);

  // Subtract matrices
  subtraction(m);

  // Return result
  return *this;
}



/***************************************************************************
 *                            Matrix multiplication                        *
 * ----------------------------------------------------------------------- *
 * Symmetric matrices must have the same dimensions to be multiplied.      *
 ***************************************************************************/
GSymMatrix& GSymMatrix::operator*= (const GSymMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows || m_rows != m.m_cols)
    throw GException::matrix_mismatch(G_OP_MAT_MUL, m_rows, m_cols, m.m_rows, m.m_cols);

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
 =                            GSymMatrix methods                           =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                       Add vector column into matrix                     *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSymMatrix::add_col(const GVector& v, int col)
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_ADD_COL, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw GException::matrix_vector_mismatch(G_ADD_COL, v.m_num, m_rows, m_cols);
  
  // Insert column into vector
  for (int row = 0; row < m_rows; ++row)
	(*this)(row,col) += v(row);
	
  // Return
  return;
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
	double diag = 0.0;
    for (int row = 0; row < m_rows; ++row) {
	  double* ptr = m_data + m_colstart[row];
      for (int col = row; col < m_cols; ++col, ++ptr) {
	    double sum = *ptr;                                // sum = M(row,col)
	    for (int k = 0; k < row; ++k) {
	      int offset = m_colstart[k] - k;            // is always positive
	      sum -= m_data[offset+row] * m_data[offset+col]; // sum -= M(row,k)*M(col,k)
	    }
	    if (row == col) {
	      if (sum <= 0.0)
		    throw GException::matrix_not_pos_definite(G_CHOL_DECOMP, row, sum);
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
	double diag = 0.0;
    for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
	  double* ptr_0 = m_data + m_colstart[*row_ptr] - *row_ptr;
      for (col = row, col_ptr = m_inx + row; col < m_num_inx; ++col, ++col_ptr) {
		double* ptr = ptr_0 + *col_ptr;
	    double  sum = *ptr;                                         // sum = M(row,col)
	    for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr) {
	      int offset = m_colstart[*k_ptr] - *k_ptr;            // is always positive
	      sum -= m_data[offset+*row_ptr] * m_data[offset+*col_ptr]; // sum -= M(row,k)*M(col,k)
	    }
	    if (*row_ptr == *col_ptr) {
	      if (sum <= 0.0)
		    throw GException::matrix_not_pos_definite(G_CHOL_DECOMP, *row_ptr, sum);
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
    throw GException::matrix_zero(G_CHOL_DECOMP);
	
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
    throw GException::matrix_vector_mismatch(G_CHOL_SOLVE, v.m_num, m_rows, m_cols);

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
    throw GException::matrix_zero(G_CHOL_SOLVE);

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
    throw GException::matrix_zero(G_CHOL_INVERT);

  // Return
  return;
}


/***************************************************************************
 *                     Extract row from matrix into vector                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GSymMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows, m_cols);
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
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GSymMatrix::extract_col(int col) const
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_EXTRACT_COL, 0, col, m_rows, m_cols);
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
 *        Extract lower triangle of symmetric matrix into full matrix      *
 * ----------------------------------------------------------------------- *
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
 * ----------------------------------------------------------------------- *
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
 *                          Determine fill of matrix                       *
 * ----------------------------------------------------------------------- *
 * The fill of a matrix is defined as the ratio of non-zero elements       *
 * versus the number of total elements.                                    *
 ***************************************************************************/
double GSymMatrix::fill() const
{
  // Determine the number of zero elements
  int zero = 0;
  for (int col = 0, i = 0; col < m_cols; ++col) {
    if (m_data[i++] == 0.0)                       // Count diag. once
	  zero++;
    for (int row = col+1; row < m_rows; ++row) {
      if (m_data[i++] == 0.0)                     // Count off-diag. twice
	    zero +=2;
	}
  }
  
  // Return the fill
  return (1.0-double(zero)/double(m_elements));
}


/***************************************************************************
 *                        Insert vector into matrix                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSymMatrix::insert_col(const GVector& v, int col)
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_INSERT_COL, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw GException::matrix_vector_mismatch(G_INSERT_COL, v.m_num, m_rows, m_cols);
  
  // Insert column into vector
  for (int row = 0; row < m_rows; ++row)
	(*this)(row,col) = v(row);
	
  // Return
  return;
}


/***************************************************************************
 *                            Sum matrix elements                          *
 * ----------------------------------------------------------------------- *
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


/*==========================================================================
 =                                                                         =
 =                      GSymMatrix private functions                       =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Constructor method                           *
 * ----------------------------------------------------------------------- *
 * This is the main constructor code, without any initialisation in it.    *
 * It is used as service function to a number of different constructors of *
 * the GSymMatrix class.                                                   *
 ***************************************************************************/
void GSymMatrix::constructor(int rows, int cols)
{
  // Throw exception if number of rows and columns is not identical
  if (rows != cols)
    throw GException::matrix_not_symmetric(G_CONSTRUCTOR, rows, cols);

  // Determine number of physical elements in matrix
  int elements = rows*(rows+1)/2;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw GException::empty(G_CONSTRUCTOR);
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new int[cols+1];
  m_inx      = new int[cols];
  if (m_data == NULL || m_colstart == NULL || m_inx == NULL)
	throw GException::mem_alloc(G_CONSTRUCTOR, elements);

  // Store matrix size (logical and physical)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;
  m_alloc    = elements;

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
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSymMatrix::init_members(void)
{
  // Initialise members
  m_num_inx = 0;
  m_inx     = NULL;

  // Return
  return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSymMatrix::copy_members(const GSymMatrix& m)
{
  // Copy index selection
  m_num_inx = m.m_num_inx;
  m_inx     = new int[m_cols];
  if (m_inx == NULL)
	throw GException::mem_alloc(G_COPY_MEMBERS, m_cols);
  for (int i = 0; i < m_cols; ++i)
	m_inx[i] = m.m_inx[i];
  
  // Return
  return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSymMatrix::free_members(void)
{
  // De-allocate only if memory has indeed been allocated by derived class
  if (m_inx != NULL) delete [] m_inx;

  // Return
  return;
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
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSymMatrix& m)
{
  // Put header in stream
  os << "=== GSymMatrix ===" << std::endl;
  if (m.m_rowsel != NULL)
    os << " Number of rows ............: " << m.m_rows << " (compressed " <<
	      m.m_num_rowsel << ")" << std::endl;
  else
    os << " Number of rows ............: " << m.m_rows << std::endl;
  if (m.m_colsel != NULL)
    os << " Number of columns .........: " << m.m_cols << " (compressed " <<
	      m.m_num_colsel << ")" << std::endl;
  else
    os << " Number of columns .........: " << m.m_cols << std::endl;
  os << " Number of elements ........: " << m.m_elements << std::endl;
  os << " Number of allocated cells .: " << m.m_alloc << std::endl;

  // Dump elements and compression schemes
  m.dump_elements(os);
  m.dump_row_comp(os);
  m.dump_col_comp(os);
  
  // Return output stream
  return os;
}


/***************************************************************************
 *                          Matrix absolute values                         *
 * ----------------------------------------------------------------------- *
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
