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
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_MUL_VEC    "GMatrix::operator* (const GVector&) const"
#define G_OP_ADD        "GMatrix::operator+= (const GMatrix&)"
#define G_OP_SUB        "GMatrix::operator-= (const GMatrix&)"
#define G_OP_MAT_MUL    "GMatrix::operator*= (const GMatrix&)"
#define G_ADD_COL       "GMatrix::add_col(const GVector&, int)"
#define G_EXTRACT_ROW   "GMatrix::extract_row(int) const"
#define G_EXTRACT_COL   "GMatrix::extract_col(int) const"
#define G_EXTRACT_LOWER "GMatrix::extract_lower_triangle() const"
#define G_EXTRACT_UPPER "GMatrix::extract_upper_triangle() const"
#define G_INSERT_COL    "GMatrix::insert_col(const GVector&, int)"
#define G_CONSTRUCTOR   "GMatrix::constructor(int, int)"


/*==========================================================================
 =                                                                         =
 =                      GMatrix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GMatrix initialisation method and finally  *
 * construct the object.                                                   *
 ***************************************************************************/
GMatrix::GMatrix(int rows, int cols) : GMatrixBase()
{
  // Initialise class members for clean destruction
  init_members();

  // Construct full matrix
  constructor(rows, cols);

  // Return
  return;
}


/***************************************************************************
 *                               Copy constructor                          *
 * ----------------------------------------------------------------------- *
 * First calls the void base class copy constructor that copys all base    *
 * class members, then call the GMatrix member copy method.                *
 ***************************************************************************/
GMatrix::GMatrix(const GMatrix& m) : GMatrixBase(m)
{
  // Initialise class members for clean destruction
  init_members();

  // Copy members
  copy_members(m);

  // Return
  return;
}


/***************************************************************************
 *               GSymMatrix -> GMatrix storage class conversion            *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GMatrix initialisation method and finally  *
 * construct the object and fill it with the GSymMatrix data.              *
 ***************************************************************************/
/*
GMatrix::GMatrix(const GSymMatrix& m) : GMatrixBase()
{ 
  // Initialise private members for clean destruction
  init_members();

  // Construct matrix
  constructor(m.rows(), m.cols());

  // Fill matrix
  int i_dst = 0;
  for (int col = 0; col < m_cols; ++col) {
    for (int row = 0; row < m_rows; ++row)
	  m_data[i_dst++] = m(row, col);
  }
  
  // Return
  return;
}
*/

/***************************************************************************
 *              GSparseMatrix -> GMatrix storage class conversion          *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GMatrix initialisation method and finally  *
 * construct the object and fill it with the GSparseMatrix data.           *
 ***************************************************************************/
/*
GMatrix::GMatrix(const GSparseMatrix& m) : GMatrixBase()
{ 
  // Initialise private members for clean destruction
  init_members();

  // Construct matrix
  constructor(m.rows(), m.cols());

  // Fill matrix
  int i_dst = 0;
  for (int col = 0; col < m_cols; ++col) {
    for (int row = 0; row < m_rows; ++row)
	  m_data[i_dst++] = m(row, col);
  }
  
  // Return
  return;
}
*/

/***************************************************************************
 *                           GMatrix destructor                            *
 * ----------------------------------------------------------------------- *
 * First destroys class members, then destroy base class members.          *
 ***************************************************************************/
GMatrix::~GMatrix()
{
  // Free members
  free_members();

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                             GMatrix operators                           =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Assignment operator                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GMatrix& GMatrix::operator= (const GMatrix& m)
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
GVector GMatrix::operator* (const GVector& v) const
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
GMatrix& GMatrix::operator+= (const GMatrix& m)
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
 *                      Unary matrix subtraction operator                  *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GMatrix& GMatrix::operator-= (const GMatrix& m)
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
 * In case of rectangular matrices the result matrix does not change and   *
 * the operations can be performed 'quasi' inplace (except of a vector).   *
 * For the general case the result matrix changes the size so for          *
 * simplicity a new matrix is allocated to hold the results.               *
 ***************************************************************************/
GMatrix& GMatrix::operator*= (const GMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows)
    throw GException::matrix_mismatch(G_OP_MAT_MUL, m_rows, m_cols, m.m_rows, m.m_cols);

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


/*==========================================================================
 =                                                                         =
 =                            GMatrix methods                              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        Add vector column into matrix                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrix::add_col(const GVector& v, int col)
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_ADD_COL, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw GException::matrix_vector_mismatch(G_ADD_COL, v.m_num, m_rows, m_cols);

  // Get start index of data in matrix
  int i = m_colstart[col];
  
  // Add column into vector
  for (int row = 0; row < m_rows; ++row)
    m_data[i++] += v.m_data[row];
	
  // Return
  return;
}


/***************************************************************************
 *               Extract row from matrix and put it into vector            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_cols);

  // Get start of data in matrix
  int i = row;
  
  // Extract row into vector
  for (int col = 0; col < m_cols; ++col, i+=m_rows)
    result.m_data[col] = m_data[i];
	
  // Return vector
  return result;

}


/***************************************************************************
 *            Extract column from matrix and put it into vector            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GMatrix::extract_col(int col) const
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_EXTRACT_COL, 0, col, m_rows, m_cols);
  #endif
  
  // Create result vector
  GVector result(m_rows);

  // Get start of data in matrix
  int i = m_colstart[col];
  
  // Extract column into vector
  for (int row = 0; row < m_rows; ++row)
    result.m_data[row] = m_data[i++];
	
  // Return vector
  return result;

}


/***************************************************************************
 *                      Extract lower triangle of matrix                   *
 * ----------------------------------------------------------------------- *
 * Triangular extraction only works for rectangular matrixes.              *
 ***************************************************************************/
GMatrix GMatrix::extract_lower_triangle() const
{
  // First check if matrix is rectangular. Raise an exception if this is not
  // the case
  if (m_rows != m_cols)
    throw GException::matrix_not_rectangular(G_EXTRACT_LOWER, m_rows, m_cols);

  // Define result matrix
  GMatrix result(m_rows, m_cols);

  // Extract all elements
  for (int col = 0; col < m_cols; ++col) {
    int i = m_colstart[col] + col;
    for (int row = col; row < m_rows; ++row, ++i)
	  result.m_data[i] = m_data[i];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                       Extract upper triangle of matrix                  *
 * ----------------------------------------------------------------------- *
 * Triangular extraction only works for rectangular matrixes.              *
 ***************************************************************************/
GMatrix GMatrix::extract_upper_triangle() const
{
  // First check if matrix is rectangular. Raise an exception if this is not
  // the case
  if (m_rows != m_cols)
    throw GException::matrix_not_rectangular(G_EXTRACT_UPPER, m_rows, m_cols);

  // Define result matrix
  GMatrix result(m_rows, m_cols);

  // Extract all elements
  for (int col = 0; col < m_cols; ++col) {
    int i = m_colstart[col];
    for (int row = 0; row <= col; ++row, ++i)
	  result.m_data[i] = m_data[i];
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
double GMatrix::fill() const
{
  // Determine the number of zero elements
  int zero = 0;
  for (int i = 0; i < m_elements; ++i) {
    if (m_data[i] == 0.0)
	  zero++;
  }
  
  // Return the fill
  return (1.0-double(zero)/double(m_elements));
}


/***************************************************************************
 *                     Insert vector column into matrix                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrix::insert_col(const GVector& v, int col)
{
  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_INSERT_COL, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw GException::matrix_vector_mismatch(G_INSERT_COL, v.m_num, m_rows, m_cols);

  // Get start index of data in matrix
  int i = m_colstart[col];
  
  // Insert column into vector
  for (int row = 0; row < m_rows; ++row)
    m_data[i++] = v.m_data[row];
	
  // Return
  return;
}


/***************************************************************************
 *                        Inplace transpose of matrix                      *
 * ----------------------------------------------------------------------- *
 * The transpose operation exchanges the number of rows against the number *
 * of columns. The element transformation is done inplace.                 *
 *                     NON RECTANGULAR NOT CODED NATIVELY !!!              *
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


/*==========================================================================
 =                                                                         =
 =                          GMatrix private methods                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Constructor method                           *
 * ----------------------------------------------------------------------- *
 * This is the main constructor code, without any initialisation in it.    *
 * It is used as service function to a number of different constructors of *
 * the GMatrix class.                                                      *
 ***************************************************************************/
void GMatrix::constructor(int rows, int cols)
{
  // Determine number of elements to store in matrix
  int elements = rows*cols;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw GException::empty(G_CONSTRUCTOR);
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new int[cols+1];
  if (m_data == NULL || m_colstart == NULL)
	throw GException::mem_alloc(G_CONSTRUCTOR, elements);
	
  // Store matrix size (logical, storage, allocated)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;
  m_alloc    = elements;
  
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


/*==========================================================================
 =                                                                         =
 =                             GMatrix friends                             =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GMatrix& m)
{
  // Put header in stream
  os << "=== GMatrix ===" << std::endl;
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
