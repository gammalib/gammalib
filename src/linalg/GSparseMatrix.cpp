/***************************************************************************
 *                  GSparseMatrix.cpp  -  sparse matrix class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2009 by Jurgen Knodlseder                           *
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
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"
#include "GSparseSymbolic.hpp"
#include "GSparseNumeric.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Method name definitions ____________________________________________ */
#define G_OP_MUL_VEC   "GSparseMatrix::operator* (const GVector&) const"
#define G_OP_ADD       "GSparseMatrix::operator+= (const GSparseMatrix&)"
#define G_OP_SUB       "GSparseMatrix::operator-= (const GSparseMatrix&)"
#define G_OP_MAT_MUL   "GSparseMatrix::operator*= (const GSparseMatrix&)"
#define G_ADD_COL      "GSparseMatrix::add_col(const GVector&, int)"
#define G_ADD_COL2     "GSparseMatrix::add_col(const double*, const int*, int, int)"
#define G_CHOL_DECOMP  "GSparseMatrix::cholesky_decompose(int)"
#define G_CHOL_SOLVE   "GSparseMatrix::cholesky_solver(const GVector&, int)"
#define G_EXTRACT_ROW  "GSparseMatrix::extract_row(int) const"
#define G_EXTRACT_COL  "GSparseMatrix::extract_col(int) const"
#define G_INSERT_COL   "GSparseMatrix::insert_col(const GVector&, int)"
#define G_INSERT_COL2  "GSparseMatrix::insert_col(const double*, const int*, int, int)"
#define G_STACK_INIT   "GSparseMatrix::stack_init(int, int)"
#define G_STACK_PUSH   "GSparseMatrix::stack_push_column(const double*, const int*, int, int)"
#define G_STACK_FLUSH  "GSparseMatrix::stack_flush()"
#define G_CONSTRUCTOR  "GSparseMatrix::constructor(int, int, int)"
#define G_COPY_MEMBERS "GSparseMatrix::copy_members(const GSparseMatrix&)"
#define G_ALLOC        "GSparseMatrix::alloc_elements(int, int)"
#define G_FREE         "GSparseMatrix::free_elements(int, int)"
#define G_REMOVE_ZERO  "GSparseMatrix::remove_zero_row_col()"
#define G_SYMPERM      "cs_symperm(GSparseMatrix*, const int*, int)"
#define G_TRANSPOSE    "cs_transpose(GSparseMatrix*, int)"

/* __ Macros _____________________________________________________________ */
#define G_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define G_MAX(a,b) (((a) > (b)) ? (a) : (b))

/* __ Coding definitions _________________________________________________ */
//#define G_USE_MEMCPY   // Do not use due to problems on 64 Bit systems
//#define G_USE_MEMMOVE  // Do not use due to problems on 64 Bit systems

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_SPARSE_PENDING                    // Analyse pending values
//#define G_DEBUG_SPARSE_INSERTION                 // Analyse value insertion
//#define G_DEBUG_SPARSE_ADDITION                   // Analyse value addition
//#define G_DEBUG_SPARSE_COMPRESSION   // Analyse zero row/column compression
//#define G_DEBUG_SPARSE_MALLOC                  // Analyse memory management
//#define G_DEBUG_SPARSE_STACK_PUSH           // Analyse stack column pushing
//#define G_DEBUG_SPARSE_STACK_FLUSH         // Analyse stack column flushing



/*==========================================================================
 =                                                                         =
 =                   GSparseMatrix constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                               Constructor                               *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GSparseMatrix initialisation method and    *
 * finally construct the object.                                           *
 ***************************************************************************/
GSparseMatrix::GSparseMatrix(int rows, int cols, int elements) : GMatrixBase()
{
  // Initialise private members for clean destruction
  init_members();

  // Construct sparse matrix
  constructor(rows, cols, elements);

  // Return
  return;
}


/***************************************************************************
 *                             Copy constructor                            *
 * ----------------------------------------------------------------------- *
 * First calls the void base class copy constructor that copys all base    *
 * class members, then call the GSparseMatrix member copy method.          *
 * NOTE: The copy constructor does NOT fill a pending element.             *
 ***************************************************************************/
GSparseMatrix::GSparseMatrix(const GSparseMatrix& m) : GMatrixBase(m)
{
  // Initialise private members for clean destruction
  init_members();

  // Copy members
  copy_members(m);

  // Return
  return;
}


/***************************************************************************
 *              GMatrix -> GSparseMatrix storage class conversion          *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GSparseMatrix initialisation method and    *
 * finally construct the object and fill it with the GMatrix data.         *
 ***************************************************************************/
/*
GSparseMatrix::GSparseMatrix(const GMatrix& m) : GMatrixBase()
{ 
  // Initialise private members for clean destruction
  init_members();

  // Construct sparse matrix
  constructor(m.rows(), m.cols());

  // Use a fill-stack for better performance
  stack_init(m_cols, G_SPARSE_MATRIX_DEFAULT_STACK_ENTRIES);
  
  // Loop over all columns
  GVector column(m_cols);
  for (int col = 0; col < m_cols; ++col) {
    column = m.extract_col(col);
	insert_col(column, col);
  }
  
  // Destroy fill stack
  stack_destroy();
  
  // Return result
  return;
}
*/

/***************************************************************************
 *          GSymMatrix -> GSparseMatrix storage class conversion           *
 * ----------------------------------------------------------------------- *
 * First calls the void base class constructor that initialises all base   *
 * class members, then call the GSparseMatrix initialisation method and    *
 * finally construct the object and fill it with the GSymMatrix data.      *
 ***************************************************************************/
/*
GSparseMatrix::GSparseMatrix(const GSymMatrix& m) : GMatrixBase()
{ 
  // Initialise private members for clean destruction
  init_members();

  // Construct sparse matrix
  constructor(m.rows(), m.cols());

  // Use a fill-stack for better performance
  stack_init(m_cols, G_SPARSE_MATRIX_DEFAULT_STACK_ENTRIES);
  
  // Loop over all columns
  GVector column(m_cols);
  for (int col = 0; col < m_cols; ++col) {
    column = m.extract_col(col);
	insert_col(column, col);
  }
  
  // Destroy fill stack
  stack_destroy();
  
  // Return result
  return;
}
*/

/***************************************************************************
 *                                Destructor                               *
 * ----------------------------------------------------------------------- *
 * First destroys class members, then destroy base class members.          *
 ***************************************************************************/
GSparseMatrix::~GSparseMatrix()
{
  // Free members
  free_members();

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                         GSparseMatrix operators                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator= (const GSparseMatrix& m)
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
    throw GException::matrix_vector_mismatch(G_OP_MUL_VEC, v.m_num, m_rows, m_cols);

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
      cout << G_OP_MUL_VEC << ": pending value " << m_fill_val << 
          " for location (" <<  m_fill_row << "," << m_fill_col <<
          ") has been used." << endl;
      #endif
    }

  } // endif: there were elements in matrix

  // Return result
  return result;
}


/***************************************************************************
 *                              Equalty operator                           *
 * ----------------------------------------------------------------------- *
 * Two matrixes are considered equal if they have the same dimensions and  *
 * identicial elements.                                                    *
 *                     !!! NOT YET IN NATIVE CODING !!!                    *
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
 *                              Equalty operator                           *
 * ----------------------------------------------------------------------- *
 * Two matrixes are considered equal if they have the same dimensions and  *
 * identicial elements.                                                    *
 *                     !!! NOT YET IN NATIVE CODING !!!                    *
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
 *                        Unary matrix addition operator                   *
 * ----------------------------------------------------------------------- *
 *                     !!! NOT YET IN NATIVE CODING !!!                    *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator+= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw GException::matrix_mismatch(G_OP_ADD, m_rows, m_cols, m.m_rows, m.m_cols);

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
 *                      Unary matrix subtraction operator                  *
 * ----------------------------------------------------------------------- *
 *                     !!! NOT YET IN NATIVE CODING !!!                    *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator-= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw GException::matrix_mismatch(G_OP_SUB, m_rows, m_cols, m.m_rows, m.m_cols);

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
 *                     !!! NOT YET IN NATIVE CODING !!!                    *
 ***************************************************************************/
GSparseMatrix& GSparseMatrix::operator*= (const GSparseMatrix& m)
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows || m_rows != m.m_cols)
    throw GException::matrix_mismatch(G_OP_MAT_MUL, m_rows, m_cols, m.m_rows, m.m_cols);

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
 =                          GSparseMatrix methods                          =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Add vector to matrix column                     *
 * ----------------------------------------------------------------------- *
 * This is the main driver routine to add data to a matrix. It handles     *
 * both normal and stack-based filled. Note that there is another instance *
 * of this function that takes a compressed array.                         *
 ***************************************************************************/
void GSparseMatrix::add_col(const GVector& v, int col)
{
  // Debug header
  #if defined(G_DEBUG_SPARSE_ADDITION)
  cout << "GSparseMatrix::add_col([" << v << "], " << col << "):" << endl;
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

  // Initialise number of non-zero elements to 0
  int non_zero = 0;

  // If we have a stack then try to push vector on stack first. Note that
  // stack_push_column does its own argument verifications, so to avoid
  // double checking we don't do anything before this call ...
  if (m_stack_data != NULL) {
    non_zero = stack_push_column(v, col);
    if (non_zero == 0)
      return;
  }

  // ... otherwise check first the arguments and determine the number of
  // non-zero elements in vector
  else {

    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
      throw GException::out_of_range(G_ADD_COL, 0, col, m_rows, m_cols);
    #endif

    // Raise an exception if the matrix and vector dimensions are incompatible
    if (m_rows != v.m_num)
      throw GException::matrix_vector_mismatch(G_ADD_COL, v.m_num, m_rows, m_cols);

    // Determine number of non-zero elements in vector
    non_zero = v.non_zeros();
  }

  // Extract vector for column, add elements, and re-insert vector (only if
  // vector to insert has non-zeros)
  if (non_zero > 0) {

    // Copy input vector
    GVector column = v;

    // Add elements to vector
    for (int i = m_colstart[col]; i < m_colstart[col+1]; ++i)
      column.m_data[m_rowinx[i]] += m_data[i];

    // If there is a pending element then put it in the vector
    if (m_fill_val != 0.0 && m_fill_col == col)
      column.m_data[m_fill_row] += m_fill_val;

    // Insert vector into matrix
    insert_col(column, col);

  }

  // Debugging: show sparse matrix after addition
  #if defined(G_DEBUG_SPARSE_ADDITION)
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
 *                   Add compressed array to matrix column                 *
 * ----------------------------------------------------------------------- *
 * This is the main driver routine to add data to a matrix. It handles     *
 * both normal and stack-based filled. Note that there is another instance *
 * of this function that takes a vector.                                   *
 ***************************************************************************/
void GSparseMatrix::add_col(const double* values, const int* rows, 
                            int number, int col)
{
  // Debug header
  #if defined(G_DEBUG_SPARSE_ADDITION)
  cout << "GSparseMatrix::add_col(v, i, n, " << col << "):" << endl;
  cout << " Matrix Data : ";
  for (int i = 0; i < m_elements; ++i)
    cout << m_data[i] << " ";
  cout << endl << " Matrix Row .: ";
  for (int i = 0; i < m_elements; ++i)
    cout << m_rowinx[i] << " ";
  cout << endl << " Matrix Col .: ";
  for (int i = 0; i < m_cols+1; ++i)
    cout << m_colstart[i] << " ";
  cout << endl;
  cout << " Array Data .: ";
  for (int i = 0; i < number; ++i)
    cout << values[i] << " ";
  cout << endl << " Array Row ..: ";
  for (int i = 0; i < number; ++i)
    cout << rows[i] << " ";
  cout << endl;
  #endif

  // If we have a stack then try to push elements on stack first. Note that
  // stack_push_column does its own argument verifications, so to avoid
  // double checking we don't do anything before this call ...
  if (m_stack_data != NULL) {
    number = stack_push_column(values, rows, number, col);
    if (number == 0)
      return;
  }

  // ... otherwise check the arguments
  else {
    // If the array is empty there is nothing to do
    if (!values || !rows || (number < 1))
      return;

    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
      throw GException::out_of_range(G_ADD_COL2, 0, col, m_rows, m_cols);
    #endif

    // Raise an exception if the index array seems incompatible with matrix 
    // dimensions
    if (rows[number-1] >= m_rows)
      throw GException::matrix_vector_mismatch(G_ADD_COL2, rows[number-1], m_rows, m_cols);
  } // endelse: there was no stack

  // Get indices of column in matrix
  int i_start = m_colstart[col];
  int i_stop  = m_colstart[col+1];

  // Case A: the column exists in the matrix, so mix new elements with existing
  // data
  if (i_start < i_stop) {

    // Allocate workspace to hold combined column
    int     wrk_size   = number + i_stop - i_start;
    double* wrk_double = new double[wrk_size];
    int*    wrk_int    = new int[wrk_size];
    if (wrk_double == NULL || wrk_int == NULL)
      throw GException::mem_alloc(G_ADD_COL2, wrk_size);

    // Mix matrix column with specified data
    int num_mix;
    mix_column(&(m_data[i_start]), &(m_rowinx[i_start]), i_stop-i_start,
               values, rows, number,
               wrk_double, wrk_int, &num_mix);

    // Insert mixed column
    insert_col(wrk_double, wrk_int, num_mix, col);

    // Free workspace
    delete [] wrk_int;
    delete [] wrk_double;

  } // endif: Case A

  // Case B: the column does not yet exist in the matrix, so just insert it
  else
    insert_col(values, rows, number, col);

  // Debugging: show sparse matrix after insertion
  #if defined(G_DEBUG_SPARSE_ADDITION)
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

  // Delete any existing symbolic and numeric analysis object and reset
  // pointers
  if (m_symbolic != NULL) delete (GSparseSymbolic*)m_symbolic;
  if (m_numeric  != NULL) delete (GSparseNumeric*)m_numeric;
  m_symbolic = NULL;
  m_numeric  = NULL;

  // Allocate symbolic analysis object
  GSparseSymbolic* symbolic = new GSparseSymbolic();
  if (symbolic == NULL)
    throw GException::mem_alloc(G_CHOL_DECOMP, 1);

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

  // Store symbolic pointer in sparse matrix object
  m_symbolic = (void*)symbolic;

  // Perform numeric Cholesky decomposition
  numeric.cholesky_numeric_analysis(*this, *symbolic);

  // Copy L matrix into this object
  free_elements(0, m_elements);
  alloc_elements(0, numeric.m_L->m_elements);
  for (int i = 0; i < m_elements; ++i) {
    m_data[i]   = numeric.m_L->m_data[i];
    m_rowinx[i] = numeric.m_L->m_rowinx[i];
  }
  for (int col = 0; col <= m_cols; ++col)
    m_colstart[col] = numeric.m_L->m_colstart[col];

  // Insert zero rows and columns if they have been removed previously.
  if (compress)
    insert_zero_row_col(matrix_rows, matrix_cols);

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
    throw GException::matrix_vector_mismatch(G_CHOL_SOLVE, v.m_num, m_rows, m_cols);

  // Raise an exception if there is no symbolic pointer
  if (!m_symbolic)
    throw GException::matrix_not_factorised(G_CHOL_SOLVE, "Cholesky decomposition");

  // Get symbolic pointer
  GSparseSymbolic* symbolic = (GSparseSymbolic*)m_symbolic;

  // Raise an exception if there is no permutation
  if (!symbolic->m_pinv)
    throw GException::matrix_not_factorised(G_CHOL_SOLVE, "Cholesky decomposition");

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
      throw GException::mem_alloc(G_CHOL_SOLVE, m_rows);
    int* col_map = new int[m_cols];
    if (col_map == NULL)
      throw GException::mem_alloc(G_CHOL_SOLVE, m_cols);

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


/***************************************************************************
 *                       Set all matrix elements to 0                      *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::clear()
{
  // Free memory
  free_elements(0, m_elements);

  // Initialise column start indices to 0
  for (int col = 0; col <= m_cols; ++col)
    m_colstart[col] = 0;

  // Return
  return;
}


/***************************************************************************
 *                     Extract row from matrix into vector                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GVector GSparseMatrix::extract_row(int row) const
{
  // Raise an exception if the row index is invalid
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows)
    throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows, m_cols);
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
    throw GException::out_of_range(G_EXTRACT_COL, 0, col, m_rows, m_cols);
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
 *                               Get matrix fill                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GSparseMatrix::fill() const
{
  // Determine number of elements in matrix
  int num = (m_fill_val == 0.0) ? m_elements : m_elements + 1;

  // Return fill
  return (double(num) / double(m_rows*m_cols));
}


/***************************************************************************
 *                      Insert vector column into matrix                   *
 * ----------------------------------------------------------------------- *
 * This is the main driver routine to insert data into a matrix. Note that *
 * there is another instance of this function that takes a compressed      *
 * array.                                                                  *
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
    throw GException::GException::out_of_range(G_INSERT_COL, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_rows != v.m_num)
    throw GException::matrix_vector_mismatch(G_INSERT_COL, v.m_num, m_rows, m_cols);

  // If there is a pending element for this column then delete it since
  // the vector overwrites this element
  if (m_fill_val != 0.0 && m_fill_col == col) {
    #if defined(G_DEBUG_SPARSE_PENDING)
    cout << G_INSERT_COL << ": pending value " << m_fill_val << 
         " for location (" << m_fill_row << "," << m_fill_col << 
         ") became obsolete" << endl;
    #endif
    m_fill_val = 0.0;
    m_fill_row = 0;
    m_fill_col = 0;
  }

  // Determine the number of non-zero elements in the vector
  int n_vector = 0;
  for (int row = 0; row < m_rows; ++row) {
    if (v.m_data[row] != 0.0)
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

  // If we need space then allocate it, if we have to much space then free it
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
      if (v.m_data[row] != 0.0) {
        m_data[i]   = v.m_data[row];
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
 *                    Insert compressed array into matrix                  *
 * ----------------------------------------------------------------------- *
 * This is the main driver routine to insert data into a matrix. Note that *
 * there is another instance of this function that takes a vector.         *
 ***************************************************************************/
void GSparseMatrix::insert_col(const double* values, const int* rows, 
                               int number, int col)
{
  // Debug header
  #if defined(G_DEBUG_SPARSE_INSERTION)
  cout << "GSparseMatrix::insert_col(v, i, n, " << col << "):" << endl;
  cout << " Matrix Data : ";
  for (int i = 0; i < m_elements; ++i)
    cout << m_data[i] << " ";
  cout << endl << " Matrix Row .: ";
  for (int i = 0; i < m_elements; ++i)
    cout << m_rowinx[i] << " ";
  cout << endl << " Matrix Col .: ";
  for (int i = 0; i < m_cols+1; ++i)
    cout << m_colstart[i] << " ";
  cout << endl;
  cout << " Array Data .: ";
  for (int i = 0; i < number; ++i)
    cout << values[i] << " ";
  cout << endl << " Array Row ..: ";
  for (int i = 0; i < number; ++i)
    cout << rows[i] << " ";
  cout << endl;
  #endif

  // Raise an exception if the column index is invalid
  #if defined(G_RANGE_CHECK)
  if (col >= m_cols)
    throw GException::out_of_range(G_INSERT_COL2, 0, col, m_rows, m_cols);
  #endif

  // Raise an exception if the index array seems incompatible with matrix 
  // dimensions
  if (rows[number-1] >= m_rows)
    throw GException::matrix_vector_mismatch(G_INSERT_COL2, rows[number-1], m_rows, m_cols);

  // If there is a pending element for this column then delete it since
  // the vector overwrites this element
  if (m_fill_val != 0.0 && m_fill_col == col) {
    #if defined(G_DEBUG_SPARSE_PENDING)
    cout << G_INSERT_COL2 << ": pending value " << m_fill_val << 
         " for location (" << m_fill_row << "," << m_fill_col << 
         ") became obsolete" << endl;
    #endif
    m_fill_val = 0.0;
    m_fill_row = 0;
    m_fill_col = 0;
  }

  // Get the start and stop indices of the actual column and compute
  // the number of exisiting elements in the column
  int i_start = m_colstart[col];
  int i_stop  = m_colstart[col+1];
  int n_exist = i_stop - i_start;

  // If the array is empty then make sure that the number of elements is 0 
  // (we then just delete the existing column)
  if (!values || !rows)
    number = 0;

  // Compute the size difference for the new matrix. It is positive if
  // the number of non-zero entries in the array is larger than the
  // number of non-zero entries in the matrix (in this case we have to
  // increase the matrix size).
  int n_diff = number - n_exist;

  // If we need space then allocate it, if we have to much space then free it
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

  // Insert the array elements into the matrix
  if (number > 0) {
    for (int row = 0, i = i_start; row < number; ++row, ++i) {
      m_data[i]   = values[row];
      m_rowinx[i] = rows[row];
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
 *                        Get minimum matrix element                       *
 * ----------------------------------------------------------------------- *
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
 * ----------------------------------------------------------------------- *
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
 *                     Initialises matrix filling stack                    *
 * ----------------------------------------------------------------------- *
 * The matrix filling stack is used to allow for a fast column-wise        *
 * filling of a sparse matrix. Columns are successively appended to a      *
 * stack which is regularily flushed when it is full. This reduces memory  *
 * copies and movements and increases filling speed.                       *
 ***************************************************************************/
void GSparseMatrix::stack_init(int size, int entries)
{
  // Free exisiting stack
  free_stack_members();

  // Initialise stack members
  init_stack_members();
  m_stack_max_entries = (entries > 0) ? entries : m_cols;
  m_stack_size        = (size    > 0) ? size    : G_SPARSE_MATRIX_DEFAULT_STACK_SIZE;

  // Allocate stack memory. Raise an exception if allocation fails
  m_stack_colinx = new int[m_stack_max_entries];
  m_stack_start  = new int[m_stack_max_entries+1];
  m_stack_data   = new double[m_stack_size];
  m_stack_rowinx = new int[m_stack_size];
  m_stack_work   = new int[m_cols];
  m_stack_buffer = new double[m_cols];
  if (m_stack_colinx == NULL)
    throw GException::mem_alloc(G_STACK_INIT, m_stack_max_entries);
  if (m_stack_start == NULL)
    throw GException::mem_alloc(G_STACK_INIT, m_stack_max_entries+1);
  if (m_stack_data == NULL || m_stack_rowinx == NULL)
    throw GException::mem_alloc(G_STACK_INIT, m_stack_size);
  if (m_stack_work == NULL || m_stack_buffer == NULL)
    throw GException::mem_alloc(G_STACK_INIT, m_cols);

  // Initialise next free stack location to the first stack element
  m_stack_start[0] = 0;

  // Return
  return;
}


/***************************************************************************
 *                Push a compressed array on the matrix stack              *
 * ----------------------------------------------------------------------- *
 * The method puts new data on the stack while assuring that each column   *
 * is present only once. If an already existing column is encountered      *
 * the data from the existing and new column are mixed and put into a new  *
 * entry one the stack; the old entry is signalled as obsolete.            *
 * On return the method indicates the number of elements in the input      *
 * array that have not been put onto the stack (due to memory limits).     *
 * This value can be either '0' (i.e. all elements have been put on the    *
 * stack) or 'number' (i.e. none of the elements have been put on the      *
 * stack). Columns are not partially put onto the stack.                   *
 ***************************************************************************/
int GSparseMatrix::stack_push_column(const double* values, const int* rows,
                                     int number, int col)
{
  // Initialise return value
  int remaining = number; 

  // Single loop for common exit point
  do {

    // If the array is empty there is nothing to do
    if (!values || !rows || (number < 1))
      continue;

    // If there is no stack or the stack can not hold the requested number of 
    // elements then report number of array elements to the caller
    if (m_stack_data == NULL || number > m_stack_size)
      continue;

    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
      throw GException::out_of_range(G_STACK_PUSH, 0, col, m_rows, m_cols);
    #endif

    // Raise an exception if the index array seems incompatible with matrix 
    // dimensions
    if (rows[number-1] >= m_rows)
      throw GException::matrix_vector_mismatch(G_STACK_PUSH, rows[number-1], m_rows, m_cols);

    // Debug header
    #if defined(G_DEBUG_SPARSE_STACK_PUSH)
    cout << "GSparseMatrix::stack_push_column(v, i, n, col=" << col << ")" << endl;
    cout << " Data to push on stack ...:";
    for (int i = 0; i < number; ++i)
      cout << " " << values[i];
    cout << endl;
    cout << " Row indices of data .....:";
    for (int i = 0; i < number; ++i)
      cout << " " << rows[i];
    cout << endl;
    #endif

    // If the stack is full then flush it before pushing the array onto it.
    // There are 2 conditions that may fill the stack: all entries are
    // occupied or there is not enough space to hold more elements
    if ((m_stack_entries >= m_stack_max_entries) ||
        (number          >= (m_stack_size - m_stack_start[m_stack_entries])))
      stack_flush();

    // If the specified column is already on the stack then mix new column with
    // old one and invalidate old column (by setting its column index to -1).
    // Since we loop here over all stack entries we don't want to have too
    // many entries in the stack. So don't set the maximum number of stack
    // entries too large ...
    for (int entry = 0; entry < m_stack_entries; ++entry) {
      if (col == m_stack_colinx[entry]) {

        // Get start and stop indices of existing column in stack
        int i_start = m_stack_start[entry];
        int i_stop  = m_stack_start[entry+1];

        // Allocate variables to hold mixing results
        int num_1;    // Number of elements that are only in stack column
        int num_2;    // Number of elements that are only in new column
        int num_mix;  // Number of elements that are in both columns

        // Determine how many elements are requested to hold the combined
        // column
        mix_column_prepare(&(m_stack_rowinx[i_start]), i_stop-i_start,
                           rows, number, &num_1, &num_2, &num_mix);
        int num_request = num_1 + num_2 + num_mix;

        // If there is not enough space in the stack to hold the combined
        // column we flush the stack and exit the loop now. In this case the
        // new column is put as the first column in the empty (flushed stack)
        if (num_request >= (m_stack_size - m_stack_start[m_stack_entries])) {
          stack_flush();
          break;
        }

        // There is enough space, so combine both columns into a fresh one
        int inx = m_stack_start[m_stack_entries];
        int num_new;
        mix_column(&(m_stack_data[i_start]), &(m_stack_rowinx[i_start]), i_stop-i_start,
                   values, rows, number,
                   &(m_stack_data[inx]), &(m_stack_rowinx[inx]), &num_new);

        // Store entry information and initialise start of next entry
        m_stack_colinx[m_stack_entries] = col;           // Store column index for entry
        m_stack_entries++;                               // Increase the # of entries
        m_stack_start[m_stack_entries] = inx + num_new;  // Set start pointer for next entry

        // Invalidate old column
        m_stack_colinx[entry] = -1;

        // Fall through to end
        remaining = 0;
        break;

      } // endif: we found an existing column
    } // endfor: looped over all columns

    // If column has already been inserted then fall through ...
    if (remaining == 0)
      continue;

    // Otherwise, push the array on the stack
    int inx = m_stack_start[m_stack_entries];
    for (int i = 0; i < number; ++i) {
      m_stack_data[inx]   = values[i];
      m_stack_rowinx[inx] = rows[i];
      inx++;
    }

    // Store entry information and initialise start of next entry
    m_stack_colinx[m_stack_entries] = col;     // Store column index for entry
    m_stack_entries++;                         // Increase the # of entries
    m_stack_start[m_stack_entries] = inx;      // Set start pointer for next entry

    // Signal success
    remaining = 0;

  } while (0); // End of main do-loop

  // Debug: show stack information  
  #if defined(G_DEBUG_SPARSE_STACK_PUSH)
  cout << " Number of stack entries .: " << m_stack_entries << endl;
  cout << " Stack entry columns .....:";
  for (int i = 0; i < m_stack_entries; ++i)
    cout << " " << m_stack_colinx[i];
  cout << endl;
  cout << " Stack entry starts ......:";
  for (int i = 0; i < m_stack_entries; ++i)
    cout << " " << m_stack_start[i];
  cout << " (next at " << m_stack_start[m_stack_entries] << ")" << endl;
  cout << " Stack data ..............:";
  for (int i = 0; i < m_stack_start[m_stack_entries]; ++i)
    cout << " " << m_stack_data[i];
  cout << endl;
  cout << " Stack rows ..............:";
  for (int i = 0; i < m_stack_start[m_stack_entries]; ++i)
    cout << " " << m_stack_rowinx[i];
  cout << endl;
  #endif

  // Return remaining number of elements
  return remaining;
}


/***************************************************************************
 *                 Push a vector column on the matrix stack                *
 * ----------------------------------------------------------------------- *
 * This method is identical to the precedent one, but it takes a full      *
 * vector instead of a compressed array. Internal working arrays are used  *
 * to convert the full column vector in a compressed array and to hand it  *
 * over to the compressed array version.                                   *
 ***************************************************************************/
int GSparseMatrix::stack_push_column(const GVector& v, int col)
{
  // Initialise number of non-zero elements
  int non_zero = 0;

  // Compress vector in the buffer and set-up row index array
  for (int i = 0; i < v.m_num; ++i) {
    if (v.m_data[i] != 0.0) {
      m_stack_buffer[non_zero] = v.m_data[i];
      m_stack_work[non_zero]   = i;
      non_zero++;
    }
  }

  // Hand over to compressed array method
  int remaining = stack_push_column(m_stack_buffer, m_stack_work, non_zero, col);

  // Return remaining number of elements
  return remaining;
}


/***************************************************************************
 *                              Flush matrix stack                         *
 * ----------------------------------------------------------------------- *
 * Adds the stack to the actual matrix. First, the total number of matrix  *
 * elements is determined. The new memory is allocated to hold all         *
 * elements. Finally, all elements are filled into a working array that is *
 * then compressed into the matrix.                                        *
 * ----------------------------------------------------------------------- *
 * NOTE: The present algorithm assumes that each column occurs only once   *
 * in the stack!                                                           *
 ***************************************************************************/
void GSparseMatrix::stack_flush(void)
{
  // Do nothing if there is no stack
  if (m_stack_data == NULL)
    return;

  // Fill pending value
  fill_pending();

  // Debug header
  #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
  cout << "GSparseMatrix::stack_flush" << endl;
  cout << " Number of stack entries .: " << m_stack_entries << endl;
  cout << " Number of stack elements : " << m_stack_start[m_stack_entries] << endl;
  cout << " Number of matrix elements: " << m_elements << endl;
  #endif

  // Use stack working array to flag all columns that exist already in the
  // matrix by 1 and all non-existing columns by 0
  for (int col = 0; col < m_cols; ++col) {
    if (m_colstart[col] < m_colstart[col+1])
      m_stack_work[col] = 1;
    else
      m_stack_work[col] = 0;
  }

  // Initialise some element counters
  int new_elements = 0;
  #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
  int num_matrix   = 0;    // Number of elements only in matrix
  int num_stack    = 0;    // Number of elements only on stack
  int num_both     = 0;    // Number of elements in matrix and on stack
  #endif

  // Loop over all stack entries and gather the number of elements that are
  // new with respect to the initial matrix. For each column set a flag with
  // the following meanings:
  //  -(entry+2): Column exists and is also found in stack=entry
  //  +(entry+2): Column exists in stack=entry only
  //           1: Column exists in matrix only
  //           0: Column does neither exist in matrix or stack  
  for (int entry = 0; entry < m_stack_entries; ++entry) {

    // Get column for actual entry
    int col = m_stack_colinx[entry];

    // Consider only valid entries
    if (col >= 0) {

      // If column exists already in matrix then determine total number
      // of additional elements
      if (m_stack_work[col] == 1) {

        // Flag that this column is a mixed column
        m_stack_work[col] = -(entry+2);

        // Setup index boundaries
        int i_start = m_colstart[col];
        int i_stop  = m_colstart[col+1];
        int k_start = m_stack_start[entry];
        int k_stop  = m_stack_start[entry+1];

        // Allocate output counters
        int num_1;
        int num_2;
        int num_mix;

        // Analyse column mixing
        mix_column_prepare(&(m_rowinx[i_start]), i_stop-i_start,
                           &(m_stack_rowinx[k_start]), k_stop-k_start,
                           &num_1, &num_2, &num_mix);

        // Update element counters
        new_elements += num_2;
        #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
        num_matrix   += num_1;
        num_stack    += num_2;
        num_both     += num_mix;
        #endif

      } // endif: column existed in the matrix

      // If column did not exists in the matrix then consider all elements
      // as new
      else {
        m_stack_work[col]  = (entry+2);
        new_elements      += (m_stack_start[entry+1] - m_stack_start[entry]);
      }
    } // endif: entry was valid
  } // endfor: looped over all entries
  int elements = m_elements + new_elements;

  // Dump number of elements in new matrix 
  #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
  cout << " New elements ............: " << new_elements << endl;
  #endif

  // Allocate memory for new matrix (always keep some elbow room)
  m_alloc = elements + m_mem_block;
  double* new_data   = new double[m_alloc];
  int*    new_rowinx = new int[m_alloc];
  if (new_data == NULL || new_rowinx == NULL)
    throw GException::mem_alloc(G_STACK_FLUSH, m_alloc);

  // Fill new matrix. For this purpose we loop over all matrix columns
  // and perform the operation that was identified in the previous scan
  int index = 0;
  for (int col = 0; col < m_cols; ++col) {

    // If column does not exist then skip
    if (m_stack_work[col] == 0) {
      m_colstart[col] = index;
      continue;
    }

    // If column exists in the matrix but not in the stack we copy the
    // existing column into the new matrix
    if (m_stack_work[col] == 1) {
      for (int i = m_colstart[col]; i < m_colstart[col+1]; ++i) {
        new_data[index]   = m_data[i];
        new_rowinx[index] = m_rowinx[i];
        index++;
        #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
        num_matrix++;
        #endif
      }
    }

    // If column is new (i.e. it did not exist in the matrix before) we 
    // copy the stack into the new matrix
    else if (m_stack_work[col] > 1) {
      int entry = m_stack_work[col] - 2;
      for (int i = m_stack_start[entry]; i < m_stack_start[entry+1]; ++i) {
        new_data[index]   = m_stack_data[i];
        new_rowinx[index] = m_stack_rowinx[i];
        index++;
        #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
        num_stack++;
        #endif
      }
    }

    // If column exists in the matrix and in the stack then we have to mix both
    else {

      // Get stack entry for column mix
      int entry = -(m_stack_work[col] + 2);

      // Setup index boundaries
      int i_start = m_colstart[col];
      int i_stop  = m_colstart[col+1];
      int k_start = m_stack_start[entry];
      int k_stop  = m_stack_start[entry+1];

      // Allocate output counter
      int num;

      // Perform column mixing
      mix_column(&(m_data[i_start]), &(m_rowinx[i_start]), i_stop-i_start,
                 &(m_stack_data[k_start]), &(m_stack_rowinx[k_start]), k_stop-k_start,
                 &(new_data[index]), &(new_rowinx[index]), &num);

      // Increment element index
      index += num;

    } // endelse: column mixing required

    // Store actual index in column start array
    m_colstart[col] = index;

  } // endfor: looped over all columns

  // Dump number of elements in new matrix after addition
  #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
  cout << " Added elements ..........: " << index << " (should be " << elements << ")" << endl;
  cout << " - Matrix only ...........: " << num_matrix << endl;
  cout << " - Stack only ............: " << num_stack << endl;
  cout << " - Matrix & Stack ........: " << num_both << endl;
  #endif

  // Correct columns start array
  for (int col = m_cols; col > 0; --col)
    m_colstart[col] = m_colstart[col-1];
  m_colstart[0] = 0;

  // Delete old matrix memory
  if (m_data   != NULL) delete [] m_data;
  if (m_rowinx != NULL) delete [] m_rowinx;

  // Update pointers to new memory and update element counter
  m_data     = new_data;
  m_rowinx   = new_rowinx;
  m_elements = elements;

  // Stack is empty now, so reset stack counters
  m_stack_entries  = 0;
  m_stack_start[0] = 0;

  // Return
  return;
}


/***************************************************************************
 *                             Destroy matrix stack                        *
 * ----------------------------------------------------------------------- *
 * Flush and destroy matrix stack                                          *
 ***************************************************************************/
void GSparseMatrix::stack_destroy(void)
{
  // Flush stack first
  stack_flush();

  // Free stack members
  free_stack_members();

  // Initialise stack (no entries)
  init_stack_members();

  // Return
  return;
}


/***************************************************************************
 *                       Compute sum of non-zero elements                  *
 * ----------------------------------------------------------------------- *
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


/*==========================================================================
 =                                                                         =
 =                     GSparseMatrix private functions                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Constructor method                           *
 * ----------------------------------------------------------------------- *
 * This is the main constructor code, without any initialisation in it.    *
 * It is used as service function to a number of different constructors of *
 * the GSparseMatrix class.                                                *
 ***************************************************************************/
void GSparseMatrix::constructor(int rows, int cols, int elements)
{
  // Allocate column start array. This is the only array that we can
  // allocate at this time. The other arrays can only be allocated during
  // filling of the matrix
  m_colstart = new int[cols+1];
  if (m_colstart == NULL)
    throw GException::mem_alloc(G_CONSTRUCTOR, cols+1);

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
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::init_members(void)
{
  // Initialise sparse matrix members
  m_rowinx    = NULL;
  m_mem_block = G_SPARSE_MATRIX_DEFAULT_MEM_BLOCK;
  m_zero      = 0.0;
  m_fill_val  = 0.0;
  m_fill_row  = 0;
  m_fill_col  = 0;
  m_symbolic  = NULL;
  m_numeric   = NULL;

  // Initialise stack members
  init_stack_members();

  // Return
  return;
}


/***************************************************************************
 *                          Initialise fill-stack                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::init_stack_members(void)
{
  // Initialise stack members
  m_stack_max_entries = 0;
  m_stack_size        = 0;
  m_stack_entries     = 0;
  m_stack_colinx      = NULL;
  m_stack_start       = NULL;
  m_stack_data        = NULL;
  m_stack_rowinx      = NULL;
  m_stack_work        = NULL;
  m_stack_buffer      = NULL;

  // Return
  return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::copy_members(const GSparseMatrix& m)
{
  // Copy GSparseMatrix members
  m_mem_block = m.m_mem_block;
  m_zero      = m.m_zero;
  m_fill_val  = m.m_fill_val;
  m_fill_row  = m.m_fill_row;
  m_fill_col  = m.m_fill_col;

  // Copy row indices
  if (m.m_rowinx != NULL && m.m_alloc > 0) {
    m_rowinx = new int[m.m_alloc];
    if (m_rowinx == NULL)
      throw GException::mem_alloc(G_COPY_MEMBERS, m.m_alloc);
    for (int i = 0; i < m.m_elements; ++i)
      m_rowinx[i] = m.m_rowinx[i];
  }

  // Copy symbolic decomposition
  if (m.m_symbolic != NULL) {
    GSparseSymbolic* symbolic = new GSparseSymbolic();
    if (symbolic == NULL)
      throw GException::mem_alloc(G_COPY_MEMBERS, 1);
    *symbolic  = *((GSparseSymbolic*)m.m_symbolic);
    m_symbolic = (void*)symbolic;
  }

  // Copy numeric decomposition
  if (m.m_numeric != NULL) {
    GSparseNumeric* numeric = new GSparseNumeric();
    if (numeric == NULL)
      throw GException::mem_alloc(G_COPY_MEMBERS, 1);
    *numeric  = *((GSparseNumeric*)m.m_numeric);
    m_numeric = (void*)numeric;
  }

  // Return
  return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::free_members(void)
{
  // De-allocate only if memory has indeed been allocated
  if (m_numeric  != NULL) delete (GSparseNumeric*)m_numeric;
  if (m_symbolic != NULL) delete (GSparseSymbolic*)m_symbolic;
  if (m_rowinx   != NULL) delete [] m_rowinx;

  // Properly mark members as free
  m_rowinx   = NULL;
  m_symbolic = NULL;
  m_numeric  = NULL;

  // Free stack members
  free_stack_members();

  // Return
  return;
}


/***************************************************************************
 *                      Delete fill-stack class members                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GSparseMatrix::free_stack_members(void)
{
  // Free stack members
  if (m_stack_colinx != NULL) delete [] m_stack_colinx;
  if (m_stack_start  != NULL) delete [] m_stack_start;
  if (m_stack_data   != NULL) delete [] m_stack_data;
  if (m_stack_rowinx != NULL) delete [] m_stack_rowinx;
  if (m_stack_work   != NULL) delete [] m_stack_work;
  if (m_stack_buffer != NULL) delete [] m_stack_buffer;

  // Properly mark members as free
  m_stack_colinx = NULL;
  m_stack_start  = NULL;
  m_stack_data   = NULL;
  m_stack_rowinx = NULL;
  m_stack_work   = NULL;
  m_stack_buffer = NULL;

  // Return
  return;
}


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
  // Dump header
  #if defined(G_DEBUG_SPARSE_MALLOC)
  cout << "GSparseMatrix::alloc_elements(start=" << start << ", num=" << 
          num << ")" << endl;
  cout << " Before allocation : " << m_elements << " " << m_alloc << endl;
  #endif

  // Continue only if we need memory
  if (num > 0) {

    // If start is after the end then append memory
    if (start > m_elements)
      start = m_elements;

    // Determine the requested new logical size of the matrix
    int new_size = m_elements + num;

    // Case A: the requested memory is already available, so just move the
    // data to make space for new elements and initialise the new cells
    if (new_size <= m_alloc) {

      // Move up all elements after index to insert
      #if defined(G_USE_MEMMOVE)
      int n_copy = m_elements - start;
      memmove(&(m_data[start+num]), &(m_data[start]), sizeof(double)*n_copy);
      memmove(&(m_rowinx[start+num]), &(m_rowinx[start]), sizeof(long)*n_copy);
      #else
      for (int i = m_elements - 1; i >= start; --i) {
        m_data[i+num]   = m_data[i];
        m_rowinx[i+num] = m_rowinx[i];
      }
      #endif

      // Clear new elements (zero content for row 0)
      for (int i = start; i < start+num; ++i) {
        m_data[i]   = 0.0;
        m_rowinx[i] = 0;
      }

      // Update element counter
      m_elements += num;

    } // endif: Case A: memory already existed

    // Case B: more memory is needed, so allocate it, copy the existing
    // content, and initialise new cells
    else {

      // Make sure that enough memory is allocated
      int new_propose = m_alloc + m_mem_block;
      m_alloc = (new_size > new_propose) ? new_size : new_propose;

      // Allocate memory for new elements
      double* new_data   = new double[m_alloc];
      int*    new_rowinx = new int[m_alloc];
      if (new_data == NULL || new_rowinx == NULL)
        throw GException::mem_alloc(G_ALLOC, m_alloc);

      // Copy all elements before index to insert
      #if defined(G_USE_MEMCPY)
      int n_copy = start;
      memcpy(new_data,   m_data,   sizeof(double)*n_copy);
      memcpy(new_rowinx, m_rowinx, sizeof(long)*n_copy);
      #else
      for (int i = 0; i < start; ++i) {
        new_data[i]   = m_data[i];
        new_rowinx[i] = m_rowinx[i];
      }
      #endif

      // Clear new elements (zero content for row 0)
      for (int i = start; i < start+num; ++i) {
        new_data[i]   = 0.0;
        new_rowinx[i] = 0;
      }

      // Copy all elements after index to insert
      #if defined(G_USE_MEMCPY)
      n_copy = m_elements - start;
      memcpy(&(new_data[start+num]),   &(m_data[start]),   sizeof(double)*n_copy);
      memcpy(&(new_rowinx[start+num]), &(m_rowinx[start]), sizeof(long)*n_copy);
      #else
      for (int i = start; i < m_elements; ++i) {
        new_data[i+num]   = m_data[i];
        new_rowinx[i+num] = m_rowinx[i];
      }
      #endif

      // Delete old memory
      if (m_data   != NULL) delete [] m_data;
      if (m_rowinx != NULL) delete [] m_rowinx;

      // Update pointers to new memory and update element counter
      m_data      = new_data;
      m_rowinx    = new_rowinx;
      m_elements += num;

    } // endelse: Case B: more memory was needed

  } // endif: needed new memory

  // Dump new memory size
  #if defined(G_DEBUG_SPARSE_MALLOC)
  cout << " After allocation .: " << m_elements << " " << m_alloc << endl;
  #endif

  // Return
  return;
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
  // Dump header
  #if defined(G_DEBUG_SPARSE_MALLOC)
  cout << "GSparseMatrix::free_elements(start=" << start << ", num=" << 
          num << ")" << endl;
  #endif

  // Continue only if we need to free memory and if start is within the
  // range
  if (num > 0 && start < m_elements) {

    // Determine the requested new logical size of the matrix
    int new_size = m_elements - num;

    // If there are no elements then simply delete all matrix elements ...
    if (new_size < 1) {
      if (m_data   != NULL) delete [] m_data;
      if (m_rowinx != NULL) delete [] m_rowinx;
        m_data     = NULL;
        m_rowinx   = NULL;
        m_elements = 0;
        m_alloc    = 0;
     }

     // ... otherwise shrink the array
     else {

       // Case A: If at least one entire memory block has been liberated then 
       // physically shrink the matrix array
       if (m_alloc - new_size > m_mem_block) {

         // Shrink, but leave a memory block for possible future filling
         m_alloc = new_size + m_mem_block;

         // Allocate new memory
         double* new_data   = new double[m_alloc];
         int*    new_rowinx = new int[m_alloc];
         if (new_data == NULL || new_rowinx == NULL)
           throw GException::mem_alloc(G_FREE, m_alloc);

         // Copy all elements before the starting index
         #if defined(G_USE_MEMCPY)
         int n_copy = start;
         memcpy(new_data,   m_data,   sizeof(double)*n_copy);
         memcpy(new_rowinx, m_rowinx, sizeof(long)*n_copy);
         #else
         for (int i = 0; i < start; ++i) {
           new_data[i]   = m_data[i];
           new_rowinx[i] = m_rowinx[i];
         }
         #endif

         // Copy all elements after the starting index
         #if defined(G_USE_MEMCPY)
         n_copy = new_size - start;
         memcpy(&(new_data[start]),   &(m_data[start+num]),   sizeof(double)*n_copy);
         memcpy(&(new_rowinx[start]), &(m_rowinx[start+num]), sizeof(long)*n_copy);
         #else
         for (int i = start; i < new_size; ++i) {
           new_data[i]   = m_data[i+num];
           new_rowinx[i] = m_rowinx[i+num];
         }
         #endif

        // Delete old memory
        if (m_data   != NULL) delete [] m_data;
        if (m_rowinx != NULL) delete [] m_rowinx;

        // Update pointers to new memory and update element counter
        m_data     = new_data;
        m_rowinx   = new_rowinx;
        m_elements = new_size;

      } // endif: Case A: memory shrinkage performed

      // Case B: we keep the memory and just move the elements
      else {

        // Move all elements after the starting index
        #if defined(G_USE_MEMMOVE)
        int n_copy = new_size - start;
        memmove(&(m_data[start]),   &(m_data[start+num]),   sizeof(double)*n_copy);
        memmove(&(m_rowinx[start]), &(m_rowinx[start+num]), sizeof(long)*n_copy);
        #else
        for (int i = start; i < new_size; ++i) {
          m_data[i]   = m_data[i+num];
          m_rowinx[i] = m_rowinx[i+num];
        }
        #endif

        // Update element counter
        m_elements = new_size;

      } // endelse: Case B

    } // endif: array shrinkage needed
  } // endif: needed new memory

  // Return
  return;
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
    throw GException::matrix_zero(G_REMOVE_ZERO);

  // Allocate row mapping array
  int* row_map = new int[m_rows];
  if (row_map == NULL)
    throw GException::mem_alloc(G_REMOVE_ZERO, m_rows);

  // Setup row mapping array that maps original matrix rows into compressed
  // matrix rows. An entry of -1 indicates that the row should be dropped
  for (int row = 0; row < m_rows; ++row)
    row_map[row] = -1;
  for (int c_row = 0; c_row < m_num_rowsel; ++c_row)
    row_map[m_rowsel[c_row]] = c_row;

  // Initialise pointers to compressed array
  double* d_data   = m_data;
  int*    d_rowinx = m_rowinx;

  // Initialise column start of first column to zero
  m_colstart[0] = 0;

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " before compression (" << m_colstart[m_cols] << "):";
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
  cout << " after compression (" << m_colstart[m_num_colsel] << ") :";
  for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
    cout << " " << m_colstart[c_col];
  cout << endl;
  #endif

  // Free row mapping array
  delete [] row_map;

  // Update matrix size
  m_rows = m_num_rowsel;
  m_cols = m_num_colsel;

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
  cout << " before restoration (" << m_colstart[m_num_colsel] << "):";
  for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
    cout << " " << m_colstart[c_col];
  cout << endl;
  #endif

  // If column selection exists then restore column counters
  if (m_colsel != NULL) {

    // Start insertion from the last original column
    int col_stop = cols - 1;

    // Loop over all compressed columns
    for (int c_col = m_num_colsel-1; c_col > 0; --c_col) {
      int col_start = m_colsel[c_col-1] + 1;
      int col_value = m_colstart[c_col];
      for (int col = col_start; col <= col_stop; ++col)
        m_colstart[col] = col_value;
      col_stop = col_start - 1;
    }

    // Set first columns if they are not yet set
    for (int col = 0; col <= col_stop; ++col)
      m_colstart[col] = 0;

    // Restore the number of elements
    m_colstart[cols] = m_elements;
  }

  // Dump column start array
  #if defined(G_DEBUG_SPARSE_COMPRESSION)
  cout << " after restoration (" << m_colstart[cols] << ") :";
  for (int col = 0; col <= cols; ++col)
    cout << " " << m_colstart[col];
  cout << endl;
  #endif

  // Update matrix size
  m_rows = rows;
  m_cols = cols;

  // Return
  return;
}


/***************************************************************************
 *                      Prepare mix of sparse columns                      *
 * ----------------------------------------------------------------------- *
 * Prepare the mix of two sparse matrix columns into a single one.         *
 * src1_row  : Row index array [0...src1_num-1] of first column            *
 * src1_num  : Number of elements in first column                          *
 * src2_row  : Row index array [0...src2_num-1] of second column           *
 * src2_num  : Number of elements in second column                         *
 * num_1     : Number of elements only found in column 1                   *
 * num_2     : Number of elements only found in column 2                   *
 * num_mix   : Number of elements found in both columns                    *
 ***************************************************************************/
void GSparseMatrix::mix_column_prepare(const int* src1_row, int src1_num,
                                       const int* src2_row, int src2_num,
                                       int* num_1, int* num_2, int* num_mix)
{
  // Initialise element counters
  *num_1   = 0;
  *num_2   = 0;
  *num_mix = 0;

  // Initialise indices and row indices of both columns
  int inx_1 = 0;                    // Column 1 element index
  int inx_2 = 0;                    // Column 2 element index
  int row_1 = src1_row[inx_1];      // Column 1 first row index
  int row_2 = src2_row[inx_2];      // Column 2 first row index

  // Mix elements of both columns while both contain still elements
  while (inx_1 < src1_num && inx_2 < src2_num) {

    // Case A: the element exist in both columns
    if (row_1 == row_2) {
      row_1 = src1_row[++inx_1];
      row_2 = src2_row[++inx_2];
      (*num_mix)++;
    }

    // Case B: the element exists only in first column
    else if (row_1 < row_2) {
      row_1 = src1_row[++inx_1];
      (*num_1)++;
    }

    // Case C: the element exists only in second column
    else {
      row_2 = src2_row[++inx_2];
      (*num_2)++;
    }

  } // endwhile: mixing

  // At this point either the first or the second column expired of elements
  // In the case that there are still elements remaining in the first column we
  // count them now ...
  if (inx_1 < src1_num)
    *num_1 += (src1_num - inx_1);

  // ... or in the case that there are still elements remaining in the second
  // column we count them now
  if (inx_2 < src2_num)
    *num_2 += (src2_num - inx_2);

  // We're done
  return;
}		


/***************************************************************************
 *                             Mix column                                  *
 * ----------------------------------------------------------------------- *
 * Mix two sparse matrix columns into one.                                 *
 * src1_data : Data array [0...src1_num-1] of first column                 *
 * src1_row  : Row index array [0...src1_num-1] of first column            *
 * src1_num  : Number of elements in first column                          *
 * src2_data : Data array [0...src2_num-1] of second column                *
 * src2_row  : Row index array [0...src2_num-1] of second column           *
 * src2_num  : Number of elements in second column                         *
 * dst_data  : Data array [0...dst_num-1] of result column                 *
 * dst_row   : Row index array [0...dst_num-1] of result column            *
 * dst_num   : Number of elements in result column                         *
 ***************************************************************************/
void GSparseMatrix::mix_column(const double* src1_data, const int* src1_row,
                               int src1_num,
                               const double* src2_data, const int* src2_row,
                               int src2_num,
                               double* dst_data, int* dst_row, int* dst_num)
{
  // Initialise indices and row indices of both columns
  int inx_1 = 0;                    // Column 1 element index
  int inx_2 = 0;                    // Column 2 element index
  int inx   = 0;                    // Result column element index
  int row_1 = src1_row[inx_1];
  int row_2 = src2_row[inx_2];

  // Mix elements of both columns while both contain still elements
  while (inx_1 < src1_num && inx_2 < src2_num) {

    // Case A: the element exists in both columns, so we add up the values
    if (row_1 == row_2) {
      dst_data[inx] = src1_data[inx_1] + src2_data[inx_2];
      dst_row[inx]  = row_1;
      row_1         = src1_row[++inx_1];
      row_2         = src2_row[++inx_2];
    }

    // Case B: the element exists only in first column, so we copy the element
    // from the first column
    else if (row_1 < row_2) {
      dst_data[inx] = src1_data[inx_1];
      dst_row[inx]  = row_1;
      row_1         = src1_row[++inx_1];
    }

    // Case C: the element exists only in second column, so we copy the element
    // from the second column
    else {
      dst_data[inx] = src2_data[inx_2];
      dst_row[inx]  = row_2;
      row_2         = src2_row[++inx_2];
    }

    // Update the destination index since we added a element
    inx++;

  } // endwhile: mixing

  // At this point either the first or the second column expired of elements
  // In the case that there are still elements remaining in the first column we
  // add them now ...
  for (int i = inx_1; i < src1_num; ++i, ++inx) {
    dst_data[inx] = src1_data[i];
    dst_row[inx]  = src1_row[i];
  }

  // ... or in the case that there are still elements remaining in the second
  // column we add them now
  for (int i = inx_2; i < src2_num; ++i, ++inx) {
    dst_data[inx] = src2_data[i];
    dst_row[inx]  = src2_row[i];
  }

  // Now store the number of columns in the second column
  *dst_num = inx;

  // We're done
  return;
}


/*==========================================================================
 =                                                                         =
 =                           GSparseMatrix friends                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
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
  if (m.m_fill_val == 0.0)
    os << " Number of non-zero elements: " << m.m_colstart[m.m_cols] << " (" <<
          m.m_elements << ")" << endl;
  else {
    os << " Number of non-zero elements: " << m.m_colstart[m.m_cols]+1 << " (" <<
          m.m_elements+1 << ")" << endl;
    os << " Pending element ...........: (" << m.m_fill_row << "," << m.m_fill_col <<
          ")=" << m.m_fill_val << endl;
  }
  os << " Number of allocated cells .: " << m.m_alloc << endl;
  os << " Memory block size .........: " << m.m_mem_block << endl;
  os << " Sparse matrix fill ........: " << m.fill() << endl;

  // Dump elements and compression schemes
  m.dump_elements(os);
  m.dump_row_comp(os);
  m.dump_col_comp(os);

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
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                                cs_symperm                               *
 * ----------------------------------------------------------------------- *
 * C = A(p,p) where A and C are symmetric the upper part stored.           *
 * ----------------------------------------------------------------------- *
 * Input:   A             Sparse matrix                                    *
 *          pinv          pinv[0..n-1]                                     *
 * Output:  C             Sparse matrix                                    *
 ***************************************************************************/
GSparseMatrix cs_symperm(const GSparseMatrix& m, const int* pinv)
{
  // Declare loop variables
  int i, j, p, q, i2, j2;

  // Assign matrix attributes
  int     n  = m.m_cols;
  int*    Ap = m.m_colstart;
  int*    Ai = m.m_rowinx; 
  double* Ax = m.m_data;

  // Allocate result matrix
  GSparseMatrix C(n, n, Ap[n]);

  // Allocate and initialise workspace
  int  wrk_size = n;
  int* wrk_int  = new int[wrk_size];
  if (wrk_int == NULL)
    throw GException::mem_alloc(G_SYMPERM, wrk_size);
  for (i = 0; i < wrk_size; ++i)
    wrk_int[i] = 0;

  // Assign result matrix attributes
  int*    Cp = C.m_colstart;
  int*    Ci = C.m_rowinx; 
  double* Cx = C.m_data;

  // Count entries in each column of C
  for (j = 0; j < n; j++) {

    // Column j of A is column j2 of C
    j2 = pinv ? pinv[j] : j;

    // Loop over entries in column j
    for (p = Ap[j]; p < Ap[j+1]; p++) {
      i = Ai [p];
      if (i > j) continue;              // skip lower triangular part of A
      i2 = pinv ? pinv[i] : i;          // row i of A is row i2 of C
      wrk_int[G_MAX(i2, j2)]++;         // column count of C
    }
  }

  // Compute column pointers of C
  cs_cumsum(Cp, wrk_int, n);

  // Loop over all columns of A
  for (j = 0 ; j < n ; j++) {

    // Column j of A is column j2 of C
    j2 = pinv ? pinv[j] : j;

    // Loop over entries in column j
    for (p = Ap[j]; p < Ap[j+1]; p++) {
      i = Ai [p] ;
      if (i > j) continue;              // skip lower triangular part of A
      i2    = pinv ? pinv[i] : i;       // row i of A is row i2 of C
      Ci[q  = wrk_int[G_MAX(i2,j2)]++] = G_MIN(i2,j2);
      if (Cx) Cx[q] = Ax[p];
    }
  }

  // Free workspace
  delete [] wrk_int;
  
  // Rectify the number of elements in matrix C
  C.free_elements(Cp[n], (C.m_elements-Cp[n]));

  // Return result
  return C;

}


/***************************************************************************
 *                              cs_transpose                               *
 * ----------------------------------------------------------------------- *
 * Transpose matrix. The flag 'values' allows to avoid copying the actual  *
 * data values. This allows to perform a logical matrix transposition, as  *
 * needed by the symbolic matrix analysis class.                           *
 * ----------------------------------------------------------------------- *
 * NOTE: This routine DOES NOT support pending elements (they have to be   *
 * filled before.                                                          *
 ***************************************************************************/
GSparseMatrix cs_transpose(const GSparseMatrix& m, int values)
{
  // Declare and allocate result matrix 
  GSparseMatrix result(m.m_cols, m.m_rows, m.m_elements);

  // Allocate and initialise workspace
  int  wrk_size = m.m_rows;
  int* wrk_int  = new int[wrk_size];
  if (wrk_int == NULL)
    throw GException::mem_alloc(G_TRANSPOSE, wrk_size);
  for (int i = 0; i < wrk_size; i++) 
    wrk_int[i] = 0;

  // Setup the number of non-zero elements in each row
  for (int p = 0; p < m.m_elements; p++) 
    wrk_int[m.m_rowinx[p]]++;

  // Set row pointers. To use a GSparseSymbolic function we have to
  // allocate and object (but this does not take memory)
  cs_cumsum(result.m_colstart, wrk_int, m.m_rows);

  // Case A: Normal transponse, including assignment of values
  if (values) {
    for (int col = 0; col < m.m_cols; col++) {
      for (int p = m.m_colstart[col] ; p < m.m_colstart[col+1] ; p++) {
        int i         = wrk_int[m.m_rowinx[p]]++;
        result.m_rowinx[i] = col;
        result.m_data[i]   = m.m_data[p] ;
      }
    }
  }

  // Case B: Logical transponse, no assignment of values is performed
  else {
    for (int col = 0; col < m.m_cols; col++) {
      for (int p = m.m_colstart[col] ; p < m.m_colstart[col+1] ; p++)
        result.m_rowinx[wrk_int[m.m_rowinx[p]]++] = col;
    }
  }

  // Return transponse matrix
  return result;

}


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
    nz2  += c[i];    // also in double to avoid int overflow
    c[i]  = p[i];    // also copy p[0..n-1] back into c[0..n-1]
  }
  p[n] = nz ;

  // Return cumulative sum of c[0..n-1]
  return nz2;
}
