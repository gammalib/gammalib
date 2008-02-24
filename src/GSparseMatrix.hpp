/***************************************************************************
 *                 GSparseMatrix.hpp  -  sparse matrix class               *
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

#ifndef GSPARSEMATRIX_HPP
#define GSPARSEMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GMatrixBase.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Definitions ________________________________________________________ */
#define G_SPARSE_MATRIX_DEFAULT_MEM_BLOCK      512  // Default Memory block size
#define G_SPARSE_MATRIX_DEFAULT_STACK_ENTRIES 1000  // Max # of stack entries
#define G_SPARSE_MATRIX_DEFAULT_STACK_SIZE     512  // Max # of stack elements


/***************************************************************************
 *                      GSparseMatrix class definition                     *
 * ----------------------------------------------------------------------- *
 * GMatrix implements a sparse matrix storage class. It derives from the   *
 * abstract base class GMatrixBase.                                        *
 ***************************************************************************/
class GSparseMatrix : public GMatrixBase {
  // Friend classes
  friend class GSparseSymbolic;
  friend class GSparseNumeric;

  // Binary operator friends
  friend GSparseMatrix operator* (const double& a,  const GSparseMatrix& b);
  friend GSparseMatrix operator* (const GSparseMatrix& a, const double& b);
  friend GSparseMatrix operator/ (const GSparseMatrix& a, const double& b);

  // I/O friends
  friend ostream& operator<< (ostream& os, const GSparseMatrix& m);

  // Friend functions
  friend GSparseMatrix transpose(const GSparseMatrix& m);
  friend GSparseMatrix fabs(const GSparseMatrix& m);
  friend GSparseMatrix cholesky_decompose(const GSparseMatrix& m, int compress = 1);
  friend GSparseMatrix cholesky_invert(const GSparseMatrix& m, int compress = 1);
  //
  friend GSparseMatrix cs_symperm(const GSparseMatrix& m, const int* pinv);
  friend GSparseMatrix cs_transpose(const GSparseMatrix& m, int values);

public:
  // Constructors and destructors
  GSparseMatrix(int rows, int cols, int elements = 0);
  GSparseMatrix(const GMatrix& m);
  GSparseMatrix(const GSymMatrix& m);
  GSparseMatrix(const GSparseMatrix& m);
  virtual ~GSparseMatrix();

  // Operators
  virtual GSparseMatrix& operator= (const GSparseMatrix& m);
  virtual double&        operator() (int row, int col);
  virtual const double&  operator() (int row, int col) const;
  virtual GSparseMatrix  operator+ (const GSparseMatrix& m) const;
  virtual GSparseMatrix  operator- (const GSparseMatrix& m) const;
  virtual GSparseMatrix  operator* (const GSparseMatrix& m) const;
  virtual GVector        operator* (const GVector& v) const;
  virtual int            operator== (const GSparseMatrix& m) const;
  virtual int            operator!= (const GSparseMatrix& m) const;
  virtual GSparseMatrix  operator- () const;
  virtual GSparseMatrix& operator+= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator-= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator*= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator*= (const double& d);
  virtual GSparseMatrix& operator/= (const double& d);

  // Matrix methods
  virtual void          add_col(const GVector& v, int col);
  virtual void          add_col(const double* values, const int* rows, int number, int col);
  virtual void          cholesky_decompose(int compress = 1);
  virtual GVector       cholesky_solver(const GVector& v, int compress = 1);
  virtual void          cholesky_invert(int compress = 1);
  virtual void          clear();
  virtual GVector       extract_row(int row) const;
  virtual GVector       extract_col(int col) const;
  //TBD virtual GSparseMatrix extract_lower_triangle() const;
  //TBD virtual GSparseMatrix extract_upper_triangle() const;
  virtual double        fill() const;
  virtual void          insert_col(const GVector& v, int col);
  virtual void          insert_col(const double* values, const int* rows, int number, int col);
  virtual double        min() const;
  virtual double        max() const;
  virtual void          set_mem_block(int block);
  virtual void          stack_init(int size = 0, int entries = 0);
  virtual int           stack_push_column(const GVector& v, int col);
  virtual int           stack_push_column(const double* values, const int* rows, int number, int col);
  virtual void          stack_flush(void);
  virtual void          stack_destroy(void);
  virtual double        sum() const;
  virtual void          transpose();
  
private:
  // Private methods
  void constructor(int rows, int cols, int elements = 0);
  void init_members(void);
  void init_stack_members(void);
  void copy_members(const GSparseMatrix& m);
  void free_members(void);
  void free_stack_members(void);
  int  get_index(int row, int col) const;
  void fill_pending(void);
  void alloc_elements(int start, int num);
  void free_elements(int start, int num);
  void remove_zero_row_col(void);
  void insert_zero_row_col(int rows, int cols);
  void mix_column_prepare(const int* src1_row, int src1_num,
						  const int* src2_row, int src2_num,
						  int* num_1, int* num_2, int* num_mix);
  void mix_column(const double* src1_data, const int* src1_row, int src1_num,
                  const double* src2_data, const int* src2_row, int src2_num,
				  double* dst_data, int* dst_row, int* dst_num);
  
  // Private data members
  int*    m_rowinx;             // Row-indices of all elements
  double  m_zero;               // The zero element (needed for data access)
  double  m_fill_val;           // Element to be filled
  int     m_fill_row;           // Row into which element needs to be filled
  int     m_fill_col;           // Column into which element needs to be filled
  int     m_mem_block;          // Memory block to be allocated at once
  void*   m_symbolic;           // Holds GSparseSymbolic object after decomposition
  void*   m_numeric;            // Holds GSparseNumeric object after decomposition

  // Fill-stack
  int     m_stack_max_entries;  // Maximum number of entries in the stack
  int     m_stack_size;         // Maximum number of elements in the stack
  int     m_stack_entries;      // Number of entries in the stack
  int*    m_stack_colinx;       // Column index for each entry [m_stack_entries]
  int*    m_stack_start;        // Start in stack for each entry [m_stack_entries+1]
  double* m_stack_data;         // Stack data [m_stack_size]
  int*    m_stack_rowinx;       // Stack row indices [m_stack_size]
  int*    m_stack_work;         // Stack integer working array [m_cols]
  double* m_stack_buffer;       // Stack double buffer [m_cols]
};


/***************************************************************************
 *              Inline members that override base class members            *
 ***************************************************************************/
// Matrix element access operator
inline
double& GSparseMatrix::operator() (int row, int col)
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw GException::out_of_range("GSparseMatrix::operator(int,int)", row, col, m_rows, m_cols);
  #endif
  fill_pending();
  int inx = get_index(row,col);
  double* value;
  if (inx < 0) {
    value      = &m_fill_val;
	m_fill_row = row;
	m_fill_col = col;
  }
  else
    value = &(m_data[inx]);
  return *value;
}

// Matrix element access operator (const version)
// We need here the zero element to return also a pointer for 0.0 entry that is
// not stored. Since we have the const version we don't have to care about
// modification of this zero value.
inline
const double& GSparseMatrix::operator() (int row, int col) const
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw GException::out_of_range("GSparseMatrix::operator(int,int)", row, col, m_rows, m_cols);
  #endif
  int inx = get_index(row,col);
  double* value;
  if (inx < 0)
    value = (double*)&m_zero;
  else if (inx == m_elements)
    value = (double*)&m_fill_val;
  else
    value = (double*)&(m_data[inx]);
  return *value;
}

// Binary matrix addition
inline
GSparseMatrix GSparseMatrix::operator+ (const GSparseMatrix& m) const
{
  GSparseMatrix result = *this;
  result += m;
  return result;
}

// Binary matrix subtraction
inline
GSparseMatrix GSparseMatrix::operator- (const GSparseMatrix& m) const
{
  GSparseMatrix result = *this;
  result -= m;
  return result;
}

// Binary matrix multiplication
inline
GSparseMatrix GSparseMatrix::operator* (const GSparseMatrix& m) const
{
  GSparseMatrix result = *this;
  result *= m;
  return result;
}

// Matrix scaling
inline
GSparseMatrix& GSparseMatrix::operator*= (const double& s)
{
  fill_pending();
  multiplication(s);
  return *this;
}

// Matrix scalar division
inline
GSparseMatrix& GSparseMatrix::operator/= (const double& s)
{
  double inverse = 1.0/s;
  fill_pending();
  multiplication(inverse);
  return *this;
}

// Negation
inline
GSparseMatrix GSparseMatrix::operator- ( ) const
{
  GSparseMatrix result = *this;
  result.fill_pending();
  result.negation();
  return result;
}


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
// Transpose
inline
void GSparseMatrix::transpose()
{
  fill_pending();
  *this = cs_transpose(*this, 1);
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline
GSparseMatrix operator* (const GSparseMatrix& a, const double& b)
{
  GSparseMatrix result = a;
  result.fill_pending();
  result *= b;
  return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GSparseMatrix operator* (const double& a, const GSparseMatrix& b)
{
  GSparseMatrix result = b;
  result.fill_pending();
  result *= a;
  return result;
}

// Binary matrix division (matrix is left operand)
inline
GSparseMatrix operator/ (const GSparseMatrix& a, const double& b)
{
  GSparseMatrix result = a;
  result.fill_pending();
  result /= b;
  return result;
}

// Matrix transpose function
inline
GSparseMatrix transpose(const GSparseMatrix& m)
{
  GSparseMatrix result = m;
  result.transpose();
  return result;
}

// Set memory block size
inline
void GSparseMatrix::set_mem_block(int block)
{
  m_mem_block = (block > 0) ? block : 1;
  return;
}

// Cholesky decomposition
inline
GSparseMatrix cholesky_decompose(const GSparseMatrix& m, int compress)
{
  GSparseMatrix result = m;
  result.cholesky_decompose(compress);
  return result;
}

// Matrix inversion using Cholesky decomposition
inline
GSparseMatrix cholesky_invert(const GSparseMatrix& m, int compress)
{
  GSparseMatrix result = m;
  result.cholesky_invert(compress);
  return result;
}


/***************************************************************************
 *                      Prototypes for other functions                     *
 ***************************************************************************/
double cs_cumsum(int* p, int* c, int n);

#endif /* GSPARSEMATRIX_HPP */
