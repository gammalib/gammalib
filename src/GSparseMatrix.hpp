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
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
double cs_cumsum(int* p, int* c, int n);


/***************************************************************************
 *                      GSparseMatrix class definition                     *
 ***************************************************************************/
class GSparseMatrix : public GMatrix {
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

public:
  // Constructors and destructors (not inherited)
  GSparseMatrix(int rows, int cols, int elements = 0);
  GSparseMatrix(const GSparseMatrix& m);
  virtual ~GSparseMatrix();

  // Matrix element access operators
  virtual double& operator() (int row, int col);
  virtual const double& operator() (int row, int col) const;

  // Matrix assignment operators (not inherited)
  virtual GSparseMatrix& operator= (const GSparseMatrix& m);

  // Binary operators
  virtual GSparseMatrix operator+ (const GSparseMatrix& m) const;
  virtual GSparseMatrix operator- (const GSparseMatrix& m) const;
  virtual GSparseMatrix operator* (const GSparseMatrix& m) const;
  virtual GVector       operator* (const GVector& v) const;
  virtual int           operator== (const GSparseMatrix& m) const;
  virtual int           operator!= (const GSparseMatrix& m) const;

  // Unary operators
  GSparseMatrix operator- () const;
  virtual GSparseMatrix& operator+= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator-= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator*= (const GSparseMatrix& m);
  virtual GSparseMatrix& operator*= (const double& d);
  virtual GSparseMatrix& operator/= (const double& d);

  // Matrix functions
  // USE_BASE: virtual int rows() const { return m_rows; }
  // USE_BASE: virtual int cols() const { return m_cols; }
  virtual double  fill() const { return (double(m_elements)/double(m_rows*m_cols)); }
  virtual double  min() const;
  virtual double  max() const;
  virtual double  sum() const;
  virtual void    transpose();
  virtual GVector extract_row(int row) const;
  virtual GVector extract_col(int col) const;
  virtual void    insert_col(const GVector& v, int col);
  virtual void    add_col(const GVector& v, int col);
  virtual GMatrix convert_to_full() const;
  //TBD virtual GMatrix extract_lower_triangle() const;
  //TBD virtual GMatrix extract_upper_triangle() const;
  virtual void    cholesky_decompose(int compress = 1);
  virtual GVector cholesky_solver(const GVector& v, int compress = 1);
  virtual void    cholesky_invert(int compress = 1);

private:
  // Private functions
  int           get_index(int row, int col) const;
  void          fill_pending(void);
  void          alloc_elements(int start, int num);
  void          free_elements(int start, int num);
  void          remove_zero_row_col(void);
  void          insert_zero_row_col(int rows, int cols);
  GSparseMatrix cs_transpose(GSparseMatrix *m, int values);
  
  // Private data members
  int*   m_rowinx;    // Row-indices of all elements
  double m_zero;      // The zero element (needed for data access)
  double m_fill_val;  // Element to be filled
  int    m_fill_row;  // Row into which element needs to be filled
  int    m_fill_col;  // Column into which element needs to be filled
  void*  m_symbolic;  // Holds GSparseSymbolic object after decomposition
  void*  m_numeric;   // Holds GSparseNumeric object after decomposition
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
    throw out_of_range("GSparseMatrix access", row, col, m_rows, m_cols);
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
    throw out_of_range("GSparseMatrix access", row, col, m_rows, m_cols);
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

// Unary scalar multiplication
inline
GSparseMatrix& GSparseMatrix::operator*= (const double& d)
{
  fill_pending();
  this->GMatrix::operator*=(d);
}

// Unary scalar division
inline
GSparseMatrix& GSparseMatrix::operator/= (const double& d)
{
  fill_pending();
  this->GMatrix::operator/=(d);
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

#endif /* GSPARSEMATRIX_HPP */
