/***************************************************************************
 *                 GSymMatrix.hpp  -  symmetric matrix class               *
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

#ifndef GSYMMATRIX_HPP
#define GSYMMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***************************************************************************
 *                       GSymMatrix class definition                       *
 ***************************************************************************/
class GSymMatrix : public GMatrix {
  // Binary operator friends
  friend GSymMatrix operator* (const double& a,  const GSymMatrix& b);
  friend GSymMatrix operator* (const GSymMatrix& a, const double& b);
  friend GSymMatrix operator/ (const GSymMatrix& a, const double& b);

  // I/O friends
  // USE_BASE: friend ostream& operator<< (ostream& os, const GSymMatrix& m);

  // Friend functions
  friend GSymMatrix transpose(const GSymMatrix& m);
  friend GSymMatrix fabs(const GSymMatrix& m);
  friend GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress = 1);
  friend GSymMatrix cholesky_invert(const GSymMatrix& m, int compress = 1);

public:
  // Constructors and destructors (not inherited)
  GSymMatrix(int rows, int cols);
  GSymMatrix(const GSymMatrix& m);
  virtual ~GSymMatrix();

  // Matrix element access operators
  virtual double& operator() (int row, int col);
  virtual const double& operator() (int row, int col) const;

  // Matrix assignment operators (not inherited)
  virtual GSymMatrix& operator= (const GSymMatrix& m);

  // Binary operators
  virtual GSymMatrix operator+ (const GSymMatrix& m) const;
  virtual GSymMatrix operator- (const GSymMatrix& m) const;
  virtual GSymMatrix operator* (const GSymMatrix& m) const;
  virtual GVector    operator* (const GVector& v) const;
  // USE_BASE: virtual int operator== (const GSymMatrix& m) const;
  // USE_BASE: virtual int operator!= (const GSymMatrix& m) const;

  // Unary operators
  GSymMatrix operator- () const;
  // USE_BASE: virtual GSymMatrix& operator+= (const GSymMatrix& m);
  // USE_BASE: virtual GSymMatrix& operator-= (const GSymMatrix& m);
  virtual GSymMatrix& operator*= (const GSymMatrix& m);
  virtual GSymMatrix& operator*= (const double& d);
  virtual GSymMatrix& operator/= (const double& d);

  // Matrix functions
  // USE_BASE: virtual int      rows() const { return m_rows; }
  // USE_BASE: virtual int      cols() const { return m_cols; }
  // USE_BASE: virtual void     clear();
  // USE_BASE: virtual double   min() const;
  // USE_BASE: virtual double   max() const;
  virtual double  sum() const;
  virtual void    transpose() { return; }
  virtual GVector extract_row(int row) const;
  virtual GVector extract_col(int col) const;
  virtual GMatrix convert_to_full() const;
  virtual GMatrix extract_lower_triangle() const;
  virtual GMatrix extract_upper_triangle() const;
  virtual void    cholesky_decompose(int compress = 1);
  virtual GVector cholesky_solver(const GVector& v, int compress = 1);
  virtual void    cholesky_invert(int compress = 1);
    
private:
  void set_inx(void);  // Set indices of non-zero rows/columns
  int  m_num_inx;      // Number of non-zero rows/columns
  int* m_inx;          // Indices of non-zero rows/columns
};


/***************************************************************************
 *                          Inline member functions                        *
 ***************************************************************************/
// Matrix element access operator
inline
double& GSymMatrix::operator() (int row, int col)
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GSymMatrix access", row, col, m_rows, m_cols);
  #endif
  int inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);
  return m_data[inx];
}

// Matrix element access operator (const version)
inline
const double& GSymMatrix::operator() (int row, int col) const
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GSymMatrix access", row, col, m_rows, m_cols);
  #endif
  int inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);
  return m_data[inx];
}

// Binary matrix addition
inline
GSymMatrix GSymMatrix::operator+ (const GSymMatrix& m) const
{
  GSymMatrix result = *this;
  result += m;
  return result;
}

// Binary matrix subtraction
inline
GSymMatrix GSymMatrix::operator- (const GSymMatrix& m) const
{
  GSymMatrix result = *this;
  result -= m;
  return result;
}

// Binary matrix multiplication
inline
GSymMatrix GSymMatrix::operator* (const GSymMatrix& m) const
{
  GSymMatrix result = *this;
  result *= m;
  return result;
}

// Unary scalar multiplication
inline
GSymMatrix& GSymMatrix::operator*= (const double& d)
{
  this->GMatrix::operator*=(d);
  return *this;
}

// Unary scalar division
inline
GSymMatrix& GSymMatrix::operator/= (const double& d)
{
  this->GMatrix::operator/=(d);
  return *this;
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline 
GSymMatrix operator* (const GSymMatrix& a, const double& b)
{
  GSymMatrix result = a;
  result *= b;
  return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GSymMatrix operator* (const double& a, const GSymMatrix& b)
{
  GSymMatrix result = b;
  result *= a;
  return result;
}

// Binary matrix division (matrix is left operand)
inline 
GSymMatrix operator/ (const GSymMatrix& a, const double& b)
{
  GSymMatrix result = a;
  result /= b;
  return result;
}

// Matrix transpose function
inline
GSymMatrix transpose(const GSymMatrix& m)
{
  return m;
}

// Cholesky decomposition
inline
GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress)
{
  GSymMatrix result = m;
  result.cholesky_decompose(compress);
  return result;
}

// Matrix inversion using Cholesky decomposition
inline
GSymMatrix cholesky_invert(const GSymMatrix& m, int compress)
{
  GSymMatrix result = m;
  result.cholesky_invert(compress);
  return result;
}

#endif /* GSYMMATRIX_HPP */
