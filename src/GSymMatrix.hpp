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
  // Operator friends
  friend GSymMatrix operator+ (const double& a, const GSymMatrix& b);
  friend GSymMatrix operator- (const double& a, const GSymMatrix& b);
  friend GSymMatrix operator* (const double& a, const GSymMatrix& b);

  // I/O friends
  // ...
  
  // Friend functions
  friend GSymMatrix transpose(const GSymMatrix& m);
  friend GSymMatrix fabs(const GSymMatrix& m);
  friend GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress = 1);
  friend GSymMatrix cholesky_invert(const GSymMatrix& m, int compress = 1);
  friend GMatrix    full(const GSymMatrix& m);
  friend GMatrix    lower_triangle(const GSymMatrix& m);
  friend GMatrix    upper_triangle(const GSymMatrix& m);
public:
  // Constructors and destructors (not inherited)
  GSymMatrix(unsigned rows, unsigned cols);
  GSymMatrix(const GSymMatrix& m);
  ~GSymMatrix();

  // Matrix element access operators
  double& operator() (unsigned row, unsigned col);
  const double& operator() (unsigned row, unsigned col) const;

  // Matrix assignment operators (not inherited)
  GSymMatrix& operator= (const GSymMatrix& m);
  GSymMatrix& operator= (const double& d);

  // Specific matrix operators (structure dependent)
  GSymMatrix operator+ (const GSymMatrix& m) const;
  GSymMatrix operator- (const GSymMatrix& m) const;
  GSymMatrix operator* (const GSymMatrix& m) const;
  GVector    operator* (const GVector& v) const;
  GSymMatrix operator+ (const double& d) const;
  GSymMatrix operator- (const double& d) const;
  GSymMatrix operator* (const double& d) const;
  GSymMatrix operator/ (const double& d) const;
  GSymMatrix operator- () const;
  
  // Specific matrix functions
  void    transpose();                                         // Inplace transpose of matrix
  void    cholesky_decompose(int compress = 1);                // Inplace Cholesky decomposition
  GVector cholesky_solver(const GVector& v, int compress = 1); // Cholesky solver
  void    cholesky_invert(int compress = 1);                   // Inplace Cholesky matrix inverter
  
  // Exception: Matrix not symmetric
  class not_sym : public GException {
  public:
    not_sym(string origin, unsigned cols, unsigned rows);
  };

  // Exception: Matrix not positive definite
  class not_pos_definite : public GException {
  public:
    not_pos_definite(string origin, unsigned row, double sum);
  };

  // Exception: All matrix elements are zero
  class all_zero : public GException {
  public:
    all_zero(string origin);
  };
private:
  void      set_inx(void);  // Set indices of non-zero rows/columns
  unsigned  m_num_inx;      // Number of non-zero rows/columns
  unsigned* m_inx;          // Indices of non-zero rows/columns
};


/***************************************************************************
 *              Inline members that override base class members            *
 ***************************************************************************/
// Scalar assignment operator: GSymMatrix = double
inline
GSymMatrix& GSymMatrix::operator= (const double& d) 
{ 
  this->GMatrix::operator=(d);
  return *this;
}

// Matrix element access operator: GSymMatrix(row,col)
inline
double& GSymMatrix::operator() (unsigned row, unsigned col)
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GSymMatrix access", row, col, m_rows, m_cols);
  #endif
  unsigned inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);
  return m_data[inx];
}

// Matrix element access operator (const version): GSymMatrix(row,col)
inline
const double& GSymMatrix::operator() (unsigned row, unsigned col) const
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GSymMatrix access", row, col, m_rows, m_cols);
  #endif
  unsigned inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);
  return m_data[inx];
}

// Matrix addition: GSymMatrix + GSymMatrix
inline
GSymMatrix GSymMatrix::operator+ (const GSymMatrix& m) const
{
  GSymMatrix result = *this;
  result += m;
  return result;
}

// Matrix subtraction: GSymMatrix - GSymMatrix
inline
GSymMatrix GSymMatrix::operator- (const GSymMatrix& m) const
{
  GSymMatrix result = *this;
  result -= m;
  return result;
}

// Matrix/Vector multiplication: GSymMatrix * GVector
inline
GVector GSymMatrix::operator* (const GVector& v) const
{
  return (this->GMatrix::operator*(v));
}

// Matrix scaling: GSymMatrix + double
inline 
GSymMatrix GSymMatrix::operator+ (const double& d) const
{
  GSymMatrix result = (*this);
  result += d;
  return result;
}

// Matrix scaling: GSymMatrix - double
inline 
GSymMatrix GSymMatrix::operator- (const double& d) const
{
  GSymMatrix result = (*this);
  result -= d;
  return result;
}

// Matrix scaling: GSymMatrix * double
inline 
GSymMatrix GSymMatrix::operator* (const double& d) const
{
  GSymMatrix result = (*this);
  result *= d;
  return result;
}

// Matrix inverse scaling: GSymMatrix * double
inline 
GSymMatrix GSymMatrix::operator/ (const double& d) const
{
  GSymMatrix result = (*this);
  result /= d;
  return result;
}

// Unary minus operator
inline
GSymMatrix GSymMatrix::operator- ( ) const
{
  GSymMatrix result = *this;
  for (unsigned i = 0; i < m_elements; ++i)
    result.m_data[i] = -result.m_data[i];
  return result;
}

// Matrix transpose
inline
void GSymMatrix::transpose()
{
  return;
}



/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Matrix addition: double + GSymMatrix
inline
GSymMatrix operator+ (const double& a, const GSymMatrix& b)
{
  GSymMatrix result = b;
  result += a;
  return result;
}

// Matrix subtracxtion: double * GSymMatrix
inline
GSymMatrix operator- (const double& a, const GSymMatrix& b)
{
  GSymMatrix result = b;
  result -= a;
  return result;
}

// Matrix scaling: double * GSymMatrix
inline
GSymMatrix operator* (const double& a, const GSymMatrix& b)
{
  GSymMatrix result = b;
  result *= a;
  return result;
}

// Matrix transpose: transpose(GSymMatrix)
inline
GSymMatrix transpose(const GSymMatrix& m)
{
  return m;
}

// Matrix absolute: fabs(GSymMatrix)
inline
GSymMatrix fabs(const GSymMatrix& m)
{
  GSymMatrix result = m;
  for (unsigned i = 0; i < result.m_elements; ++i)
    result.m_data[i] = fabs(result.m_data[i]);
  return result;
}

// Cholesky decomposition: cholesky_decompose(GSymMatrix)
inline
GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress)
{
  GSymMatrix result = m;
  result.cholesky_decompose(compress);
  return result;
}

// Inversion using Cholesky decomposition: cholesky_invert(GSymMatrix)
inline
GSymMatrix cholesky_invert(const GSymMatrix& m, int compress)
{
  GSymMatrix result = m;
  result.cholesky_invert(compress);
  return result;
}

// Convert symmetric matrix into full matrix: full(GSymMatrix)
inline
GMatrix full(const GSymMatrix& m)
{
  GMatrix result(m.m_rows,m.m_cols);
  for (unsigned row = 0; row < m.m_rows; ++row) {
    for (unsigned col = 0; col < m.m_cols; ++col)
	  result(row,col) = m(row,col);
  }
  return result;
}

// Convert symmetric matrix into lower triangle full matrix: lower_triangle(GSymMatrix)
inline
GMatrix lower_triangle(const GSymMatrix& m)
{
  GMatrix result(m.m_rows,m.m_cols);
  for (unsigned row = 0; row < m.m_rows; ++row) {
    for (unsigned col = 0; col <= row; ++col)
	  result(row,col) = m.m_data[m.m_colstart[col]+(row-col)];
  }
  return result;
}

// Convert symmetric matrix into upper triangle full matrix: lower_triangle(GSymMatrix)
inline
GMatrix upper_triangle(const GSymMatrix& m)
{
  GMatrix result(m.m_rows,m.m_cols);
  for (unsigned row = 0; row < m.m_rows; ++row) {
    for (unsigned col = row; col < m.m_cols; ++col)
	  result(row,col) = m.m_data[m.m_colstart[row]+(col-row)];
  }
  return result;
}

#endif /* GSYMMATRIX_HPP */
