/***************************************************************************
 *                         GMatrix.hpp  -  matrix class                    *
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

#ifndef GMATRIX_HPP
#define GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***************************************************************************
 *                         GMatrix class definition                        *
 ***************************************************************************/
class GMatrix {
  // Operator friends
  friend GMatrix operator+ (const double& a, const GMatrix& b);
  friend GMatrix operator- (const double& a, const GMatrix& b);
  friend GMatrix operator* (const double& a, const GMatrix& b);

  // I/O friends
  friend ostream& operator<< (ostream& os, const GMatrix& m); 

  // Friend functions
  friend GMatrix transpose(const GMatrix& m);
  friend GMatrix fabs(const GMatrix& m);
public:
  // Constructors and destructors
  GMatrix(unsigned rows, unsigned cols);
  GMatrix(const GMatrix& m);
  virtual ~GMatrix();

  // Matrix element access operators
  virtual double& operator() (unsigned row, unsigned col);
  virtual const double& operator() (unsigned row, unsigned col) const;

  // Matrix assignment operators
  GMatrix& operator= (const GMatrix& m);
  GMatrix& operator= (const double& d);

  // Generic matrix operators (structure independent)
  virtual GMatrix& operator+= (const GMatrix& m);
  virtual GMatrix& operator-= (const GMatrix& m);
  virtual GMatrix& operator+= (const double& d);
  virtual GMatrix& operator-= (const double& d);
  virtual GMatrix& operator*= (const double& d);
  virtual GMatrix& operator/= (const double& d);
  virtual int      operator== (const GMatrix& m) const;
  virtual int      operator!= (const GMatrix& m) const;
  
  // Specific matrix operators (structure dependent)
  GMatrix operator+ (const GMatrix& m) const;
  GMatrix operator- (const GMatrix& m) const;
  GMatrix operator* (const GMatrix& m) const;
  virtual GVector operator* (const GVector& v) const;
  GMatrix operator+ (const double& d) const;
  GMatrix operator- (const double& d) const;
  GMatrix operator* (const double& d) const;
  GMatrix operator/ (const double& d) const;
  GMatrix operator- () const;

  // Generic matrix functions (structure independent)
  virtual double   min() const;      // Minimum physical element
  virtual double   max() const;      // Maximum physical element
  virtual double   sum() const;      // Sum of matrix elements
  virtual unsigned rows() const;     // Return number of rows
  virtual unsigned cols() const;     // Return number of columns

  // Specific matrix functions (structure dependent)
  virtual void     transpose();      // Inplace transpose of matrix
  
  // Exception: Matrix indices out of range
  class out_of_range : public GException {
  public:
    out_of_range(string origin, unsigned row, unsigned col, 
	             unsigned rows, unsigned cols);
  };
  
  // Exception: Vector - Matrix mismatch
  class vec_mat_mismatch : public GException {
  public:
    vec_mat_mismatch(string origin, unsigned num, unsigned rows, unsigned cols);
  };

  // Exception: Matrix dimensions mismatch (addition, etc.)
  class dim_add_mismatch : public GException {
  public:
    dim_add_mismatch(string origin, unsigned rows1, unsigned rows2,
									unsigned cols1, unsigned cols2);
  };

  // Exception: Matrix dimensions mismatch (multiplication)
  class dim_mult_mismatch : public GException {
  public:
    dim_mult_mismatch(string origin, unsigned cols, unsigned rows);
  };

protected:
  // Void constructor to be used by derived classes which have own constructor
  GMatrix() { }

  // Private data area
  unsigned  m_rows;     // Number of rows
  unsigned  m_cols;     // Number of columns
  unsigned  m_elements; // Number of (physical) elements in matrix
  unsigned* m_colstart; // Column start indices (m_cols+1)
  double*   m_data;     // Matrix data
};


/***************************************************************************
 *                  Generic inline members (to be inherited)               *
 ***************************************************************************/
// Matrix unary addition operator: GMatrix += GMatrix
inline
GMatrix& GMatrix::operator+= (const GMatrix& m)
{
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw dim_add_mismatch("GMatrix += operator", m_rows, m.m_rows, m_cols, m.m_cols);
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] += m.m_data[i];
  return *this;
}

// Matrix unary subtraction operator: GMatrix -= GMatrix
inline
GMatrix& GMatrix::operator-= (const GMatrix& m)
{
  if (m_rows != m.m_rows || m_cols != m.m_cols)
    throw dim_add_mismatch("GMatrix -= operator", m_rows, m.m_rows, m_cols, m.m_cols);
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] -= m.m_data[i];
  return *this;
}

// Scalar unary addition operator: GMatrix += double
inline
GMatrix& GMatrix::operator+= (const double& d)
{
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] += d;
  return *this;
}

// Scalar unary subtraction operator: GMatrix -= double
inline
GMatrix& GMatrix::operator-= (const double& d)
{
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] -= d;
  return *this;
}

// Scalar unary multiplication operator: GMatrix *= double
inline
GMatrix& GMatrix::operator*= (const double& d)
{
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] *= d;
  return *this;
}

// Scalar unary division operator: GMatrix /= double
inline
GMatrix& GMatrix::operator/= (const double& d)
{
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] /= d;
  return *this;
}

// Equality operator: GMatrix == GMatrix
inline
int GMatrix::operator== (const GMatrix &m) const
{
  int result = 1;
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (unsigned i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 0;
		break;
	  }
	}
  }
  else
    result = 0;
  return result;
}

// Non equality operator: GMatrix != GMatrix
inline
int GMatrix::operator!= (const GMatrix &m) const
{
  int result = 0;
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (unsigned i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 1;
		break;
	  }
	}
  }
  else
    result = 1;
  return result;
}

// Matrix minimum: GMatrix.min()
inline
double GMatrix::min() const
{
  double result = m_data[0];
  for (unsigned i = 1; i < m_elements; ++i) {
    if (m_data[i] < result)
	  result = m_data[i];
  }
  return result;
}

// Matrix maximum: GMatrix.max()
inline
double GMatrix::max() const
{
  double result = m_data[0];
  for (unsigned i = 1; i < m_elements; ++i) {
    if (m_data[i] > result)
	  result = m_data[i];
  }
  return result;
}

// Matrix sum: GMatrix.sum()
inline
double GMatrix::sum() const
{
  double result = 0.0;
  for (unsigned row = 0; row < m_rows; ++row) {
    for (unsigned col = 0; col < m_cols; ++col)
      result += (*this)(row,col);
  }
  return result;
}

// Return number of rows: GMatrix.rows()
inline
unsigned GMatrix::rows() const
{
  return m_rows;
}

// Return number of columns: GMatrix.cols()
inline
unsigned GMatrix::cols() const
{
  return m_cols;
}


/***************************************************************************
 * Specific inline members (may need specific function for derived class)  *
 ***************************************************************************/
// Scalar assignment operator: GMatrix = double
inline
GMatrix& GMatrix::operator= (const double& d)
{
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] = d;
  return *this;
}

// Matrix element access operator: GMatrix(row,col)
inline
double& GMatrix::operator() (unsigned row, unsigned col)
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GMatrix access", row, col, m_rows, m_cols);
  #endif
  return m_data[m_colstart[col]+row];
}

// Matrix element access operator (const version): GMatrix(row,col)
inline
const double& GMatrix::operator() (unsigned row, unsigned col) const
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GMatrix access", row, col, m_rows, m_cols);
  #endif
  return m_data[m_colstart[col]+row];
}

// Matrix addition: GMatrix + GMatrix
inline
GMatrix GMatrix::operator+ (const GMatrix& m) const
{
  GMatrix result = *this;
  result += m;
  return result;
}

// Matrix subtraction: GMatrix - GMatrix
inline
GMatrix GMatrix::operator- (const GMatrix& m) const
{
  GMatrix result = *this;
  result -= m;
  return result;
}

// Matrix scaling: GMatrix + double
inline 
GMatrix GMatrix::operator+ (const double& d) const
{
  GMatrix result = (*this);
  result += d;
  return result;
}

// Matrix scaling: GMatrix - double
inline 
GMatrix GMatrix::operator- (const double& d) const
{
  GMatrix result = (*this);
  result -= d;
  return result;
}

// Matrix scaling: GMatrix * double
inline 
GMatrix GMatrix::operator* (const double& d) const
{
  GMatrix result = (*this);
  result *= d;
  return result;
}

// Matrix inverse scaling: GMatrix * double
inline 
GMatrix GMatrix::operator/ (const double& d) const
{
  GMatrix result = (*this);
  result /= d;
  return result;
}

// Unary minus operator: -GMatrix
inline
GMatrix GMatrix::operator- ( ) const
{
  GMatrix result = *this;
  for (unsigned i = 0; i < m_elements; ++i)
    result.m_data[i] = -result.m_data[i];
  return result;
}

// Matrix transpose: GMatrix.transpose()
inline
void GMatrix::transpose()
{
  for (unsigned row = 0; row < m_rows; ++row) {
	double* ptr0 = m_data + row;
	double* ptr2 = m_data + m_colstart[row]+row+1;
    for (unsigned col = row+1; col < m_cols; ++col) {
	  double* ptr1 = ptr0 + m_colstart[col];
	  double  swap = *ptr1;
	  *ptr1        = *ptr2;
	  *ptr2++      = swap;
    }
  }
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Matrix addition: double + GMatrix
inline
GMatrix operator+ (const double& a, const GMatrix& b)
{
  GMatrix result = b;
  result += a;
  return result;
}

// Matrix subtraction: double - GMatrix
inline
GMatrix operator- (const double& a, const GMatrix& b)
{
  GMatrix result = b;
  result -= a;
  return result;
}

// Matrix scaling: double * GMatrix
inline
GMatrix operator* (const double& a, const GMatrix& b)
{
  GMatrix result = b;
  result *= a;
  return result;
}

// Matrix transpose: transpose(GMatrix)
inline
GMatrix transpose(const GMatrix& m)
{
  GMatrix result = m;
  result.transpose();
  return result;
}

// Matrix absolute: fabs(GMatrix)
inline
GMatrix fabs(const GMatrix& m)
{
  GMatrix result = m;
  for (unsigned i = 0; i < result.m_elements; ++i)
    result.m_data[i] = fabs(result.m_data[i]);
  return result;
}

#endif /* GMATRIX_HPP */
