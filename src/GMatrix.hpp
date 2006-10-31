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
  // Binary operator friends
  friend GMatrix operator* (const double& a,  const GMatrix& b);
  friend GMatrix operator* (const GMatrix& a, const double& b);
  friend GMatrix operator/ (const GMatrix& a, const double& b);

  // I/O friends
  friend ostream& operator<< (ostream& os, const GMatrix& m);

  // Friend functions
  friend GMatrix transpose(const GMatrix& m);
  friend GMatrix fabs(const GMatrix& m);

public:
  // Constructors and destructors
  GMatrix(int rows, int cols);
  GMatrix(const GMatrix& m);
  virtual ~GMatrix();

  // Matrix element access operators
  virtual double& operator() (int row, int col);
  virtual const double& operator() (int row, int col) const;

  // Matrix assignment operators
  virtual GMatrix& operator= (const GMatrix& m);

  // Binary operators
  virtual GMatrix operator+ (const GMatrix& m) const;
  virtual GMatrix operator- (const GMatrix& m) const;
  virtual GMatrix operator* (const GMatrix& m) const;
  virtual GVector operator* (const GVector& v) const;
  virtual int     operator== (const GMatrix& m) const;
  virtual int     operator!= (const GMatrix& m) const;

  // Unary operators
  GMatrix operator- () const;
  virtual GMatrix& operator+= (const GMatrix& m);
  virtual GMatrix& operator-= (const GMatrix& m);
  virtual GMatrix& operator*= (const GMatrix& m);
  virtual GMatrix& operator*= (const double& d);
  virtual GMatrix& operator/= (const double& d);

  // Matrix functions
  virtual int     rows() const { return m_rows; }
  virtual int     cols() const { return m_cols; }
  virtual void    clear();
  virtual double  min() const;
  virtual double  max() const;
  virtual double  sum() const;
  virtual void    transpose();
  virtual GVector extract_row(int row) const;
  virtual GVector extract_col(int col) const;
  
  // Exception: Matrix indices out of range
  class out_of_range : public GException {
  public:
    out_of_range(string origin, int row, int col, int rows, int cols);
  };
  
  // Exception: Vector - Matrix mismatch
  class matrix_vector_mismatch : public GException {
  public:
    matrix_vector_mismatch(string origin, int num, int rows, int cols);
  };

  // Exception: Matrix dimensions mismatch
  class matrix_mismatch : public GException {
  public:
    matrix_mismatch(string origin, int rows1, int cols1, int cols1, int cols2);
  };

  // Exception: Matrix not rectangular
  class matrix_not_rect : public GException {
  public:
    matrix_not_rect(string origin, int rows1, int cols1, int cols1, int cols2);
  };

  // Exception: Matrix not positive definite
  class matrix_not_pos_definite : public GException {
  public:
    matrix_not_pos_definite(string origin, int row, double sum);
  };

  // Exception: Matrix not symmetric
  class matrix_not_symmetric : public GException {
  public:
    matrix_not_symmetric(string origin, int cols, int rows);
  };

  // Exception: Matrix not factorised
  class matrix_not_factorised : public GException {
  public:
    matrix_not_factorised(string origin, string type);
  };

  // Exception: All matrix elements are zero
  class matrix_zero : public GException {
  public:
    matrix_zero(string origin);
  };

protected:
  // Void constructor to be used by derived classes which have own constructor
  GMatrix();

  // Private member functions
  void select_non_zero(void);

  // Private data area
  int     m_rows;         // Number of rows
  int     m_cols;         // Number of columns
  int     m_elements;     // Number of elements stored in matrix
  int     m_alloc;        // Size of allocated matrix memory
  int     m_num_rowsel;   // Number of selected rows (for compressed decomposition)
  int     m_num_colsel;   // Number of selected columns (for compressed decomposition)
  int*    m_colstart;     // Column start indices (m_cols+1)
  int*    m_rowsel;       // Row selection (for compressed decomposition)
  int*    m_colsel;       // Column selection (for compressed decomposition)
  double* m_data;         // Matrix data
};


/***************************************************************************
 *                          Inline member functions                        *
 ***************************************************************************/
// Matrix element access operator
inline
double& GMatrix::operator() (int row, int col)
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GMatrix access", row, col, m_rows, m_cols);
  #endif
  return m_data[m_colstart[col]+row];
}

// Matrix element access operator (const version)
inline
const double& GMatrix::operator() (int row, int col) const
{
  #if defined(G_RANGE_CHECK)
  if (row >= m_rows || col >= m_cols)
    throw out_of_range("GMatrix access", row, col, m_rows, m_cols);
  #endif
  return m_data[m_colstart[col]+row];
}

// Binary matrix addition
inline
GMatrix GMatrix::operator+ (const GMatrix& m) const
{
  GMatrix result = *this;
  result += m;
  return result;
}

// Binary matrix subtraction
inline
GMatrix GMatrix::operator- (const GMatrix& m) const
{
  GMatrix result = *this;
  result -= m;
  return result;
}

// Binary matrix multiplication
inline
GMatrix GMatrix::operator* (const GMatrix& m) const
{
  GMatrix result = *this;
  result *= m;
  return result;
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline 
GMatrix operator* (const GMatrix& a, const double& b)
{
  GMatrix result = a;
  result *= b;
  return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GMatrix operator* (const double& a, const GMatrix& b)
{
  GMatrix result = b;
  result *= a;
  return result;
}

// Binary matrix division (matrix is left operand)
inline 
GMatrix operator/ (const GMatrix& a, const double& b)
{
  GMatrix result = a;
  result /= b;
  return result;
}

// Matrix transpose function
inline
GMatrix transpose(const GMatrix& m)
{
  GMatrix result = m;
  result.transpose();
  return result;
}

#endif /* GMATRIX_HPP */
