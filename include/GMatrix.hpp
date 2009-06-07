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
#include "GMatrixBase.hpp"


/***************************************************************************
 *                         GMatrix class definition                        *
 * ----------------------------------------------------------------------- *
 * GMatrix implements a full matrix storage class. It derives from the     *
 * abstract base class GMatrixBase.                                        *
 ***************************************************************************/
class GMatrix : public GMatrixBase {
  // Binary operator friends
  friend GMatrix operator* (const double& a,  const GMatrix& b);
  friend GMatrix operator* (const GMatrix& a, const double& b);
  friend GMatrix operator/ (const GMatrix& a, const double& b);

  // I/O friends
  friend std::ostream& operator<< (std::ostream& os, const GMatrix& m);

  // Friend functions
  friend GMatrix transpose(const GMatrix& m);
  friend GMatrix fabs(const GMatrix& m);

public:
  // Constructors and destructors
  GMatrix(int rows, int cols);
  GMatrix(const GMatrix& m);
  GMatrix(const GSymMatrix& m);
  GMatrix(const GSparseMatrix& m);
  virtual ~GMatrix();

  // Operators
  virtual GMatrix&      operator= (const GMatrix& m);
  virtual double&       operator() (int row, int col);
  virtual const double& operator() (int row, int col) const;
  virtual GMatrix       operator+ (const GMatrix& m) const;
  virtual GMatrix       operator- (const GMatrix& m) const;
  virtual GMatrix       operator* (const GMatrix& m) const;
  virtual GVector       operator* (const GVector& v) const;
  //_USE_BASE virtual int           operator== (const GMatrix& m) const;
  //_USE_BASE virtual int           operator!= (const GMatrix& m) const;
  virtual GMatrix       operator- () const;
  virtual GMatrix&      operator+= (const GMatrix& m);
  virtual GMatrix&      operator-= (const GMatrix& m);
  virtual GMatrix&      operator*= (const GMatrix& m);
  virtual GMatrix&      operator*= (const double& s);
  virtual GMatrix&      operator/= (const double& s);

  // Methods
  virtual void    add_col(const GVector& v, int col);
  //TBD virtual void    cholesky_decompose(int compress = 1);
  //TBD virtual GVector cholesky_solver(const GVector& v, int compress = 1);
  //TBD virtual void    cholesky_invert(int compress = 1);
  virtual void    clear();
  virtual GVector extract_row(int row) const;
  virtual GVector extract_col(int col) const;
  virtual GMatrix extract_lower_triangle() const;
  virtual GMatrix extract_upper_triangle() const;
  virtual double  fill() const;
  virtual void    insert_col(const GVector& v, int col);
  virtual double  min() const;
  virtual double  max() const;
  virtual double  sum() const;
  virtual void    transpose();
  
private:
  // Private methods
  void constructor(int rows, int cols);
  void init_members(void) { ; }
  void copy_members(const GMatrix& m) { ; }
  void free_members(void) { ; }
};


/***************************************************************************
 *                            Inline operators                             *
 ***************************************************************************/
// Matrix element access operator
inline
double& GMatrix::operator() (int row, int col)
{
  #if defined(G_RANGE_CHECK)
  if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
    throw GException::out_of_range("GMatrix::operator(int,int)", row, col, m_rows, m_cols);
  #endif
  return m_data[m_colstart[col]+row];
}

// Matrix element access operator (const version)
inline
const double& GMatrix::operator() (int row, int col) const
{
  #if defined(G_RANGE_CHECK)
  if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
    throw GException::out_of_range("GMatrix::operator(int,int) const", row, col, m_rows, m_cols);
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

// Matrix scaling
inline
GMatrix& GMatrix::operator*= (const double& s)
{
  multiplication(s);
  return *this;
}

// Matrix scalar division
inline
GMatrix& GMatrix::operator/= (const double& s)
{
  double inverse = 1.0/s;
  multiplication(inverse);
  return *this;
}

// Negation
inline
GMatrix GMatrix::operator- ( ) const
{
  GMatrix result = *this;
  result.negation();
  return result;
}


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline void   GMatrix::clear() { set_all_elements(0.0); }
inline double GMatrix::min() const { return get_min_element(); }
inline double GMatrix::max() const { return get_max_element(); }
inline double GMatrix::sum() const { return get_element_sum(); }


/***************************************************************************
 *                              Inline friends                             *
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
