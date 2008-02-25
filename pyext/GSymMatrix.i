/***************************************************************************
 *            GSymMatrix.i  -  Symmetric Matrix class SWIG file            *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSymMatrix.hpp"
%}

%feature("notabstract") GSymMatrix;

/***************************************************************************
 *                       GSymMatrix class definition                       *
 * ----------------------------------------------------------------------- *
 * GSymMatrix implements a symmetric matrix storage class. It derives from *
 * the abstract base class GMatrixBase.                                    *
 ***************************************************************************/
class GSymMatrix : public GMatrixBase {
public:
  // Constructors and destructors (not inherited)
  GSymMatrix(int rows, int cols);
  GSymMatrix(const GMatrix& m);
  GSymMatrix(const GSymMatrix& m);
  GSymMatrix(const GSparseMatrix& m);
  virtual ~GSymMatrix();

  // Operators
  virtual GSymMatrix    operator+ (const GSymMatrix& m) const;
  virtual GSymMatrix    operator- (const GSymMatrix& m) const;
  virtual GSymMatrix    operator* (const GSymMatrix& m) const;
  virtual GVector       operator* (const GVector& v) const;
  //_USE_BASE virtual int           operator== (const GSymMatrix& m) const;
  //_USE_BASE virtual int           operator!= (const GSymMatrix& m) const;
  virtual GSymMatrix    operator- () const;
  virtual GSymMatrix&   operator+= (const GSymMatrix& m);
  virtual GSymMatrix&   operator-= (const GSymMatrix& m);
  virtual GSymMatrix&   operator*= (const GSymMatrix& m);
  virtual GSymMatrix&   operator*= (const double& d);
  virtual GSymMatrix&   operator/= (const double& d);

  // Methods
  virtual void    add_col(const GVector& v, int col);
  virtual void    cholesky_decompose(int compress = 1);
  virtual GVector cholesky_solver(const GVector& v, int compress = 1);
  virtual void    cholesky_invert(int compress = 1);
  virtual void    clear();
  virtual GVector extract_row(int row) const;
  virtual GVector extract_col(int col) const;
  virtual GMatrix extract_lower_triangle() const;
  virtual GMatrix extract_upper_triangle() const;
  virtual void    insert_col(const GVector& v, int col);
  virtual double  fill() const;
  virtual double  min() const;
  virtual double  max() const;
  virtual double  sum() const;
  virtual void    transpose() { return; }
};


/***************************************************************************
 *                    GSymMatrix class SWIG extension                      *
 ***************************************************************************/
%extend GSymMatrix {
  char *__str__() {
    static char str_buffer[100001];
    std::ostringstream buffer;
    buffer << *self;
	std::string str = buffer.str();
    strncpy(str_buffer, (char*)str.c_str(), 100001);
	str_buffer[100000] = '\0';
	return str_buffer;
  }
  double __getslice__(int row, int col) {
    if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
	  return (*self)(row,col);
    else
	  throw GException::out_of_range("__getitem__(int,int)", row, col, (int)self->rows(), (int)self->cols());
  }
  void __setslice__(int row, int col, const double val) {
    if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
	  (*self)(row,col) = val;
    else
	  throw GException::out_of_range("__setitem__(int,int)", row, col, (int)self->rows(), (int)self->cols());
  }
  GSymMatrix __mul__(const double &a) {
    return (*self) * a;
  }
  GSymMatrix __div__(const double &a) {
    return (*self) / a;
  }
  GSymMatrix copy() {
    return (*self);
  }
};
