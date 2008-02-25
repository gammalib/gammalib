/***************************************************************************
 *                   GMatrix.i  -  Matrix class SWIG file                  *
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
#include "GMatrix.hpp"
%}

%feature("notabstract") GMatrix;

/***************************************************************************
 *                         GMatrix class definition                        *
 * ----------------------------------------------------------------------- *
 * GMatrix implements a full matrix storage class. It derives from the     *
 * abstract base class GMatrixBase.                                        *
 ***************************************************************************/
class GMatrix : public GMatrixBase {
public:
  // Constructors and destructors
  GMatrix(int rows, int cols);
  GMatrix(const GMatrix& m);
  GMatrix(const GSymMatrix& m);
  GMatrix(const GSparseMatrix& m);
  virtual ~GMatrix();

  // Operators
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
};


/***************************************************************************
 *                      GMatrix class SWIG extension                       *
 ***************************************************************************/
%extend GMatrix {
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
  GMatrix __mul__(const double &a) {
    return (*self) * a;
  }
  GMatrix __div__(const double &a) {
    return (*self) / a;
  }
  GMatrix copy() {
    return (*self);
  }
};
