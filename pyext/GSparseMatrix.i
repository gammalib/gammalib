/***************************************************************************
 *            GSparseMatrix.i  -  Sparse Matrix class SWIG file            *
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
#include "GSparseMatrix.hpp"
%}

%feature("notabstract") GSparseMatrix;

/***************************************************************************
 *                      GSparseMatrix class definition                     *
 * ----------------------------------------------------------------------- *
 * GMatrix implements a sparse matrix storage class. It derives from the   *
 * abstract base class GMatrixBase.                                        *
 ***************************************************************************/
class GSparseMatrix : public GMatrixBase {
public:
  // Constructors and destructors
  GSparseMatrix(int rows, int cols, int elements = 0);
  GSparseMatrix(const GMatrix& m);
  GSparseMatrix(const GSymMatrix& m);
  GSparseMatrix(const GSparseMatrix& m);
  virtual ~GSparseMatrix();

  // Operators
  //virtual GSparseMatrix& operator= (const GSparseMatrix& m);
  //virtual double&        operator() (int row, int col);
  //virtual const double&  operator() (int row, int col) const;
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
};


/***************************************************************************
 *                  GSparseMatrix class SWIG extension                     *
 ***************************************************************************/
%extend GSparseMatrix {
  char *__str__() {
    std::ostringstream buffer;
	buffer << *self;
	static std::string str = buffer.str();
	char* ptr = (char*)str.c_str();
	return ptr;
  }
  void __setslice__(int row, int col, const double val) throw (std::out_of_range) {
    if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
	  (*self)(row,col) = val;
    else
      throw std::out_of_range("row,column");
  }
  double __getslice__(int row, int col) throw (std::out_of_range) {
    if (row >=0 && row < (int)self->rows() && col >= 0 || col < (int)self->cols())
	  return (*self)(row,col);
    else
	  throw std::out_of_range("row,column");
  }
  GSparseMatrix __mul__(const double &a) {
    return (*self) * a;
  }
  GSparseMatrix __div__(const double &a) {
    return (*self) / a;
  }
  GSparseMatrix copy() {
    return (*self);
  }
};
