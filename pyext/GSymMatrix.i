/***************************************************************************
 *            GSymMatrix.i  -  Symmetric Matrix class SWIG file            *
 * ----------------------------------------------------------------------- *
 *  copyright - (C) 2008-2011 by Jurgen Knodlseder                         *
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
#include "GTools.hpp"
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
    GSymMatrix(const GSymMatrix& m);
    virtual ~GSymMatrix();

    // Operators
    GSymMatrix    operator+ (const GSymMatrix& m) const;
    GSymMatrix    operator- (const GSymMatrix& m) const;
    GSymMatrix    operator* (const GSymMatrix& m) const;
    GVector       operator* (const GVector& v) const;
    //_USE_BASE int           operator== (const GSymMatrix& m) const;
    //_USE_BASE int           operator!= (const GSymMatrix& m) const;
    GSymMatrix    operator- () const;
    GSymMatrix&   operator+= (const GSymMatrix& m);
    GSymMatrix&   operator-= (const GSymMatrix& m);
    GSymMatrix&   operator*= (const GSymMatrix& m);
    GSymMatrix&   operator*= (const double& d);
    GSymMatrix&   operator/= (const double& d);

    // Methods
    void    add_col(const GVector& v, int col);
    void    cholesky_decompose(int compress = 1);
    GVector cholesky_solver(const GVector& v, int compress = 1);
    void    cholesky_invert(int compress = 1);
    void    clear();
    GVector extract_row(int row) const;
    GVector extract_col(int col) const;
    GMatrix extract_lower_triangle() const;
    GMatrix extract_upper_triangle() const;
    void    insert_col(const GVector& v, int col);
    double  fill() const;
    double  min() const;
    double  max() const;
    double  sum() const;
    void    transpose() { return; }
};


/***************************************************************************
 *                    GSymMatrix class SWIG extension                      *
 ***************************************************************************/
%extend GSymMatrix {
    char *__str__() {
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        return tochar(str);
    }
    double __getslice__(int row, int col) {
        if (row >=0 && row < (int)self->rows() && col >= 0 && col < (int)self->cols())
            return (*self)(row,col);
        else
            throw GException::out_of_range("__getitem__(int,int)", row, col,
                                           (int)self->rows(), (int)self->cols());
    }
    void __setslice__(int row, int col, const double val) {
        if (row >=0 && row < (int)self->rows() && col >= 0 && col < (int)self->cols())
            (*self)(row,col) = val;
        else
            throw GException::out_of_range("__setitem__(int,int)", row, col,
                                           (int)self->rows(), (int)self->cols());
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
