/***************************************************************************
 *                   GMatrix.i  -  Matrix class SWIG file                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
#include "GTools.hpp"
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
    virtual ~GMatrix();

    // Operators
    GMatrix  operator+ (const GMatrix& m) const;
    GMatrix  operator- (const GMatrix& m) const;
    GMatrix  operator* (const GMatrix& m) const;
    GVector  operator* (const GVector& v) const;
    GMatrix  operator- () const;
    GMatrix& operator+= (const GMatrix& m);
    GMatrix& operator-= (const GMatrix& m);
    GMatrix& operator*= (const GMatrix& m);
    GMatrix& operator*= (const double& s);
    GMatrix& operator/= (const double& s);

    // Methods
    void    clear(void);
    void    transpose(void);
    void    invert(void);
    void    add_col(const GVector& v, int col);
    void    insert_col(const GVector& v, int col);
    GVector extract_row(int row) const;
    GVector extract_col(int col) const;
    GMatrix extract_lower_triangle(void) const;
    GMatrix extract_upper_triangle(void) const;
    double  fill(void) const;
    double  min(void) const;
    double  max(void) const;
    double  sum(void) const;
    void    eulerx(const double& angle);
    void    eulery(const double& angle);
    void    eulerz(const double& angle);
};


/***************************************************************************
 *                      GMatrix class SWIG extension                       *
 ***************************************************************************/
%extend GMatrix {
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
