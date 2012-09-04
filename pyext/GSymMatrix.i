/***************************************************************************
 *                GSymMatrix.i  -  Symmetric Matrix class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSymMatrix.hpp"
#include "GTools.hpp"
%}

%feature("notabstract") GSymMatrix;

/***********************************************************************//**
 * @class GSymMatrix
 *
 * @brief Symmetric matrix class definition
 *
 * GSymMatrix implements a symmetric matrix storage class. It derives from
 * the abstract base class GMatrixBase.
 ***************************************************************************/
class GSymMatrix : public GMatrixBase {
public:
    // Constructors and destructors
    GSymMatrix(void);
    GSymMatrix(int rows, int cols);
    GSymMatrix(const GMatrix& matrix);
    GSymMatrix(const GSymMatrix& matrix);
    GSymMatrix(const GSparseMatrix& matrix);
    virtual ~GSymMatrix(void);

    // Operators
    GSymMatrix    operator+ (const GSymMatrix& m) const;
    GSymMatrix    operator- (const GSymMatrix& m) const;
    GSymMatrix    operator* (const GSymMatrix& m) const;
    GVector       operator* (const GVector& v) const;
    GSymMatrix    operator- () const;
    GSymMatrix&   operator+= (const GSymMatrix& m);
    GSymMatrix&   operator-= (const GSymMatrix& m);
    GSymMatrix&   operator*= (const GSymMatrix& m);
    GSymMatrix&   operator*= (const double& d);
    GSymMatrix&   operator/= (const double& d);

    // Methods
    virtual void    add_col(const GVector& v, int col);
    virtual void    cholesky_decompose(int compress = 1);
    virtual GVector cholesky_solver(const GVector& v, int compress = 1);
    virtual void    cholesky_invert(int compress = 1);
    virtual void    clear(void);
    virtual GVector extract_row(int row) const;
    virtual GVector extract_col(int col) const;
    virtual GMatrix extract_lower_triangle() const;
    virtual GMatrix extract_upper_triangle() const;
    virtual void    insert_col(const GVector& v, int col);
    virtual double  fill(void) const;
    virtual double  min(void) const;
    virtual double  max(void) const;
    virtual double  sum(void) const;
    virtual void    transpose(void);
};


/***********************************************************************//**
 * @brief GSymMatrix class extension
 *
 * @todo Use typemap to allow for a [row,col] matrix access
 ***************************************************************************/
%extend GSymMatrix {
    char *__str__() {
        return tochar(self->print());
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
