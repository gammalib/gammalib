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
/**
 * @file GSymMatrix.i
 * @brief SYmmetric matrix class definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSymMatrix.hpp"
#include "GTools.hpp"
%}


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
    GSymMatrix(const int& rows, const int& cols);
    GSymMatrix(const GMatrix& matrix);
    GSymMatrix(const GSymMatrix& matrix);
    GSymMatrix(const GSparseMatrix& matrix);
    virtual ~GSymMatrix(void);

    // Implemented pure virtual base class operators
    virtual GVector       operator*(const GVector& v) const;

    // Other operators
    virtual GSymMatrix    operator+(const GSymMatrix& matrix) const;
    virtual GSymMatrix    operator-(const GSymMatrix& matrix) const;
    virtual GMatrix       operator*(const GSymMatrix& matrix) const;
    virtual GSymMatrix    operator-(void) const;
    virtual GSymMatrix&   operator+=(const GSymMatrix& matrix);
    virtual GSymMatrix&   operator-=(const GSymMatrix& matrix);
    virtual GSymMatrix&   operator*=(const double& scaler);
    virtual GSymMatrix&   operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual void        transpose(void);
    virtual void        invert(void);
    virtual void        add_col(const GVector& vector, const int& col);
    virtual void        insert_col(const GVector& vector, const int& col);
    virtual GVector     extract_row(const int& row) const;
    virtual GVector     extract_col(const int& col) const;
    virtual double      fill(void) const;
    virtual double      min(void) const;
    virtual double      max(void) const;
    virtual double      sum(void) const;

    // Other methods
    virtual GMatrix     extract_lower_triangle(void) const;
    virtual GMatrix     extract_upper_triangle(void) const;
    virtual void        cholesky_decompose(bool compress = true);
    virtual GVector     cholesky_solver(const GVector& vector, bool compress = true);
    virtual void        cholesky_invert(bool compress = true);
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
