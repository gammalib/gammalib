/***************************************************************************
 *                   GMatrix.i - General Matrix class                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GMatrix.i
 * @brief Generic matrix class definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMatrix.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GMatrix
 *
 * @brief General matrix class definition
 *
 * GMatrix implements a general matrix storage class. It derives from the
 * abstract base class GMatrixBase.
 ***************************************************************************/
class GMatrix : public GMatrixBase {
public:
    // Constructors and destructors
    GMatrix(void);
    explicit GMatrix(const int& rows, const int& columns);
    GMatrix(const GMatrix& matrix);
    GMatrix(const GMatrixSymmetric& matrix);
    GMatrix(const GMatrixSparse& matrix);
    virtual ~GMatrix(void);

    // Implemented pure virtual base class operators
    virtual GVector       operator*(const GVector& vector) const;

    // Other operators
    virtual GMatrix       operator+(const GMatrix& matrix) const;
    virtual GMatrix       operator-(const GMatrix& matrix) const;
    virtual GMatrix       operator*(const GMatrix& matrix) const;
    virtual GMatrix       operator-(void) const;
    virtual GMatrix&      operator+=(const GMatrix& matrix);
    virtual GMatrix&      operator-=(const GMatrix& matrix);
    virtual GMatrix&      operator*=(const GMatrix& matrix);
    virtual GMatrix&      operator*=(const double& scalar);
    virtual GMatrix&      operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GMatrix*    clone(void) const;
    virtual GVector     row(const int& row) const;
    virtual void        row(const int& row, const GVector& vector);
    virtual GVector     column(const int& column) const;
    virtual void        column(const int& column, const GVector& vector);
    virtual void        add_to_row(const int& row, const GVector& vector);
    virtual void        add_to_column(const int& column, const GVector& vector);
    virtual void        transpose(void);
    virtual void        invert(void);
    virtual void        negate(void);
    virtual void        abs(void);
    virtual double      fill(void) const;
    virtual double      min(void) const;
    virtual double      max(void) const;
    virtual double      sum(void) const;

    // Other methods
    virtual GMatrix extract_lower_triangle(void) const;
    virtual GMatrix extract_upper_triangle(void) const;
    virtual void    eulerx(const double& angle);
    virtual void    eulery(const double& angle);
    virtual void    eulerz(const double& angle);
};


/***********************************************************************//**
 * @brief GMatrix class extension
 ***************************************************************************/
%extend GMatrix {
    double __getitem__(int GTuple[2]) {
        return (*self)(GTuple[0], GTuple[1]);
    }
    void __setitem__(int GTuple[2], double value) {
        (*self)(GTuple[0], GTuple[1]) = value;
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
