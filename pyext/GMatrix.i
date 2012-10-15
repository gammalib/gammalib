/***************************************************************************
 *                  GMatrix.i  -  General Matrix class                     *
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
    GMatrix(const int& rows, const int& cols);
    GMatrix(const GMatrix& matrix);
    GMatrix(const GSymMatrix& matrix);
    GMatrix(const GSparseMatrix& matrix);
    virtual ~GMatrix(void);

    // Implemented pure virtual base class operators
    GVector       operator*(const GVector& vector) const;

    // Other operators
    GMatrix       operator+(const GMatrix& matrix) const;
    GMatrix       operator-(const GMatrix& matrix) const;
    GMatrix       operator*(const GMatrix& matrix) const;
    GMatrix       operator-(void) const;
    GMatrix&      operator+=(const GMatrix& matrix);
    GMatrix&      operator-=(const GMatrix& matrix);
    GMatrix&      operator*=(const GMatrix& matrix);
    GMatrix&      operator*=(const double& scalar);
    GMatrix&      operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    void        clear(void);
    GMatrix*    clone(void) const;
    void        transpose(void);
    void        invert(void);
    void        add_col(const GVector& vector, const int& col);
    void        insert_col(const GVector& vector, const int& col);
    GVector     extract_row(const int& row) const;
    GVector     extract_col(const int& col) const;
    double      fill(void) const;
    double      min(void) const;
    double      max(void) const;
    double      sum(void) const;

    // Other methods
    GMatrix     extract_lower_triangle(void) const;
    GMatrix     extract_upper_triangle(void) const;
    void        eulerx(const double& angle);
    void        eulery(const double& angle);
    void        eulerz(const double& angle);
};


/***********************************************************************//**
 * @brief GMatrix class extension
 ***************************************************************************/
%extend GMatrix {
    char *__str__() {
        return tochar(self->print());
    }
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
