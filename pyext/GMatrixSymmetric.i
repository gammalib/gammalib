/***************************************************************************
 *              GMatrixSymmetric.i - Symmetric matrix class                *
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
 * @file GMatrixSymmetric.i
 * @brief Symmetric matrix class definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMatrixSymmetric.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GMatrixSymmetric
 *
 * @brief Symmetric matrix class definition
 *
 * GMatrixSymmetric implements a symmetric matrix storage class. It derives from
 * the abstract base class GMatrixBase.
 ***************************************************************************/
class GMatrixSymmetric : public GMatrixBase {
public:
    // Constructors and destructors
    GMatrixSymmetric(void);
    explicit GMatrixSymmetric(const int& rows, const int& columns);
    GMatrixSymmetric(const GMatrix& matrix);
    GMatrixSymmetric(const GMatrixSparse& matrix);
    GMatrixSymmetric(const GMatrixSymmetric& matrix);
    virtual ~GMatrixSymmetric(void);

    // Other operators
    virtual GMatrixSymmetric  operator+(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric  operator-(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric  operator-(void) const;
    virtual GMatrixSymmetric& operator+=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric& operator-=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric& operator*=(const double& scaler);
    virtual GMatrixSymmetric& operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void              clear(void);
    virtual GMatrixSymmetric* clone(void) const;
    virtual double&           at(const int& row, const int& column);
    virtual GVector           row(const int& row) const;
    virtual void              row(const int& row, const GVector& vector);
    virtual GVector           column(const int& column) const;
    virtual void              column(const int& column, const GVector& vector);
    virtual void              add_to_row(const int& row, const GVector& vector);
    virtual void              add_to_column(const int& column, const GVector& vector);
    virtual void              transpose(void);
    virtual void              invert(void);
    virtual void              negate(void);
    virtual void              abs(void);
    virtual double            fill(void) const;
    virtual double            min(void) const;
    virtual double            max(void) const;
    virtual double            sum(void) const;

    // Other methods
    virtual GMatrix extract_lower_triangle(void) const;
    virtual GMatrix extract_upper_triangle(void) const;
    virtual void    cholesky_decompose(bool compress = true);
    virtual GVector cholesky_solver(const GVector& vector, bool compress = true);
    virtual void    cholesky_invert(bool compress = true);
};


/***********************************************************************//**
 * @brief GMatrixSymmetric class extension
 ***************************************************************************/
%extend GMatrixSymmetric {
    double __getitem__(int GTuple[2]) {
        return (*self)(GTuple[0], GTuple[1]);
    }
    void __setitem__(int GTuple[2], double value) {
        (*self)(GTuple[0], GTuple[1]) = value;
    }
    GVector __mul__(const GVector& vector) {
        return ((*self) * vector);
    }
    GMatrixSymmetric __mul__(const GMatrixSymmetric& matrix) {
        return ((*self) * matrix);
    }
    GMatrixSymmetric __mul__(const double &scalar) {
        return ((*self) * scalar);
    }
    GMatrixSymmetric __div__(const double &scalar) {
        return ((*self) / scalar);
    }
    GMatrixSymmetric copy() {
        return (*self);
    }
    GMatrixSymmetric set(const double& value) {
        (*self) = value;
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GMatrixSymmetric friends
 ***************************************************************************/
GMatrixSymmetric transpose(const GMatrixSymmetric& matrix);
GMatrixSymmetric invert(const GMatrixSymmetric& matrix);
GMatrixSymmetric abs(const GMatrixSymmetric& matrix);
