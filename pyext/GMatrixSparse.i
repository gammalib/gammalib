/***************************************************************************
 *                   GMatrixSparse.i - Sparse matrix class                 *
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
 * @file GMatrixSparse.i
 * @brief Sparse matrix class definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMatrixSparse.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GMatrixSparse
 *
 * @brief Sparse matrix class definition
 *
 * GMatrixSparse implements a sparse matrix storage class. It derives from
 * the abstract base class GMatrixBase.
 ***************************************************************************/
class GMatrixSparse : public GMatrixBase {
public:
    // Constructors and destructors
    GMatrixSparse(void);
    explicit GMatrixSparse(const int& rows, const int& columns,
                           const int& elements = 0);
    GMatrixSparse(const GMatrix& matrix);
    GMatrixSparse(const GMatrixSparse& matrix);
    GMatrixSparse(const GMatrixSymmetric& matrix);
    virtual ~GMatrixSparse(void);

    // Implemented pure virtual base class operators
    virtual GVector        operator*(const GVector& vector) const;

    // Overloaded virtual base class operators
    virtual bool           operator==(const GMatrixSparse& matrix) const;
    virtual bool           operator!=(const GMatrixSparse& matrix) const;

    // Other operators
    virtual GMatrixSparse  operator+(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator-(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator*(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator-(void) const;
    virtual GMatrixSparse& operator+=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator-=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator*=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator*=(const double& scalar);
    virtual GMatrixSparse& operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GMatrixSparse* clone(void) const;
    virtual double&        at(const int& row, const int& column);
    virtual GVector        row(const int& row) const;
    virtual void           row(const int& row, const GVector& vector);
    virtual GVector        column(const int& column) const;
    virtual void           column(const int& column, const GVector& vector);
    virtual void           add_to_row(const int& row, const GVector& vector);
    virtual void           add_to_column(const int& column, const GVector& vector);
    virtual void           transpose(void);
    virtual void           invert(void);
    virtual void           negate(void);
    virtual void           abs(void);
    virtual double         fill(void) const;
    virtual double         min(void) const;
    virtual double         max(void) const;
    virtual double         sum(void) const;

    // Other methods
    void    column(const int& column, const double* values,
                   const int* rows, int number);
    void    add_to_column(const int& column, const double* values,
                          const int* rows, int number);
    void    cholesky_decompose(bool compress = true);
    GVector cholesky_solver(const GVector& vector, bool compress = true);
    void    cholesky_invert(bool compress = true);
    void    set_mem_block(const int& block);
    void    stack_init(const int& size = 0, const int& entries = 0);
    int     stack_push_column(const GVector& vector, const int& col);
    int     stack_push_column(const double* values, const int* rows,
                              const int& number, const int& col);
    void    stack_flush(void);
    void    stack_destroy(void);
};


/***********************************************************************//**
 * @brief GMatrixSparse class extension
 ***************************************************************************/
%extend GMatrixSparse {
    double __getitem__(int GTuple[2]) {
        return (*self)(GTuple[0], GTuple[1]);
    }
    void __setitem__(int GTuple[2], double value) {
        (*self)(GTuple[0], GTuple[1]) = value;
    }
    GMatrixSparse __mul__(const double &a) {
        return (*self) * a;
    }
    GMatrixSparse __div__(const double &a) {
        return (*self) / a;
    }
    GMatrixSparse copy() {
        return (*self);
    }
    GMatrixSparse set(const double& value) {
        (*self) = value;
        return (*self);
    }
};
