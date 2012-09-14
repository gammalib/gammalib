/***************************************************************************
 *            GSparseMatrix.i  -  Sparse Matrix class SWIG file            *
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
 * @file GSparseMatrix.i
 * @brief Sparse matrix class definition.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSparseMatrix.hpp"
#include "GTools.hpp"
%}

/* __ Includes ___________________________________________________________ */
%include "GTypemaps.i"


/***********************************************************************//**
 * @class GSparseMatrix
 *
 * @brief Sparse matrix class definition
 *
 * GSparseMatrix implements a sparse matrix storage class. It derives from
 * the abstract base class GMatrixBase.
 ***************************************************************************/
class GSparseMatrix : public GMatrixBase {
public:
    // Constructors and destructors
    GSparseMatrix(void);
    GSparseMatrix(const int& rows, const int& cols, const int& elements = 0);
    GSparseMatrix(const GMatrix& matrix);
    GSparseMatrix(const GSymMatrix& matrix);
    GSparseMatrix(const GSparseMatrix& matrix);
    virtual ~GSparseMatrix(void);

    // Implemented pure virtual base class operators
    GVector        operator*(const GVector& vector) const;

    // Overloaded virtual base class operators
    bool           operator==(const GSparseMatrix& matrix) const;
    bool           operator!=(const GSparseMatrix& matrix) const;

    // Other operators
    GSparseMatrix  operator+(const GSparseMatrix& matrix) const;
    GSparseMatrix  operator-(const GSparseMatrix& matrix) const;
    GSparseMatrix  operator*(const GSparseMatrix& matrix) const;
    GSparseMatrix  operator-(void) const;
    GSparseMatrix& operator+=(const GSparseMatrix& matrix);
    GSparseMatrix& operator-=(const GSparseMatrix& matrix);
    GSparseMatrix& operator*=(const GSparseMatrix& matrix);
    GSparseMatrix& operator*=(const double& scalar);
    GSparseMatrix& operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    void        clear(void);
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
    void        add_col(const double* values, const int* rows,
                        int number, const int& col);
    void        insert_col(const double* values, const int* rows,
                           int number, const int& col);
    void        cholesky_decompose(bool compress = true);
    GVector     cholesky_solver(const GVector& vector, bool compress = true);
    void        cholesky_invert(bool compress = true);
    void        set_mem_block(const int& block);
    void        stack_init(const int& size = 0, const int& entries = 0);
    int         stack_push_column(const GVector& vector, const int& col);
    int         stack_push_column(const double* values, const int* rows,
                                  const int& number, const int& col);
    void        stack_flush(void);
    void        stack_destroy(void);
};


/***********************************************************************//**
 * @brief GSparseMatrix class extension
 ***************************************************************************/
%extend GSparseMatrix {
    char *__str__() {
        return tochar(self->print());
    }
    double __getitem__(int GTuple[2]) {
        return (*self)(GTuple[0], GTuple[1]);
    }
    void __setitem__(int GTuple[2], double value) {
        (*self)(GTuple[0], GTuple[1]) = value;
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
