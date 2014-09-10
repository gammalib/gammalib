/***************************************************************************
 *                GMatrixBase.hpp - Abstract matrix base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2014 by Juergen Knoedlseder                         *
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
 * @file GMatrixBase.hpp
 * @brief Abstract matrix base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIXBASE_HPP
#define GMATRIXBASE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GMatrixBase
 *
 * @brief Abstract matrix base class interface definition
 *
 * This class provides the abstract interface for all matrix classes. A
 * GammaLib matrix is a two dimensional array of double precision values
 * that is commonly used for linear algebra computations.
 *
 * The matrix classes were designed in a manner so that matrices can be used
 * in a natural syntax that is comparable to the syntax used for scalar
 * values. In particular, the following matrix operations are available:
 *
 *     assigned_matrix     = matrix;
 *     assigned_matrix     = value;
 *     added_matrix        = matrix1 + matrix2;
 *     added_matrix       += matrix;
 *     subtracted_matrix   = matrix1 - matrix2;
 *     subtracted_matrix  -= matrix;
 *     multiplied_matrix   = matrix1 * matrix2;
 *     multiplied_matrix  *= matrix;
 *     scaled_matrix       = matrix * scale;
 *     scaled_matrix       = scale * matrix;
 *     scaled_matrix      *= scale;
 *     divided_matrix      = matrix / factor;
 *     divided_matrix     /= factor;
 *     negated_matrix      = -matrix;
 *
 * where @p matrix, @p matrix1 and @p matrix2 are matrices, and
 * @p value, @p scale and @p factor are double precision values. For the
 * value assignment operator, each matrix element is assigned the specified
 * value. For the scaling and division operators, every matrix element is
 * multiplied or divided by the specified value.
 *
 * Matrices can also be compared using
 *
 *     if (matrix == another_matrix) ...
 *     if (matrix != another_matrix) ...
 *
 * where identify is defined as having the same number of rows and columns
 * and every element having the same value.
 *
 * Matrix elements are accessed using the access operators or the at() methods
 *
 *     double element = matrix(row, column);
 *     double element = matrix.at(row, column);
 *
 * the difference being that the access operators do NOT but the at()
 * methods do perform validity checking of the @p row and @p column
 * arguments. Access operators and at() methods exist for non-const and
 * const references.
 *
 * For more efficient matrix access, methods operating on matrix rows and
 * columns are available. The
 *
 *     GVector row    = matrix.row(index);
 *     GVector column = matrix.column(index);
 *
 * methods return a row or a column in a GVector object, while the
 *
 *     GVector row;
 *     GVector column;
 *     ...
 *     matrix.row(index, row);
 *     matrix.column(index, column);
 *
 * store a row or a column provided by a GVector object into the matrix.
 * Obviously, the vector needs to be of the same size as the number of rows
 * or columns available in the matrix. Adding of matrix elements row-by-row
 * or column-by-column is done using the
 *
 *     matrix.add_to_row(index, row);
 *     matrix.add_to_column(index, column);
 *
 * methods.
 *
 * Matrix multiplication with a vector is done using
 *
 *     GVector vector;
 *     ...
 *     GVector product = matrix * vector;
 *
 * where the size of @p vector needs to be identical to the number of
 * columns in the matrix.
 *
 * Furthermore, the following methods exist for matrix transformations:
 *
 *     new_matrix = matrix.transpose();    // Transposes matrix
 *     new_matrix = matrix.invert();       // Inverts matrix
 *     new_vector = matrix.solve(vector);  // Solves linear matrix equation
 *     new_matrix = matrix.abs();          // Returns matrix with all elements replaced by their absolute values
 *
 * The following methods allow to access matrix attributes:
 *
 *     matrix.fill();     // Percentage of non-zero elements in matrix [0,1]
 *     matrix.min();      // Minimum matrix element
 *     matrix.max();      // Maximum matrix element
 *     matrix.sum();      // Sum of all matrix elements
 *     matrix.rows();     // Number of rows in matrix
 *     matrix.columns();  // Number of columns in matrix
 *     matrix.size();     // Number of elements in matrix
 ***************************************************************************/
class GMatrixBase : public GBase {

public:
    // Constructors and destructors
    GMatrixBase(void);
    GMatrixBase(const GMatrixBase& matrix);
    virtual ~GMatrixBase(void);

    // Pure virtual operators
    virtual double&       operator()(const int& row, const int& column) = 0;
    virtual const double& operator()(const int& row, const int& column) const = 0;
    virtual GVector       operator*(const GVector& vector) const = 0;

    // Base class operators
    virtual GMatrixBase&  operator=(const GMatrixBase& matrix);
    virtual bool          operator==(const GMatrixBase& matrix) const;
    virtual bool          operator!=(const GMatrixBase& matrix) const;

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GMatrixBase*  clone(void) const = 0;
    virtual std::string   classname(void) const = 0;
    virtual double&       at(const int& row, const int& column) = 0;
    virtual const double& at(const int& row, const int& column) const = 0;
    virtual GVector       row(const int& row) const = 0;
    virtual void          row(const int& row, const GVector& vector) = 0;
    virtual GVector       column(const int& column) const = 0;
    virtual void          column(const int& column, const GVector& vector) = 0;
    virtual void          add_to_row(const int& row, const GVector& vector) = 0;
    virtual void          add_to_column(const int& column, const GVector& vector) = 0;
    virtual double        fill(void) const = 0;
    virtual double        min(void) const = 0;
    virtual double        max(void) const = 0;
    virtual double        sum(void) const = 0;
    virtual std::string   print(const GChatter& chatter = NORMAL) const = 0;

    // Base class methods
    const int& size(void) const;
    const int& columns(void) const;
    const int& rows(void) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GMatrixBase& matrix);
    void        free_members(void);
    void        select_non_zero(void);
    void        scale_elements(const double& scalar);
    void        set_all_elements(const double& value);
    double      get_min_element(void) const;
    double      get_max_element(void) const;
    double      get_element_sum(void) const;
    std::string print_elements(const GChatter& chatter = NORMAL,
                               const int&      num = 10) const;
    std::string print_row_compression(const GChatter& chatter = NORMAL) const;
    std::string print_col_compression(const GChatter& chatter = NORMAL) const;

    // Protected data area
    int     m_rows;       //!< Number of rows
    int     m_cols;       //!< Number of columns
    int     m_elements;   //!< Number of elements stored in matrix
    int     m_alloc;      //!< Size of allocated matrix memory
    int     m_num_rowsel; //!< Number of selected rows (for compressed decomposition)
    int     m_num_colsel; //!< Number of selected columns (for compressed decomposition)
    int*    m_colstart;   //!< Column start indices (m_cols+1)
    int*    m_rowsel;     //!< Row selection (for compressed decomposition)
    int*    m_colsel;     //!< Column selection (for compressed decomposition)
    double* m_data;       //!< Matrix data
};


/***********************************************************************//**
 * @brief Return number of matrix elements
 *
 * @return Number of matrix elements.
 *
 * Returns the number of matrix elements that are stored in the matrix.
 ***************************************************************************/
inline
const int& GMatrixBase::size(void) const
{
    return (m_elements);
}


/***********************************************************************//**
 * @brief Return number of matrix columns
 *
 * @return Number of matrix columns.
 *
 * Returns the number of matrix columns.
 ***************************************************************************/
inline
const int& GMatrixBase::columns(void) const
{
    return (m_cols);
}


/***********************************************************************//**
 * @brief Return number of matrix rows
 *
 * @return Number of matrix rows.
 *
 * Returns the number of matrix rows.
 ***************************************************************************/
inline
const int& GMatrixBase::rows(void) const
{
    return (m_rows);
}

#endif /* GMATRIXBASE_HPP */
