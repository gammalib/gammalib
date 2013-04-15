/***************************************************************************
 *                      GMatrix.hpp - General matrix class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
 * @file GMatrix.hpp
 * @brief Generic matrix class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIX_HPP
#define GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GMatrixBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GMatrixSymmetric;
class GMatrixSparse;
class GVector;


/***********************************************************************//**
 * @class GMatrix
 *
 * @brief Generic matrix class defintion
 *
 * This class implements a generic matrix. The class is a non-spezialized
 * representation of a matrix, and all other matrix storage classes can be
 * converted into that class.
 *
 * For a description of the common matrix methods, please refer to the
 * GMatrixBase class.
 *
 * Matrix allocation is done using the constructors
 *
 *     GMatrix matrix(rows, columns);
 *     GMatrix matrix(matrix);
 *     GMatrix matrix(sparsematrix);
 *     GMatrix matrix(symmetricmatrix);
 *
 * where @p rows and @p columns specify the number of rows and columns of
 * the matrix. Storage conversion constructors exist that allow allocating
 * a generic matrix by copying from a sparse matrix of type GMatrixSparse
 * and a symmetrix matrix of type GMatrixSymmetric.
 *
 * To support using the GMatrix for coordinate transformations, methods
 * are available to compute matrices for the three Euler angles:
 *
 *     matrix.eulerx(angle);
 *     matrix.eulery(angle);
 *     matrix.eulerz(angle);
 *
 * Methods are also available to extract the lower or the upper triangle of
 * a matrix:
 *
 *     matrix.extract_lower_triangle();
 *     matrix.extract_upper_triangle();
 *
 * Matrix elements are stored column-wise by the class.
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
    virtual double&       operator()(const int& row, const int& column);
    virtual const double& operator()(const int& row, const int& column) const;
    virtual GVector       operator*(const GVector& vector) const;

    // Other operators
    virtual GMatrix&      operator=(const GMatrix& matrix);
    virtual GMatrix&      operator=(const double& value);
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
    virtual void          clear(void);
    virtual GMatrix*      clone(void) const;
    virtual double&       at(const int& row, const int& column);
    virtual const double& at(const int& row, const int& column) const;
    virtual GVector       row(const int& row) const;
    virtual void          row(const int& row, const GVector& vector);
    virtual GVector       column(const int& column) const;
    virtual void          column(const int& column, const GVector& vector);
    virtual void          add_to_row(const int& row, const GVector& vector);
    virtual void          add_to_column(const int& column, const GVector& vector);
    virtual double        fill(void) const;
    virtual double        min(void) const;
    virtual double        max(void) const;
    virtual double        sum(void) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GMatrix transpose(void) const;
    GMatrix invert(void) const;
    GVector solve(const GVector& vector) const;
    GMatrix negate(void) const;
    GMatrix abs(void) const;
    GMatrix extract_lower_triangle(void) const;
    GMatrix extract_upper_triangle(void) const;
    void    eulerx(const double& angle);
    void    eulery(const double& angle);
    void    eulerz(const double& angle);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrix& matrix);
    void free_members(void);
    void alloc_members(const int& rows, const int& columns);
};


/***********************************************************************//**
 * @brief Return reference to matrix element
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 * @return Reference to matrix element.
 *
 * Returns a reference to the matrix element at @p row and @p column.
 ***************************************************************************/
inline
double& GMatrix::operator()(const int& row, const int& column)
{
    return (m_data[m_colstart[column]+row]);
}


/***********************************************************************//**
 * @brief Return reference to matrix element (const version)
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 * @return Const reference to matrix element.
 *
 * Returns a const reference to the matrix element at @p row and @p column.
 ***************************************************************************/
inline
const double& GMatrix::operator()(const int& row, const int& column) const
{
    return (m_data[m_colstart[column]+row]);
}


/***********************************************************************//**
 * @brief Binary matrix addition
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix addition.
 *
 * Returns the sum of two matrices. The method makes use of the unary
 * addition operator.
 ***************************************************************************/
inline
GMatrix GMatrix::operator+(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result += matrix;
    return result;
}


/***********************************************************************//**
 * @brief Binary matrix subtraction
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix subtraction.
 *
 * Returns the difference between two matrices. The method makes use of the
 * unary subtraction operator.
 ***************************************************************************/
inline
GMatrix GMatrix::operator-(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result -= matrix;
    return result;
}


/***********************************************************************//**
 * @brief Binary matrix multiplication
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix multiplication.
 *
 * Returns the product of two matrices. The method makes use of the unary
 * multiplication operator.
 ***************************************************************************/
inline
GMatrix GMatrix::operator*(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result *= matrix;
    return result;
}


/***********************************************************************//**
 * @brief Scale matrix elements
 *
 * @param[in] scalar Scale factor.
 * @return Matrix with elements multiplied by @p scalar.
 *
 * Returns a matrix where all elements have been multiplied by the specified
 * @p scalar value.
 ***************************************************************************/
inline
GMatrix& GMatrix::operator*=(const double& scalar)
{
    scale_elements(scalar);
    return *this;
}


/***********************************************************************//**
 * @brief Divide matrix elements
 *
 * @param[in] scalar Scalar.
 * @return Matrix with elements divided by @p scalar.
 *
 * Returns a matrix where all elements have been divided by the specified
 * @p scalar value.
 ***************************************************************************/
inline
GMatrix& GMatrix::operator/=(const double& scalar)
{
    *this *= 1.0/scalar;
    return *this;
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 *
 * @return Minimum element in matrix.
 *
 * Returns the smallest element in the matrix.
 ***************************************************************************/
inline
double GMatrix::min(void) const
{
    return get_min_element();
}


/***********************************************************************//**
 * @brief Return maximum matrix element
 *
 * @return Maximum element in matrix.
 *
 * Returns the largest element in the matrix.
 ***************************************************************************/
inline
double GMatrix::max(void) const
{
    return get_max_element();
}


/***********************************************************************//**
 * @brief Return matrix element sum
 *
 * @return Sum of all matrix elements.
 *
 * Returns the sum of all matrix elements.
 ***************************************************************************/
inline
double GMatrix::sum(void) const
{
    return get_element_sum();
}


/***********************************************************************//**
 * @brief Multiply matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 *
 * Returns a matrix where each element has been multiplied by a @p scalar.
 ***************************************************************************/
inline
GMatrix operator*(const GMatrix& matrix, const double& scalar)
{
    GMatrix result = matrix;
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Multiply matrix by scalar
 *
 * @param[in] scalar Scalar.
 * @param[in] matrix Matrix.
 * @return Matrix divided by @p scalar.
 *
 * Returns a matrix where each element has been multiplied by a @p scalar.
 ***************************************************************************/
inline
GMatrix operator*(const double& scalar, const GMatrix& matrix)
{
    GMatrix result = matrix;
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Divide matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 *
 * Returns a matrix where each element has been divided by a @p scalar.
 ***************************************************************************/
inline 
GMatrix operator/(const GMatrix& matrix, const double& scalar)
{
    GMatrix result = matrix;
    result /= scalar;
    return result;
}

#endif /* GMATRIX_HPP */
