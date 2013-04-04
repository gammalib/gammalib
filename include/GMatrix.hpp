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


/***********************************************************************//**
 * @class GMatrix
 *
 * @brief Generic matrix class defintion
 *
 * This class implements a generic matrix class. This class is a
 * non-spezialized representation of a matrix, and all other matrix storage
 * classes can be converted into that class.
 *
 * The matrix stores the elements column-wise.
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
    virtual std::string print(void) const;

    // Other methods
    virtual GMatrix extract_lower_triangle(void) const;
    virtual GMatrix extract_upper_triangle(void) const;
    virtual void    eulerx(const double& angle);
    virtual void    eulery(const double& angle);
    virtual void    eulerz(const double& angle);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrix& matrix);
    void free_members(void);
    void alloc_members(const int& rows, const int& columns);
};


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
 * @brief Negate matrix elements
 *
 * @return Matrix with negated elements.
 ***************************************************************************/
inline
GMatrix GMatrix::operator-(void) const
{
    GMatrix result = *this;
    result.negate();
    return result;
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 *
 * @return Minimum element in matrix.
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
 ***************************************************************************/
inline 
GMatrix operator/(const GMatrix& matrix, const double& scalar)
{
    GMatrix result = matrix;
    result /= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Return transpose of matrix
 *
 * @param[in] matrix Matrix.
 * @return Transpose of matrix.
 ***************************************************************************/
inline
GMatrix transpose(const GMatrix& matrix)
{
    GMatrix result = matrix;
    result.transpose();
    return result;
}


/***********************************************************************//**
 * @brief Return inverse of matrix
 *
 * @param[in] matrix Matrix.
 * @return Inverse of matrix.
 ***************************************************************************/
inline
GMatrix invert(const GMatrix& matrix)
{
    GMatrix result = matrix;
    result.invert();
    return result;
}


/***********************************************************************//**
 * @brief Return matrix with absolute values of all elements
 *
 * @param[in] matrix Matrix.
 * @return Matrix with elements being the absolute elements of the input
 *         matrix.
 ***************************************************************************/
inline
GMatrix abs(const GMatrix& matrix)
{
    GMatrix result = matrix;
    result.abs();
    return result;
}

#endif /* GMATRIX_HPP */
