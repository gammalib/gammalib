/***************************************************************************
 *                GMatrixSymmetric.hpp - Symmetric matrix class            *
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
 * @file GMatrixSymmetric.hpp
 * @brief Symmetric matrix class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIXSYMMETRIC_HPP
#define GMATRIXSYMMETRIC_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GMatrixBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GMatrix;
class GMatrixSparse;


/***********************************************************************//**
 * @class GMatrixSymmetric
 *
 * @brief Symmetric matrix class interface defintion
 *
 * This class implements a symmetric matrix class. Only one triangle of the
 * matrix is physically stored, reducing the memory requirements and
 * imposing strict matrix symmetry.
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

    // Implemented pure virtual base class operators
    virtual double&           operator()(const int& row, const int& column);
    virtual const double&     operator()(const int& row, const int& column) const;
    virtual GVector           operator*(const GVector& v) const;

    // Other operators
    virtual GMatrixSymmetric& operator=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric  operator+(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric  operator-(const GMatrixSymmetric& matrix) const;
    virtual GMatrix           operator*(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric  operator-(void) const;
    virtual GMatrixSymmetric& operator+=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric& operator-=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric& operator*=(const double& scaler);
    virtual GMatrixSymmetric& operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void              clear(void);
    virtual GMatrixSymmetric* clone(void) const;
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
    virtual std::string       print(void) const;

    // Other methods
    virtual GMatrix extract_lower_triangle(void) const;
    virtual GMatrix extract_upper_triangle(void) const;
    virtual void    cholesky_decompose(bool compress = true);
    virtual GVector cholesky_solver(const GVector& vector, bool compress = true);
    virtual void    cholesky_invert(bool compress = true);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrixSymmetric& matrix);
    void free_members(void);
    void alloc_members(const int& rows, const int& columns);
    void set_inx(void);

    // Private data area
    int  m_num_inx;          //!< Number of indices in array
    int* m_inx;              //!< Index array of non-zero rows/columns
};


/***********************************************************************//**
 * @brief Transpose matrix
 *
 * Transpose matrix. As the transposed matrix of a symmetric matrix is
 * identical to the original matrix, this method simply returns without
 * taking any action.
 ***************************************************************************/
inline
void GMatrixSymmetric::transpose(void)
{
    return;
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 *
 * @return Minimum matrix element
 ***************************************************************************/
inline
double GMatrixSymmetric::min(void) const
{
    return (get_min_element());
}


/***********************************************************************//**
 * @brief Return maximum matrix element
 *
 * @return Maximum matrix element
 ***************************************************************************/
inline
double GMatrixSymmetric::max(void) const
{
    return (get_max_element());
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
GMatrixSymmetric GMatrixSymmetric::operator+(const GMatrixSymmetric& matrix) const
{
    GMatrixSymmetric result = *this;
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
GMatrixSymmetric GMatrixSymmetric::operator-(const GMatrixSymmetric& matrix) const
{
    GMatrixSymmetric result = *this;
    result -= matrix;
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
GMatrixSymmetric& GMatrixSymmetric::operator*=(const double& scalar)
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
GMatrixSymmetric& GMatrixSymmetric::operator/=(const double& scalar)
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
GMatrixSymmetric GMatrixSymmetric::operator-(void) const
{
    GMatrixSymmetric result = *this;
    result.negate();
    return result;
}


/***********************************************************************//**
 * @brief Multiply matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 ***************************************************************************/
inline 
GMatrixSymmetric operator*(const GMatrixSymmetric& matrix, const double& scalar)
{
    GMatrixSymmetric result = matrix;
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
GMatrixSymmetric operator*(const double& scalar, const GMatrixSymmetric& matrix)
{
    GMatrixSymmetric result = matrix;
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
GMatrixSymmetric operator/(const GMatrixSymmetric& matrix, const double& scalar)
{
    GMatrixSymmetric result = matrix;
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
GMatrixSymmetric transpose(const GMatrixSymmetric& matrix)
{
    return matrix;
}


/***********************************************************************//**
 * @brief Return inverse of matrix
 *
 * @param[in] matrix Matrix.
 * @return Inverse of matrix.
 ***************************************************************************/
inline
GMatrixSymmetric invert(const GMatrixSymmetric& matrix)
{
    GMatrixSymmetric result = matrix;
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
GMatrixSymmetric abs(const GMatrixSymmetric& matrix)
{
    GMatrixSymmetric result = matrix;
    result.abs();
    return result;
}


/***********************************************************************//**
 * @brief Return Cholesky decomposition of matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] compress Use matrix compression (defaults to true).
 * @return Cholesky decomposition of matrix.
 ***************************************************************************/
inline
GMatrixSymmetric cholesky_decompose(const GMatrixSymmetric& matrix, bool compress)
{
    GMatrixSymmetric result = matrix;
    result.cholesky_decompose(compress);
    return result;
}


/***********************************************************************//**
 * @brief Return inverse matrix using Cholesky decomposition
 *
 * @param[in] matrix Matrix.
 * @param[in] compress Use matrix compression (defaults to true).
 * @return Inverse of matrix.
 ***************************************************************************/
inline
GMatrixSymmetric cholesky_invert(const GMatrixSymmetric& matrix, bool compress)
{
    GMatrixSymmetric result = matrix;
    result.cholesky_invert(compress);
    return result;
}

#endif /* GMATRIXSYMMETRIC_HPP */
