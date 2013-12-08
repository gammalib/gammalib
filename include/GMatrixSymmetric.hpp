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
class GVector;


/***********************************************************************//**
 * @class GMatrixSymmetric
 *
 * @brief Symmetric matrix class interface defintion
 *
 * This class implements a symmetric matrix. The class stores only a triangle
 * physically, imposing thus strict matrix symmetry.
 *
 * For a description of the common matrix methods, please refer to the
 * GMatrixBase class.
 *
 * Matrix allocation is done using the constructors
 *
 *     GMatrixSymmetric symmetricmatrix(rows, columns);
 *     GMatrixSymmetric symmetricmatrix(matrix);
 *     GMatrixSymmetric symmetricmatrix(sparsematrix);
 *     GMatrixSymmetric symmetricmatrix(symmetricmatrix);
 *
 * where @p rows and @p columns specify the number of rows and columns of
 * the matrix. Storage conversion constructors exist that allow allocating
 * a symmetric matrix by copying from a sparse matrix of type GMatrixSparse
 * and a general matrix of type GMatrix. Exceptions will be thrown in case
 * that the matrix from which the object should be allocated is not
 * symmetric.
 *
 * Methods are available to extract the lower or the upper triangle of the
 * matrix:
 *
 *     matrix.extract_lower_triangle();
 *     matrix.extract_upper_triangle();
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
    virtual GMatrixSymmetric& operator=(const double& value);
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
    virtual double&           at(const int& row, const int& column);
    virtual const double&     at(const int& row, const int& column) const;
    virtual GVector           row(const int& row) const;
    virtual void              row(const int& row, const GVector& vector);
    virtual GVector           column(const int& column) const;
    virtual void              column(const int& column, const GVector& vector);
    virtual void              add_to_row(const int& row, const GVector& vector);
    virtual void              add_to_column(const int& column, const GVector& vector);
    virtual double            fill(void) const;
    virtual double            min(void) const;
    virtual double            max(void) const;
    virtual double            sum(void) const;
    virtual std::string       print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GMatrixSymmetric transpose(void) const;
    GMatrixSymmetric invert(void) const;
    GVector          solve(const GVector& vector) const;
    GMatrixSymmetric abs(void) const;
    GMatrix          extract_lower_triangle(void) const;
    GMatrix          extract_upper_triangle(void) const;
    GMatrixSymmetric cholesky_decompose(const bool& compress = true) const;
    GVector          cholesky_solver(const GVector& vector, const bool& compress = true) const;
    GMatrixSymmetric cholesky_invert(const bool& compress = true) const;

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
 * @brief Return transposed matrix
 *
 * @return Transposed matrix.
 *
 * Returns transposed matrix of the matrix. As the transposed matrix of a
 * symmetric matrix is identical to the original matrix, this method simply
 * returns the actual matrix.
 ***************************************************************************/
inline
GMatrixSymmetric GMatrixSymmetric::transpose(void) const
{
    return (*this);
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 *
 * @return Minimum matrix element
 *
 * Returns the smallest element in the matrix.
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
 *
 * Returns the largest element in the matrix.
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
 * @brief Multiply matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 *
 * Returns a matrix where each element is multiplied by @p scalar.
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
 *
 * Returns a matrix where each element is multiplied by @p scalar.
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
 *
 * Returns a matrix where each element is divided by @p scalar.
 ***************************************************************************/
inline 
GMatrixSymmetric operator/(const GMatrixSymmetric& matrix, const double& scalar)
{
    GMatrixSymmetric result = matrix;
    result /= scalar;
    return result;
}

#endif /* GMATRIXSYMMETRIC_HPP */
