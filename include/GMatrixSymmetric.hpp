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

    // Binary operator friends
    friend GMatrixSymmetric operator*(const double& a,  const GMatrixSymmetric& b);
    friend GMatrixSymmetric operator*(const GMatrixSymmetric& a, const double& b);
    friend GMatrixSymmetric operator/(const GMatrixSymmetric& a, const double& b);

    // Friend functions
    friend GMatrixSymmetric transpose(const GMatrixSymmetric& matrix);
    friend GMatrixSymmetric abs(const GMatrixSymmetric& matrix);
    friend GMatrixSymmetric cholesky_decompose(const GMatrixSymmetric& matrix, bool compress = true);
    friend GMatrixSymmetric cholesky_invert(const GMatrixSymmetric& matrix, bool compress = true);

public:
    // Constructors and destructors
    GMatrixSymmetric(void);
    GMatrixSymmetric(const int& rows, const int& cols);
    GMatrixSymmetric(const GMatrix& matrix);
    GMatrixSymmetric(const GMatrixSparse& matrix);
    GMatrixSymmetric(const GMatrixSymmetric& matrix);
    virtual ~GMatrixSymmetric(void);

    // Implemented pure virtual base class operators
    virtual double&       operator()(const int& row, const int& col);
    virtual const double& operator()(const int& row, const int& col) const;
    virtual GVector       operator*(const GVector& v) const;

    // Other operators
    virtual GMatrixSymmetric&   operator=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric    operator+(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric    operator-(const GMatrixSymmetric& matrix) const;
    virtual GMatrix             operator*(const GMatrixSymmetric& matrix) const;
    virtual GMatrixSymmetric    operator-(void) const;
    virtual GMatrixSymmetric&   operator+=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric&   operator-=(const GMatrixSymmetric& matrix);
    virtual GMatrixSymmetric&   operator*=(const double& scaler);
    virtual GMatrixSymmetric&   operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void              clear(void);
    virtual GMatrixSymmetric* clone(void) const;
    virtual void              transpose(void) { return; }
    virtual void              invert(void);
    virtual void              add_col(const GVector& vector, const int& col);
    virtual void              insert_col(const GVector& vector, const int& col);
    virtual GVector           extract_row(const int& row) const;
    virtual GVector           extract_col(const int& col) const;
    virtual double            fill(void) const;
    virtual double            min(void) const;
    virtual double            max(void) const;
    virtual double            sum(void) const;
    virtual std::string       print(void) const;

    // Other methods
    virtual GMatrix     extract_lower_triangle(void) const;
    virtual GMatrix     extract_upper_triangle(void) const;
    virtual void        cholesky_decompose(bool compress = true);
    virtual GVector     cholesky_solver(const GVector& vector, bool compress = true);
    virtual void        cholesky_invert(bool compress = true);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrixSymmetric& matrix);
    void free_members(void);
    void alloc_members(const int& rows, const int& cols);
    void set_inx(void);

    // Private data area
    int  m_num_inx;          //!< Number of indices in array
    int* m_inx;              //!< Index array of non-zero rows/columns
};


/***************************************************************************
 *                            Inline operators                             *
 ***************************************************************************/
// Binary matrix addition
inline
GMatrixSymmetric GMatrixSymmetric::operator+(const GMatrixSymmetric& matrix) const
{
    GMatrixSymmetric result = *this;
    result += matrix;
    return result;
}

// Binary matrix subtraction
inline
GMatrixSymmetric GMatrixSymmetric::operator-(const GMatrixSymmetric& matrix) const
{
    GMatrixSymmetric result = *this;
    result -= matrix;
    return result;
}

// Matrix scaling
inline
GMatrixSymmetric& GMatrixSymmetric::operator*=(const double& scalar)
{
    multiplication(scalar);
    return *this;
}

// Matrix scalar division
inline
GMatrixSymmetric& GMatrixSymmetric::operator/=(const double& scalar)
{
    double inverse = 1.0/scalar;
    multiplication(inverse);
    return *this;
}

// Negation
inline
GMatrixSymmetric GMatrixSymmetric::operator-(void) const
{
    GMatrixSymmetric result = *this;
    result.negation();
    return result;
}


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline double GMatrixSymmetric::min(void) const { return get_min_element(); }
inline double GMatrixSymmetric::max(void) const { return get_max_element(); }


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline 
GMatrixSymmetric operator*(const GMatrixSymmetric& a, const double& b)
{
    GMatrixSymmetric result = a;
    result *= b;
    return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GMatrixSymmetric operator*(const double& a, const GMatrixSymmetric& b)
{
    GMatrixSymmetric result = b;
    result *= a;
    return result;
}

// Binary matrix division (matrix is left operand)
inline 
GMatrixSymmetric operator/(const GMatrixSymmetric& a, const double& b)
{
    GMatrixSymmetric result = a;
    result /= b;
    return result;
}

// Matrix transpose function
inline
GMatrixSymmetric transpose(const GMatrixSymmetric& matrix)
{
    return matrix;
}

// Cholesky decomposition
inline
GMatrixSymmetric cholesky_decompose(const GMatrixSymmetric& matrix, bool compress)
{
    GMatrixSymmetric result = matrix;
    result.cholesky_decompose(compress);
    return result;
}

// Matrix inversion using Cholesky decomposition
inline
GMatrixSymmetric cholesky_invert(const GMatrixSymmetric& matrix, bool compress)
{
    GMatrixSymmetric result = matrix;
    result.cholesky_invert(compress);
    return result;
}

#endif /* GMATRIXSYMMETRIC_HPP */
