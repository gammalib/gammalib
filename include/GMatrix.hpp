/***************************************************************************
 *                         GMatrix.hpp  -  matrix class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2012 by Juergen Knoedlseder                         *
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
 * @brief Generic matrix class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIX_HPP
#define GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GMatrixBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GSymMatrix;
class GSparseMatrix;


/***********************************************************************//**
 * @class GMatrix
 *
 * @brief Generic matrix class defintion
 *
 * This class implements a generic matrix class. This class is a
 * non-spezialized representation of a matrix, and all other matrix storage
 * classes can be converted into that class.
 ***************************************************************************/
class GMatrix : public GMatrixBase {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GMatrix& matrix);
    friend GLog&         operator<<(GLog& log,        const GMatrix& matrix);

    // Binary operator friends
    friend GMatrix operator*(const double& a,  const GMatrix& b);
    friend GMatrix operator*(const GMatrix& a, const double& b);
    friend GMatrix operator/(const GMatrix& a, const double& b);

    // Friend functions
    friend GMatrix transpose(const GMatrix& matrix);
    friend GMatrix invert(const GMatrix& matrix);
    friend GMatrix abs(const GMatrix& matrix);

public:
    // Constructors and destructors
    GMatrix(void);
    GMatrix(const int& rows, const int& cols);
    GMatrix(const GMatrix& matrix);
    GMatrix(const GSymMatrix& matrix);
    GMatrix(const GSparseMatrix& matrix);
    virtual ~GMatrix(void);

    // Implemented pure virtual base class operators
    double&       operator()(const int& row, const int& col);
    const double& operator()(const int& row, const int& col) const;
    GVector       operator*(const GVector& vector) const;

    // Other operators
    GMatrix&      operator=(const GMatrix& matrix);
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
    std::string print(void) const;

    // Other methods
    GMatrix     extract_lower_triangle(void) const;
    GMatrix     extract_upper_triangle(void) const;
    void        eulerx(const double& angle);
    void        eulery(const double& angle);
    void        eulerz(const double& angle);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrix& matrix);
    void free_members(void);
    void alloc_members(const int& rows, const int& cols);
};


/***************************************************************************
 *                            Inline operators                             *
 ***************************************************************************/
// Binary matrix addition
inline
GMatrix GMatrix::operator+(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result += matrix;
    return result;
}

// Binary matrix subtraction
inline
GMatrix GMatrix::operator-(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result -= matrix;
    return result;
}

// Binary matrix multiplication
inline
GMatrix GMatrix::operator*(const GMatrix& matrix) const
{
    GMatrix result = *this;
    result *= matrix;
    return result;
}

// Matrix scaling
inline
GMatrix& GMatrix::operator*=(const double& scalar)
{
    multiplication(scalar);
    return *this;
}

// Matrix scalar division
inline
GMatrix& GMatrix::operator/=(const double& scalar)
{
    double inverse = 1.0/scalar;
    multiplication(inverse);
    return *this;
}

// Negation
inline
GMatrix GMatrix::operator-(void) const
{
    GMatrix result = *this;
    result.negation();
    return result;
}


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline void   GMatrix::clear(void) { set_all_elements(0.0); }
inline double GMatrix::min(void) const { return get_min_element(); }
inline double GMatrix::max(void) const { return get_max_element(); }
inline double GMatrix::sum(void) const { return get_element_sum(); }


/***************************************************************************
 *                              Inline friends                             *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline
GMatrix operator* (const GMatrix& a, const double& b)
{
    GMatrix result = a;
    result *= b;
    return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GMatrix operator* (const double& a, const GMatrix& b)
{
    GMatrix result = b;
    result *= a;
    return result;
}

// Binary matrix division (matrix is left operand)
inline 
GMatrix operator/ (const GMatrix& a, const double& b)
{
    GMatrix result = a;
    result /= b;
    return result;
}

// Matrix transpose function
inline
GMatrix transpose(const GMatrix& matrix)
{
    GMatrix result = matrix;
    result.transpose();
    return result;
}

// Matrix inversion
inline
GMatrix invert(const GMatrix& matrix)
{
    GMatrix result = matrix;
    result.invert();
    return result;
}

#endif /* GMATRIX_HPP */
