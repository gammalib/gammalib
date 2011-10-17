/***************************************************************************
 *                 GSymMatrix.hpp  -  symmetric matrix class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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
 * @file GSymMatrix.hpp
 * @brief GSymMatrix class definition.
 * @author J. Knodlseder
 */

#ifndef GSYMMATRIX_HPP
#define GSYMMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GMatrixBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GMatrix;
class GSparseMatrix;


/***********************************************************************//**
 * @class GSymMatrix
 *
 * @brief GSymMatrix class interface defintion
 *
 * Implements a sparse matrix storage class. It derives from the abstract
 * base class GMatrixBase.
 ***************************************************************************/
class GSymMatrix : public GMatrixBase {

    // Binary operator friends
    friend GSymMatrix operator* (const double& a,  const GSymMatrix& b);
    friend GSymMatrix operator* (const GSymMatrix& a, const double& b);
    friend GSymMatrix operator/ (const GSymMatrix& a, const double& b);

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GSymMatrix& m);

    // Friend functions
    friend GSymMatrix transpose(const GSymMatrix& m);
    friend GSymMatrix abs(const GSymMatrix& m);
    friend GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress = 1);
    friend GSymMatrix cholesky_invert(const GSymMatrix& m, int compress = 1);

public:
    // Constructors and destructors (not inherited)
    GSymMatrix(int rows, int cols);
    GSymMatrix(const GMatrix& m);
    GSymMatrix(const GSymMatrix& m);
    GSymMatrix(const GSparseMatrix& m);
    virtual ~GSymMatrix();

    // Operators
    virtual GSymMatrix&   operator= (const GSymMatrix& m);
    virtual double&       operator() (int row, int col);
    virtual const double& operator() (int row, int col) const;
    virtual GSymMatrix    operator+ (const GSymMatrix& m) const;
    virtual GSymMatrix    operator- (const GSymMatrix& m) const;
    virtual GSymMatrix    operator* (const GSymMatrix& m) const;
    virtual GVector       operator* (const GVector& v) const;
    //_USE_BASE virtual int           operator== (const GSymMatrix& m) const;
    //_USE_BASE virtual int           operator!= (const GSymMatrix& m) const;
    virtual GSymMatrix    operator- () const;
    virtual GSymMatrix&   operator+= (const GSymMatrix& m);
    virtual GSymMatrix&   operator-= (const GSymMatrix& m);
    virtual GSymMatrix&   operator*= (const GSymMatrix& m);
    virtual GSymMatrix&   operator*= (const double& d);
    virtual GSymMatrix&   operator/= (const double& d);

    // Methods
    virtual void    add_col(const GVector& v, int col);
    virtual void    cholesky_decompose(int compress = 1);
    virtual GVector cholesky_solver(const GVector& v, int compress = 1);
    virtual void    cholesky_invert(int compress = 1);
    virtual void    clear();
    virtual GVector extract_row(int row) const;
    virtual GVector extract_col(int col) const;
    virtual GMatrix extract_lower_triangle() const;
    virtual GMatrix extract_upper_triangle() const;
    virtual void    insert_col(const GVector& v, int col);
    virtual double  fill() const;
    virtual double  min() const;
    virtual double  max() const;
    virtual double  sum() const;
    virtual void    transpose() { return; }

private:
    // Private methods
    void constructor(int rows, int cols);
    void init_members(void);
    void copy_members(const GSymMatrix& m);
    void free_members(void);
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
GSymMatrix GSymMatrix::operator+ (const GSymMatrix& m) const
{
    GSymMatrix result = *this;
    result += m;
    return result;
}

// Binary matrix subtraction
inline
GSymMatrix GSymMatrix::operator- (const GSymMatrix& m) const
{
    GSymMatrix result = *this;
    result -= m;
    return result;
}

// Binary matrix multiplication
inline
GSymMatrix GSymMatrix::operator* (const GSymMatrix& m) const
{
    GSymMatrix result = *this;
    result *= m;
    return result;
}

// Matrix scaling
inline
GSymMatrix& GSymMatrix::operator*= (const double& s)
{
    multiplication(s);
    return *this;
}

// Matrix scalar division
inline
GSymMatrix& GSymMatrix::operator/= (const double& s)
{
    double inverse = 1.0/s;
    multiplication(inverse);
    return *this;
}

// Negation
inline
GSymMatrix GSymMatrix::operator- ( ) const
{
    GSymMatrix result = *this;
    result.negation();
    return result;
}


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline void   GSymMatrix::clear() { set_all_elements(0.0); }
inline double GSymMatrix::min() const { return get_min_element(); }
inline double GSymMatrix::max() const { return get_max_element(); }


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Binary matrix scaling (matrix is left operand)
inline 
GSymMatrix operator* (const GSymMatrix& a, const double& b)
{
    GSymMatrix result = a;
    result *= b;
    return result;
}

// Binary matrix scaling (matrix is right operand)
inline
GSymMatrix operator* (const double& a, const GSymMatrix& b)
{
    GSymMatrix result = b;
    result *= a;
    return result;
}

// Binary matrix division (matrix is left operand)
inline 
GSymMatrix operator/ (const GSymMatrix& a, const double& b)
{
    GSymMatrix result = a;
    result /= b;
    return result;
}

// Matrix transpose function
inline
GSymMatrix transpose(const GSymMatrix& m)
{
    return m;
}

// Cholesky decomposition
inline
GSymMatrix cholesky_decompose(const GSymMatrix& m, int compress)
{
    GSymMatrix result = m;
    result.cholesky_decompose(compress);
    return result;
}

// Matrix inversion using Cholesky decomposition
inline
GSymMatrix cholesky_invert(const GSymMatrix& m, int compress)
{
    GSymMatrix result = m;
    result.cholesky_invert(compress);
    return result;
}

#endif /* GSYMMATRIX_HPP */
