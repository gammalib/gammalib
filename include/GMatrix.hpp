/***************************************************************************
 *                         GMatrix.hpp  -  matrix class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMatrix.hpp
 * @brief GMatrix class definition.
 * @author J. Knodlseder
 */

#ifndef GMATRIX_HPP
#define GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GMatrixBase.hpp"


/***********************************************************************//**
 * @class GMatrix
 *
 * @brief GMatrix class interface defintion
 ***************************************************************************/
class GMatrix : public GMatrixBase {

    // Binary operator friends
    friend GMatrix operator* (const double& a,  const GMatrix& b);
    friend GMatrix operator* (const GMatrix& a, const double& b);
    friend GMatrix operator/ (const GMatrix& a, const double& b);

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GMatrix& m);

    // Friend functions
    friend GMatrix transpose(const GMatrix& m);
    friend GMatrix invert(const GMatrix& m);
    friend GMatrix fabs(const GMatrix& m);

public:
    // Constructors and destructors
    GMatrix(void);
    GMatrix(int rows, int cols);
    GMatrix(const GMatrix& m);
    ~GMatrix(void);

    // Operators
    GMatrix&      operator= (const GMatrix& m);
    double&       operator() (int row, int col);
    const double& operator() (int row, int col) const;
    GMatrix       operator+ (const GMatrix& m) const;
    GMatrix       operator- (const GMatrix& m) const;
    GMatrix       operator* (const GMatrix& m) const;
    GVector       operator* (const GVector& v) const;
    GMatrix       operator- () const;
    GMatrix&      operator+= (const GMatrix& m);
    GMatrix&      operator-= (const GMatrix& m);
    GMatrix&      operator*= (const GMatrix& m);
    GMatrix&      operator*= (const double& s);
    GMatrix&      operator/= (const double& s);

    // Methods
    void    clear(void);
    void    transpose(void);
    void    invert(void);
    void    add_col(const GVector& v, int col);
    void    insert_col(const GVector& v, int col);
    GVector extract_row(int row) const;
    GVector extract_col(int col) const;
    GMatrix extract_lower_triangle(void) const;
    GMatrix extract_upper_triangle(void) const;
    double  fill(void) const;
    double  min(void) const;
    double  max(void) const;
    double  sum(void) const;
    void    eulerx(const double& angle);
    void    eulery(const double& angle);
    void    eulerz(const double& angle);

private:
    // Private methods
    void constructor(int rows, int cols);
    void init_members(void) { ; }
    void copy_members(const GMatrix& m) { ; }
    void free_members(void) { ; }
};


/***************************************************************************
 *                            Inline operators                             *
 ***************************************************************************/
// Binary matrix addition
inline
GMatrix GMatrix::operator+ (const GMatrix& m) const
{
    GMatrix result = *this;
    result += m;
    return result;
}

// Binary matrix subtraction
inline
GMatrix GMatrix::operator- (const GMatrix& m) const
{
    GMatrix result = *this;
    result -= m;
    return result;
}

// Binary matrix multiplication
inline
GMatrix GMatrix::operator* (const GMatrix& m) const
{
    GMatrix result = *this;
    result *= m;
    return result;
}

// Matrix scaling
inline
GMatrix& GMatrix::operator*= (const double& s)
{
    multiplication(s);
    return *this;
}

// Matrix scalar division
inline
GMatrix& GMatrix::operator/= (const double& s)
{
    double inverse = 1.0/s;
    multiplication(inverse);
    return *this;
}

// Negation
inline
GMatrix GMatrix::operator- ( ) const
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
GMatrix transpose(const GMatrix& m)
{
    GMatrix result = m;
    result.transpose();
    return result;
}

// Matrix inversion
inline
GMatrix invert(const GMatrix& m)
{
    GMatrix result = m;
    result.invert();
    return result;
}

#endif /* GMATRIX_HPP */
