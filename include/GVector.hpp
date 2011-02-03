/***************************************************************************
 *                         GVector.hpp  -  vector class                    *
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
 * @file GVector.hpp
 * @brief GVector class definition.
 * @author J. Knodlseder
 */

#ifndef GVECTOR_HPP
#define GVECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include <cmath>
#include "GException.hpp"
#include "GLog.hpp"


/***********************************************************************//**
 * @class GVector
 *
 * @brief GVector class interface defintion
 ***************************************************************************/
class GVector {

    // Friend classes
    friend class GMatrix;
    friend class GSymMatrix;
    friend class GSparseMatrix;

    // Operator friends
    friend GVector  operator+ (const GVector &a, const GVector &b);
    friend GVector  operator+ (const GVector &a, const double &b);
    friend GVector  operator+ (const double &a, const GVector &b);
    friend GVector  operator- (const GVector &a, const GVector &b);
    friend GVector  operator- (const GVector &a, const double &b);
    friend GVector  operator- (const double &a, const GVector &b);
    friend double   operator* (const GVector &a, const GVector &b);
    friend GVector  operator* (const GVector &a, const double &b);
    friend GVector  operator* (const double &a, const GVector &b);
    friend GVector  operator/ (const GVector &a, const double &b);
    friend int      operator== (const GVector &a, const GVector &b);
    friend int      operator!= (const GVector &a, const GVector &b);

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GVector& v);
    friend GLog&         operator<< (GLog& log, const GVector& v);

    // Friend functions
    friend GVector  cross(const GVector &a, const GVector &b);
    friend double   norm(const GVector &v);
    friend double   min(const GVector &v);
    friend double   max(const GVector &v);
    friend double   sum(const GVector &v);
    friend GVector  perm(const GVector &v, const int *p);
    friend GVector  iperm(const GVector &v, const int *p);
    friend GVector  acos(const GVector &v);
    friend GVector  acosh(const GVector &v);
    friend GVector  asin(const GVector &v);
    friend GVector  asinh(const GVector &v);
    friend GVector  atan(const GVector &v);
    friend GVector  atanh(const GVector &v);
    friend GVector  cos(const GVector &v);
    friend GVector  cosh(const GVector &v);
    friend GVector  exp(const GVector &v);
    friend GVector  abs(const GVector &v);
    friend GVector  log(const GVector &v);
    friend GVector  log10(const GVector &v);
    friend GVector  sin(const GVector &v);
    friend GVector  sinh(const GVector &v);
    friend GVector  sqrt(const GVector &v);
    friend GVector  tan(const GVector &v);
    friend GVector  tanh(const GVector &v);
    friend GVector  pow(const GVector &v, const double &q);

public:
    // Constructors and destructors
    GVector(void);
    explicit GVector(const int& num);
    explicit GVector(const double& a);
    explicit GVector(const double& a, const double& b);
    explicit GVector(const double& a, const double& b, const double& c);
    GVector(const GVector& v);
    virtual ~GVector(void);

    // Vector element access operators
    double& operator() (int inx);
    const double& operator() (int inx) const;

    // Vector operators
    GVector& operator= (const GVector& v);
    GVector& operator+= (const GVector& v);
    GVector& operator-= (const GVector& v);
    GVector& operator= (const double& v);
    GVector& operator+= (const double& v);
    GVector& operator-= (const double& v);
    GVector& operator*= (const double& v);
    GVector& operator/= (const double& v);
    GVector  operator- () const;

    // Vector methods
    void        clear(void);
    int         size(void) const;
    int         non_zeros(void) const;
    std::string print(void) const;

private:
    // Private methods
    void init_members(void);
    void alloc_members(void);
    void copy_members(const GVector& v);
    void free_members(void);

    // Private data area
    int     m_num;    //!< Number of elements in vector
    double* m_data;   //!< Vector array
};


/***************************************************************************
 *                               Inline members                            *
 ***************************************************************************/
// Vector unary addition operator
inline
GVector& GVector::operator+= (const GVector& v)
{
    if (m_num != v.m_num)
        throw GException::vector_mismatch("GVector::operator+=(GVector)",
                                          m_num, v.m_num);
    for (int i = 0; i < m_num; ++i)
        m_data[i] += v.m_data[i];
    return *this;
}

// Vector unary subtraction operator
inline
GVector& GVector::operator-= (const GVector& v)
{
    if (m_num != v.m_num)
        throw GException::vector_mismatch("GVector::operator-=(GVector)",
                                          m_num, v.m_num);
    for (int i = 0; i < m_num; ++i)
        m_data[i] -= v.m_data[i];
    return *this;
}

// Scalar assignment operator
inline
GVector& GVector::operator= (const double& v)
{
    for (int i = 0; i < m_num; ++i)
        m_data[i] = v;
    return *this;
}

// Scalar unary addition operator
inline
GVector& GVector::operator+= (const double& v)
{
    for (int i = 0; i < m_num; ++i)
        m_data[i] += v;
    return *this;
}

// Scalar unary subtraction operator
inline
GVector& GVector::operator-= (const double& v)
{
    for (int i = 0; i < m_num; ++i)
        m_data[i] -= v;
    return *this;
}

// Scalar unary multiplication operator
inline
GVector& GVector::operator*= (const double& v)
{
    for (int i = 0; i < m_num; ++i)
        m_data[i] *= v;
    return *this;
}

// Scalar unary division operator
inline
GVector& GVector::operator/= (const double& v)
{
    for (int i = 0; i < m_num; ++i)
        m_data[i] /= v;
    return *this;
}

// Unary minus operator
inline
GVector GVector::operator- ( ) const
{
    GVector result = *this;
    for (int i = 0; i < m_num; ++i)
        result.m_data[i] = -result.m_data[i];
    return result;
}

// Return size of vector
inline
int GVector::size() const
{
    return m_num;
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
// Vector + Vector
inline
GVector operator+ (const GVector& a, const GVector& b)
{
    GVector result = a;
    result += b;
    return result;
}

// Vector + Scalar
inline
GVector operator+ (const GVector& a, const double& b)
{
    GVector result = a;
    result += b;
    return result;
}

// Scalar + Vector
inline
GVector operator+ (const double& a, const GVector& b)
{
    GVector result = b;
    result += a;
    return result;
}

// Vector - Vector
inline
GVector operator- (const GVector& a, const GVector& b)
{
    GVector result = a;
    result -= b;
    return result;
}

// Vector - Scalar
inline
GVector operator- (const GVector& a, const double& b)
{
    GVector result = a;
    result -= b;
    return result;
}

// Scalar - Vector
inline
GVector operator- (const double& a, const GVector& b)
{
    GVector result = -b;
    result += a;
    return result;
}

// Scalar product
inline
double operator* (const GVector& a, const GVector& b)
{
    if (a.m_num != b.m_num)
        throw GException::vector_mismatch("operator*(GVector, GVector)",
                                          a.m_num, b.m_num);
    double result = 0.0;
    for (int i = 0; i < a.m_num; ++i)
        result += (a.m_data[i] * b.m_data[i]);
    return result;
}

// Vector * double
inline
GVector operator* (const GVector& a, const double& b)
{
    GVector result = a;
    result *= b;
    return result;
}

// double * Vector
inline
GVector operator* (const double& a, const GVector& b)
{
    GVector result = b;
    result *= a;
    return result;
}

// Vector / double
inline
GVector operator/ (const GVector& a, const double& b)
{
    GVector result = a;
    result /= b;
    return result;
}

// Vector norm
inline
double norm(const GVector &v)
{
    double result = 0.0;
    for (int i = 0; i < v.m_num; ++i)
        result += (v.m_data[i] * v.m_data[i]);
    result = (result > 0.0) ? std::sqrt(result) : 0.0;
    return result;
}

// Vector minimum
inline
double min(const GVector &v)
{
    double result = v.m_data[0];
    for (int i = 1; i < v.m_num; ++i) {
        if (v.m_data[i] < result)
            result = v.m_data[i];
  }
  return result;
}

// Vector maximum
inline
double max(const GVector &v)
{
    double result = v.m_data[0];
    for (int i = 1; i < v.m_num; ++i) {
        if (v.m_data[i] > result)
            result = v.m_data[i];
    }
    return result;
}

// Vector sum
inline
double sum(const GVector &v)
{
    double result = 0.0;
    for (int i = 0; i < v.m_num; ++i)
        result += v.m_data[i];
    return result;
}

// Vector permutation
inline
GVector perm(const GVector &v, const int *p)
{
    GVector result(v.m_num);
    if (p == NULL)
        result = v;
    else {
        for (int i = 0; i < v.m_num; ++i)
            result.m_data[i] = v.m_data[p[i]];
    }
    return result;
}

// Inverse vector permutation
inline
GVector iperm(const GVector &v, const int *p)
{
    GVector result(v.m_num);
    if (p == NULL)
        result = v;
    else {
        for (int i = 0; i < v.m_num; ++i)
            result.m_data[p[i]] = v.m_data[i];
    }
    return result;
}

// Vector acos
inline
GVector acos(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::acos(v.m_data[i]);
    return result;
}

// Vector acosh
inline
GVector acosh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = acosh(v.m_data[i]);
    return result;
}

// Vector asin
inline
GVector asin(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::asin(v.m_data[i]);
    return result;
}

// Vector asinh
inline
GVector asinh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = asinh(v.m_data[i]);
    return result;
}

// Vector atan
inline
GVector atan(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::atan(v.m_data[i]);
    return result;
}

// Vector atanh
inline
GVector atanh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = atanh(v.m_data[i]);
    return result;
}

// Vector cos
inline
GVector cos(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::cos(v.m_data[i]);
    return result;
}

// Vector cosh
inline
GVector cosh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::cosh(v.m_data[i]);
    return result;
}

// Vector exp
inline
GVector exp(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::exp(v.m_data[i]);
    return result;
}

// Vector abs
inline
GVector abs(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::abs(v.m_data[i]);
    return result;
}

// Vector log
inline
GVector log(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::log(v.m_data[i]);
    return result;
}

// Vector log10
inline
GVector log10(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::log10(v.m_data[i]);
    return result;
}

// Vector sin
inline
GVector sin(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::sin(v.m_data[i]);
    return result;
}

// Vector sinh
inline
GVector sinh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::sinh(v.m_data[i]);
    return result;
}

// Vector sqrt
inline
GVector sqrt(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::sqrt(v.m_data[i]);
    return result;
}

// Vector tan
inline
GVector tan(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::tan(v.m_data[i]);
    return result;
}

// Vector tanh
inline
GVector tanh(const GVector &v)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::tanh(v.m_data[i]);
    return result;
}

// Vector pow
inline
GVector pow(const GVector &v, const double &a)
{
    GVector result(v.m_num);
    for (int i = 0; i < v.m_num; ++i)
        result.m_data[i] = std::pow(v.m_data[i],a);
    return result;
}

// Equality operator
inline
int operator== (const GVector &a, const GVector &b)
{
    int result = 1;
    if (a.m_num == b.m_num) {
        for (int i = 0; i < a.m_num; ++i) {
            if (a.m_data[i] != b.m_data[i]) {
                result = 0;
                break;
            }
        }
    }
    else
      result = 0;
    return result;
}

// Not equality operator
inline
int operator!= (const GVector &a, const GVector &b)
{
    int result = 0;
    if (a.m_num == b.m_num) {
        for (int i = 0; i < a.m_num; ++i) {
            if (a.m_data[i] != b.m_data[i]) {
                result = 1;
                break;
            }
        }
    }
    else
        result = 1;
    return result;
}

#endif /* GVECTOR_HPP */
