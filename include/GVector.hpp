/***************************************************************************
 *                         GVector.hpp - Vector class                      *
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
 * @file GVector.hpp
 * @brief Vector class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GVECTOR_HPP
#define GVECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <string>
#include "GBase.hpp"
#include "GException.hpp"


/***********************************************************************//**
 * @class GVector
 *
 * @brief Vector class
 *
 * This class implement a double precision floating point vector class that
 * is intended to be used for numerical computation (it is not meant to
 * replace the std::vector template class).
 ***************************************************************************/
class GVector : public GBase {

    // Friend functions
    friend GVector cross(const GVector& a, const GVector& b);
    friend double  operator*(const GVector& a, const GVector& b);
    friend double  norm(const GVector& vector);
    friend double  min(const GVector& vector);
    friend double  max(const GVector& vector);
    friend double  sum(const GVector& vector);
    friend GVector perm(const GVector& vector, const int *p);
    friend GVector iperm(const GVector& vector, const int *p);
    friend GVector acos(const GVector& vector);
    friend GVector acosh(const GVector& vector);
    friend GVector asin(const GVector& vector);
    friend GVector asinh(const GVector& vector);
    friend GVector atan(const GVector& vector);
    friend GVector atanh(const GVector& vector);
    friend GVector cos(const GVector& vector);
    friend GVector cosh(const GVector& vector);
    friend GVector exp(const GVector& vector);
    friend GVector abs(const GVector& vector);
    friend GVector log(const GVector& vector);
    friend GVector log10(const GVector& vector);
    friend GVector sin(const GVector& vector);
    friend GVector sinh(const GVector& vector);
    friend GVector sqrt(const GVector& vector);
    friend GVector tan(const GVector& vector);
    friend GVector tanh(const GVector& vector);
    friend GVector pow(const GVector& vector, const double& power);

public:
    // Constructors and destructors
    GVector(void);
    explicit GVector(const int& num);
    explicit GVector(const double& a);
    explicit GVector(const double& a, const double& b);
    explicit GVector(const double& a, const double& b, const double& c);
    GVector(const GVector& vector);
    virtual ~GVector(void);

    // Vector element access operators
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Vector operators
    bool     operator==(const GVector& vector) const;
    bool     operator!=(const GVector& vector) const;
    GVector& operator=(const GVector& vector);
    GVector& operator+=(const GVector& vector);
    GVector& operator-=(const GVector& vector);
    GVector& operator=(const double& scalar);
    GVector& operator+=(const double& scalar);
    GVector& operator-=(const double& scalar);
    GVector& operator*=(const double& scalar);
    GVector& operator/=(const double& scalar);
    GVector  operator-(void) const;

    // Vector methods
    void          clear(void);
    GVector*      clone(void) const;
    const int&    size(void) const;
    double&       at(const int& index);
    const double& at(const int& index) const;
    int           non_zeros(void) const;
    int           first_nonzero(void) const;
    int           last_nonzero(void) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void alloc_members(void);
    void copy_members(const GVector& vector);
    void free_members(void);

    // Private data area
    int     m_num;    //!< Number of elements in vector
    double* m_data;   //!< Vector array
};


/***********************************************************************//**
 * @brief Vector element access operator
 *
 * @param[in] index Element index [0,...,size()-1]
 * @return Reference to vector element.
 ***************************************************************************/
inline
double& GVector::operator[](const int& index)
{
    // Return vector element
    return m_data[index];
}


/***********************************************************************//**
 * @brief Vector element access operator (const variant)
 *
 * @param[in] index Element index [0,...,size()-1]
 * @return Reference to vector element.
 ***************************************************************************/
inline
const double& GVector::operator[](const int& index) const
{
    // Return vector element
    return m_data[index];
}


/***********************************************************************//**
 * @brief Return size of vector
 *
 * @return Size of vector
 *
 * Returns the number of elements in the vector.
 ***************************************************************************/
inline
const int& GVector::size() const
{
    return m_num;
}


/***********************************************************************//**
 * @brief Add two vectors
 *
 * @param[in] a Vector.
 * @param[in] b Vector.
 * @return Sum of vectors @p a and @p b.
 *
 * Returns the sum of vectors @p a and @p b.
 ***************************************************************************/
inline
GVector operator+(const GVector& a, const GVector& b)
{
    GVector result = a;
    result += b;
    return result;
}


/***********************************************************************//**
 * @brief Add scalar to vector
 *
 * @param[in] vector Vector.
 * @param[in] scalar Scalar.
 * @return Vector with @p scalar added to all elements.
 *
 * Returns a vector for which the @p scalar has been added to all elements.
 ***************************************************************************/
inline
GVector operator+(const GVector& vector, const double& scalar)
{
    GVector result = vector;
    result += scalar;
    return result;
}


/***********************************************************************//**
 * @brief Add scalar to vector
 *
 * @param[in] scalar Scalar.
 * @param[in] vector Vector.
 * @return Vector with @p scalar added to all elements.
 *
 * Returns a vector for which the @p scalar has been added to all elements.
 ***************************************************************************/
inline
GVector operator+(const double& scalar, const GVector& vector)
{
    GVector result = vector;
    result += scalar;
    return result;
}


/***********************************************************************//**
 * @brief Subtract vector from vector
 *
 * @param[in] a Vector.
 * @param[in] b Vector.
 * @return Difference between vector @p a and @p b.
 *
 * Returns the difference between vector @p a and @p b.
 ***************************************************************************/
inline
GVector operator- (const GVector& a, const GVector& b)
{
    GVector result = a;
    result -= b;
    return result;
}


/***********************************************************************//**
 * @brief Subtract scalar from vector
 *
 * @param[in] vector Vector.
 * @param[in] scalar Scalar.
 * @return Vector with @p scalar subtracted from all elements.
 *
 * Returns a vector for which the @p scalar has been subtracted from all
 * elements.
 ***************************************************************************/
inline
GVector operator-(const GVector& vector, const double& scalar)
{
    GVector result = vector;
    result -= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Subtract vector from scalar
 *
 * @param[in] scalar Scalar.
 * @param[in] vector Vector.
 * @return Vector with all elements subtracted from @p scalar.
 *
 * Returns a vector for which all elements have been subtracted from the 
 * @p scalar.
 ***************************************************************************/
inline
GVector operator-(const double& scalar, const GVector& vector)
{
    GVector result = -vector;
    result += scalar;
    return result;
}


/***********************************************************************//**
 * @brief Multiply vector by scalar
 *
 * @param[in] vector Vector.
 * @param[in] scalar Scalar.
 * @return Vector for which all elements have be multiplied by @p scalar.
 *
 * Returns a vector for which all elements have be multiplied by @p scalar.
 ***************************************************************************/
inline
GVector operator*(const GVector& vector, const double& scalar)
{
    GVector result = vector;
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Multiply vector by scalar
 *
 * @param[in] scalar Scalar.
 * @param[in] vector Vector.
 * @return Vector for which all elements have be multiplied by @p scalar.
 *
 * Returns a vector for which all elements have be multiplied by @p scalar.
 ***************************************************************************/
inline
GVector operator*(const double& scalar, const GVector& vector)
{
    GVector result = vector;
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Divide vector by scalar
 *
 * @param[in] vector Vector.
 * @param[in] scalar Scalar.
 * @return Vector for which all elements have be divided by @p scalar.
 *
 * Returns a vector for which all elements have be divided by @p scalar.
 ***************************************************************************/
inline
GVector operator/(const GVector& vector, const double& scalar)
{
    GVector result = vector;
    result /= scalar;
    return result;
}

#endif /* GVECTOR_HPP */
