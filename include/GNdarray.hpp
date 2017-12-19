/***************************************************************************
 *                 GNdarray.hpp - N-dimensional array class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file GNdarray.hpp
 * @brief N-dimensional array class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GNDARRAY_HPP
#define GNDARRAY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GNdarray
 *
 * @brief N-dimensional array class
 *
 * This class implements a n-dimensional double precision floating point
 * array.
 ***************************************************************************/
class GNdarray : public GBase {

    // Friend functions
    friend double   min(const GNdarray& array);
    friend double   max(const GNdarray& array);
    friend double   sum(const GNdarray& array);
    friend GNdarray acos(const GNdarray& array);
    friend GNdarray acosh(const GNdarray& array);
    friend GNdarray asin(const GNdarray& array);
    friend GNdarray asinh(const GNdarray& array);
    friend GNdarray atan(const GNdarray& array);
    friend GNdarray atanh(const GNdarray& array);
    friend GNdarray cos(const GNdarray& array);
    friend GNdarray cosh(const GNdarray& array);
    friend GNdarray exp(const GNdarray& array);
    friend GNdarray abs(const GNdarray& array);
    friend GNdarray log(const GNdarray& array);
    friend GNdarray log10(const GNdarray& array);
    friend GNdarray sign(const GNdarray& array);
    friend GNdarray sin(const GNdarray& array);
    friend GNdarray sinh(const GNdarray& array);
    friend GNdarray sqrt(const GNdarray& array);
    friend GNdarray tan(const GNdarray& array);
    friend GNdarray tanh(const GNdarray& array);
    friend GNdarray pow(const GNdarray& array, const double& power);

public:
    // Constructors and destructors
    GNdarray(void);
    explicit GNdarray(const int& nx);
    GNdarray(const int& nx, const int& ny);
    GNdarray(const int& nx, const int& ny, const int& nz);
    GNdarray(const std::vector<int>& n);
    GNdarray(const GNdarray& array);
    virtual ~GNdarray(void);

    // Element access operators
    double&       operator()(const int& ix);
    double&       operator()(const int& ix, const int& iy);
    double&       operator()(const int& ix, const int& iy, const int& iz);
    double&       operator()(const std::vector<int>& i);
    const double& operator()(const int& ix) const;
    const double& operator()(const int& ix, const int& iy) const;
    const double& operator()(const int& ix, const int& iy, const int& iz) const;
    const double& operator()(const std::vector<int>& i) const;

    // Operators
    GNdarray& operator=(const GNdarray& array);
    bool      operator==(const GNdarray& array) const;
    bool      operator!=(const GNdarray& array) const;
    GNdarray& operator+=(const GNdarray& array);
    GNdarray& operator-=(const GNdarray& array);
    GNdarray& operator*=(const GNdarray& array);
    GNdarray& operator/=(const GNdarray& array);
    GNdarray& operator+=(const double& value);
    GNdarray& operator-=(const double& value);
    GNdarray& operator*=(const double& value);
    GNdarray& operator/=(const double& value);
    GNdarray  operator-(void) const;

    // Methods
    void                    clear(void);
    GNdarray*               clone(void) const;
    std::string             classname(void) const;
    int                     dim(void) const;
    int                     size(void) const;
    const std::vector<int>& shape(void) const;
    const std::vector<int>& strides(void) const;
    void                    shape(const std::vector<int>& shape);
    double&                 at(const int& ix);
    double&                 at(const int& ix, const int& iy);
    double&                 at(const int& ix, const int& iy, const int& iz);
    double&                 at(const std::vector<int>& i);
    const double&           at(const int& ix) const;
    const double&           at(const int& ix, const int& iy) const;
    const double&           at(const int& ix, const int& iy, const int& iz) const;
    const double&           at(const std::vector<int>& i) const;
    const double*           data(void) const;
    double*                 data(void);
    std::string             print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GNdarray& array);
    void free_members(void);
    int  index(const std::vector<int>& i) const;
    bool has_same_shape(const GNdarray& array) const;
    void require_same_shape(const std::string& method, const GNdarray& array) const;

    // Protected members
    std::vector<int>    m_shape;   //!< Array dimensions
    std::vector<int>    m_strides; //!< Steps in each dimension when traversing array
    std::vector<double> m_data;    //!< Array data
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GNdarray").
 ***************************************************************************/
inline
std::string GNdarray::classname(void) const
{
    return ("GNdarray");
}


/***********************************************************************//**
 * @brief 1-dimensional array element access operator
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::operator()(const int& ix)
{
    // Return array element
    return m_data[ix];
}


/***********************************************************************//**
 * @brief 2-dimensional array element access operator
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::operator()(const int& ix, const int& iy)
{
    // Return array element
    return m_data[ix+m_strides[1]*iy];
}


/***********************************************************************//**
 * @brief 3-dimensional array element access operator
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @param[in] iz Index in third dimension [0,...,shape(2)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::operator()(const int& ix, const int& iy, const int& iz)
{
    // Return array element
    return m_data[ix+m_strides[1]*iy+m_strides[2]*iz];
}


/***********************************************************************//**
 * @brief n-dimensional array element access operator
 *
 * @param[in] i Index vector.
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::operator()(const std::vector<int>& i)
{
    // Return array element
    return m_data[index(i)];
}


/***********************************************************************//**
 * @brief 1-dimensional array element access operator (const variant)
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Const reference to array element.
 ***************************************************************************/
inline
const double& GNdarray::operator()(const int& ix) const
{
    // Return array element
    return m_data[ix];
}


/***********************************************************************//**
 * @brief 2-dimensional array element access operator (const variant)
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Const reference to array element.
 ***************************************************************************/
inline
const double& GNdarray::operator()(const int& ix, const int& iy) const
{
    // Return array element
    return m_data[ix+m_strides[1]*iy];
}


/***********************************************************************//**
 * @brief 3-dimensional array element access operator (const variant)
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @param[in] iz Index in third dimension [0,...,shape(2)-1].
 * @return Const reference to array element.
 ***************************************************************************/
inline
const double& GNdarray::operator()(const int& ix, const int& iy, const int& iz) const
{
    // Return array element
    return m_data[ix+m_strides[1]*iy+m_strides[2]*iz];
}


/***********************************************************************//**
 * @brief n-dimensional array element access operator (const variant)
 *
 * @param[in] i Index vector.
 * @return Const reference to array element.
 ***************************************************************************/
inline
const double& GNdarray::operator()(const std::vector<int>& i) const
{
    // Return array element
    return m_data[index(i)];
}


/***********************************************************************//**
 * @brief Return dimension of array
 *
 * @return Dimension of array
 *
 * Returns the dimension of the array.
 ***************************************************************************/
inline
int GNdarray::dim(void) const
{
    return int(m_shape.size());
}


/***********************************************************************//**
 * @brief Return number of elements in array
 *
 * @return Number of elements in array
 *
 * Returns the number of elements in the array.
 ***************************************************************************/
inline
int GNdarray::size(void) const
{
    return int(m_data.size());
}


/***********************************************************************//**
 * @brief Return shape of array
 *
 * @return Shape of array
 *
 * Returns the shape of the array.
 ***************************************************************************/
inline
const std::vector<int>& GNdarray::shape(void) const
{
    return m_shape;
}


/***********************************************************************//**
 * @brief Return strides of array
 *
 * @return Strides of array
 *
 * Returns the strides of the array.
 ***************************************************************************/
inline
const std::vector<int>& GNdarray::strides(void) const
{
    return m_strides;
}


/***********************************************************************//**
 * @brief 1-dimensional array element access with range checking
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::at(const int& ix)
{
    return const_cast<double &>(static_cast<const GNdarray &>(*this).at(ix));
}


/***********************************************************************//**
 * @brief 2-dimensional array element access with range checking
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::at(const int& ix, const int& iy)
{
    return const_cast<double &>(static_cast<const GNdarray &>(*this).at(ix,iy));
}


/***********************************************************************//**
 * @brief 3-dimensional array element access with range checking
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @param[in] iz Index in third dimension [0,...,shape(2)-1].
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::at(const int& ix, const int& iy, const int& iz)
{
    return const_cast<double &>(static_cast<const GNdarray &>(*this).at(ix,iy,iz));
}


/***********************************************************************//**
 * @brief n-dimensional array element access with range checking
 *
 * @param[in] i Index vector.
 * @return Reference to array element.
 ***************************************************************************/
inline
double& GNdarray::at(const std::vector<int>& i)
{
    return const_cast<double &>(static_cast<const GNdarray &>(*this).at(i));
}


/***********************************************************************//**
 * @brief Data access method (const version)
 *
 * @return Const reference to array data.
 ***************************************************************************/
inline
const double* GNdarray::data(void) const
{
    return (&(m_data[0]));
}


/***********************************************************************//**
 * @brief Data access method
 *
 * @return Reference to array data.
 ***************************************************************************/
inline
double* GNdarray::data(void)
{
    return (&(m_data[0]));
}


/***********************************************************************//**
 * @brief Return sum of two arrays
 *
 * @param[in] a First array.
 * @param[in] b Second array.
 * @return Sum of arrays @p a and @p b.
 *
 * Returns the sum of arrays @p a and @p b.
 ***************************************************************************/
inline
GNdarray operator+(const GNdarray& a, const GNdarray& b)
{
    GNdarray result = a;
    result += b;
    return result;
}


/***********************************************************************//**
 * @brief Add value to array (right addition)
 *
 * @param[in] array Array.
 * @param[in] value Value.
 * @return Array with @p value added to all elements.
 *
 * Returns an array for which the @p value has been added to all elements.
 ***************************************************************************/
inline
GNdarray operator+(const GNdarray& array, const double& value)
{
    GNdarray result = array;
    result += value;
    return result;
}


/***********************************************************************//**
 * @brief Add value to array (left addition)
 *
 * @param[in] value Value.
 * @param[in] array Array.
 * @return Array with @p value added to all elements.
 *
 * Returns an array for which the @p value has been added to all elements.
 ***************************************************************************/
inline
GNdarray operator+(const double& value, const GNdarray& array)
{
    GNdarray result = array;
    result += value;
    return result;
}


/***********************************************************************//**
 * @brief Return difference of arrays
 *
 * @param[in] a First array.
 * @param[in] b Second array.
 * @return Difference between array @p a and @p b.
 *
 * Returns the difference between array @p a and @p b.
 ***************************************************************************/
inline
GNdarray operator-(const GNdarray& a, const GNdarray& b)
{
    GNdarray result = a;
    result -= b;
    return result;
}


/***********************************************************************//**
 * @brief Subtract value from array
 *
 * @param[in] array Array.
 * @param[in] value Value.
 * @return Array with @p value subtracted from all elements.
 *
 * Returns an array for which the @p value has been subtracted from all
 * elements. For example
 *
 *     double   value  = 5.0;
 *     GNdarray result = array - value;
 ***************************************************************************/
inline
GNdarray operator-(const GNdarray& array, const double& value)
{
    GNdarray result = array;
    result -= value;
    return result;
}


/***********************************************************************//**
 * @brief Subtract array from value
 *
 * @param[in] value Value.
 * @param[in] array Array.
 * @return Array with @p value subtracted from all elements.
 *
 * Returns an array for which all elements have been subtracted from
 * @p value. For example
 *
 *     double   value  = 5.0;
 *     GNdarray result = value - array;
 ***************************************************************************/
inline
GNdarray operator-(const double& value, const GNdarray& array)
{
    GNdarray result = -array;
    result += value;
    return result;
}


/***********************************************************************//**
 * @brief Multiply array by value (right multiplication)
 *
 * @param[in] array Array.
 * @param[in] value Value.
 * @return Array for which all elements have be multiplied by @p value.
 *
 * Returns an array for which all elements have be multiplied by @p value.
 ***************************************************************************/
inline
GNdarray operator*(const GNdarray& array, const double& value)
{
    GNdarray result = array;
    result *= value;
    return result;
}


/***********************************************************************//**
 * @brief Multiply array by value (left multiplication)
 *
 * @param[in] value Value.
 * @param[in] array Array.
 * @return Array for which all elements have be multiplied by @p value.
 *
 * Returns an array for which all elements have be multiplied by @p value.
 ***************************************************************************/
inline
GNdarray operator*(const double& value, const GNdarray& array)
{
    GNdarray result = array;
    result *= value;
    return result;
}


/***********************************************************************//**
 * @brief Divide array by value
 *
 * @param[in] array Array.
 * @param[in] value Value.
 * @return Array for which all elements have be divided by @p value.
 *
 * Returns an array for which all elements have be divided by @p value.
 ***************************************************************************/
inline
GNdarray operator/(const GNdarray& array, const double& value)
{
    GNdarray result = array;
    result /= value;
    return result;
}

#endif /* GNDARRAY_HPP */
