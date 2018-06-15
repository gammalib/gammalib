/***************************************************************************
 *                GFft.hpp - Fast Fourier transformation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2018 by Juergen Knoedlseder                         *
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
 * @file GFft.hpp
 * @brief Fast Fourier transformation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFFT_HPP
#define GFFT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <complex>
#include "GBase.hpp"
#include "GFftWavetable.hpp"

/* __ Forward declarations _______________________________________________ */
class GNdarray;


/***********************************************************************//**
 * @class GFft
 *
 * @brief Fast Fourier Transformation class
 *
 * This class implements a Fast Fourier Transformation of an n-dimensional
 * double precision floating point array.
 *
 * The class implementation is based on the Fast Fourier Transformation
 * functions that are provided by the GNU Scientific Library (GSL)
 * (version 2.2.1). The mixed-radix routines for complex numbers have been
 * adopted that work for any array lengths. The mixed-radix routines are a
 * reimplementation of the FFTPACK library of Paul Swarztrauber.
 *
 * The GSL is documented at https://www.gnu.org/software/gsl/manual
 ***************************************************************************/
class GFft : public GBase {

public:
    // Constructors and destructors
    GFft(void);
    explicit GFft(const GNdarray& array);
    GFft(const GFft& fft);
    virtual ~GFft(void);

    // FFT access operators
    std::complex<double>&       operator()(const int& ix);
    std::complex<double>&       operator()(const int& ix, const int& iy);
    const std::complex<double>& operator()(const int& ix) const;
    const std::complex<double>& operator()(const int& ix, const int& iy) const;

    // Operators
    GFft& operator=(const GFft& fft);
    GFft& operator+=(const GFft& fft);
    GFft& operator-=(const GFft& fft);
    GFft& operator*=(const GFft& fft);
    GFft& operator/=(const GFft& fft);
    GFft  operator-(void) const;

    // Methods
    void                    clear(void);
    GFft*                   clone(void) const;
    std::string             classname(void) const;
    int                     dim(void) const;
    int                     size(void) const;
    const std::vector<int>& shape(void) const;
    const std::vector<int>& strides(void) const;
    void                    forward(const GNdarray& array);
    GNdarray                backward(void) const;
    std::string             print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFft& fft);
    void free_members(void);
    void set_data(const GNdarray& array);
    bool has_same_shape(const GFft& fft) const;
    void require_same_shape(const std::string& method, const GFft& fft) const;

    // Low-level FFT methods
    void transform(std::complex<double>* data,
                   const int&            stride,
                   const int&            n,
                   const GFftWavetable&  wavetable,
                   const bool&           forward = true);
    void factor2(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factor3(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factor4(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factor5(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factor6(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factor7(const std::complex<double>* in,
                 const int&                  istride,
                 std::complex<double>*       out,
                 const int&                  ostride,
                 const GFftWavetable&        wavetable,
                 const int&                  sign,
                 const int&                  product,
                 const int&                  n,
                 const int&                  index);
    void factorn(std::complex<double>* in,
                 const int&            istride,
                 std::complex<double>* out,
                 const int&            ostride,
                 const GFftWavetable&  wavetable,
                 const int&            sign,
                 const int&            factor,
                 const int&            product,
                 const int&            n,
                 const int&            index);
    std::vector<std::complex<double> > get_w(const GFftWavetable& wavetable,
                                             const int&           index,
                                             const int&           k,
                                             const int&           q,
                                             const int&           n,
                                             const int&           sign) const;
    std::complex<double> timesi(const std::complex<double>& value) const;

    // Protected members
    int                        m_size;      //!< Size of data array
    std::complex<double>*      m_data;      //!< Pointer on array data
    std::vector<int>           m_shape;     //!< Array dimensions
    std::vector<int>           m_strides;   //!< Steps in each dimension
    std::vector<GFftWavetable> m_wavetable; //!< Trigonometric coefficients
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GFft").
 ***************************************************************************/
inline
std::string GFft::classname(void) const
{
    return ("GFft");
}


/***********************************************************************//**
 * @brief 1-dimensional FFT element access operator
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Reference to FFT element.
 ***************************************************************************/
inline
std::complex<double>& GFft::operator()(const int& ix)
{
	// Return array element
    return *(m_data + ix);
}


/***********************************************************************//**
 * @brief 2-dimensional FFT element access operator
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Reference to FFT element.
 ***************************************************************************/
inline
std::complex<double>& GFft::operator()(const int& ix, const int& iy)
{
	// Return array element
    return *(m_data + ix + m_strides[1]*iy);
}


/***********************************************************************//**
 * @brief 1-dimensional FFT element access operator (const variant)
 *
 * @param[in] ix Element index [0,...,shape(0)-1].
 * @return Const reference to FFT element.
 ***************************************************************************/
inline
const std::complex<double>& GFft::operator()(const int& ix) const
{
    // Return array element
    return *(m_data + ix);
}


/***********************************************************************//**
 * @brief 2-dimensional FFT element access operator (const variant)
 *
 * @param[in] ix Index in first dimension [0,...,shape(0)-1].
 * @param[in] iy Index in second dimension [0,...,shape(1)-1].
 * @return Const reference to FFT element.
 ***************************************************************************/
inline
const std::complex<double>& GFft::operator()(const int& ix, const int& iy) const
{
	// Return array element
    return *(m_data + ix + m_strides[1]*iy);
}




/***********************************************************************//**
 * @brief Return dimension of Fast Fourier Transformation
 *
 * @return Dimension of Fast Fourier Transformation.
 *
 * Returns the dimension of the Fast Fourier Transformation.
 ***************************************************************************/
inline
int GFft::dim(void) const
{
    return ((int)m_shape.size());
}


/***********************************************************************//**
 * @brief Return number of elements in array
 *
 * @return Number of elements in array
 *
 * Returns the number of elements in the array.
 ***************************************************************************/
inline
int GFft::size(void) const
{
    return m_size;
}


/***********************************************************************//**
 * @brief Return shape of array
 *
 * @return Shape of array
 *
 * Returns the shape of the array.
 ***************************************************************************/
inline
const std::vector<int>& GFft::shape(void) const
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
const std::vector<int>& GFft::strides(void) const
{
    return m_strides;
}


/***********************************************************************//**
 * @brief Return sum of two Fast Fourier Transformations
 *
 * @param[in] a First Fast Fourier Transformation.
 * @param[in] b Second Fast Fourier Transformation.
 * @return Sum of the Fast Fourier Transformations @p a and @p b.
 *
 * Returns the sum of the Fast Fourier Transformations @p a and @p b.
 ***************************************************************************/
inline
GFft operator+(const GFft& a, const GFft& b)
{
    GFft result = a;
    result += b;
    return result;
}


/***********************************************************************//**
 * @brief Return difference of two Fast Fourier Transformations
 *
 * @param[in] a First Fast Fourier Transformation.
 * @param[in] b Second Fast Fourier Transformation.
 * @return Difference of the Fast Fourier Transformations @p a and @p b.
 *
 * Returns the difference of the Fast Fourier Transformations @p a and @p b.
 ***************************************************************************/
inline
GFft operator-(const GFft& a, const GFft& b)
{
    GFft result = a;
    result -= b;
    return result;
}


/***********************************************************************//**
 * @brief Return product of two Fast Fourier Transformations
 *
 * @param[in] a First Fast Fourier Transformation.
 * @param[in] b Second Fast Fourier Transformation.
 * @return Product of the Fast Fourier Transformations @p a and @p b.
 *
 * Returns the product of the Fast Fourier Transformations @p a and @p b.
 ***************************************************************************/
inline
GFft operator*(const GFft& a, const GFft& b)
{
    GFft result = a;
    result *= b;
    return result;
}


/***********************************************************************//**
 * @brief Return quotient of two Fast Fourier Transformations
 *
 * @param[in] a First Fast Fourier Transformation.
 * @param[in] b Second Fast Fourier Transformation.
 * @return Quotient of the Fast Fourier Transformations @p a and @p b.
 *
 * Returns the quotient of the Fast Fourier Transformations @p a and @p b.
 ***************************************************************************/
inline
GFft operator/(const GFft& a, const GFft& b)
{
    GFft result = a;
    result /= b;
    return result;
}


/***********************************************************************//**
 * @brief Return complex value times i
 *
 * @param[in] value Complex value.
 * @return Complex value times i.
 *
 * Returns complex value times i.
 ***************************************************************************/
inline
std::complex<double> GFft::timesi(const std::complex<double>& value) const
{
    std::complex<double> result(-value.imag(), value.real());
    return result;
}

#endif /* GFFT_HPP */
