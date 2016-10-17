/***************************************************************************
 *                GFft.hpp - Fast Fourier transformation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GFft.cpp
 * @brief Fast Fourier transformation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GFft.hpp"
#include "GNdarray.hpp"

/* __ Method name definitions ____________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFft::GFft(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief N-dimensional array constructor
 *
 * @param[in] array N-dimensional array.
 *
 * Constructs a Fast Fourier Transformation from a n-dimensional array.
 ***************************************************************************/
GFft::GFft(const GNdarray& array)
{
    // Initialise class members
    init_members();

    // Perform foward transformation
    forward(array);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] fft Fast Fourier Transform.
 ***************************************************************************/
GFft::GFft(const GFft& fft)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(fft);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFft::~GFft(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] fft Fast Fourier Transform.
 * @return Fast Fourier Transform.
 ***************************************************************************/
GFft& GFft::operator=(const GFft& fft)
{
    // Execute only if object is not identical
    if (this != &fft) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(fft);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Unary addition operator
 *
 * @param[in] fft Fast Fourier Transform.
 * @return Fast Fourier Transform.
 *
 * Adds a Fast Fourier Transform to another Fast Fourier Transform.
 *
 * @todo Implement operator
 ***************************************************************************/
GFft& GFft::operator+=(const GFft& fft)
{
    // TODO: Implement operator

    // Return Fast Fourier Transform
    return *this;
}


/***********************************************************************//**
 * @brief Unary subtraction operator
 *
 * @param[in] fft Fast Fourier transform.
 * @return Fast Fourier transform.
 *
 * Subtracts a Fast Fourier Transform from another Fast Fourier Transform.
 *
 * @todo Implement operator
 ***************************************************************************/
GFft& GFft::operator-=(const GFft& fft)
{
    // TODO: Implement operator

    // Return Fast Fourier Transform
    return *this;
}


/***********************************************************************//**
 * @brief Unary multiplication operator
 *
 * @param[in] fft Fast Fourier transform.
 * @return Fast Fourier transform.
 *
 * Multiplies a Fast Fourier Transform to another Fast Fourier Transform.
 *
 * @todo Implement operator
 ***************************************************************************/
GFft& GFft::operator*=(const GFft& fft)
{
    // TODO: Implement operator

    // Return Fast Fourier Transform
    return *this;
}


/***********************************************************************//**
 * @brief Unary division operator
 *
 * @param[in] fft Fast Fourier transform.
 * @return Fast Fourier transform.
 *
 * Divides a Fast Fourier Transform by another Fast Fourier Transform.
 *
 * @todo Implement operator
 ***************************************************************************/
GFft& GFft::operator/=(const GFft& fft)
{
    // TODO: Implement operator

    // Return Fast Fourier Transform
    return *this;
}


/***********************************************************************//**
 * @brief Unary minus operator
 *
 * @return Fast Fourier transformation.
 *
 * Negates all elements of the Fast Fourier Transformation.
 ***************************************************************************/
GFft GFft::operator-(void) const
{
    // Copy FFT
    GFft fft = *this;

    // Negate all elements
    for (int i = 0; i < fft.size(); ++i) {
        fft(i) = -fft(i);
    }

    // Return FFT
    return fft;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Fast Fourier Transform
 ***************************************************************************/
void GFft::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone Fast Fourier Transform
 *
 * @return Pointer to deep copy of Fast Fourier Transform.
 ***************************************************************************/
GFft* GFft::clone(void) const
{
    // Clone FFT
    return new GFft(*this);
}


/***********************************************************************//**
 * @brief Forward Fast Fourier Transform
 *
 * @param[in] array N-dimensional array.
 *
 * @todo Method should support N-dim arrays, only supports 1-dim so far
 ***************************************************************************/
void GFft::forward(const GNdarray& array)
{
    // Set data
    set_data(array);

    // Continue only if there are elements in the data member
    if (m_size > 0) {
/*
        // Loop over all dimensions of array
        for (int d = 0; d < dim(); ++d) {

            // Loop over all dimensions of array
            for (int i = 0; i < dim(); ++i) {

                // Skip
                if (d == i) {
                    continue;
                }

                // Loop over all
                int n      = shape()[d];
                int stride = strides()[d];
                int istart =
*/

        // Perform FFT (is only correct for 1D for now)
        transform(m_data, 1, m_shape[0], m_wavetable[0], true);

    } // endif: there were elements in the FFT data member

    // Return
    return;
}


/***********************************************************************//**
 * @brief Backward Fast Fourier Transform
 *
 * @return N-dimensional array.
 *
 * @todo Method should support N-dim arrays, only supports 1-dim so far
 ***************************************************************************/
GNdarray GFft::backward(void) const
{
    // Initialise n-dimensional array
    GNdarray array(shape());

    // Continue only if there are elements in the array
    if (array.size() > 0) {

        // Create copy of complex array
        std::complex<double>* data = new std::complex<double>[array.size()];
        for (int i = 0; i < array.size(); ++i) {
            *(data+i) = *(m_data+i);
        }

        // Perform FFT (circumvent const correctness)
        const_cast<GFft*>(this)->transform(data, 1, m_shape[0], m_wavetable[0],
                                           false);

        // Get normalisation factor
        const double norm  = 1.0 / double(array.size());

        // Extract N-dimensional array and normalise with 1/n
        for (int i = 0; i < array.size(); ++i) {
            array(i) = (data+i)->real() * norm;
        }

        // Delete copy of complex array
        delete [] data;

    } // endif: there were elements in the array

    // Return N-dimensional array
    return array;
}


/***********************************************************************//**
 * @brief Print Fast Fourier Transform information
 *
 * @param[in] chatter Chattiness.
 * @return String containing Fast Fourier Transform information.
 ***************************************************************************/
std::string GFft::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFft ===");

        // Append array information
        result.append("\n"+gammalib::parformat("Dimension"));
        result.append(gammalib::str(dim()));
        result.append("\n"+gammalib::parformat("Shape"));
        result.append("(");
        for (int i = 0; i < dim(); ++i) {
            if (i > 0) {
                result.append(",");
            }
            result.append(gammalib::str(m_shape[i]));
        }
        result.append(")");
        result.append("\n"+gammalib::parformat("Size"));
        result.append(gammalib::str(m_size));

        // VERBOSE: Put all FFT elements in string
        if (chatter >= VERBOSE) {
            result.append("\n"+gammalib::parformat("Elements"));
            for (int i = 0; i < size(); ++i) {
                if (i > 0) {
                    result.append(",");
                }
                result.append(gammalib::str((*this)(i)));
            }
            result.append(")");
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFft::init_members(void)
{
    // Initialise members
    m_size = 0;
    m_data = NULL;
    m_shape.clear();
    m_strides.clear();
    m_wavetable.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fft Fast Fourier Transform.
 ***************************************************************************/
void GFft::copy_members(const GFft& fft)
{
    // Copy members
    m_size      = fft.m_size;
    m_shape     = fft.m_shape;
    m_strides   = fft.m_strides;
    m_wavetable = fft.m_wavetable;

    // Copy data
    if (m_size > 0) {
        m_data = new std::complex<double>[m_size];
        for (int i = 0; i < m_size; ++i) {
            *(m_data+i) = *(fft.m_data+i);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFft::free_members(void)
{
    // Delete data
    if (m_data != NULL) {
        delete [] m_data;
    }

    // Reset size and pointer
    m_size = 0;
    m_data = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set data from n-dimensional array
 *
 * @param[in] array N-dimensional array
 ***************************************************************************/
void GFft::set_data(const GNdarray& array)
{
    // Clear information
    free_members();
    init_members();

    // Continue only if there are elements in the array
    if (array.size() > 0) {

        // Set size of array
        m_size = array.size();

        // Allocate data and initialise data to zero
        m_data = new std::complex<double>[m_size];

        // Set real part of data from n-dimensional array
        for (int i = 0; i < m_size; ++i) {
            *(m_data+i) = std::complex<double>(array(i), 0.0);
        }

        // Save shape
        m_shape = array.shape();

        // Set strides
        int stride = 1;
        for (int i = 0; i < m_shape.size(); ++i) {
            m_strides.push_back(stride);
            stride *= m_shape[i];
        }

        // Set trigonometric coefficients
        for (int i = 0; i < m_shape.size(); ++i) {
            GFftWavetable wavetable(m_shape[i]);
            m_wavetable.push_back(wavetable);
        }

    } // endif: there were elements in the array

    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform Fast Fourier Transform
 *
 * @param[in] data Pointer to complex array to be transformed.
 * @param[in] stride Step size when traversing complex array.
 * @param[in] n Number of elements in complex array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] forward Forward transform if true, otherwise backward transform.
 ***************************************************************************/
void GFft::transform(std::complex<double>* data,
                     const int&            stride,
                     const int&            n,
                     const GFftWavetable&  wavetable,
                     const bool&           forward)
{
    // Allocate scratch space and initialize to zero
    std::complex<double>* scratch = new std::complex<double>[n];
    std::complex<double>  zero(0.0, 0.0);
    for (int i = 0; i < n; ++i) {
        *(scratch+i) = zero;
    }

    // Initialise transformation state
    int state = 0;

    // Initialise factorisation product
    int product = 1;

    // Set sign
    int sign = (forward) ? -1 : +1;

    // Loop over all factors
    for (int i = 0; i < wavetable.factors(); ++i) {

        // Get factorisation factor and start index
        int factor = wavetable.factor(i);
        int index  = wavetable.index(i);

        // Update factorisation product
        product *= factor;

        // Compute number of coefficients per factor
        int q = n / product;

        // Set state dependent stuff
        std::complex<double>* in;
        std::complex<double>* out;
        int                   istride;
        int                   ostride;
        if (state == 0) {
            in      = data;
            istride = stride;
            out     = scratch;
            ostride = 1;
            state   = 1;
        }
        else {
            in      = scratch;
            istride = 1;
            out     = data;
            ostride = stride;
            state   = 0;
        }

        // Call factor dependent method
        if (factor == 2) {
            factor2(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 3) {
            factor3(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 4) {
            factor4(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 5) {
            factor5(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 6) {
            factor6(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 7) {
            factor7(in, istride, out, ostride, wavetable, sign, product, n, index);
        }
        else {
            factorn(in, istride, out, ostride, wavetable, sign, factor, product, n, index);
        }

    } // endfor: looped over all factors

    // In case the loop was exited in state 1 the results are in the scratch
    // array and we need to copy the results back from the scratch array to
    // the data array
    if (state == 1) {
        for (int i = 0; i < n; ++i) {
            *(data+stride*i) = *(scratch+i);
        }
    }

    // Free scratch space
    delete [] scratch;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 2
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor2(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 2;
    const int m      = n / factor;
    const int q      = n / product;
    const int p_1    = product / factor;
    const int jump   = (factor - 1) * p_1;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 1, sign);

        // Compute x = W(2) z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));

            // Compute x
            const std::complex<double> x0 = z0 + z1;
            const std::complex<double> x1 = z0 - z1;

            // out = w * x
            *(out+ostride*j)       = x0;
            *(out+ostride*(j+p_1)) = w[0] * x1;

        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 3
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor3(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 3;
    const int m      = n / factor;
    const int q      = n / product;
    const int p_1    = product / factor;
    const int jump   = (factor - 1) * p_1;


    // Precompute some factors
    const double tau = std::sqrt(3.0) / 2.0;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 2, sign);

        // Compute x = W(3) z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));
            const std::complex<double> z2 = *(in+istride*(i+2*m));

            // Compute t
            const std::complex<double> t1 = z1 + z2;
            const std::complex<double> t2 = z0 - t1/2.0;
            const std::complex<double> t3 = double((int)sign) * tau * (z1 - z2);

            // Compute x
            const std::complex<double> x0 = z0 + t1;
            const std::complex<double> x1 = t2 + timesi(t3);
            const std::complex<double> x2 = t2 - timesi(t3);

            // out = w * x
            *(out+ostride*j)         = x0;
            *(out+ostride*(j+p_1))   = w[0] * x1;
            *(out+ostride*(j+2*p_1)) = w[1] * x2;

        } // endfor: k1
        
    } // endfor: k

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 4
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor4(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 4;
    const int m = n / factor;
    const int q = n / product;
    const int p_1 = product / factor;
    const int jump = (factor - 1) * p_1;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 3, sign);

        // Compute x = W(4) z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));
            const std::complex<double> z2 = *(in+istride*(i+2*m));
            const std::complex<double> z3 = *(in+istride*(i+3*m));
            
            // Compute t
            const std::complex<double> t1 = z0 + z2;
            const std::complex<double> t2 = z1 + z3;
            const std::complex<double> t3 = z0 - z2;
            const std::complex<double> t4 = double((int)sign) * (z1 - z3);

            // Compute x
            const std::complex<double> x0 = t1 + t2;
            const std::complex<double> x1 = t3 + timesi(t4);
            const std::complex<double> x2 = t1 - t2;
            const std::complex<double> x3 = t3 - timesi(t4);

            // out = w * x
            *(out+ostride*j)         = x0;
            *(out+ostride*(j+p_1))   = w[0] * x1;
            *(out+ostride*(j+2*p_1)) = w[1] * x2;
            *(out+ostride*(j+3*p_1)) = w[2] * x3;
          
        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 5
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor5(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 5;
    const int m      = n / factor;
    const int q      = n / product;
    const int p_1    = product / factor;
    const int jump   = (factor - 1) * p_1;

    // Precompute some factors
    const double sin_2pi_by_5  = std::sin(gammalib::twopi / 5.0);
    const double sin_2pi_by_10 = std::sin(gammalib::twopi / 10.0);
    const double tau           = std::sqrt(5.0) / 4.0;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 4, sign);

        // Compute x = W(5) z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));
            const std::complex<double> z2 = *(in+istride*(i+2*m));
            const std::complex<double> z3 = *(in+istride*(i+3*m));
            const std::complex<double> z4 = *(in+istride*(i+4*m));

            // Compute t
            const std::complex<double> t1 = z1 + z4;
            const std::complex<double> t2 = z2 + z3;
            const std::complex<double> t3 = z1 - z4;
            const std::complex<double> t4 = z2 - z3;
            const std::complex<double> t5 = t1 + t2;
            const std::complex<double> t6 = tau * (t1 - t2);
            const std::complex<double> t7 = z0 - t5/4.0;
            const std::complex<double> t8 = t7 + t6;
            const std::complex<double> t9 = t7 - t6;
            const std::complex<double> t10 = double((int)sign) *
                                             (sin_2pi_by_5*t3 + sin_2pi_by_10*t4);
            const std::complex<double> t11 = double((int)sign) *
                                             (sin_2pi_by_10*t3 - sin_2pi_by_5*t4);
          
            // Compute x
            const std::complex<double> x0 = z0 + t5;
            const std::complex<double> x1 = t8 + timesi(t10);
            const std::complex<double> x2 = t9 + timesi(t11);
            const std::complex<double> x3 = t9 - timesi(t11);
            const std::complex<double> x4 = t8 - timesi(t10);

            // Compute out = w * x
            *(out+ostride*j)         = x0;
            *(out+ostride*(j+p_1))   = w[0] * x1;
            *(out+ostride*(j+2*p_1)) = w[1] * x2;
            *(out+ostride*(j+3*p_1)) = w[2] * x3;
            *(out+ostride*(j+4*p_1)) = w[3] * x4;
          
        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 6
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor6(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 6;
    const int m = n / factor;
    const int q = n / product;
    const int p_1 = product / factor;
    const int jump = (factor - 1) * p_1;

    // Precompute some factors
    const double tau = std::sqrt(3.0) / 2.0;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 5, sign);

        // Compute x = W(6) z. W(6) is a combination of sums and differences of
        // W(3) acting on the even and odd elements of z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));
            const std::complex<double> z2 = *(in+istride*(i+2*m));
            const std::complex<double> z3 = *(in+istride*(i+3*m));
            const std::complex<double> z4 = *(in+istride*(i+4*m));
            const std::complex<double> z5 = *(in+istride*(i+5*m));

            // Compute ta
            const std::complex<double> ta1 = z2 + z4;
            const std::complex<double> ta2 = z0 - ta1/2.0;
            const std::complex<double> ta3 = double((int)sign) * tau * (z2 - z4);
          
            // Compute a
            const std::complex<double> a0 = z0 + ta1;
            const std::complex<double> a1 = ta2 + timesi(ta3);
            const std::complex<double> a2 = ta2 - timesi(ta3);
          
            // Compute tb
            const std::complex<double> tb1 = z5 + z1;
            const std::complex<double> tb2 = z3 - tb1/2.0;
            const std::complex<double> tb3 = double((int)sign) * tau * (z5 - z1);
          
            // Compute b
            const std::complex<double> b0 = z3 + tb1;
            const std::complex<double> b1 = tb2 + timesi(tb3);
            const std::complex<double> b2 = tb2 - timesi(tb3);
          
            // Compute x
            const std::complex<double> x0 = a0 + b0;
            const std::complex<double> x1 = a1 - b1;
            const std::complex<double> x2 = a2 + b2;
            const std::complex<double> x3 = a0 - b0;
            const std::complex<double> x4 = a1 + b1;
            const std::complex<double> x5 = a2 - b2;
            
            // Compute out = w * x
            *(out+ostride*j)         = x0;
            *(out+ostride*(j+p_1))   = w[0] * x1;
            *(out+ostride*(j+2*p_1)) = w[1] * x2;
            *(out+ostride*(j+3*p_1)) = w[2] * x3;
            *(out+ostride*(j+4*p_1)) = w[3] * x4;
            *(out+ostride*(j+5*p_1)) = w[4] * x5;

        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 7
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor7(const std::complex<double>* in,
                   const int&                  istride,
                   std::complex<double>*       out,
                   const int&                  ostride,
                   const GFftWavetable&        wavetable,
                   const int&                  sign,
                   const int&                  product,
                   const int&                  n,
                   const int&                  index)
{
    // Compute ...
    const int factor = 7;
    const int m      = n / factor;
    const int q      = n / product;
    const int p_1    = product / factor;
    const int jump   = (factor - 1) * p_1;

    // Precompute some factors
    static const double twopi7 = gammalib::twopi / 7.0;
    static const double c1     = std::cos(1.0 * twopi7);
    static const double c2     = std::cos(2.0 * twopi7);
    static const double c3     = std::cos(3.0 * twopi7);
    static const double s1     = std::sin(1.0 * twopi7);
    static const double s2     = std::sin(2.0 * twopi7);
    static const double s3     = std::sin(3.0 * twopi7);
    static const double tau1   = (c1 + c2 + c3) / 3.0 - 1.0;
    static const double tau2   = (2.0 * c1 - c2 - c3) / 3.0;
    static const double tau3   = (c1 - 2.0*c2 + c3) / 3.0;
    static const double tau4   = (c1 + c2 - 2.0 * c3) / 3.0;
    static const double tau5   = (s1 + s2 - s3) / 3.0;
    static const double tau6   = (2.0 * s1 - s2 + s3) / 3.0;
    static const double tau7   = (s1 - 2.0 * s2 - s3) / 3.0;
    static const double tau8   = (s1 + s2 + 2.0 * s3) / 3.0;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Extract coefficients from wavetable
        std::vector<std::complex<double> > w = get_w(wavetable, index, k, q, 6, sign);

        // Compute x = W(7) z
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const std::complex<double> z0 = *(in+istride*i);
            const std::complex<double> z1 = *(in+istride*(i+m));
            const std::complex<double> z2 = *(in+istride*(i+2*m));
            const std::complex<double> z3 = *(in+istride*(i+3*m));
            const std::complex<double> z4 = *(in+istride*(i+4*m));
            const std::complex<double> z5 = *(in+istride*(i+5*m));
            const std::complex<double> z6 = *(in+istride*(i+6*m));

            // Compute t
            const std::complex<double> t0 = z1 + z6;
            const std::complex<double> t1 = z1 - z6;
            const std::complex<double> t2 = z2 + z5;
            const std::complex<double> t3 = z2 - z5;
            const std::complex<double> t4 = z4 + z3;
            const std::complex<double> t5 = z4 - z3;
            const std::complex<double> t6 = t2 + t0;
            const std::complex<double> t7 = t5 + t3;
          
            // Compute b
            const std::complex<double> b0 = z0 + t6 + t4;
            const std::complex<double> b1 = tau1 * (t6 + t4);
            const std::complex<double> b2 = tau2 * (t0 - t4);
            const std::complex<double> b3 = tau3 * (t4 - t2);
            const std::complex<double> b4 = tau4 * (t2 - t0);
            const std::complex<double> b5 = double(-(int)sign) * tau5 * (t7 + t1);
            const std::complex<double> b6 = double(-(int)sign) * tau6 * (t1 - t5);
            const std::complex<double> b7 = double(-(int)sign) * tau7 * (t5 - t3);
            const std::complex<double> b8 = double(-(int)sign) * tau8 * (t3 - t1);
          
            // Compute T
            const std::complex<double> T0  = b0 + b1;
            const std::complex<double> T1  = b2 + b3;
            const std::complex<double> T2  = b4 - b3;
            const std::complex<double> T3  = -b2 - b4;
            const std::complex<double> T4  = b6 + b7;
            const std::complex<double> T5  = b8 - b7;
            const std::complex<double> T6  = -b8 - b6;
            const std::complex<double> T7  = T0 + T1;
            const std::complex<double> T8  = T0 + T2;
            const std::complex<double> T9  = T0 + T3;
            const std::complex<double> T10 = T4 + b5;
            const std::complex<double> T11 = T5 + b5;
            const std::complex<double> T12 = T6 + b5;
          
            // Compute x
            const std::complex<double> x0 = b0;
            const std::complex<double> x1 = T7 - timesi(T10);
            const std::complex<double> x2 = T9 - timesi(T12);
            const std::complex<double> x3 = T8 + timesi(T11);
            const std::complex<double> x4 = T8 - timesi(T11);
            const std::complex<double> x5 = T9 + timesi(T12);
            const std::complex<double> x6 = T7 + timesi(T10);
          
            // Compute out = w * x
            *(out+ostride*j)         = x0;
            *(out+ostride*(j+p_1))   = w[0] * x1;
            *(out+ostride*(j+2*p_1)) = w[1] * x2;
            *(out+ostride*(j+3*p_1)) = w[2] * x3;
            *(out+ostride*(j+4*p_1)) = w[3] * x4;
            *(out+ostride*(j+5*p_1)) = w[4] * x5;
            *(out+ostride*(j+6*p_1)) = w[5] * x6;

        } // endfor: k1

    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for arbitrary factor
 *
 * @param[in] in Pointer to input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Pointer to output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] factor Factor.
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factorn(std::complex<double>* in,
                   const int&            istride,
                   std::complex<double>* out,
                   const int&            ostride,
                   const GFftWavetable&  wavetable,
                   const int&            sign,
                   const int&            factor,
                   const int&            product,
                   const int&            n,
                   const int&            index)
{
    // Compute ...
    const int m    = n / factor;
    const int q    = n / product;
    const int p_1  = product / factor;
    const int jump = (factor - 1) * p_1;

    // ...
    for (int i = 0; i < m; ++i) {
        *(out+ostride*i) = *(in+istride*i);
    }

    // ...
    for (int e = 1; e < (factor - 1) / 2 + 1; ++e) {
        for (int i = 0; i < m; ++i) {
            const int idx  = i + e * m;
            const int idxc = i + (factor - e) * m;
            *(out+ostride*idx)  = *(in+istride*idx) + *(in+istride*idxc);
            *(out+ostride*idxc) = *(in+istride*idx) - *(in+istride*idxc);
        }
    }

    // e = 0
    for (int i = 0; i < m; ++i) {
        *(in+istride*i) = *(out+ostride*i);
    }
    for (int e = 1; e < (factor - 1) / 2 + 1; ++e) {
        for (int i = 0; i < m; ++i) {
            const int idx = i + e * m;
            *(in+istride*i) += *(out+ostride*idx);
        }
    }

    // ...
    for (int e = 1; e < (factor-1)/2 + 1; ++e) {

        //
        int       idx      = e * q;
        const int idx_step = e * q;

        //
        double w_real, w_imag;

        //
        const int em  = e * m;
        const int ecm = (factor - e) * m;

        // ...
        for (int i = 0; i < m; ++i) {
            *(in+istride*(i+em))  = *(out+ostride*i);
            *(in+istride*(i+ecm)) = *(out+ostride*i);
        }

        // ...
        for (int e1 = 1; e1 < (factor - 1) / 2 + 1; ++e1) {

            //
            if (idx == 0) {
                w_real = 1.0;
                w_imag = 0.0;
            }
            else {

                // Compute indices
                int twiddle = index + idx - 1;

                // Set trigonometric coefficients for forward transform
                if (sign == -1) {
                    w_real = wavetable[twiddle].real();
                    w_imag = wavetable[twiddle].imag();
                }
                // ... otherwise set trigonometric coefficients for backward
                // tranform: w -> conjugate(w)
                else {
                    w_real =  wavetable[twiddle].real();
                    w_imag = -wavetable[twiddle].imag();
                }
            }

            // Loop over ...
            for (int i = 0; i < m; ++i) {

                //
                const double xp_real = (out+ostride*(i + e1 * m))->real();
                const double xp_imag = (out+ostride*(i + e1 * m))->imag();
                const double xm_real = (out+ostride*(i + (factor-e1)*m))->real();
                const double xm_imag = (out+ostride*(i + (factor-e1)*m))->imag();

                //
                const double ap = w_real * xp_real;
                const double am = w_imag * xm_imag;

                //
                double sum_real  = ap - am;
                double sumc_real = ap + am;

                //
                const double bp = w_real * xp_imag;
                const double bm = w_imag * xm_real;

                //
                double sum_imag  = bp + bm;
                double sumc_imag = bp - bm;

                //
                *(in+istride*(i+em))  = std::complex<double>((in+istride*(i+em))->real() + sum_real,
                                                             (in+istride*(i+em))->imag() + sum_imag);
                *(in+istride*(i+ecm)) = std::complex<double>((in+istride*(i+ecm))->real() + sumc_real,
                                                             (in+istride*(i+ecm))->imag() + sumc_imag);
                

            } // endfor: i
            
            // Increment
            idx += idx_step ;
            idx %= factor * q ;
            
        } // endfor: e1
        
    } // endfor: e

    // k = 0
    for (int k1 = 0; k1 < p_1; ++k1) {
        *(out+ostride*k1) = *(in+istride*k1);
    }
    
    // k > 0
    for (int e1 = 1; e1 < factor; ++e1) {
        for (int k1 = 0; k1 < p_1; ++k1) {
            *(out+ostride*(k1 + e1 * p_1)) = *(in+istride*(k1 + e1 * m));
        }
    }

    // e = 0
    for (int k = 1, i = p_1, j = product; k < q; ++k, j += jump) {
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {
            *(out+ostride*j) = *(in+istride*i);
        }
    }

    // e > 0
    for (int k = 1, i = p_1, j = product; k < q; ++k, j += jump) {

        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            for (int e1 = 1; e1 < factor; ++e1) {
 
                // Get x
                double x_real = (in+istride*(i + e1 * m))->real();
                double x_imag = (in+istride*(i + e1 * m))->imag();

                // Get w
                double w_real, w_imag;
             
                // Compute indices
                int twiddle = index + (e1-1)*q + k-1;

                // Set trigonometric coefficients for forward transform
                if (sign == -1) {
                    w_real = wavetable[twiddle].real();
                    w_imag = wavetable[twiddle].imag();
                }
                // ... otherwise set trigonometric coefficients for backward
                // tranform: w -> conjugate(w)
                else {
                    w_real =  wavetable[twiddle].real();
                    w_imag = -wavetable[twiddle].imag();
                }

                //
                *(out+ostride*(j + e1 * p_1)) = std::complex<double>(w_real * x_real - w_imag * x_imag,
                                                                     w_real * x_imag + w_imag * x_real);

            } // endfor: e1

        } // endfor: k1
        
    } // endfor: k

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract coefficients from wavetable
 *
 * @param[in] index Start index of trigonometric coefficients.
 * @param[in] k ...
 * @param[in] q ...
 * @param[in] n Number of coefficients.
 * @param[in] sign Forward (-1) or backward (+1) transformation.
 * @param[in] wavetable Trigonometric coefficients.
 ***************************************************************************/
std::vector<std::complex<double> > GFft::get_w(const GFftWavetable& wavetable,
                                               const int&           index,
                                               const int&           k,
                                               const int&           q,
                                               const int&           n,
                                               const int&           sign) const
{
    // Allocate w vectors
    std::vector<std::complex<double> > w(n);

    // Set trigonometric coefficients for k=0 since they are not stored in the
    // wavetable object
    if (k == 0) {
        w.assign(n, std::complex<double>(1.0, 0.0));
    }
 
    // ... otherwise use coefficients stored in wavetable
    else {

        // Compute indices
        int twiddle = index + k - 1;

        // If sign==-1 then set trigonometric coefficients for forward transform
        if (sign == -1) {
            for (int i = 0; i < n; ++i, twiddle += q) {
                w[i] = wavetable[twiddle];
            }
        }

        // ... otherwise set trigonometric coefficients for backward tranform
        // w -> conjugate(w)
        else {
            for (int i = 0; i < n; ++i, twiddle += q) {
                w[i] = std::conj(wavetable[twiddle]);
            }
        }

    } // endelse: k != 0

    // Return w vectors
    return w;
}
