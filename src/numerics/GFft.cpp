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
    for (int i = 0; i < m_data.size(); ++i) {
        fft.m_data[i] = -fft.m_data[i];
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
    if (m_data.size() > 0) {

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
        std::vector<std::complex<double> > data = m_data;

        // Perform FFT (circumvent const correctness)
        const_cast<GFft*>(this)->transform(data, 1, m_shape[0], m_wavetable[0],
                                           false);

        // Get normalisation factor
        const double norm  = 1.0 / double(array.size());

        // Extract N-dimensional array and normalise with 1/n
        for (int i = 0; i < array.size(); ++i) {
            array(i) = data[i].real() * norm;
        }

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
        result.append(gammalib::str(size()));

        // VERBOSE: Put all FFT elements in string
        if (chatter >= VERBOSE) {
            result.append("\n"+gammalib::parformat("Elements"));
            for (int i = 0; i < m_data.size(); ++i) {
                if (i > 0) {
                    result.append(",");
                }
                result.append(gammalib::str(m_data[i]));
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
    m_shape.clear();
    m_strides.clear();
    m_data.clear();
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
    m_shape     = fft.m_shape;
    m_strides   = fft.m_strides;
    m_data      = fft.m_data;
    m_wavetable = fft.m_wavetable;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFft::free_members(void)
{
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
    m_shape.clear();
    m_strides.clear();
    m_data.clear();
    m_wavetable.clear();

    // Continue only if there are elements in the array
    if (array.size() > 0) {

        // Allocate data and initialise data to zero
        std::complex<double> zero(0.0, 0.0);
        m_data.assign(array.size(), zero);

        // Set real part of data from n-dimensional array
        for (int i = 0; i < array.size(); ++i) {
            m_data[i] = std::complex<double>(array(i), 0.0);
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
 * @param[in] data FFT array.
 * @param[in] stride Step size when traversing FFT array in one dimension.
 * @param[in] n Number of elements in the dimension.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] forward Forward transform if true, otherwise backward transform.
 ***************************************************************************/
void GFft::transform(std::vector<std::complex<double> >& data,
                     const int&                          stride,
                     const int&                          n,
                     const GFftWavetable&                wavetable,
                     const bool&                         forward)
{
    // Allocate scratch space
    std::complex<double> zero(0.0, 0.0);
    std::vector<std::complex<double> > scratch(n, zero);

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
        std::vector<std::complex<double> >* in  = &data;
        std::vector<std::complex<double> >* out = &scratch;
        int                                 istride;
        int                                 ostride;
        if (state == 0) {
            in      = &data;
            istride = stride;
            out     = &scratch;
            ostride = 1;
            state   = 1;
        }
        else {
            in      = &scratch;
            istride = 1;
            out     = &data;
            ostride = stride;
            state   = 0;
        }

        // Call factor dependent method
        if (factor == 2) {
            factor2(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 3) {
            factor3(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 4) {
            factor4(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 5) {
            factor5(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 6) {
            factor6(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else if (factor == 7) {
            factor7(*in, istride, *out, ostride, wavetable, sign, product, n, index);
        }
        else {
            factorn(*in, istride, *out, ostride, wavetable, sign, factor, product, n, index);
        }

    } // endfor: looped over all factors

    // In case the loop was exited in state 1 the results are in the scratch
    // array and we need to copy the results back from the scratch array to
    // the data array
    if (state == 1) {
        for (int i = 0; i < n; ++i) {
            data[stride*i] = scratch[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 2
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor2(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
{
    // Compute ...
    const int factor = 2;
    const int m      = n / factor;
    const int q      = n / product;
    const int p_1    = product / factor;
    const int jump   = (factor - 1) * p_1;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {

        // Allocate variables
        double w1_real, w1_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
        }
        else {

            // Compute indices
            int twiddle = index + k - 1;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle].real();
                w1_imag = wavetable[twiddle].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle].real();
                w1_imag = -wavetable[twiddle].imag();
            }

        } // endelse: k != 0

        // Loop over ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();

            // compute x = W(2) z

            // x0 = z0 + z1
            const double x0_real = z0_real + z1_real;
            const double x0_imag = z0_imag + z1_imag;

            // x1 = z0 - z1
            const double x1_real = z0_real - z1_real;
            const double x1_imag = z0_imag - z1_imag;

            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 3
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor3(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
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
    
        // Allocate variables
        double w1_real, w1_imag, w2_real, w2_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
            w2_real = 1.0;
            w2_imag = 0.0;
        }
        else {

            // Compute indices
            int twiddle1 = index + k - 1;
            int twiddle2 = twiddle1 + q;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle1].real();
                w1_imag = wavetable[twiddle1].imag();
                w2_real = wavetable[twiddle2].real();
                w2_imag = wavetable[twiddle2].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle1].real();
                w1_imag = -wavetable[twiddle1].imag();
                w2_real =  wavetable[twiddle2].real();
                w2_imag = -wavetable[twiddle2].imag();
            }
        } // endelse: k != 0

        // ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();
            const double z2_real = in[istride*(i+2*m)].real();
            const double z2_imag = in[istride*(i+2*m)].imag();

            // compute x = W(3) z

            // t1 = z1 + z2
            const double t1_real = z1_real + z2_real;
            const double t1_imag = z1_imag + z2_imag;
          
            // t2 = z0 - t1/2
            const double t2_real = z0_real - t1_real / 2.0;
            const double t2_imag = z0_imag - t1_imag / 2.0;
          
            // t3 = (+/-) sin(pi/3)*(z1 - z2)
            const double t3_real = ((int) sign) * tau * (z1_real - z2_real);
            const double t3_imag = ((int) sign) * tau * (z1_imag - z2_imag);
          
            // x0 = z0 + t1
            const double x0_real = z0_real + t1_real;
            const double x0_imag = z0_imag + t1_imag;
          
            // x1 = t2 + i t3
            const double x1_real = t2_real - t3_imag;
            const double x1_imag = t2_imag + t3_real;
          
            // x2 = t2 - i t3
            const double x2_real = t2_real + t3_imag;
            const double x2_imag = t2_imag - t3_real;

            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

            // out2 = w2 * x2
            //out[ostride*(j+2*p_1)].real(w2_real * x2_real - w2_imag * x2_imag);
            //out[ostride*(j+2*p_1)].imag(w2_real * x2_imag + w2_imag * x2_real);
            out[ostride*(j+2*p_1)] = std::complex<double>(w2_real * x2_real - w2_imag * x2_imag,
                                                          w2_real * x2_imag + w2_imag * x2_real);

        } // endfor: k1
        
    } // endfor: k

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 4
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor4(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
{
    // Compute ...
    const int factor = 4;
    const int m = n / factor;
    const int q = n / product;
    const int p_1 = product / factor;
    const int jump = (factor - 1) * p_1;

    // Loop over ...
    for (int k = 0, i = 0, j = 0; k < q; ++k, j += jump) {
    
        // Allocate w
        double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
            w2_real = 1.0;
            w2_imag = 0.0;
            w3_real = 1.0;
            w3_imag = 0.0;
        }
        else
        {

            // Compute indices
            int twiddle1 = index + k - 1;
            int twiddle2 = twiddle1 + q;
            int twiddle3 = twiddle2 + q;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle1].real();
                w1_imag = wavetable[twiddle1].imag();
                w2_real = wavetable[twiddle2].real();
                w2_imag = wavetable[twiddle2].imag();
                w3_real = wavetable[twiddle3].real();
                w3_imag = wavetable[twiddle3].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle1].real();
                w1_imag = -wavetable[twiddle1].imag();
                w2_real =  wavetable[twiddle2].real();
                w2_imag = -wavetable[twiddle2].imag();
                w3_real =  wavetable[twiddle3].real();
                w3_imag = -wavetable[twiddle3].imag();
            }

        } // endelse: k != 0

        // Loop over ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {


            // Get z
            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();
            const double z2_real = in[istride*(i+2*m)].real();
            const double z2_imag = in[istride*(i+2*m)].imag();
            const double z3_real = in[istride*(i+3*m)].real();
            const double z3_imag = in[istride*(i+3*m)].imag();

            // compute x = W(4) z
            
            // t1 = z0 + z2
            const double t1_real = z0_real + z2_real;
            const double t1_imag = z0_imag + z2_imag;
          
            // t2 = z1 + z3
            const double t2_real = z1_real + z3_real;
            const double t2_imag = z1_imag + z3_imag;
          
            // t3 = z0 - z2
            const double t3_real = z0_real - z2_real;
            const double t3_imag = z0_imag - z2_imag;
          
            // t4 = (+/-) (z1 - z3)
            const double t4_real = ((int) sign) * (z1_real - z3_real);
            const double t4_imag = ((int) sign) * (z1_imag - z3_imag);

            // x0 = t1 + t2
            const double x0_real = t1_real + t2_real;
            const double x0_imag = t1_imag + t2_imag;

            // x1 = t3 + i t4
            const double x1_real = t3_real - t4_imag;
            const double x1_imag = t3_imag + t4_real;

            // x2 = t1 - t2
            const double x2_real = t1_real - t2_real;
            const double x2_imag = t1_imag - t2_imag;

            // x3 = t3 - i t4
            const double x3_real = t3_real + t4_imag;
            const double x3_imag = t3_imag - t4_real;

            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

            // out2 = w2 * x2
            //out[ostride*(j+2*p_1)].real(w2_real * x2_real - w2_imag * x2_imag);
            //out[ostride*(j+2*p_1)].imag(w2_real * x2_imag + w2_imag * x2_real);
            out[ostride*(j+2*p_1)] = std::complex<double>(w2_real * x2_real - w2_imag * x2_imag,
                                                          w2_real * x2_imag + w2_imag * x2_real);

            // out3 = w3 * x3
            //out[ostride*(j+3*p_1)].real(w3_real * x3_real - w3_imag * x3_imag);
            //out[ostride*(j+3*p_1)].imag(w3_real * x3_imag + w3_imag * x3_real);
            out[ostride*(j+3*p_1)] = std::complex<double>(w3_real * x3_real - w3_imag * x3_imag,
                                                          w3_real * x3_imag + w3_imag * x3_real);
          
        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 5
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor5(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
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

        // Allocate w
        double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
            w2_real = 1.0;
            w2_imag = 0.0;
            w3_real = 1.0;
            w3_imag = 0.0;
            w4_real = 1.0;
            w4_imag = 0.0;
        }
        else {

            // Compute indices
            int twiddle1 = index + k - 1;
            int twiddle2 = twiddle1 + q;
            int twiddle3 = twiddle2 + q;
            int twiddle4 = twiddle3 + q;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle1].real();
                w1_imag = wavetable[twiddle1].imag();
                w2_real = wavetable[twiddle2].real();
                w2_imag = wavetable[twiddle2].imag();
                w3_real = wavetable[twiddle3].real();
                w3_imag = wavetable[twiddle3].imag();
                w4_real = wavetable[twiddle4].real();
                w4_imag = wavetable[twiddle4].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle1].real();
                w1_imag = -wavetable[twiddle1].imag();
                w2_real =  wavetable[twiddle2].real();
                w2_imag = -wavetable[twiddle2].imag();
                w3_real =  wavetable[twiddle3].real();
                w3_imag = -wavetable[twiddle3].imag();
                w4_real =  wavetable[twiddle4].real();
                w4_imag = -wavetable[twiddle4].imag();
            }

        } // endelse: k != 0

        // Loop over ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();
            const double z2_real = in[istride*(i+2*m)].real();
            const double z2_imag = in[istride*(i+2*m)].imag();
            const double z3_real = in[istride*(i+3*m)].real();
            const double z3_imag = in[istride*(i+3*m)].imag();
            const double z4_real = in[istride*(i+4*m)].real();
            const double z4_imag = in[istride*(i+4*m)].imag();

            // compute x = W(5) z

            // t1 = z1 + z4
            const double t1_real = z1_real + z4_real;
            const double t1_imag = z1_imag + z4_imag;
          
            // t2 = z2 + z3
            const double t2_real = z2_real + z3_real;
            const double t2_imag = z2_imag + z3_imag;
          
            // t3 = z1 - z4
            const double t3_real = z1_real - z4_real;
            const double t3_imag = z1_imag - z4_imag;
          
            // t4 = z2 - z3
            const double t4_real = z2_real - z3_real;
            const double t4_imag = z2_imag - z3_imag;
          
            // t5 = t1 + t2
            const double t5_real = t1_real + t2_real;
            const double t5_imag = t1_imag + t2_imag;
          
            // t6 = (sqrt(5)/4)(t1 - t2)
            const double t6_real = (tau) * (t1_real - t2_real);
            const double t6_imag = (tau) * (t1_imag - t2_imag);
          
            // t7 = z0 - ((t5)/4)
            const double t7_real = z0_real - t5_real / 4.0;
            const double t7_imag = z0_imag - t5_imag / 4.0;
          
            // t8 = t7 + t6
            const double t8_real = t7_real + t6_real;
            const double t8_imag = t7_imag + t6_imag;
          
            // t9 = t7 - t6
            const double t9_real = t7_real - t6_real;
            const double t9_imag = t7_imag - t6_imag;
          
            // t10 = sin(2 pi/5) t3 + sin(2 pi/10) t4 */
            const double t10_real = ((int) sign) * (sin_2pi_by_5 * t3_real +
                                                    sin_2pi_by_10 * t4_real);
            const double t10_imag = ((int) sign) * (sin_2pi_by_5 * t3_imag +
                                                    sin_2pi_by_10 * t4_imag);
          
            // t11 = sin(2 pi/10) t3 - sin(2 pi/5) t4
            const double t11_real = ((int) sign) * (sin_2pi_by_10 * t3_real -
                                                    sin_2pi_by_5 * t4_real);
            const double t11_imag = ((int) sign) * (sin_2pi_by_10 * t3_imag -
                                                    sin_2pi_by_5 * t4_imag);
          
            // x0 = z0 + t5
            const double x0_real = z0_real + t5_real;
            const double x0_imag = z0_imag + t5_imag;
          
            // x1 = t8 + i t10
            const double x1_real = t8_real - t10_imag;
            const double x1_imag = t8_imag + t10_real;
          
            // x2 = t9 + i t11
            const double x2_real = t9_real - t11_imag;
            const double x2_imag = t9_imag + t11_real;
          
            // x3 = t9 - i t11
            const double x3_real = t9_real + t11_imag;
            const double x3_imag = t9_imag - t11_real;
          
            // x4 = t8 - i t10
            const double x4_real = t8_real + t10_imag;
            const double x4_imag = t8_imag - t10_real;

            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

            // out2 = w2 * x2
            //out[ostride*(j+2*p_1)].real(w2_real * x2_real - w2_imag * x2_imag);
            //out[ostride*(j+2*p_1)].imag(w2_real * x2_imag + w2_imag * x2_real);
            out[ostride*(j+2*p_1)] = std::complex<double>(w2_real * x2_real - w2_imag * x2_imag,
                                                          w2_real * x2_imag + w2_imag * x2_real);

            // out3 = w3 * x3
            //out[ostride*(j+3*p_1)].real(w3_real * x3_real - w3_imag * x3_imag);
            //out[ostride*(j+3*p_1)].imag(w3_real * x3_imag + w3_imag * x3_real);
            out[ostride*(j+3*p_1)] = std::complex<double>(w3_real * x3_real - w3_imag * x3_imag,
                                                          w3_real * x3_imag + w3_imag * x3_real);

            // out4 = w4 * x4
            //out[ostride*(j+4*p_1)].real(w4_real * x4_real - w4_imag * x4_imag);
            //out[ostride*(j+4*p_1)].imag(w4_real * x4_imag + w4_imag * x4_real);
            out[ostride*(j+4*p_1)] = std::complex<double>(w4_real * x4_real - w4_imag * x4_imag,
                                                          w4_real * x4_imag + w4_imag * x4_real);
          
        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 6
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor6(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
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
    
        // Allocate w
        double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag, w5_real, w5_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
            w2_real = 1.0;
            w2_imag = 0.0;
            w3_real = 1.0;
            w3_imag = 0.0;
            w4_real = 1.0;
            w4_imag = 0.0;
            w5_real = 1.0;
            w5_imag = 0.0;
        }
        else {

            // Compute indices
            int twiddle1 = index + k - 1;
            int twiddle2 = twiddle1 + q;
            int twiddle3 = twiddle2 + q;
            int twiddle4 = twiddle3 + q;
            int twiddle5 = twiddle4 + q;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle1].real();
                w1_imag = wavetable[twiddle1].imag();
                w2_real = wavetable[twiddle2].real();
                w2_imag = wavetable[twiddle2].imag();
                w3_real = wavetable[twiddle3].real();
                w3_imag = wavetable[twiddle3].imag();
                w4_real = wavetable[twiddle4].real();
                w4_imag = wavetable[twiddle4].imag();
                w5_real = wavetable[twiddle5].real();
                w5_imag = wavetable[twiddle5].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle1].real();
                w1_imag = -wavetable[twiddle1].imag();
                w2_real =  wavetable[twiddle2].real();
                w2_imag = -wavetable[twiddle2].imag();
                w3_real =  wavetable[twiddle3].real();
                w3_imag = -wavetable[twiddle3].imag();
                w4_real =  wavetable[twiddle4].real();
                w4_imag = -wavetable[twiddle4].imag();
                w5_real =  wavetable[twiddle5].real();
                w5_imag = -wavetable[twiddle5].imag();
            }

        } // endelse: k != 0

        // ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();
            const double z2_real = in[istride*(i+2*m)].real();
            const double z2_imag = in[istride*(i+2*m)].imag();
            const double z3_real = in[istride*(i+3*m)].real();
            const double z3_imag = in[istride*(i+3*m)].imag();
            const double z4_real = in[istride*(i+4*m)].real();
            const double z4_imag = in[istride*(i+4*m)].imag();
            const double z5_real = in[istride*(i+5*m)].real();
            const double z5_imag = in[istride*(i+5*m)].imag();

            // compute x = W(6) z

            // W(6) is a combination of sums and differences of W(3) acting
            // on the even and odd elements of z
          
            // ta1 = z2 + z4
            const double ta1_real = z2_real + z4_real;
            const double ta1_imag = z2_imag + z4_imag;
          
            // ta2 = z0 - ta1/2
            const double ta2_real = z0_real - ta1_real / 2.0;
            const double ta2_imag = z0_imag - ta1_imag / 2.0;
          
            // ta3 = (+/-) sin(pi/3)*(z2 - z4)
            const double ta3_real = ((int) sign) * tau * (z2_real - z4_real);
            const double ta3_imag = ((int) sign) * tau * (z2_imag - z4_imag);
          
            // a0 = z0 + ta1
            const double a0_real = z0_real + ta1_real;
            const double a0_imag = z0_imag + ta1_imag;
          
            // a1 = ta2 + i ta3
            const double a1_real = ta2_real - ta3_imag;
            const double a1_imag = ta2_imag + ta3_real;
          
            // a2 = ta2 - i ta3
            const double a2_real = ta2_real + ta3_imag;
            const double a2_imag = ta2_imag - ta3_real;
          
            // tb1 = z5 + z1
            const double tb1_real = z5_real + z1_real;
            const double tb1_imag = z5_imag + z1_imag;
          
            // tb2 = z3 - tb1/2
            const double tb2_real = z3_real - tb1_real / 2.0;
            const double tb2_imag = z3_imag - tb1_imag / 2.0;
          
            // tb3 = (+/-) sin(pi/3)*(z5 - z1)
            const double tb3_real = ((int) sign) * tau * (z5_real - z1_real);
            const double tb3_imag = ((int) sign) * tau * (z5_imag - z1_imag);
          
            // b0 = z3 + tb1
            const double b0_real = z3_real + tb1_real;
            const double b0_imag = z3_imag + tb1_imag;
          
            // b1 = tb2 + i tb3
            const double b1_real = tb2_real - tb3_imag;
            const double b1_imag = tb2_imag + tb3_real;
          
            // b2 = tb2 - i tb3
            const double b2_real = tb2_real + tb3_imag;
            const double b2_imag = tb2_imag - tb3_real;
          
            // x0 = a0 + b0
            const double x0_real = a0_real + b0_real;
            const double x0_imag = a0_imag + b0_imag;
          
            // x4 = a1 + b1
            const double x4_real = a1_real + b1_real;
            const double x4_imag = a1_imag + b1_imag;
          
            // x2 = a2 + b2
            const double x2_real = a2_real + b2_real;
            const double x2_imag = a2_imag + b2_imag;
          
            // x3 = a0 - b0
            const double x3_real = a0_real - b0_real;
            const double x3_imag = a0_imag - b0_imag;
          
            // x1 = a1 - b1
            const double x1_real = a1_real - b1_real;
            const double x1_imag = a1_imag - b1_imag;
          
            // x5 = a2 - b2
            const double x5_real = a2_real - b2_real;
            const double x5_imag = a2_imag - b2_imag;

            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

            // out2 = w2 * x2
            //out[ostride*(j+2*p_1)].real(w2_real * x2_real - w2_imag * x2_imag);
            //out[ostride*(j+2*p_1)].imag(w2_real * x2_imag + w2_imag * x2_real);
            out[ostride*(j+2*p_1)] = std::complex<double>(w2_real * x2_real - w2_imag * x2_imag,
                                                          w2_real * x2_imag + w2_imag * x2_real);

            // out3 = w3 * x3
            //out[ostride*(j+3*p_1)].real(w3_real * x3_real - w3_imag * x3_imag);
            //out[ostride*(j+3*p_1)].imag(w3_real * x3_imag + w3_imag * x3_real);
            out[ostride*(j+3*p_1)] = std::complex<double>(w3_real * x3_real - w3_imag * x3_imag,
                                                          w3_real * x3_imag + w3_imag * x3_real);

            // out4 = w4 * x4
            //out[ostride*(j+4*p_1)].real(w4_real * x4_real - w4_imag * x4_imag);
            //out[ostride*(j+4*p_1)].imag(w4_real * x4_imag + w4_imag * x4_real);
            out[ostride*(j+4*p_1)] = std::complex<double>(w4_real * x4_real - w4_imag * x4_imag,
                                                          w4_real * x4_imag + w4_imag * x4_real);
            
            // out5 = w5 * x5
            //out[ostride*(j+5*p_1)].real(w5_real * x5_real - w5_imag * x5_imag);
            //out[ostride*(j+5*p_1)].imag(w5_real * x5_imag + w5_imag * x5_real);
            out[ostride*(j+5*p_1)] = std::complex<double>(w5_real * x5_real - w5_imag * x5_imag,
                                                          w5_real * x5_imag + w5_imag * x5_real);

        } // endfor: k1
        
    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for factor 7
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factor7(const std::vector<std::complex<double> >& in,
                   const int&                                istride,
                   std::vector<std::complex<double> >&       out,
                   const int&                                ostride,
                   const GFftWavetable&                      wavetable,
                   const int&                                sign,
                   const int&                                product,
                   const int&                                n,
                   const int&                                index)
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
    
        //
        double w1_real, w1_imag, w2_real, w2_imag, w3_real, w3_imag, w4_real,
        w4_imag, w5_real, w5_imag, w6_real, w6_imag;

        // Set trigonometric coefficients for k=0 since they are not stored
        // in the wavetable object
        if (k == 0) {
            w1_real = 1.0;
            w1_imag = 0.0;
            w2_real = 1.0;
            w2_imag = 0.0;
            w3_real = 1.0;
            w3_imag = 0.0;
            w4_real = 1.0;
            w4_imag = 0.0;
            w5_real = 1.0;
            w5_imag = 0.0;
            w6_real = 1.0;
            w6_imag = 0.0;
        }
        else {

            // Compute indices
            int twiddle1 = index + k - 1;
            int twiddle2 = twiddle1 + q;
            int twiddle3 = twiddle2 + q;
            int twiddle4 = twiddle3 + q;
            int twiddle5 = twiddle4 + q;
            int twiddle6 = twiddle5 + q;

            // Set trigonometric coefficients for forward transform
            if (sign == -1) {
                w1_real = wavetable[twiddle1].real();
                w1_imag = wavetable[twiddle1].imag();
                w2_real = wavetable[twiddle2].real();
                w2_imag = wavetable[twiddle2].imag();
                w3_real = wavetable[twiddle3].real();
                w3_imag = wavetable[twiddle3].imag();
                w4_real = wavetable[twiddle4].real();
                w4_imag = wavetable[twiddle4].imag();
                w5_real = wavetable[twiddle5].real();
                w5_imag = wavetable[twiddle5].imag();
                w6_real = wavetable[twiddle6].real();
                w6_imag = wavetable[twiddle6].imag();
            }

            // ... otherwise set trigonometric coefficients for backward
            // tranform: w -> conjugate(w)
            else {
                w1_real =  wavetable[twiddle1].real();
                w1_imag = -wavetable[twiddle1].imag();
                w2_real =  wavetable[twiddle2].real();
                w2_imag = -wavetable[twiddle2].imag();
                w3_real =  wavetable[twiddle3].real();
                w3_imag = -wavetable[twiddle3].imag();
                w4_real =  wavetable[twiddle4].real();
                w4_imag = -wavetable[twiddle4].imag();
                w5_real =  wavetable[twiddle5].real();
                w5_imag = -wavetable[twiddle5].imag();
                w6_real =  wavetable[twiddle6].real();
                w6_imag = -wavetable[twiddle6].imag();
            }

        } // endelse: k != 0

        // ...
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            // Get z
            const double z0_real = in[istride*i].real();
            const double z0_imag = in[istride*i].imag();
            const double z1_real = in[istride*(i+m)].real();
            const double z1_imag = in[istride*(i+m)].imag();
            const double z2_real = in[istride*(i+2*m)].real();
            const double z2_imag = in[istride*(i+2*m)].imag();
            const double z3_real = in[istride*(i+3*m)].real();
            const double z3_imag = in[istride*(i+3*m)].imag();
            const double z4_real = in[istride*(i+4*m)].real();
            const double z4_imag = in[istride*(i+4*m)].imag();
            const double z5_real = in[istride*(i+5*m)].real();
            const double z5_imag = in[istride*(i+5*m)].imag();
            const double z6_real = in[istride*(i+6*m)].real();
            const double z6_imag = in[istride*(i+6*m)].imag();

            // Compute x = W(7) z
          
            // t0 = z1 + z6
            const double t0_real = z1_real + z6_real;
            const double t0_imag = z1_imag + z6_imag;
          
            // t1 = z1 - z6
            const double t1_real = z1_real - z6_real;
            const double t1_imag = z1_imag - z6_imag;
          
            // t2 = z2 + z5
            const double t2_real = z2_real + z5_real;
            const double t2_imag = z2_imag + z5_imag;
          
            // t3 = z2 - z5
            const double t3_real = z2_real - z5_real;
            const double t3_imag = z2_imag - z5_imag;
          
            // t4 = z4 + z3
            const double t4_real = z4_real + z3_real;
            const double t4_imag = z4_imag + z3_imag;
          
            // t5 = z4 - z3
            const double t5_real = z4_real - z3_real;
            const double t5_imag = z4_imag - z3_imag;
          
            // t6 = t2 + t0
            const double t6_real = t2_real + t0_real;
            const double t6_imag = t2_imag + t0_imag;
          
            // t7 = t5 + t3
            const double t7_real = t5_real + t3_real;
            const double t7_imag = t5_imag + t3_imag;
          
            // b0 = z0 + t6 + t4
            const double b0_real = z0_real + t6_real + t4_real;
            const double b0_imag = z0_imag + t6_imag + t4_imag;

            // b1 = ((cos(2pi/7) + cos(4pi/7) + cos(6pi/7))/3-1) (t6 + t4)
            const double b1_real = tau1 * (t6_real + t4_real);
            const double b1_imag = tau1 * (t6_imag + t4_imag);
          
            // b2 = ((2*cos(2pi/7) - cos(4pi/7) - cos(6pi/7))/3) (t0 - t4)
            const double b2_real = tau2 * (t0_real - t4_real);
            const double b2_imag = tau2 * (t0_imag - t4_imag);
          
            // b3 = ((cos(2pi/7) - 2*cos(4pi/7) + cos(6pi/7))/3) (t4 - t2)
            const double b3_real = tau3 * (t4_real - t2_real);
            const double b3_imag = tau3 * (t4_imag - t2_imag);
          
            // b4 = ((cos(2pi/7) + cos(4pi/7) - 2*cos(6pi/7))/3) (t2 - t0)
            const double b4_real = tau4 * (t2_real - t0_real);
            const double b4_imag = tau4 * (t2_imag - t0_imag);
          
            // b5 = sign * ((sin(2pi/7) + sin(4pi/7) - sin(6pi/7))/3) (t7 + t1)
            const double b5_real = (-(int)sign) * tau5 * (t7_real + t1_real);
            const double b5_imag = (-(int)sign) * tau5 * (t7_imag + t1_imag);
          
            // b6 = sign * ((2sin(2pi/7) - sin(4pi/7) + sin(6pi/7))/3) (t1 - t5)
            const double b6_real = (-(int)sign) * tau6 * (t1_real - t5_real);
            const double b6_imag = (-(int)sign) * tau6 * (t1_imag - t5_imag);
          
            // b7 = sign * ((sin(2pi/7) - 2sin(4pi/7) - sin(6pi/7))/3) (t5 - t3)
            const double b7_real = (-(int)sign) * tau7 * (t5_real - t3_real);
            const double b7_imag = (-(int)sign) * tau7 * (t5_imag - t3_imag);
          
            // b8 = sign * ((sin(2pi/7) + sin(4pi/7) + 2sin(6pi/7))/3) (t3 - t1) */
            const double b8_real = (-(int)sign) * tau8 * (t3_real - t1_real);
            const double b8_imag = (-(int)sign) * tau8 * (t3_imag - t1_imag);
          
            // T0 = b0 + b1
            const double T0_real = b0_real + b1_real;
            const double T0_imag = b0_imag + b1_imag;
          
            // T1 = b2 + b3
            const double T1_real = b2_real + b3_real;
            const double T1_imag = b2_imag + b3_imag;
          
            // T2 = b4 - b3
            const double T2_real = b4_real - b3_real;
            const double T2_imag = b4_imag - b3_imag;
          
            // T3 = -b2 - b4
            const double T3_real = -b2_real - b4_real;
            const double T3_imag = -b2_imag - b4_imag;
          
            // T4 = b6 + b7
            const double T4_real = b6_real + b7_real;
            const double T4_imag = b6_imag + b7_imag;
          
            // T5 = b8 - b7
            const double T5_real = b8_real - b7_real;
            const double T5_imag = b8_imag - b7_imag;
          
            // T6 = -b8 - b6
            const double T6_real = -b8_real - b6_real;
            const double T6_imag = -b8_imag - b6_imag;
          
            // T7 = T0 + T1
            const double T7_real = T0_real + T1_real;
            const double T7_imag = T0_imag + T1_imag;
          
            // T8 = T0 + T2
            const double T8_real = T0_real + T2_real;
            const double T8_imag = T0_imag + T2_imag;
          
            // T9 = T0 + T3
            const double T9_real = T0_real + T3_real;
            const double T9_imag = T0_imag + T3_imag;
          
            // T10 = T4 + b5
            const double T10_real = T4_real + b5_real;
            const double T10_imag = T4_imag + b5_imag;
          
            // T11 = T5 + b5
            const double T11_real = T5_real + b5_real;
            const double T11_imag = T5_imag + b5_imag;
          
            // T12 = T6 + b5
            const double T12_real = T6_real + b5_real;
            const double T12_imag = T6_imag + b5_imag;
          
            // x0 = b0
            const double x0_real = b0_real;
            const double x0_imag = b0_imag;
          
            // x1 = T7 - i T10
            const double x1_real = T7_real + T10_imag;
            const double x1_imag = T7_imag - T10_real;
          
            // x2 = T9 - i T12
            const double x2_real = T9_real + T12_imag;
            const double x2_imag = T9_imag - T12_real;
          
            // x3 = T8 + i T11
            const double x3_real = T8_real - T11_imag;
            const double x3_imag = T8_imag + T11_real;
          
            // x4 = T8 - i T11
            const double x4_real = T8_real + T11_imag;
            const double x4_imag = T8_imag - T11_real;
          
            // x5 = T9 + i T12
            const double x5_real = T9_real - T12_imag;
            const double x5_imag = T9_imag + T12_real;
          
            // x6 = T7 + i T10
            const double x6_real = T7_real - T10_imag;
            const double x6_imag = T7_imag + T10_real;
          
            // out0 = 1 * x0
            //out[ostride*j].real(x0_real);
            //out[ostride*j].imag(x0_imag);
            out[ostride*j] = std::complex<double>(x0_real, x0_imag);
          
            // out1 = w1 * x1
            //out[ostride*(j+p_1)].real(w1_real * x1_real - w1_imag * x1_imag);
            //out[ostride*(j+p_1)].imag(w1_real * x1_imag + w1_imag * x1_real);
            out[ostride*(j+p_1)] = std::complex<double>(w1_real * x1_real - w1_imag * x1_imag,
                                                        w1_real * x1_imag + w1_imag * x1_real);

            // out2 = w2 * x2
            //out[ostride*(j+2*p_1)].real(w2_real * x2_real - w2_imag * x2_imag);
            //out[ostride*(j+2*p_1)].imag(w2_real * x2_imag + w2_imag * x2_real);
            out[ostride*(j+2*p_1)] = std::complex<double>(w2_real * x2_real - w2_imag * x2_imag,
                                                          w2_real * x2_imag + w2_imag * x2_real);

            // out3 = w3 * x3
            //out[ostride*(j+3*p_1)].real(w3_real * x3_real - w3_imag * x3_imag);
            //out[ostride*(j+3*p_1)].imag(w3_real * x3_imag + w3_imag * x3_real);
            out[ostride*(j+3*p_1)] = std::complex<double>(w3_real * x3_real - w3_imag * x3_imag,
                                                          w3_real * x3_imag + w3_imag * x3_real);

            // out4 = w4 * x4
            //out[ostride*(j+4*p_1)].real(w4_real * x4_real - w4_imag * x4_imag);
            //out[ostride*(j+4*p_1)].imag(w4_real * x4_imag + w4_imag * x4_real);
            out[ostride*(j+4*p_1)] = std::complex<double>(w4_real * x4_real - w4_imag * x4_imag,
                                                          w4_real * x4_imag + w4_imag * x4_real);

            // out5 = w5 * x5
            //out[ostride*(j+5*p_1)].real(w5_real * x5_real - w5_imag * x5_imag);
            //out[ostride*(j+5*p_1)].imag(w5_real * x5_imag + w5_imag * x5_real);
            out[ostride*(j+5*p_1)] = std::complex<double>(w5_real * x5_real - w5_imag * x5_imag,
                                                          w5_real * x5_imag + w5_imag * x5_real);
          
            // out6 = w6 * x6
            //out[ostride*(j+6*p_1)].real(w6_real * x6_real - w6_imag * x6_imag);
            //out[ostride*(j+6*p_1)].imag(w6_real * x6_imag + w6_imag * x6_real);
            out[ostride*(j+6*p_1)] = std::complex<double>(w6_real * x6_real - w6_imag * x6_imag,
                                                          w6_real * x6_imag + w6_imag * x6_real);
       
        } // endfor: k1

    } // endfor: k
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT for arbitrary factor
 *
 * @param[in] in Input array.
 * @param[in] istride Step size when traversing input array.
 * @param[in] out Output array.
 * @param[in] ostride Step size when traversing output array.
 * @param[in] wavetable Trigonometric coefficients.
 * @param[in] sign Forward (-1) or backward (+1).
 * @param[in] factor Factor.
 * @param[in] product ...
 * @param[in] n Logical array length.
 * @param[in] index Start index of trigonometric coefficients.
 ***************************************************************************/
void GFft::factorn(std::vector<std::complex<double> >& in,
                   const int&                          istride,
                   std::vector<std::complex<double> >& out,
                   const int&                          ostride,
                   const GFftWavetable&                wavetable,
                   const int&                          sign,
                   const int&                          factor,
                   const int&                          product,
                   const int&                          n,
                   const int&                          index)
{
    // Compute ...
    const int m    = n / factor;
    const int q    = n / product;
    const int p_1  = product / factor;
    const int jump = (factor - 1) * p_1;

    // ...
    for (int i = 0; i < m; ++i) {
        out[ostride*i] = in[istride*i];
    }

    // ...
    for (int e = 1; e < (factor - 1) / 2 + 1; ++e) {
        for (int i = 0; i < m; ++i) {
            const int idx  = i + e * m;
            const int idxc = i + (factor - e) * m;
            out[ostride*idx]  = in[istride*idx] + in[istride*idxc];
            out[ostride*idxc] = in[istride*idx] - in[istride*idxc];
        }
    }

    // e = 0
    for (int i = 0; i < m; ++i) {
        in[istride*i] = out[ostride*i];
    }
    for (int e = 1; e < (factor - 1) / 2 + 1; ++e) {
        for (int i = 0; i < m; ++i) {
            const int idx  = i + e * m;
            in[istride*i] += out[ostride*idx];
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
            in[istride*(i+em)]  = out[ostride*i];
            in[istride*(i+ecm)] = out[ostride*i];
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
                const double xp_real = out[ostride*(i + e1 * m)].real();
                const double xp_imag = out[ostride*(i + e1 * m)].imag();
                const double xm_real = out[ostride*(i + (factor-e1)*m)].real();
                const double xm_imag = out[ostride*(i + (factor-e1)*m)].imag();

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
                in[istride*(i+em)].real(in[istride*(i+em)].real() + sum_real);
                in[istride*(i+em)].imag(in[istride*(i+em)].imag() + sum_imag);
                in[istride*(i+ecm)].real(in[istride*(i+ecm)].real() + sumc_real);
                in[istride*(i+ecm)].imag(in[istride*(i+ecm)].imag() + sumc_imag);

            } // endfor: i
            
            // Increment
            idx += idx_step ;
            idx %= factor * q ;
            
        } // endfor: e1
        
    } // endfor: e

    // k = 0
    for (int k1 = 0; k1 < p_1; ++k1) {
        out[ostride*k1] = in[istride*k1];
    }
    
    // k > 0
    for (int e1 = 1; e1 < factor; ++e1) {
        for (int k1 = 0; k1 < p_1; ++k1) {
            out[ostride*(k1 + e1 * p_1)] = in[istride*(k1 + e1 * m)];
        }
    }

    // e = 0
    for (int k = 1, i = p_1, j = product; k < q; ++k, j += jump) {
        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {
            out[ostride*j] = in[istride*i];
        }
    }

    // e > 0
    for (int k = 1, i = p_1, j = product; k < q; ++k, j += jump) {

        for (int k1 = 0; k1 < p_1; ++k1, ++i, ++j) {

            for (int e1 = 1; e1 < factor; ++e1) {
 
                // Get x
                double x_real = in[istride*(i + e1 * m)].real();
                double x_imag = in[istride*(i + e1 * m)].imag();

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
                //out[ostride*(j + e1 * p_1)].real(w_real * x_real - w_imag * x_imag);
                //out[ostride*(j + e1 * p_1)].imag(w_real * x_imag + w_imag * x_real);
                out[ostride*(j + e1 * p_1)] = std::complex<double>(w_real * x_real - w_imag * x_imag,
                                                                   w_real * x_imag + w_imag * x_real);

            } // endfor: e1

        } // endfor: k1
        
    } // endfor: k

    // Return
    return;
}
