/***************************************************************************
 *  GFftWavetable.cpp - Lookup table class for Fast Fourier transformation *
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
 * @file GFftWavetable.cpp
 * @brief Lookup table class implementation for Fast Fourier transformation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GFftWavetable.hpp"


/* __ Method name definitions ____________________________________________ */
#define G_FACTOR                                "GFftWavetable::factor(int&)"
#define G_SET_MEMBERS                      "GFftWavetable::set_members(int&)"
#define G_SET_FACTORS                      "GFftWavetable::set_factors(int&)"


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFftWavetable::GFftWavetable(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Wave table constructor
 *
 * @param[in] size Array size.
 *
 * Constructs the lookup table for Fast Fourier Transformation for a given
 * array @p size.
 ***************************************************************************/
GFftWavetable::GFftWavetable(const int& size)
{
    // Initialise class members
    init_members();

    // Set members
    set_members(size);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] wavetable Lookup table for Fast Fourier Transform.
 ***************************************************************************/
GFftWavetable::GFftWavetable(const GFftWavetable& wavetable)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(wavetable);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFftWavetable::~GFftWavetable(void)
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
 * @param[in] wavetable Lookup table for Fast Fourier Transform.
 * @return Lookup table for Fast Fourier Transform.
 ***************************************************************************/
GFftWavetable& GFftWavetable::operator=(const GFftWavetable& wavetable)
{
    // Execute only if object is not identical
    if (this != &wavetable) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(wavetable);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear lookup table for Fast Fourier Transform
 ***************************************************************************/
void GFftWavetable::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone lookup table for Fast Fourier Transform
 *
 * @return Pointer to deep copy of lookup table for Fast Fourier Transform.
 ***************************************************************************/
GFftWavetable* GFftWavetable::clone(void) const
{
    // Clone FFT
    return new GFftWavetable(*this);
}


/***********************************************************************//**
 * @brief Return factorisation factor
 *
 * @param[in] index Factorisation index [0,...,factors()-1].
 * @return Factorisation factor.
 *
 * @exception GException::out_of_range
 *            Index is outside valid range.
 *
 * Returns factorisation factor.
 ***************************************************************************/
int GFftWavetable::factor(const int& index) const
{
    // Throw an exception if index is out of range
    if (index < 0 || index >= factors()) {
        throw GException::out_of_range(G_FACTOR, "Factorisation factor index",
                                       index, factors());
    }

    // Return factorisation factor
    return (m_factors[index]);
}


/***********************************************************************//**
 * @brief Print lookup table for Fast Fourier Transform information
 *
 * @param[in] chatter Chattiness.
 * @return String containing lookup table for Fast Fourier Transform
 *         information.
 ***************************************************************************/
std::string GFftWavetable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFftWavetable ===");

        // Append information
        result.append("\n"+gammalib::parformat("Coefficients"));
        result.append(gammalib::str(size()));
        result.append("\n"+gammalib::parformat("Factorisation"));
        if (factors() > 0) {
            for (int i = 0; i < factors(); ++i) {
                if (i > 0) {
                    result.append("*");
                }
                result.append(gammalib::str(factor(i)));
            }
        }
        else {
            result.append("none");
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
void GFftWavetable::init_members(void)
{
    // Initialise members
    m_factors.clear();
    m_twiddle.clear();
    m_trig.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wavetable Lookup table for Fast Fourier Transform.
 ***************************************************************************/
void GFftWavetable::copy_members(const GFftWavetable& wavetable)
{
    // Copy members
    m_factors = wavetable.m_factors;
    m_twiddle = wavetable.m_twiddle;
    m_trig    = wavetable.m_trig;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFftWavetable::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set wavetable
 *
 * @param[in] n Length.
 *
 * Computes the coefficients
 *
 * \f[
 *    \frac{-i 2 \pi j k}{n}
 * \f]
 *
 * of the discrete Fourier transformation
 *
 * \f[
 *    x_j = \sum_{k=0}^{n-1} z_k \exp \left( \frac{-i 2 \pi j k}{n} \right)
 * \f]
 *
 * The method is a C++ implementation of the gsl_fft_complex_wavetable_alloc()
 * function of the GSL library (version 2.2.1), defined in the file
 * fft/c_init.c.
 ***************************************************************************/
void GFftWavetable::set_members(const int& n)
{
    // Initialise trigonometric lookup table
    init_members();

    // Make sure that length is positive
    if (n < 1) {
        std::string msg = "Invalid array length "+gammalib::str(n)+". Array "
                          "needs to be of positive length.";
        throw GException::invalid_argument(G_SET_MEMBERS, msg);
    }

    // Compute factorisation
    set_factors(n);

    // Compute delta theta angle (2 pi / N)
    double d_theta = -gammalib::twopi / double(n);

    // Initialise factorisation product
    int product = 1;

    // Loop over all factors in the factorisation
    for (int i = 0; i < m_factors.size(); ++i) {

        // Set start index of the current factor
        m_twiddle.push_back(m_trig.size());

        // Store current factorisation product
        int product_1 = product;

        // Compute next factorisation product
        product *= m_factors[i];

        // Compute number of coefficients for current factor
        int q = n / product;
/*
std::cout << "i=" << i << " q=" << q << " factor=" << m_factors[i];
std::cout << " product=" << product << " n=" << n << std::endl;
*/
        // Loop over all output coefficients
        for (int j = 1; j < m_factors[i]; ++j) {

            // Loop over all coefficients
            for (int k = 0, m = 0; k < q; ++k) {
            
                // Compute theta value
                m = m + j * product_1;
                m = m % n;
                double theta = d_theta * m;      /*  d_theta*j*k*p_(i-1) */
/*
std::cout << " j=" << j << " k=" << k << " theta=" << theta << std::endl;
*/
                // Compute trigonometric lookup table value
                std::complex<double> value(std::cos(theta), std::sin(theta));

                // Push back value on lookup table
                m_trig.push_back(value);

            } // endfor: loop over all elements for a given factor

        } // endfor

    } // endfor: computed ...

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute FFT factorisation
 *
 * @param[in] n Length of array to be transformed.
 *
 * @exception GException::invalid_argument
 *            Non positive array length specified
 * @expection GException::invalid_value
 *            Product of factorisation factors are not equal to array length
 *
 * Computes the factorisation
 *
 * \f[
 *    n = \prod_i p_i
 * \f]
 *
 * of an array of length @p n.
 *
 * The method is a C++ implementation of the fft_complex_factorize() and the
 * fft_factorize() functions of the GSL library (version 2.2.1), defined in
 * the file fft/factorize.c.
 ***************************************************************************/
void GFftWavetable::set_factors(const int& n)
{
    // Make sure that length is positive
    if (n < 1) {
        std::string msg = "Invalid array length "+gammalib::str(n)+". Array "
                          "needs to be of positive length.";
        throw GException::invalid_argument(G_SET_FACTORS, msg);
    }

    // Set FFT factorisations. Other factors can be added here if their
    // transform modules are implemented.
    std::vector<int> subtransforms;
    subtransforms.push_back(7);
    subtransforms.push_back(6);
    subtransforms.push_back(5);
    subtransforms.push_back(4);
    subtransforms.push_back(3);
    subtransforms.push_back(2);

    // If array has unit length then set the first and only factor to 1 ...
    if (n == 1) {
        m_factors.push_back(1);
    }

    // ... otherwise
    else {

        // Save array length
        int ntest = n;

        // Deal with the implemented factors first
        for (int i = 0; i < subtransforms.size() && ntest > 1; ++i) {
            int factor = subtransforms[i];
            while ((ntest % factor) == 0) {
                ntest /= factor;
                m_factors.push_back(factor);
            }
        }
        /*
        int i = 0;
        while (subtransforms[i] && ntest != 1) {
            int factor = subtransforms[i];
            while ((ntest % factor) == 0) {
                ntest /= factor;
                m_factors.push_back(factor);
            }
            i++;
        }
        */

        // Deal with any other even prime factors (there is only one)
        int factor = 2;
        while ((ntest % factor) == 0 && (ntest != 1)) {
            ntest /= factor;
            m_factors.push_back(factor);
        }

        // Deal with any other odd prime factors
        factor = 3;
        while (ntest != 1) {
            while ((ntest % factor) != 0) {
                factor += 2;
            }
            ntest /= factor;
            m_factors.push_back(factor);
        }

        // Check that the factorization is correct
        int product = 1;
        for (int i = 0; i < m_factors.size(); ++i) {
            product *= m_factors[i];
        }
        if (product != n) {
            std::string msg = "Product "+gammalib::str(product)+" of "
                              "factorisation factors differs from array "
                              "length "+gammalib::str(n)+".";
            throw GException::invalid_value(G_SET_FACTORS, msg);
        }

    } // endelse: array length is larger than 1

    // Return
    return;
}
