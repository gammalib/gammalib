/***************************************************************************
 *                     GMath.cpp - Mathematical functions                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GMath.cpp
 * @brief Mathematical function implementations
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GMath.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_ASIN_FOR_ACOS             //!< Use asin for acos computations

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                               Functions                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Computes acos by avoiding NaN due to rounding errors
 *
 * @param[in] arg Argument.
 * @return Arc cosine of argument.
 *
 * Returns the arc cosine by restricting the argument to [-1,1].
 *
 * If the compile option G_USE_ASIN_FOR_ACOS is defined, the function will
 * compute
 *
 * \f[
 *    acos(x) = \frac{\pi}{2} - asin(x)
 * \f]
 *
 * which happens to be faster on most systems.
 ***************************************************************************/
double gammalib::acos(const double& arg)
{
    // Allocate result
    double arccos;

    // Compute acos
    if (arg >= 1) {
        arccos = 0.0;
    }
    else if (arg <= -1.0) {
        arccos = gammalib::pi;
    }
    else {
        #if defined(G_USE_ASIN_FOR_ACOS)
        arccos = gammalib::pihalf - std::asin(arg);
        #else
        arccos = std::acos(arg);
        #endif
    }

    // Return result
    return arccos;
}


/***********************************************************************//**
 * @brief Computes logarithm of gamma function
 *
 * @param[in] x Argument.
 * @return Logarithm of gamma function.
 ***************************************************************************/
double gammalib::gammln(const double& arg) {

    // Define static constants
    static const double cof[6] = { 76.18009172947146,
                                  -86.50532032941677,  
                                   24.01409824083091,
                                   -1.231739572450155,
                                    0.1208650973866179e-2,
                                   -0.5395239384953e-5};

    // Evaluate logarithm of gamma function
    double a = arg;
    double b = arg;
    double c = a + 5.5;
    c -= (a + 0.5) * std::log(c);
    double d = 1.000000000190015;
    for (int i = 0; i < 6; ++i) {
        d += cof[i]/++b;
    }
    double result = std::log(2.5066282746310005 * d/a) - c;
	
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns the remainder of the division \a v1/v2.
 *
 * @param[in] v1 Argument 1.
 * @param[in] v2 Argument 2.
 *
 * Returns the remainder of the division \a v1/v2.
 * The result is non-negative.
 * \a v1 can be positive or negative; \a v2 must be positive.
 ***************************************************************************/
double gammalib::modulo(const double& v1, const double& v2)
{
    // Return
    return (v1 >= 0) ? ((v1 < v2) ? v1 : std::fmod(v1,v2)) : (std::fmod(v1,v2)+v2);
}

