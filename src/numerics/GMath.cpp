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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
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
 * @brief Compute cosine of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::cosd().
 ***************************************************************************/
double gammalib::cosd(const double& angle)
{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 + 0.5)) % 4;
        switch (i) {
        case 0:
            return 1.0;
        case 1:
            return 0.0;
        case 2:
            return -1.0;
        case 3:
            return 0.0;
        }
    }

    // Return cosine
    return std::cos(angle * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Compute sine of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::sind().
 ***************************************************************************/
double gammalib::sind(const double& angle)
{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 - 0.5))%4;
        switch (i) {
        case 0:
            return 1.0;
        case 1:
            return 0.0;
        case 2:
            return -1.0;
        case 3:
            return 0.0;
        }
    }

    // Return sine
    return std::sin(angle * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Compute tangens of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::tand().
 ***************************************************************************/
double gammalib::tand(const double& angle)
{
    // Check for rounding errors
    double resid = fmod(angle, 360.0);
    if (resid == 0.0 || std::abs(resid) == 180.0) {
        return 0.0;
    }
    else if (resid == 45.0 || resid == 225.0) {
        return 1.0;
    }
    else if (resid == -135.0 || resid == -315.0) {
        return -1.0;
    }

    // Return tangens
    return std::tan(angle * gammalib::deg2rad);
}


/***********************************************************************//**
 * @brief Compute arc cosine in degrees
 *
 * @param[in] value Value
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::acosd().
 ***************************************************************************/
double gammalib::acosd(const double& value)
{
    // Domain tolerance
    const double wcstrig_tol = 1.0e-10;
    
    // Check for rounding errors
    if (value >= 1.0) {
        if (value-1.0 <  wcstrig_tol) {
            return 0.0;
        }
    } 
    else if (value == 0.0) {
        return 90.0;
    }
    else if (value <= -1.0) {
        if (value+1.0 > -wcstrig_tol) {
            return 180.0;
        }
    }

    // Return arc cosine
    return std::acos(value) * gammalib::rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc sine in degrees
 *
 * @param[in] value Value
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::asind().
 ***************************************************************************/
double gammalib::asind(const double& value)
{
    // Domain tolerance
    const double wcstrig_tol = 1.0e-10;
    
    // Check for rounding errors
    if (value <= -1.0) {
        if (value+1.0 > -wcstrig_tol) {
            return -90.0;
        }
    } 
    else if (value == 0.0) {
        return 0.0;
    }
    else if (value >= 1.0) {
        if (value-1.0 <  wcstrig_tol) {
            return 90.0;
        }
    }

    // Return arc sine
    return std::asin(value) * gammalib::rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc tangens in degrees
 *
 * @param[in] value Value
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::atand().
 ***************************************************************************/
double gammalib::atand(const double& value)
{
    // Check for rounding errors
    if (value == -1.0) {
        return -45.0;
    }
    else if (value == 0.0) {
        return 0.0;
    }
    else if (value == 1.0) {
        return 45.0;
    }

    // Return arc sine
    return std::atan(value) * gammalib::rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc tangens in degrees
 *
 * @param[in] y Nominator
 * @param[in] x Denominator
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::atan2d().
 ***************************************************************************/
double gammalib::atan2d(const double& y, const double& x)
{
    // Check for rounding errors
    if (y == 0.0) {
        if (x >= 0.0) {
            return 0.0;
        }
        else if (x < 0.0) {
            return 180.0;
        }
    }
    else if (x == 0.0) {
        if (y > 0.0) {
            return 90.0;
        }
        else if (y < 0.0) {
            return -90.0;
        }
    }

    // Return arc sine
    return std::atan2(y,x) * gammalib::rad2deg;
}


/***********************************************************************//**
 * @brief Compute sine and cosine of angle in degrees
 *
 * @param[in] angle Angle [degrees].
 * @param[out] s Sine of angle.
 * @param[out] c Cosine of angle.
 *
 * This code has been adapted from the WCSLIB function wcstrig.c::sincosd().
 ***************************************************************************/
void gammalib::sincosd(const double& angle, double *s, double *c)

{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 + 0.5)) % 4;
        switch (i) {
        case 0:
            *s = 0.0;
            *c = 1.0;
            return;
        case 1:
            *s = (angle > 0.0) ? 1.0 : -1.0;
            *c = 0.0;
            return;
        case 2:
            *s =  0.0;
            *c = -1.0;
            return;
        case 3:
            *s = (angle > 0.0) ? -1.0 : 1.0;
            *c = 0.0;
            return;
        }
    }
  
    // Compute sine and cosine
    *s = std::sin(angle * gammalib::deg2rad);
    *c = std::cos(angle * gammalib::deg2rad);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Computes logarithm of gamma function
 *
 * @param[in] arg Argument.
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
