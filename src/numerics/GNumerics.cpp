/***************************************************************************
 *                    GNumerics.cpp  -  Numerical functions                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GNumerics.cpp
 * @brief Numerical function implementations
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GNumerics.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                               Functions                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Computes logarithm of gamma function
 *
 * @param[in] x Argument.
 ***************************************************************************/
double gammln(const double& x) {

    // Define static constants
    static const double cof[6] = { 76.18009172947146,
                                  -86.50532032941677,  
                                   24.01409824083091,
                                   -1.231739572450155,
                                    0.1208650973866179e-2,
                                   -0.5395239384953e-5};

    // Evaluate logarithm of gamma function
    double a = x;
    double b = x;
    double c = a + 5.5;
    c -= (a + 0.5) * std::log(c);
    double d = 1.000000000190015;
    for (int i = 0; i < 6; ++i) {
        d += cof[i]/++b;
    }
    double result = log(2.5066282746310005 * d/a) - c;
	
    // Return result
    return result;
}
