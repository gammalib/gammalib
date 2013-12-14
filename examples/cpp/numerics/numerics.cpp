/***************************************************************************
 *            numerics.cpp - Illustrates numerical class usage             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file numerics.cpp
 * @brief Illustrates numerical class usage
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Gaussian function
 *
 * This code illustrates the definition of a Gaussian function.
 ***************************************************************************/
class function : public GFunction {
public:
    function(const double& sigma) : m_a(1.0/(sigma*std::sqrt(gammalib::twopi))),
                                    m_sigma(sigma) {}
    double eval(const double& x) { return m_a*std::exp(-0.5*x*x/(m_sigma*m_sigma)); }
protected:
    double m_a;     //!< Amplitude parameter
    double m_sigma; //!< Width parameter
};


/***********************************************************************//**
 * @brief Create XML file
 *
 * This code illustrates the creation of a XML file.
 ***************************************************************************/
int main(void) {

    // Create instance of function
    function fct(3.0);

    // Create an integral object based on the function
    GIntegral integral(&fct);

    // Set relative integration precision
    integral.eps(1.0e-8);

    // Integrate over the interval [-10,10]
    double result = integral.romb(-15.0, +15.0);

    // Print result (should be basically 1)
    std::cout << "Integral:       " << result << std::endl;

    // Create a derivative object based on the function
    GDerivative derivative(&fct);

    // Print derivatives
    std::cout << "Derivative(0):  " << derivative.value(0.0) << std::endl;
    std::cout << "Derivative(3):  " << derivative.value(3.0) << std::endl;
    std::cout << "Derivative(-3): " << derivative.value(-3.0) << std::endl;

    // Exit
    return 0;
}
