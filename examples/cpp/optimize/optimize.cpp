/***************************************************************************
 *            optimize.cpp - Illustrates optimizer class usage             *
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
 * @file optimize.cpp
 * @brief Illustrates optimizer class usage
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @brief Parabola function
 *
 * This code illustrates the definition of a parabola function.
 ***************************************************************************/
class parabola : public GOptimizerFunction {
public:
    parabola(void) : m_value(0), m_gradient(1), m_curvature(1,1) {}
    void           eval(const GOptimizerPars& pars);
    double         value(void) { return m_value; }
    GVector*       gradient(void) { return &m_gradient; }
    GMatrixSparse* curvature(void) { return &m_curvature; }
protected:
    double        m_value;     //!< Function value
    GVector       m_gradient;  //!< Function gradient vector
    GMatrixSparse m_curvature; //!< Curvature matrix
};
void parabola::eval(const GOptimizerPars& pars)
{
    const double a =  2.0;
    const double b = -4.0;
    const double c =  2.0;
    double x       = pars[0]->value();
    m_value        = a*x*x + b*x + c;
    m_gradient[0]  = 2.0*a*x + b;
    m_curvature(0,0)   = m_gradient[0] * m_gradient[0];
};


/***********************************************************************//**
 * @brief Optimize function
 *
 * This code illustrates the optimization of a function.
 ***************************************************************************/
int main(void) {

    // Allocate logger
    GLog log;
    log.cout(true);
    
    // Allocate optimizer
    //GOptimizerLM opt(log);
    GOptimizerLM opt;

    // Allocate function
    parabola fct;

    // Allocate parameters and set initial value
    GOptimizerPars pars(1);
    pars[0]->value(1.5);

    // Optimize parameters
    opt.optimize(fct, pars);

    // Print derivatives
    std::cout << "Function value .....: " << fct.value() << std::endl;
    std::cout << "Parameter value ....: " << pars[0]->value() << std::endl;

    // Exit
    return 0;
}
