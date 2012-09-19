/***************************************************************************
 *                test_GNumerics.hpp  -  test numerics modules             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Jurgen Knodlseder                           *
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

#ifndef TEST_NUMERICS_HPP
#define TEST_NUMERICS_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @class Gauss
 *
 * @brief Gaussian function.
 ***************************************************************************/
class Gauss : public GIntegrand {
public:
    Gauss(const double& sigma) : m_sigma(sigma) { return; }
    virtual ~Gauss(void) { return; }
    double eval(double x) {
        double arg = -0.5*x*x/m_sigma/m_sigma;
        double val = 1.0/sqrt(twopi)/m_sigma * exp(arg);
        return val;
    }
protected:
    double m_sigma;
};

class TestGNumerics : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGNumerics(void) : GTestSuite(), m_sigma(1.0) { return; }
        virtual ~TestGNumerics(void){ return; }

        // Methods
        virtual void set(void);
        void test_integral(void);
        void test_romberg_integration(void);

    // Private attributes
    private:
        double m_sigma;

};

#endif /* TEST_NUMERICS_HPP */
