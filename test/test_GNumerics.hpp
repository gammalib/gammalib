/***************************************************************************
 *                 test_GNumerics.hpp - test numerics modules              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file test_GNumerics.hpp
 * @brief Definition of unit tests for numerics module
 * @author Juergen Knoedlseder
 */

#ifndef TEST_NUMERICS_HPP
#define TEST_NUMERICS_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GammaLib.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @class Gauss
 *
 * @brief Gaussian function
 ***************************************************************************/
class Gauss : public GFunction {
public:
    Gauss(const double& sigma) : m_sigma(sigma) { return; }
    virtual ~Gauss(void) { return; }
    double eval(const double& x) {
        double arg = -0.5*x*x/m_sigma/m_sigma;
        double val = 1.0/std::sqrt(gammalib::twopi)/m_sigma * std::exp(arg);
        return val;
    }
protected:
    double m_sigma;
};


/***********************************************************************//**
 * @class GaussArray
 *
 * @brief Gaussian function array
 ***************************************************************************/
class GaussArray : public GFunctions {
public:
    GaussArray(const GNdarray& sigma) : m_sigma(sigma) { return; }
    virtual ~GaussArray(void) { return; }
    GNdarray eval(const double& x) {
        GNdarray val(m_sigma.shape());
        for (int i = 0; i < m_sigma.size(); ++i) {
            double arg = -0.5*x*x/m_sigma(i)/m_sigma(i);
            val(i)     = 1.0/std::sqrt(gammalib::twopi)/m_sigma(i) * std::exp(arg);
        }
        return val;
    }
protected:
    GNdarray m_sigma;
};


/***********************************************************************//**
 * @class TestGNumerics
 *
 * @brief Test suite for numerical functions
 ***************************************************************************/
class TestGNumerics : public GTestSuite
{
public:
    // Constructors and destructors
    TestGNumerics(void) : GTestSuite(), m_sigma(1.0) { }
    virtual ~TestGNumerics(void) { }

    // Methods
    virtual void           set(void);
    virtual TestGNumerics* clone(void) const;
    virtual std::string    classname(void) const { return "TestGNumerics"; }
    void                   test_math(void);
    void                   test_ndarray(void);
    void                   test_fft(void);
    void                   check_fft1(const GNdarray& array,
                                      const GFft&     fft,
                                      const GNdarray& back);
    void                   check_fft2(const GNdarray& array,
                                      const GFft&     fft,
                                      const GNdarray& back);
    void                   test_function(void);
    void                   test_functions(void);
    void                   test_integral(void);
    void                   test_romberg_integration(void);
    void                   test_adaptive_simpson_integration(void);
    void                   test_gauss_kronrod_integration(void);

private:
    // Private members
    double m_sigma;
};

#endif /* TEST_NUMERICS_HPP */
