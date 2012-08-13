/***************************************************************************
 *                test_GNumerics.cpp  -  test numerics modules             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <ostream>
#include <stdexcept>
#include <stdlib.h>
#include "test_GNumerics.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***********************************************************************//**
 * @brief Test model parameter handling.
 ***************************************************************************/
void test_GIntegral(GTestSuite& testsuite)
{
    // Set sigma
    double sigma = 2.5;

    // Test integral and integrand allocation
    testsuite.test_try("Test integral and integrand allocation");
    try {
        Gauss     integrand(sigma);
        GIntegral integral(&integrand);
        
        testsuite.success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test Romberg integration
    testsuite.test_try("Test Romberg integration");
    try {
        Gauss     integrand(sigma);
        GIntegral integral(&integrand);
        double    result = integral.romb(-10.0*sigma, 10.0*sigma);
        if (fabs(result-1.0) > 1.0e-6) {
            throw testsuite.exception_failure("Gaussian integral is not 1.0 (integral="+str(result)+")");
        }
        result = integral.romb(-sigma, sigma);
        if (fabs(result-0.68268948130801355) > 1.0e-6) {
            throw testsuite.exception_failure("Gaussian integral is not 0.682689 (difference="+str((result-0.68268948130801355))+")");
        }
        result = integral.romb(0.0, sigma);
        if (fabs(result-0.3413447460687748) > 1.0e-6) {
            throw testsuite.exception_failure("Gaussian integral is not 0.341345 (difference="+str((result-0.3413447460687748))+")");
        }

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    //Create a test suites container
    GTestSuites testsuites("GNumerics testing");

    // Create a test suite
    GTestSuite testsuite("Test GIntegral");

    // Append it
    testsuites.append(testsuite);

    // Execute Healpix tests
    test_GIntegral(testsuite);

    //End of the tests
    testsuite.end_test();

    //save xml report
    testsuites.save("reports/GNumerics.xml");

    // Return
    return 0;
}
