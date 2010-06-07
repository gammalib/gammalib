/***************************************************************************
 *                test_GNumerics.cpp  -  test numerics modules             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
void test_GIntegral(void)
{
    // Write header
    std::cout << "Test GIntegral: ";

    // Set sigma
    double sigma = 2.5;

    // Test integral and integrand allocation
    try {
        Gauss     integrand(sigma);
        GIntegral integral(&integrand);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to allocated integral and integrand."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test Romberg integration
    try {
        Gauss     integrand(sigma);
        GIntegral integral(&integrand);
        double    result = integral.romb(-10.0*sigma, 10.0*sigma);
        if (fabs(result-1.0) > 1.0e-6) {
            std::cout << std::endl
                      << "TEST ERROR: Gaussian integral is not 1.0 "
                      << " (integral=" << result << ")"
                      << std::endl;
            throw;
        }
        result = integral.romb(-sigma, sigma);
        if (fabs(result-0.68268948130801355) > 1.0e-6) {
            std::cout << std::endl
                      << "TEST ERROR: Gaussian integral is not 0.682689 "
                      << " (difference=" << (result-0.68268948130801355) << ")"
                      << std::endl;
            throw;
        }
        result = integral.romb(0.0, sigma);
        if (fabs(result-0.3413447460687748) > 1.0e-6) {
            std::cout << std::endl
                      << "TEST ERROR: Gaussian integral is not 0.341345 "
                      << " (difference=" << (result-0.3413447460687748) << ")"
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to allocated integral and integrand."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "********************" << std::endl;
    std::cout << "* Numerics testing *" << std::endl;
    std::cout << "********************" << std::endl;

    // Execute Healpix tests
    test_GIntegral();

    // Return
    return 0;
}
