/***************************************************************************
 *                  GDerivative.cpp  -  Derivative class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GDerivative.cpp
 * @brief GDerivative class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cstdlib>          // Definition of NULL
#include <cmath>            // for abs()
#include <cfloat>           // for DBL_MAX
#include "GDerivative.hpp"
#include "GMatrix.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_RIDDER             "GDerivative::ridder(double&, double&, double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_VALUE_SIGNAL_BAD_TOLERANCE 1              //!< Signal bad tolerance

/* __ Debug definitions __________________________________________________ */
#define G_VALUE_DEBUG 0                               //!< Debug value method


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GDerivative::GDerivative(void)
{
    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Function constructor
 *
 * @param[in] func Function pointer.
 *
 * Construct class instance from single parameter function. Note that the
 * function is not copied but just a pointer to it is stored. The function
 * should thus not be destroyed before computation of derivatives is
 * finished.
 ***************************************************************************/
GDerivative::GDerivative(GFunction* func)
{
    // Initialise private members
    init_members();

    // Set function
    this->function(func);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dx Derivative.
 ***************************************************************************/
GDerivative::GDerivative(const GDerivative& dx)
{ 
    // Initialise members
    init_members();

    // Copy members
    copy_members(dx);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GDerivative::~GDerivative(void)
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
 * @param[in] dx Derivative.
 ***************************************************************************/
GDerivative& GDerivative::operator= (const GDerivative& dx)
{
    // Execute only if object is not identical
    if (this != &dx) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(dx);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns derivative
 *
 * @param[in] x Function value.
 *
 * Computes derivative using the Ridders' method. The method starts with a
 * step size h=0.02*|x| (with h >= 1e-6) and iteratively decreases the
 * step size until the error becomes less than 1e-6. A maximum number of
 * 5 iterations is done.
 ***************************************************************************/
double GDerivative::value(const double& x)
{
    // Set constants
    const int    max_iter = 5;
    const double min_h    = 1.0e-6;
    const double max_err  = 1.0e-6;

    // Initialise result
    double result = 0.0;

    // Set initial step size
    double h = 0.02*fabs(x);
    if (h < min_h)
        h = min_h;

    // Initialise tolerance
    double err = 1.0;

    // Loop over Ridder's method until we have an acceptable tolerance
    int iter = 0;
    for (; iter < max_iter; ++iter) {

        // Compute derivative using Ridder's method
        result = ridder(x, h, err);

        // Debug option: Show actual results
        #if G_VALUE_DEBUG
        std::cout << "GDerivative::value(";
        std::cout << "iter=" << iter;
        std::cout << ", x=" << x;
        std::cout << ", dy/dx=" << result;
        std::cout << ", h=" << h;
        std::cout << ", err=" << err;
        std::cout << ")" << std::endl;
        #endif

        // If uncertainty is below tolerance then exit now
        if (err < max_err)
            break;

        // ... otherwise reduce the step size
        h /= 5.0;

    } // endfor: tolerance matching loop

    // Compile option: signal if we exceed the tolerance
    #if G_VALUE_SIGNAL_BAD_TOLERANCE
    if (err >= max_err) {
        std::cout << "WARNING: GDerivative::value(";
        std::cout << "iter=" << iter;
        std::cout << ", x=" << x;
        std::cout << ", dy/dx=" << result;
        std::cout << ", h=" << h;
        std::cout << ", err=" << err;
        std::cout << "): error exceeds tolerance of ";
        std::cout << max_err << std::endl;
    }
    #endif
    

    // Return derivative
    return result;
}


/***********************************************************************//**
 * @brief Returns derivative by Ridders' method
 *
 * @param[in] x Function value.
 * @param[in] h Estimated initial step size.
 * @param[out] err Estimated error in the derivative.
 *
 * @exception GException::invalid_argument
 *            Step size of 0 has been specified.
 *
 * Returns the derivative of a function at point x by Ridders's method of
 * polynomial extrapolation. The value h is input as an estimated initial
 * stepsize; it needs not to be small, but rather should be an increment
 * in x over which the function changes substantially.
 *
 * This method has been adopted from Numerical Recipes, 3rd edition, page
 * 321f.
 ***************************************************************************/
double GDerivative::ridder(const double& x, const double& h, double& err)
{
    // Set constants
    const int    ntab = 10;         // Maximum table size
    const double con  = 1.4;        // Stepsize decreases by con at each iteration
    const double con2 = con * con;
    const double big  = DBL_MAX;
    const double safe = 2.0;        // Return when error is safe worse than best so far

    // Initialise result
    double dx = 0.0;
    err       = big;

    // Continue only if we have a valid function
    if (m_func != NULL) {

        // Allocate Neville table
        GMatrix a(ntab, ntab);

        // Throw exception if step size is zero
        if (h <= 0)
            throw GException::invalid_argument(G_RIDDER, "Step size must be positive.");

        // Initialise adaptive step size
        double hh = h;

        // Perform initial function evaluation
        a(0,0) = (m_func->eval(x+hh) - m_func->eval(x-hh))/(2.0*hh);

        // Perform successive computations
        for (int i = 1; i < ntab; ++i) {

            // Reduce step size
            hh /= con;

            // Try new smaller step size
            a(0,i) = (m_func->eval(x+hh) - m_func->eval(x-hh))/(2.0*hh);

            // Compute extrapolations of various orders, requiring no new
            // function evaluations
            double fac = con2;
            for (int j = 1; j <= i; ++j) {

                // Compute extrapolation
                a(j,i) = (  a(j-1,i) * fac - a(j-1,i-1) ) / (fac - 1.0);
                fac   *= con2;

                // If the error is reduced then save the improved answer
                double err1 = std::abs(a(j,i) - a(j-1,i));
                double err2 = std::abs(a(j,i) - a(j-1,i-1));
                double errt = (err1 > err2) ? err1 : err2;
                if (errt <= err) {
                    err = errt;
                    dx  = a(j,i);
                }

            } // endfor: looped over tableau

            // If higher order is worse by a significant factor, the quit early
            if (std::abs( a(i,i) - a(i-1,i-1) ) >= safe*err)
                break;

        } // endfor: looped over computations

    } // endif: we had a valid function

    // Return derivative
    return dx;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GDerivative::init_members(void)
{
    // Initialise members
    m_func = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dx Derivative.
 ***************************************************************************/
void GDerivative::copy_members(const GDerivative& dx)
{
    // Copy attributes
    m_func = dx.m_func;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GDerivative::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
