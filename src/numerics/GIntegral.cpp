/***************************************************************************
 *                  GIntegral.cpp  -  Integration class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GIntegral.cpp
 * @brief Integration class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>            // For std::abs()
#include <vector>
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GIntegral::GIntegral(void)
{
    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Integrand constructor
 *
 * @param[in] integrand Pointer to integrand that should be assigned.
 *
 * The Integrand constructor assigns the integrand pointer in constructing
 * the object.
 ***************************************************************************/
GIntegral::GIntegral(GIntegrand* integrand)
{
    // Initialise private members
    init_members();

    // Set integrand
    m_integrand = integrand;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] integral Object from which the instance should be built.
 ***************************************************************************/
GIntegral::GIntegral(const GIntegral& integral)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(integral);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GIntegral::~GIntegral(void)
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
 * @param[in] integral Object to be assigned.
 ***************************************************************************/
GIntegral& GIntegral::operator= (const GIntegral& integral)
{
    // Execute only if object is not identical
    if (this != &integral) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(integral);

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
 * @brief Perform Romberg integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] k Integration order (default: k=5)
 *
 * Returns the integral of the integrand from a to b. Integration is
 * performed by Romberg's method of order 2K, where e.g. K=2 in Simpson's
 * rule.
 * The number of iterations is limited by m_max_iter. m_eps specifies the
 * requested fractional accuracy. By default it is set to 1e-6.
 *
 * @todo Check that k is smaller than m_max_iter
 ***************************************************************************/
double GIntegral::romb(double a, double b, int k)
{
    // Initialise result
    bool   converged = false;
    double result    = 0.0;
    double ss        = 0.0;
    double dss       = 0.0;

    // Allocate temporal storage
    double* s = new double[m_max_iter+2];
    double* h = new double[m_max_iter+2];

    // Initialise step size
    h[1] = 1.0;
    s[0] = 0.0;

    // Loop
    for (m_iter = 1; m_iter <= m_max_iter; ++m_iter) {

        // Integration using Trapezoid rule
        s[m_iter] = trapzd(a, b, m_iter, s[m_iter-1]);

        // Starting from iteration k on, use polynomial interpolation
        if (m_iter >= k) {
            ss = polint(&h[m_iter-k], &s[m_iter-k], k, 0.0, &dss);
            if (std::abs(dss) <= m_eps * std::abs(ss)) {
                converged = true;
                result    = ss;
                break;
            }
        }

        // Reduce step size
        h[m_iter+1]= 0.25 * h[m_iter];

    }

    // Free temporal storage
    delete [] s;
    delete [] h;

    // Dump warning
    if (!converged) {
        std::cout << "*** ERROR: GIntegral::romb: ";
        std::cout << "Integration did not converge (result=" << ss << ")";
        std::cout << std::endl;
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GIntegral::init_members(void)
{
    // Initialise members
    m_integrand = NULL;
    m_eps       = 1.0e-6;
    m_max_iter  = 20;
    m_iter      = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] integral Object from which members are to be copied.
 ***************************************************************************/
void GIntegral::copy_members(const GIntegral& integral)
{
    // Copy attributes
    m_integrand = integral.m_integrand;
    m_eps       = integral.m_eps;
    m_max_iter  = integral.m_max_iter;
    m_iter      = integral.m_iter;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GIntegral::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform Polynomial interpolation
 *
 * @param[in] xa Pointer to array of X values.
 * @param[in] ya Pointer to array of Y values.
 * @param[in] n Number of elements in arrays.
 * @param[in] x X value at which interpolations should be performed.
 * @param[out] dy Error estimate for interpolated values.
 *
 * Given arrays xa[1,..,n] and ya[1,..,n], and given a value x, this
 * method returns a value y, and an error estimate dy. If P(x) is the
 * polynomial of degree n-1, then the returned value y=P(x).
 *
 * @todo Implement exceptions instead of screen dump.
 * @todo Use std::vector for xa and ya and start at 0
 ***************************************************************************/
double GIntegral::polint(double* xa, double* ya, int n, double x, double *dy)
{
    // Initialise result
    double y = 0.0;

    // Allocate temporary memory
    std::vector<double> c(n, 0.0);
    std::vector<double> d(n, 0.0);

    // Compute initial distance to first node
    double dif = std::abs(x-xa[1]);

    // Find index ns of the closest table entry
    int ns = 0;
    for (int i = 0; i < n; ++i) {
        double dift = std::abs(x-xa[i+1]);
        if (dift < dif) {
            ns  = i;
            dif = dift;
        }
        c[i] = ya[i+1];
        d[i] = ya[i+1];
    }

    // Get initial approximation to y
    y = ya[ns+1];
    ns--;

    // Loop over each column of the tableau
    for (int m = 1; m < n; ++m) {

        // Update current c's and d's
        for (int i = 0; i < n-m; ++i) {
            double ho  = xa[i+1]   - x;
            double hp  = xa[i+m+1] - x;
            double w   = c[i+1] - d[i];
            double den = ho - hp;
            if (den == 0.0) {
                std::cout << "*** ERROR: GIntegral::polint: ";
                std::cout << "This error can only occur if two input xa's are identical.";
                std::cout << std::endl;
            }
            den  = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }

        // Compute y correction
        *dy = (2*(ns+1) < (n-m)) ? c[ns+1] : d[ns--];

        // Update y
        y += *dy;

    } // endfor: looped over columns of tableau

    // Return
    return y;
}


/***********************************************************************//**
 * @brief Perform Trapezoidal integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] n Number of steps.
 * @param[in] result Result from a previous trapezoidal integration step.
 *
 * The original Numerical Recipes function had result declared as a static
 * variable, yet this led to some untrackable integration problems. For this
 * reason, previous results are now passed using an argument.
 * Result initialisation is done if n=1.
 ***************************************************************************/
double GIntegral::trapzd(double a, double b, int n, double result)
{
    // Case A: Only a single step is requested
    if (n == 1)
        result = 0.5*(b-a)*(m_integrand->eval(a) + m_integrand->eval(b));

    // Case B: More than a single step is requested
    else {
        //
        int it = 1;
        for (int j = 1; j < n-1; ++j)
            it <<= 1;

        // Set step size
        double tnm = double(it);
        double del = (b-a)/tnm;

        // Sum up values
        double x   = a + 0.5*del;
        double sum = 0.0;
        for (int j = 0; j < it; ++j, x+=del)
            sum += m_integrand->eval(x);

        // Set result
        result = 0.5*(result + (b-a)*sum/tnm);
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] integral Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GIntegral& integral)
{
    // Put object in stream
    os << "=== GIntegral ===" << std::endl;
    os << " Integration precision .....: " << integral.m_eps << std::endl;
    os << " Max. number of iterations .: " << integral.m_max_iter;

    // Return output stream
    return os;
}
