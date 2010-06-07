/***************************************************************************
 *                  GIntegral.cpp  -  Integration class                    *
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
/**
 * @file GIntegral.cpp
 * @brief GIntegral class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <math.h>
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
    double* s = new double[m_max_iter];
    double* h = new double[m_max_iter+1];

    // Initialise step size
    h[0] = 1.0;

    // Loop
    for (int iter = 1; iter <= m_max_iter; ++iter) {

        // Integration using Trapezoid rule (j=0 initialises method)
        s[iter-1] = trapzd(a, b, iter);

        // Starting from iteration k on, use polynomial interpolation
        if (iter >= k) {
            polint(&h[iter-k], &s[iter-k], k, 0.0, &ss, &dss);
            if (fabs(dss) <= m_eps * fabs(ss)) {
                converged = true;
                result    = ss;
                break;
            }
        }

        // Reduce step size
        h[iter]= 0.25 * h[iter-1];
    }

    // Free temporal storage
    delete [] s;
    delete [] h;

    // Dump warning
    if (!converged) {
        std::cout << "GIntegral::romb: Integration did not converge (result="
                  << ss << ")" << std::endl;
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
 * @param[out] y Interpolated Y value.
 * @param[out] dy Error estimate for interpolated values.
 *
 * Given arrays xa[0,..,n-1] and ya[0,..,n-1], and given a value x, this
 * method returns a value y, and an error estimate dy. If P(x) is the
 * polynomial of degree n-1, then the returned value y=P(x).
 *
 * @todo Implement exceptions instead of screen dump.
 ***************************************************************************/
void GIntegral::polint(double* xa, double* ya, int n, double x, double *y,
                       double *dy)
{
    // Allocate temporary memory
    double* c = new double[n];
    double* d = new double[n];

    // Compute initial distance to first node
    double dif = fabs(x-xa[0]);

    // Find index ns of the closest table entry
    int ns = 0;
    for (int i = 0; i < n; ++i) {
        double dift = fabs(x-xa[i]);
        if (dift < dif) {
            ns  = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    // Get initial approximation to y
    *y = ya[ns--];

    // Loop over each column of the tableau
    for (int m = 1; m < n; ++m) {

        // Update current c's and d's
        for (int i = 0; i < n-m; ++i) {
            double ho  = xa[i]   - x;
            double hp  = xa[i+m] - x;
            double w   = c[i+1] - d[i];
            double den = ho - hp;
            if (den == 0.0) {
                std::cout << "GIntegral::polint: an error occured. "
                          << "This error can only occur if two input xa's are identical."
                          << std::endl;
            }
            den  = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }

        // Compute y correction
        *dy = (2*(ns+1) < (n-m)) ? c[ns+1] : d[ns--];

        // Update y
        *y += *dy;

    } // endfor: looped over columns of tableau

    // Delete temporary memory
    delete [] d;
    delete [] c;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform Trapezoidal integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] n Number of steps.
 *
 * Note that this method keeps track about the last result using a static
 * variable. Initialisation is done by setting n=1.
 ***************************************************************************/
double GIntegral::trapzd(double a, double b, int n)
{
    // Declare result as static
    static double result;

    // Case A: Only a single step is requested
    if (n == 1)
        result = 0.5*(b-a)*(m_integrand->eval(a) + m_integrand->eval(b));

    // Case B: Mone than a single step is requested
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
        for (int j = 1; j <= it; ++j, x+=del)
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
