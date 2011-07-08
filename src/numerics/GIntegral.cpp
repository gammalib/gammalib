/***************************************************************************
 *                  GIntegral.cpp  -  Integration class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GIntegral.cpp
 * @brief Integration class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>            // For std::abs()
#include <vector>
#include "GIntegral.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_CHECK_NAN                                     //!< Check for NaN


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GIntegral::GIntegral(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Integrand constructor
 *
 * @param[in] integrand Pointer to integrand.
 *
 * The Integrand constructor assigns the integrand pointer in constructing
 * the object.
 ***************************************************************************/
GIntegral::GIntegral(GIntegrand* integrand)
{
    // Initialise members
    init_members();

    // Set integrand
    m_integrand = integrand;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] integral Integral.
 ***************************************************************************/
GIntegral::GIntegral(const GIntegral& integral)
{ 
    // Initialise members
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
 * @param[in] integral Integral.
 ***************************************************************************/
GIntegral& GIntegral::operator= (const GIntegral& integral)
{
    // Execute only if object is not identical
    if (this != &integral) {

        // Free members
        free_members();

        // Initialise integral
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
 * @brief Clear instance
 ***************************************************************************/
void GIntegral::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GIntegral* GIntegral::clone(void) const
{
    return new GIntegral(*this);
}


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
 * @todo Make use of std::vector class
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

        // Debug: Check for NaN
        #if defined(G_CHECK_NAN)
        if (std::isnan(s[m_iter]) || std::isinf(s[m_iter])) {
            std::cout << "*** ERROR: GIntegral::romb";
            std::cout << "(a=" << a << ", b=" << b << ", k=" << k << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (s[" << m_iter << "]=" << s[m_iter] << ")";
            std::cout << std::endl;
        }
        #endif

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
    if (!m_silent) {
        if (!converged) {
            std::cout << "*** WARNING: GIntegral::romb: ";
            std::cout << "Integration did not converge ";
            std::cout << "(iter=" << m_iter;
            std::cout << ", result=" << ss;
            std::cout << ", d=" << std::abs(dss);
            std::cout << " > " << m_eps * std::abs(ss) << ")";
            std::cout << std::endl;
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print integral information
 ***************************************************************************/
std::string GIntegral::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GIntegral ===");

    // Append information
    result.append("\n"+parformat("Relative precision")+str(eps()));
    result.append("\n"+parformat("Max. number of iterations")+str(max_iter()));
    if (silent())
        result.append("\n"+parformat("Warnings")+"suppressed");
    else
        result.append("\n"+parformat("Warnings")+"in standard output");

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
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
    m_silent    = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] integral Integral.
 ***************************************************************************/
void GIntegral::copy_members(const GIntegral& integral)
{
    // Copy attributes
    m_integrand = integral.m_integrand;
    m_eps       = integral.m_eps;
    m_max_iter  = integral.m_max_iter;
    m_iter      = integral.m_iter;
    m_silent    = integral.m_silent;

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
    // Handle case of identical boundaries
    if (a == b) {
        result = 0.0;
    }
    
    // ... otherwise use trapeziodal rule
    else {
    
        // Case A: Only a single step is requested
        if (n == 1) {
        
            // Evaluate integrand at boundaries
            double y_a = m_integrand->eval(a);
            double y_b = m_integrand->eval(b);
            
            // Compute result
            result = 0.5*(b-a)*(m_integrand->eval(a) + m_integrand->eval(b));
            
        } // endif: only a single step was requested

        // Case B: More than a single step is requested
        else {
            // Compute step level 2^(n-1)
            int it = 1;
            for (int j = 1; j < n-1; ++j)
                it <<= 1;

            // Verify that step level is valid
            if (it == 0) {
                std::cout << "*** ERROR: GIntegral::trapzd(";
                std::cout << "a=" << a << ", b=" << b << ", n=" << n;
                std::cout << ", result=" << result << "): it=" << it << ". ";
                std::cout << "Looks like an overflow?";
                std::cout << std::endl;
            }

            // Set step size
            double tnm = double(it);
            double del = (b-a)/tnm;

            // Verify that step is >0
            if (it == 0) {
                std::cout << "*** ERROR: GIntegral::trapzd(";
                std::cout << "a=" << a << ", b=" << b << ", n=" << n;
                std::cout << ", result=" << result << "): del=" << del << ". ";
                std::cout << "Step is too small to make sense.";
                std::cout << std::endl;
            }

            // Sum up values
            double x   = a + 0.5*del;
            double sum = 0.0;
            for (int j = 0; j < it; ++j, x+=del) {
                
                // Evaluate integrand
                double y = m_integrand->eval(x);

                // Add integrand
                sum += y;
                
            } // endfor: looped over steps

            // Set result
            result = 0.5*(result + (b-a)*sum/tnm);
        }
        
    } // endelse: trapeziodal rule was applied

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] integral Integral.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GIntegral& integral)
{
    // Write integral in output stream
    os << integral.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] integral Integral.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GIntegral& integral)
{
    // Write integral into logger
    log << integral.print();

    // Return logger
    return log;
}
