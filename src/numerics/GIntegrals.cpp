/***************************************************************************
 *         GIntegrals.cpp - Integration class for set of functions         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GIntegrals.cpp
 * @brief Integration class for set of functions implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>            // std::abs()
#include <vector>
#include <algorithm>        // std::sort
#include "GIntegrals.hpp"
#include "GException.hpp"
#include "GTools.hpp"
#include "GFunctions.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ROMBERG1           "GIntegrals::romberg(std::vector<double>, int&)"
#define G_ROMBERG2              "GIntegrals::romberg(double&, double&, int&)"
#define G_GAUSS_KRONROD         "GIntegrals::gauss_kronrod(double&, double&)"
#define G_TRAPZD        "GIntegrals::trapzd(double&, double&, int&, GVector)"
#define G_POLINT   "GIntegrals::polint(double*, GVector*, int, int, double, "\
                                                                  "GVector*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GIntegrals::GIntegrals(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector function kernels constructor
 *
 * @param[in] kernels Pointer to function kernels.
 *
 * The vector function kernels constructor assigns the vector function
 * kernels pointer in constructing the object.
 ***************************************************************************/
GIntegrals::GIntegrals(GFunctions* kernels)
{
    // Initialise members
    init_members();

    // Set function kernel
    m_kernels = kernels;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] integrals Integrals.
 ***************************************************************************/
GIntegrals::GIntegrals(const GIntegrals& integrals)
{ 
    // Initialise members
    init_members();

    // Copy members
    copy_members(integrals);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GIntegrals::~GIntegrals(void)
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
 * @param[in] integrals Integrals.
 * @return Integrals.
 ***************************************************************************/
GIntegrals& GIntegrals::operator=(const GIntegrals& integrals)
{
    // Execute only if object is not identical
    if (this != &integrals) {

        // Free members
        free_members();

        // Initialise integral
        init_members();

        // Copy members
        copy_members(integrals);

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
 * @brief Clear integral
 ***************************************************************************/
void GIntegrals::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone integral
 *
 * @return Pointer to deep copy of integrals.
 ***************************************************************************/
GIntegrals* GIntegrals::clone(void) const
{
    return new GIntegrals(*this);
}


/***********************************************************************//**
 * @brief Perform Romberg integration
 *
 * @param[in] bounds Integration boundaries.
 * @param[in] order Integration order (default: 5)
 * @return Vector of integration results.
 *
 * Returns the integral of the integrand, computed over a number of
 * intervals [a0,a1], [a1,a2], ... that are given as an unordered vector
 * by the @p bounds argument.
 *
 * Integration is performed by Romberg's method of order 2*order, where
 *
 *     order=1 is equivalent to the trapezoidal rule,
 *     order=2 is equivalent to Simpson's rule, and
 *     order=3 is equivalent to Boole's rule.
 *
 * The number of iterations is limited by m_max_iter. m_eps specifies the
 * requested fractional accuracy. By default it is set to 1e-6.
 ***************************************************************************/
GVector GIntegrals::romberg(std::vector<double> bounds, const int& order)
{
    // Throw an exception if the instance has no kernels
    if (m_kernels == NULL) {
        std::string msg = "Function kernels not set. Please set function "
                          "kernels before calling the method.";
        throw GException::invalid_value(G_ROMBERG1, msg);
    }

    // Sort integration boundaries in ascending order
    std::sort(bounds.begin(), bounds.end());

    // Initialise result
    GVector result(m_kernels->size());

    // Initialise integration status information
    int calls = 0;

    // Add integral of all intervals
    for (int i = 0; i < bounds.size()-1; ++i) {
        result += romberg(bounds[i], bounds[i+1], order);
        calls  += m_calls;
    }

    // Set integration status information
    m_calls = calls;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Perform Romberg integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] order Integration order (default: 5)
 * @return Vector of integration results.
 *
 * @exception GException::invalid_value
 *            Function kernels not set.
 * @exception GException::invalid_argument
 *            Integration order incompatible with number of iterations.
 *
 * Returns the integral of the integrand from a to b. Integration is
 * performed by Romberg's method of order 2*order, where
 *
 *     order=1 is equivalent to the trapezoidal rule,
 *     order=2 is equivalent to Simpson's rule, and
 *     order=3 is equivalent to Boole's rule.
 *
 * The number of iterations is limited by m_max_iter. m_eps specifies the
 * requested fractional accuracy. By default it is set to 1e-6.
 ***************************************************************************/
GVector GIntegrals::romberg(const double& a, const double& b, const int& order)
{
    // Throw an exception if the instance has no kernels
    if (m_kernels == NULL) {
        std::string msg = "Function kernels not set. Please set function "
                          "kernels before calling the method.";
        throw GException::invalid_value(G_ROMBERG2, msg);
    }

    // Initialise result
    GVector result(m_kernels->size());

    // Initialise integration status information
    m_isvalid    = true;
    m_calls      = 0;
    m_has_abserr = false;
    m_has_relerr = false;

    // Continue only if integration range is valid
    if (b > a) {

        // Initialise variables
        bool    converged = false;
        GVector dss(m_kernels->size());

        // Determine (maximum) number of iterations
        int max_iter = (m_fix_iter > 0) ? m_fix_iter : m_max_iter;

        // Check whether maximum number of iterations is compliant with
        // order
        if (order > max_iter) {
            std::string msg = "Requested integration order "+
                              gammalib::str(order)+" is larger than the "
                              "maximum number of iterations "+
                              gammalib::str(max_iter)+". Either reduce the "
                              "integration order or increase the (maximum) "
                              "number of iterations.";
            throw GException::invalid_argument(G_ROMBERG2, msg);
        }

        // Allocate temporal storage
        double*  h = new double[max_iter+2];
        GVector* s = new GVector[max_iter+2];

        // Initialise step size and initial result
        h[1] = 1.0;
        s[0] = GVector(m_kernels->size());

        // Iterative loop
        for (m_iter = 1; m_iter <= max_iter; ++m_iter) {

            // Integration using Trapezoid rule
            s[m_iter] = trapzd(a, b, m_iter, s[m_iter-1]);

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (is_notanumber(s[m_iter]) || is_infinite(s[m_iter])) {
                m_message = "*** ERROR: GIntegrals::romberg"
                            "(a="+gammalib::str(a)+", b="+gammalib::str(b)+""
                            ", k="+gammalib::str(k)+"): NaN/Inf encountered"
                            " (s["+gammalib::str(m_iter)+"]="
                            ""+s[m_iter]->print()+")";
                std::cout << m_message << std::endl;
                m_isvalid = false;
            }
            #endif

            // Starting from iteration order on, use polynomial interpolation
            if (m_iter >= order) {

                // Compute result using polynom interpolation
                result = polint(&h[m_iter-order], &s[m_iter-order], order, 0.0, &dss);

                // If a fixed number of iterations has been requested and if
                // the final iteration was reached then set the converged
                // flag
                if (m_fix_iter > 0) {
                    if (m_iter == max_iter) {
                        converged = true;
                    }
                }

                // ... otherwise if the requested precision was reached for
                // all functions then set the converged flag
                else {
                    converged = true;
                    for (int i = 0; i < result.size(); ++i) {
                        if (std::abs(dss[i]) > m_eps * std::abs(result[i])) {
                            converged = false;
                            break;
                        }
                    }
                }

                // If the integration has converged then compute the absolute
                // and the relative integration errors. If an intergation result
                // is zero, the relative error is set to zero. The relative
                // error flag is set to true if at least one of the relative
                // integration errors is valid.
                if (converged) {
                    m_has_abserr = true;
                    m_abserr     = abs(dss);
                    m_relerr     = m_abserr;
                    for (int i = 0; i < result.size(); ++i) {
                        double abs_result = std::abs(result[i]);
                        if (abs_result > 0.0) {
                            m_relerr[i] /= abs_result;
                            m_has_relerr = true;
                        }
                        else {
                            m_relerr[i] = 0.0;
                        }
                    }
                    break;
                }

            } // endif: polynomial interpolation performed

            // Reduce step size
            h[m_iter+1] = 0.25 * h[m_iter];

        } // endfor: iterative loop

        // Delete temporal storage
        delete [] s;
        delete [] h;

        // Set status and optionally dump warning
        if (!converged) {
            m_isvalid = false;
            m_message = "Integration uncertainty "+(abs(dss)).print()+
                        " exceeds absolute tolerance of "+
                        (m_eps * abs(result)).print()+
                        " after "+gammalib::str(m_iter)+" iterations. Result "+
                        result.print()+" is inaccurate.";
            if (!m_silent) {
                std::string origin = "GIntegrals::romberg("+
                                     gammalib::str(a)+", "+
                                     gammalib::str(b)+", "+
                                     gammalib::str(order)+")";
                gammalib::warning(origin, m_message);
            }
        }

    } // endif: integration range was valid

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Perform Trapezoidal integration for a set of functions
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] n Number of steps.
 * @param[in] previous_result Result vector from a previous trapezoidal
 *                            integration step.
 * @return Vector of integration results.
 *
 * @exception GException::invalid_value
 *            Function kernels not set.
 *
 * The method performs a trapezoidal integration of a set of functions for
 * the integration interval [a,b].
 *
 * If @p n = 1 the integral for each function \f$j\f$ is computed using
 *
 * \f[
 *    \int_a^b f_j(x) \, dx = \frac{1}{2} (b-a) (f_j(a) + f_j(b))
 * \f]
 *
 * For @p n > 1 the integral is computed using
 *
 * \f[
 *    \int_a^b f_j(x) \, dx = \frac{1}{2} \left[{\tt previous\_result}_j +
 *    \frac{b-a}{2^{n-1}}
 *    \sum_{i=0}^{2^{n-1}-1} f_j \left( a + (0.5+i) \frac{b-a}{2^{n-1}}
 *    \right) \right]
 * \f]
 *
 * where \f${\tt previous\_result}_j\f$ is the integration result for
 * function \f$j\f$ from a previous call to the method with @p n = @p n - 1.
 ***************************************************************************/
GVector GIntegrals::trapzd(const double&  a,
                           const double&  b,
                           const int&     n,
                           const GVector& previous_result)
{
    // Throw an exception if the instance has no kernel
    if (m_kernels == NULL) {
        std::string msg = "Function kernels not set. Please set function "
                          "kernels before calling the method.";
        throw GException::invalid_value(G_TRAPZD, msg);
    }

    // If lower boundary is smaller than upper boundary than use trapezoidal
    // rule
    if (a < b) {

        // Case A: Only a single step is requested
        if (n == 1) {

            // Evaluate integrand at boundaries
            GVector result = 0.5 * (b-a) * (m_kernels->eval(a) +
                                            m_kernels->eval(b));

            // Bookkeeping of function calls
            m_calls += 2;

            // Return result here to avoid an extra vector copy
            return result;

        } // endif: only a single step was requested

        // Case B: More than a single step is requested
        else {

            // Compute step level 2^(n-1)
            int it = 1;
            for (int j = 1; j < n-1; ++j) {
                it <<= 1;
            }

            // Verify that step level is valid
            if (it == 0) {
                m_isvalid = false;
                m_message = "Invalid step level "+gammalib::str(it)+
                            " encountered for"
                            " a="+gammalib::str(a)+
                            ", b="+gammalib::str(b)+
                            ", n="+gammalib::str(n)+
                            ", result="+previous_result.print()+
                            ". Looks like n is too large.";
                if (!m_silent) {
                    gammalib::warning(G_TRAPZD, m_message);
                }
            }

            // Set step size
            double tnm = double(it);
            double del = (b-a)/tnm;

            // Verify that step is >0
            if (del == 0) {
                m_isvalid = false;
                m_message = "Invalid step size "+gammalib::str(del)+
                            " encountered for"
                            " a="+gammalib::str(a)+
                            ", b="+gammalib::str(b)+
                            ", n="+gammalib::str(n)+
                            ", result="+previous_result.print()+
                            ". Step is too small to make sense.";
                if (!m_silent) {
                    gammalib::warning(G_TRAPZD, m_message);
                }
            }

            // Initialise result
            GVector result(previous_result.size());
            
            // Initialise argument
            double x = a + 0.5 * del;

            // Sum up values
            for (int j = 0; j < it; ++j, x += del) {

                // Evaluate and add integrand
                result += m_kernels->eval(x);

                // Bookkeeping of function calls
                m_calls++;

            } // endfor: looped over steps

            // Compute result benefiting from previous result
            result *= del;
            result += previous_result;
            result *= 0.5;

            // Return result here to avoid an extra vector copy
            return result;

        } // endelse: Case B

    } // endif: trapezoidal rule was applied

    // ... otherwise handle case a >= b
    else {

        // Set empty vector
        GVector result(previous_result.size());

        // Return result here to avoid an extra vector copy
        return result;
    }

}


/***********************************************************************//**
 * @brief Print integral information
 *
 * @param[in] chatter Chattiness.
 * @return String containing integral information.
 ***************************************************************************/
std::string GIntegrals::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GIntegrals ===");

        // Append information
        result.append("\n"+gammalib::parformat("Relative precision"));
        result.append(gammalib::str(eps()));
        if (m_has_abserr) {
            result.append("\n"+gammalib::parformat("Absolute errors"));
            result.append(m_abserr.print());
        }
        if (m_has_relerr) {
            result.append("\n"+gammalib::parformat("Relative errors"));
            result.append(m_relerr.print());
        }
        result.append("\n"+gammalib::parformat("Function calls"));
        result.append(gammalib::str(calls()));
        result.append("\n"+gammalib::parformat("Iterations"));
        result.append(gammalib::str(iter()));
        if (m_fix_iter > 0) {
            result.append(" (fixed: ");
            result.append(gammalib::str(fixed_iter()));
            result.append(")");
        }
        else {
            result.append(" (maximum: ");
            result.append(gammalib::str(max_iter()));
            result.append(")");
        }

        // Append status information
        result.append("\n"+gammalib::parformat("Status"));
        if (is_valid()) {
            result.append("Result accurate.");
        }
        else {
            result.append(message());
        }
        if (silent()) {
            result.append("\n"+gammalib::parformat("Warnings")+"suppressed");
        }
        else {
            result.append("\n"+gammalib::parformat("Warnings"));
            result.append("in standard output");
        }

    } // endif: chatter was not silent

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
void GIntegrals::init_members(void)
{
    // Initialise members
    m_kernels   = NULL;
    m_eps       = 1.0e-6;
    m_max_iter  = 20;
    m_fix_iter  = 0;
    m_message.clear();
    m_silent    = false;

    // Initialise results
    m_iter       = 0;
    m_calls      = 0;
    m_isvalid    = true;
    m_has_abserr = false;
    m_has_relerr = false;
    m_abserr.clear();
    m_relerr.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] integral Integral.
 ***************************************************************************/
void GIntegrals::copy_members(const GIntegrals& integral)
{
    // Copy attributes
    m_kernels  = integral.m_kernels;
    m_eps      = integral.m_eps;
    m_max_iter = integral.m_max_iter;
    m_fix_iter = integral.m_fix_iter;
    m_message  = integral.m_message;
    m_silent   = integral.m_silent;

    // Copy results
    m_iter       = integral.m_iter;
    m_calls      = integral.m_calls;
    m_isvalid    = integral.m_isvalid;
    m_has_abserr = integral.m_has_abserr;
    m_has_relerr = integral.m_has_relerr;
    m_abserr     = integral.m_abserr;
    m_relerr     = integral.m_relerr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GIntegrals::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Perform Polynomial interpolation
 *
 * @param[in] xa Pointer to array of X values.
 * @param[in] ya Pointer to GVector of Y values.
 * @param[in] order Number of elements in arrays.
 * @param[in] x X value for which interpolations should be performed.
 * @param[out] dy Pointer to vector of error estimates for interpolated values.
 * @return Vector of interpolated values.
 *
 * Given arrays @p xa[1,..,n] and @p ya[1,..,n], and given a value @p x, this
 * method returns a value y, and an error estimate @p dy. If P(x) is the
 * polynomial of degree n-1, then the returned value y=P(x).
 ***************************************************************************/
GVector GIntegrals::polint(const double*  xa,
                           const GVector* ya,
                           const int&     order,
                           const double&  x,
                           GVector*       dy)
{
    // Get vector dimension
    int nsize = ya->size();

    // Initialise result
    GVector y(nsize);

    // Allocate temporary memory
    double* c = new double[order];
    double* d = new double[order];

    // Find index nclosest of the closest table entry
    int    nclosest = 0;
    double dif      = std::abs(x-xa[1]);
    for (int i = 1; i < order; ++i) {
        double dift = std::abs(x-xa[i+1]);
        if (dift < dif) {
            nclosest  = i;
            dif       = dift;
        }
    }

    // Loop over vector dimension
    for (int index = 0; index < nsize; ++index) {

        // Set closest table entry
        int ns = nclosest;

        // Initialise table values
        for (int i = 0; i < order; ++i) {
            double value = ya[i+1][index];
            c[i]         = value;
            d[i]         = value;
        }

        // Get initial approximation to y
        y[index] = ya[ns+1][index];
        ns--;

        // Loop over each column of the tableau
        for (int m = 1; m < order; ++m) {

            // Update current c's and d's
            for (int i = 0; i < order-m; ++i) {
                double ho  = xa[i+1]   - x;
                double hp  = xa[i+m+1] - x;
                double w   = c[i+1] - d[i];
                double den = ho - hp;
                if (den == 0.0) {
                    m_isvalid = false;
                    m_message = "Invalid step size "+gammalib::str(den)+
                                " encountered. Two values in xa array are"
                                " identical.";
                    if (!m_silent) {
                        gammalib::warning(G_POLINT, m_message);
                    }
                }
                den  = w/den;
                d[i] = hp*den;
                c[i] = ho*den;
            }

            // Compute y correction
            double delta_y = (2*(ns+1) < (order-m)) ? c[ns+1] : d[ns--];
        
            // Update y
            y[index] += delta_y;

            // Store result
            (*dy)[index] = delta_y;

        } // endfor: looped over columns of tableau

    } // endfor: looped over vector dimension

    // Free temporary memory
    delete [] c;
    delete [] d;

    // Return
    return y;
}
