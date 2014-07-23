/***************************************************************************
 *                   GIntegral.cpp - Integration class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>            // For std::abs()
#include <vector>
#include "GIntegral.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ROMB                      "GIntegral::romb(double&, double&, int&)"
#define G_TRAPZD          "GIntegral::trapzd(double&, double&, int&, double)"
#define G_POLINT  "GIntegral::polint(double*, double*, int, double, double*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


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
 * @brief Function kernel constructor
 *
 * @param[in] kernel Pointer to function kernel.
 *
 * The function kernel constructor assigns the function kernel pointer in
 * constructing the object.
 ***************************************************************************/
GIntegral::GIntegral(GFunction* kernel)
{
    // Initialise members
    init_members();

    // Set function kernel
    m_kernel = kernel;

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
 * @return Integral.
 ***************************************************************************/
GIntegral& GIntegral::operator=(const GIntegral& integral)
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
 * @brief Clear integral
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
 * @brief Clone integral
 *
 * @return Pointer to deep copy of integral.
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
 * performed by Romberg's method of order 2k, where
 * k=1 is equivalent to the trapezoidal rule,
 * k=2 is equivalent to Simpson's rule, and
 * k=3 is equivalent to Boole's rule.
 * The number of iterations is limited by m_max_iter. m_eps specifies the
 * requested fractional accuracy. By default it is set to 1e-6.
 *
 * @todo Check that k is smaller than m_max_iter
 ***************************************************************************/
double GIntegral::romb(const double& a, const double& b, const int& k)
{
    // Initialise result and status
    double result = 0.0;

    // Initialise integration status information
    m_isvalid = true;
    m_calls   = 0;
    
    // Continue only if integration range is valid
    if (b > a) {

        // Initialise variables
        bool   converged = false;
        double ss        = 0.0;
        double dss       = 0.0;

        // Allocate temporal storage
        double* s = new double[m_max_iter+2];
        double* h = new double[m_max_iter+2];

        // Initialise step size
        h[1] = 1.0;
        s[0] = 0.0;

        // Iterative loop
        for (m_iter = 1; m_iter <= m_max_iter; ++m_iter) {

            // Integration using Trapezoid rule
            s[m_iter] = trapzd(a, b, m_iter, s[m_iter-1]);

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (is_notanumber(s[m_iter]) || is_infinite(s[m_iter])) {
                m_message = "*** ERROR: GIntegral::romb"
                            "(a="+gammalib::str(a)+", b="+gammalib::str(b)+""
                            ", k="+gammalib::str(k)+"): NaN/Inf encountered"
                            " (s["+gammalib::str(m_iter)+"]="
                            ""+gammalib::str(s[m_iter])+")";
                std::cout << m_message << std::endl;
                m_isvalid = false;
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

        } // endfor: iterative loop

        // Free temporal storage
        delete [] s;
        delete [] h;

        // Set status and optionally dump warning
        if (!converged) {
            m_isvalid = false;
            m_message = "Integration uncertainty "+
                        gammalib::str(std::abs(dss))+
                        " exceeds absolute tolerance of "+
                        gammalib::str(m_eps * std::abs(ss))+
                        " after "+gammalib::str(m_iter)+
                        " iterations. Result "+
                        gammalib::str(ss)+
                        " is inaccurate.";
            if (!m_silent) {
                std::string origin = "GIntegral::romb("+
                                     gammalib::str(a)+", "+
                                     gammalib::str(b)+", "+
                                     gammalib::str(k)+")";
                gammalib::warning(origin, m_message);
            }
        }
    
    } // endif: integration range was valid

    // Return result
    return result;
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
double GIntegral::trapzd(const double& a, const double& b, const int& n,
                         double result)
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
            double y_a = m_kernel->eval(a);
            double y_b = m_kernel->eval(b);
            m_calls += 2;
            
            // Compute result
            result = 0.5*(b-a)*(y_a + y_b);
            
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
                            ", result="+gammalib::str(result)+
                            ". Looks like n is too large.";
                gammalib::warning(G_TRAPZD, m_message);
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
                            ", result="+gammalib::str(result)+
                            ". Step is too small to make sense.";
                gammalib::warning(G_TRAPZD, m_message);
            }

            // Sum up values
            double x   = a + 0.5*del;
            double sum = 0.0;
            for (int j = 0; j < it; ++j, x+=del) {
                
                // Evaluate integrand
                double y = m_kernel->eval(x);
                m_calls++;

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


/***********************************************************************//**
 * @brief Adaptive Simpson's integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 *
 * Integrates the function using an adaptive Simpson's rule. The initial
 * interval [a,b] is split into two sub-intervals [a,c] and [c,b] for which
 * the integral is computed using
 *
 * \f[
 *    \frac{b-a}{6} f(a) + 4f(c) + f(b)
 * \f]
 *
 * where \f$c=(a+b)/2\f$ is the mid-point of interval [a,b]. Each
 * sub-interval is then recursively divided into sub-interval and the process
 * is repeated. Dividing of sub-intervals is stopped when the difference
 * between subsequent intervals falls below the relative tolerance specified
 * by eps(). The maximum recursion depth is set by the max_iter() method.
 *
 * I almost do not dare to confess, but the code has been taken from
 * http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
 * It's really pretty simple ...
 ***************************************************************************/
double GIntegral::adaptive_simpson(const double& a, const double& b) const
{
    // Initialise integration status information
    m_isvalid = true;
    m_calls   = 0;
    m_iter    = m_max_iter;

    // Compute mid-point c
    double c = 0.5*(a + b);    //!< Mid-point of interval [a,b]
    double h = b - a;          //!< Length of interval [a,b]

    // Evaluate function at boundaries and mid-point c
    double fa = m_kernel->eval(a);
    double fb = m_kernel->eval(b);
    double fc = m_kernel->eval(c);
    m_calls  += 3;
    
    // Compute integral using Simpson's rule
    double S = (h/6.0) * (fa + 4.0*fc + fb);

    // Call recursive auxiliary function
    double value = adaptive_simpson_aux(a, b, m_eps, S, fa, fb, fc, m_max_iter);

    // Deduce the number of iterations from the iteration counter
    m_iter = m_max_iter - m_iter;

    // If result is not valid, set and output status message
    if (!m_isvalid) {
        m_message = "Integration uncertainty exceeds relative tolerance "
                    "of "+gammalib::str(m_eps)+" after "+gammalib::str(m_iter)+
                    " iterations. Result "+gammalib::str(value)+" inaccurate.";
        if (!m_silent) {
            std::string origin = "GIntegral::adaptive_simpson("+
                                 gammalib::str(a)+", "+
                                 gammalib::str(b)+")";
            gammalib::warning(origin, m_message);
        }
    }

    // Return result
    return value;
}


/***********************************************************************//**
 * @brief Print integral information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing integral information.
 ***************************************************************************/
std::string GIntegral::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GIntegral ===");

        // Append information
        result.append("\n"+gammalib::parformat("Relative precision"));
        result.append(gammalib::str(eps()));
        result.append("\n"+gammalib::parformat("Function calls"));
        result.append(gammalib::str(calls()));
        result.append("\n"+gammalib::parformat("Iterations"));
        result.append(gammalib::str(iter()));
        result.append(" (maximum: ");
        result.append(gammalib::str(max_iter()));
        result.append(")");

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
void GIntegral::init_members(void)
{
    // Initialise members
    m_kernel    = NULL;
    m_eps       = 1.0e-6;
    m_max_iter  = 20;
    m_iter      = 0;
    m_isvalid   = true;
    m_message.clear();
    m_silent    = false;
    m_calls     = 0;

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
    m_kernel   = integral.m_kernel;
    m_eps      = integral.m_eps;
    m_max_iter = integral.m_max_iter;
    m_iter     = integral.m_iter;
    m_isvalid  = integral.m_isvalid;
    m_message  = integral.m_message;
    m_silent   = integral.m_silent;
    m_calls    = integral.m_calls;

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
double GIntegral::polint(double* xa, double* ya, int n, double x, double* dy)
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
                m_isvalid = false;
                m_message = "Invalid step size "+gammalib::str(den)+
                            " encountered. Two values in xa array are"
                            " identical.";
                gammalib::warning(G_POLINT, m_message);
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
 * @brief Auxiliary function for adaptive Simpson's method.
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] S Integral of last computation.
 * @param[in] fa Function value at left integration boundary.
 * @param[in] fb Function value at right integration boundary.
 * @param[in] fc Function value at mid-point of interval [a,b]
 * @param[in] bottom Iteration counter (stop when 0)
 *
 * Implements a recursive auxiliary method for the adative_simpson()
 * integrator.
 ***************************************************************************/
double GIntegral::adaptive_simpson_aux(const double& a, const double& b,
                                       const double& eps, const double& S,
                                       const double& fa, const double& fb,
                                       const double& fc,
                                       const int& bottom) const
{
    // Store the iteration counter
    m_iter = bottom;

    // Compute mid-point c bet
    double c = 0.5*(a + b);    //!< Mid-point of interval [a,b]
    double h = b - a;          //!< Length of interval [a,b]
    double d = 0.5*(a + c);    //!< Mid-point of interval [a,c]
    double e = 0.5*(c + b);    //!< Mid-point of interval [c,b]

    // Evaluate function at mid-points d and e
    double fd = m_kernel->eval(d);
    double fe = m_kernel->eval(e);
    m_calls  += 2;
    
    // Compute integral using Simpson's rule for the left and right interval
    double h12    = h / 12.0;
    double Sleft  = h12 * (fa + 4.0*fd + fc);
    double Sright = h12 * (fc + 4.0*fe + fb);
    double S2     = Sleft + Sright;

    // Allocate result
    double value;
 
    // If converged then compute the result ...
//    if (std::abs(S2 - S) <= 15.0 * eps * std::abs(S2)) {
    if (std::abs(S2 - S) <= 15.0 * m_eps * std::abs(S2)) {
        value = S2 + (S2 - S)/15.0;
    }
    
    // ... else if the maximum recursion depth was reached then compute the
    // result and signal result invalidity
    else if (bottom <= 0) {
        value     = S2 + (S2 - S)/15.0;
        m_isvalid = false;
    }
    
    // ... otherwise call this method recursively
    else {
        value = adaptive_simpson_aux(a, c, 0.5*eps, Sleft,  fa, fc, fd, bottom-1) +
                adaptive_simpson_aux(c, b, 0.5*eps, Sright, fc, fb, fe, bottom-1);
    }

    // Return result
    return value;
}
