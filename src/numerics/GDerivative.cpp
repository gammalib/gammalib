/***************************************************************************
 *                   GDerivative.cpp - Derivative class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GDerivative.cpp
 * @brief GDerivative class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cstdlib>          // Definition of NULL
#include <cmath>            // for abs()
#include <cfloat>           // for DBL_MAX
#include "GDerivative.hpp"
#include "GMatrix.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_VALUE                        "GDerivative::value(double&, double&)"
#define G_RIDDER             "GDerivative::ridder(double&, double&, double*)"
#define G_MINUIT2                    "GDerivative::minuit2(double&, double*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_VALUE_DEBUG  0                             //!< Debug value method
#define G_MINUIT_DEBUG 0                             //!< Debug minuit method


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
 * @return Derivative.
 ***************************************************************************/
GDerivative& GDerivative::operator=(const GDerivative& dx)
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
 * @brief Clear derivative
 ***************************************************************************/
void GDerivative::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone derivative
 *
 * @return Pointer to deep copy of derivative.
 ***************************************************************************/
GDerivative* GDerivative::clone(void) const
{
    return new GDerivative(*this);
}


/***********************************************************************//**
 * @brief Returns derivative
 *
 * @param[in] x Function value.
 * @param[in] step Initial step size (default: automatic)
 *
 * Computes derivative using the Ridders' method. If the specified initial
 * step size is non-zero, the method starts with this step size and
 * iteratively decreases the step size until the error becomes less than
 * a threshold defined with the eps() method (by default the threshold is
 * set to 1e-6). If the initial step size is zero, the method takes
 * \f$h={\tt m\_step\_frac} |x|\f$ as initial step size.
 * The maximum number of iterations is controlled by the max_iter()
 * method (by default, the maximum is set to 5).
 ***************************************************************************/
double GDerivative::value(const double& x, const double& step)
{
    // Set constants
    const double min_h = 1.0e-6;

    // Initialise result
    double result = 0.0;

    // Set initial step size. Use either the specified initial step size,
    // or if the initial step size is 0, use a fixed fraction of the function
    // value
    double h = step;
    if (h == 0.0) {
        h = m_step_frac * std::abs(x);
        if (h < min_h) {
            h = min_h;
        }
    }

    // Initialise tolerance
    double err = DBL_MAX;

    // Loop over Ridder's method until we have an acceptable tolerance
    m_iter = 0;
    for (; m_iter < m_max_iter; ++m_iter) {

        // Compute derivative using Ridder's method
        result = ridder(x, h, &err);

        // Debug option: Show actual results
        #if G_VALUE_DEBUG
        std::cout << "GDerivative::value(";
        std::cout << "iter=" << m_iter;
        std::cout << ", x=" << x;
        std::cout << ", dy/dx=" << result;
        std::cout << ", h=" << h;
        std::cout << ", err=" << err;
        std::cout << ")" << std::endl;
        #endif

        // If uncertainty is below tolerance then exit now
        if (err < m_eps) {
            break;
        }

        // Reduce the step size
        h /= 5.0;

    } // endfor: tolerance matching loop

    // Compile option: signal if we exceed the tolerance
    if (!silent()) {
        if (err >= m_eps) {
            std::string msg = "Derivative uncertainty "+gammalib::str(err)+
                              " exceeds tolerance of "+gammalib::str(m_eps)+
                              " at function value "+gammalib::str(x)+
                              " (df/dx="+gammalib::str(result)+
                              ", iter="+gammalib::str(m_iter)+
                              ", h="+gammalib::str(h)+").";
            gammalib::warning(G_VALUE, msg);
        }
    }

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
double GDerivative::ridder(const double& x, const double& h, double* err)
{
    // Set constants
    const int    ntab = 10;         // Maximum table size
    const double con  = 1.4;        // Stepsize decreases by con at each iteration
    const double con2 = con * con;
    const double big  = DBL_MAX;
    const double safe = 2.0;        // Return when error is safe worse than best so far

    // Initialise result
    double dx = 0.0;
    *err      = big;

    // Continue only if we have a valid function
    if (m_func != NULL) {

        // Allocate Neville table
        GMatrix a(ntab, ntab);

        // Throw exception if step size is zero
        if (h <= 0) {
            throw GException::invalid_argument(G_RIDDER, "Step size must be positive.");
        }

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
                if (errt <= *err) {
                    *err = errt;
                    dx   = a(j,i);
                }

            } // endfor: looped over tableau

            // If higher order is worse by a significant factor, then quit early
            if (std::abs( a(i,i) - a(i-1,i-1) ) >= safe*(*err)) {
                break;
            }

        } // endfor: looped over computations

    } // endif: we had a valid function

    // Return derivative
    return dx;
}


/***********************************************************************//**
 * @brief Returns derivative using Minuit2 algorithm
 *
 * @param[in] x Function value.
 * @param[out] err Estimated error in the derivative.
 *
 * This method has been inspired by corresponding code in the Minuit2
 * package. Please check the following files in Minuit2 for more information:
 *
 *    Numerical2PGradientCalculator.cxx for the main code
 *    MnStrategy.cxx for the definition of algorithm parameters
 *    MnMachinePrecision.h for the definition of eps and eps2
 *    InitialGradientCalculator.cxx for the Minuit2 way of estimating initial
 *    gradients.
 *
 * The exact value of fcnup has little impact on the results.
 ***************************************************************************/
double GDerivative::minuit2(const double& x, double* err)
{
    // Evaluate function at x
    double fcnmin = m_func->eval(x);

    // Set precision
//    double eps  = m_eps;       // Minuit2: smallest number that gives 1.+eps>1
//    double eps  = m_tiny;         // Bad
//    double eps  = 1.0e-30;        // Bad
//    double eps  = 1.0e-10;        // Poor
//    double eps  = 1.0e-6;        // Quite good
//    double eps  = 1.0e-8;        // Worse than 1.0e-6
    double eps  = 1.0e-4;        // Pretty good !
//    double eps  = 1.0e-2;        // Wrong
    double eps2 = 2.0 * std::sqrt(eps);

    // Set ...
    //double fcnup  = 1.0;     // Minuit2: Fcn().Up() - function delta for 1 sigma errors
    double fcnup  = 0.01 * fcnmin;
    double dfmin  = 8.0  * eps2 * (std::abs(fcnmin) + fcnup);
    double vrysml = 8.0  * eps * eps;

    // Set high strategy (see MnStrategy)
    int    ncycles  = 5;
    double step_tol = 0.1;
    double grad_tol = 0.02;

    // Initial gradient estimation. Note that this is brute force method
    // that does not correspond to the Minuit2 approach. For the Minuit2
    // approach, check InitialGradientCalculator.cxx
    double gstep = 0.001;    // Initial gradient step
    double fs1   = m_func->eval(x + gstep);
    double fs2   = m_func->eval(x - gstep);
    double grd   = 0.5 * (fs1 - fs2) / gstep;
    double g2    = (fs1 + fs2 - 2.0*fcnmin) / gstep / gstep;

    // Debug option: Show actual results
    #if G_MINUIT_DEBUG
    std::cout << "GDerivative::minuit(";
    std::cout << "iter=0";
    std::cout << ", x=" << x;
    std::cout << ", grd=" << grd;
    std::cout << ", step=" << gstep;
    std::cout << ")" << std::endl;
    #endif

    // Initialise values
    double epspri      = eps2 + std::abs(grd*eps2);
    double step        = 0.0;
    double stepb4      = 0.0;
    double step_change = 0.0;
    *err               = 0.0;

    // Loop over cycles
    for (m_iter = 1; m_iter <= ncycles; ++m_iter) {

        // Compute optimum step size
        double optstp = std::sqrt(dfmin/(std::abs(g2)+epspri));

        // Compute step size
        step = std::max(optstp, std::abs(0.1*gstep));
        double stpmax = 10.0 * std::abs(gstep);
        double stpmin = std::max(vrysml, 8.0 * std::abs(eps2*x));
        if (step > stpmax) step = stpmax;
        if (step < stpmin) step = stpmin;

        // Break is we are below the step tolerance
        step_change = std::abs((step-stepb4)/step);
        if (step_change < step_tol) {
            break;
        }

        // Store step size information
        gstep  = step;
        stepb4 = step;

        // Bookkeeping of last gradient
        double grdb4 = grd;

        // Evaluate gradient
        double fs1 = m_func->eval(x + step);
        double fs2 = m_func->eval(x - step);
        grd = 0.5 * (fs1 - fs2) / step;
        g2  = (fs1 + fs2 - 2.0*fcnmin) / step / step;

        // Compute error
        *err = std::abs(grdb4-grd) / (std::abs(grd)+dfmin/step);

        // Debug option: Show actual results
        #if G_MINUIT_DEBUG
        std::cout << "GDerivative::minuit(";
        std::cout << "iter=" << m_iter;
        std::cout << ", x=" << x;
        std::cout << ", grd=" << grd;
        std::cout << ", step=" << step;
        std::cout << ", err=" << *err;
        std::cout << ")" << std::endl;
        #endif

        // Break if gradient is accurate enough
        if (*err < grad_tol) {
            break;
        }

    } // endfor: looped over cycles

    // Compile option: signal if we exceed the tolerance
    if (!silent()) {
        if (*err >= grad_tol) {
            std::string msg = "Derivative uncertainty "+gammalib::str(*err)+
                              " exceeds tolerance of "+gammalib::str(grad_tol)+
                              " at function value "+gammalib::str(x)+
                              " (df/dx="+gammalib::str(grd)+
                              ", step="+gammalib::str(step)+
                              ", step_change="+gammalib::str(step_change)+").";
            gammalib::warning(G_MINUIT2, msg);
        }
    }

    // Return gradient
    return grd;
}


/***********************************************************************//**
 * @brief Returns gradient computed from function difference
 *
 * @param[in] x Function value.
 * @param[in] h Step size.
 *
 * This is the most simple and dumb gradient computation method we can think
 * of. It does a reasonable good job if the problem is well controlled, and
 * in particular, if the step size is known in advance.
 ***************************************************************************/
double GDerivative::difference(const double& x, const double& h)
{
    // Compute gradient
    double fs1 = m_func->eval(x + h);
    double fs2 = m_func->eval(x - h);
    double grd = 0.5 * (fs1 - fs2) / h;

    // Return gradient
    return grd;
}


/***********************************************************************//**
 * @brief Print derivative information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing derivative information.
 ***************************************************************************/
std::string GDerivative::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GDerivative ===");

        // Append information
        result.append("\n"+gammalib::parformat("Relative precision"));
        result.append(gammalib::str(eps()));
        result.append("\n"+gammalib::parformat("Max. number of iterations"));
        result.append(gammalib::str(max_iter()));
        result.append("\n"+gammalib::parformat("Initial step fraction"));
        result.append(gammalib::str(step_frac()));
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
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GDerivative::init_members(void)
{
    // Initialise members
    m_func      = NULL;
    m_eps       = 1.0e-6;
    m_step_frac = 0.02;
    m_tiny      = 1.0e-7;
    m_max_iter  = 5;
    m_iter      = 0;
    m_silent    = false;

    // Compute machine precision
    //set_tiny();

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
    m_func      = dx.m_func;
    m_eps       = dx.m_eps;
    m_step_frac = dx.m_step_frac;
    m_tiny      = dx.m_tiny;
    m_max_iter  = dx.m_max_iter;
    m_iter      = dx.m_iter;
    m_silent    = dx.m_silent;

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


/***********************************************************************//**
 * @brief Compute tiny number for Minuit2
 *
 * On fermi gives m_tiny = 8.88178e-16 for i=51.
 ***************************************************************************/
void GDerivative::set_tiny(void)
{
    // Allocate tiny instance
    GDerivative::tiny tiny;

    // Calculate machine precision
    double epstry = 0.5;
    double epsbak = 0.0;
    double epsp1  = 0.0;
    double one    = 1.0;

    // Loop until we found
    for (int i = 0; i < 100; ++i) {
        epstry *= 0.5;
        epsp1   = one + epstry;
        epsbak  = tiny(epsp1);
        if (epsbak < epstry) {
            m_tiny = 8.0 * epstry;
            break;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return one
 ***************************************************************************/
double GDerivative::tiny::one(void) const
{
    // Return
    return m_one;
}


/***********************************************************************//**
 * @brief Evaluate minimal difference between two floating points
 ***************************************************************************/
double GDerivative::tiny::operator()(double eps) const
{
    // Compute difference
    double result = eps - one();

    // Return result
    return result;
}
