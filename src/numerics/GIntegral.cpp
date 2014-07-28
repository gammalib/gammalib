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
#include <cmath>            // std::abs()
#include <vector>
#include <algorithm>        // std::sort
#include "GIntegral.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ROMBERG                "GIntegral::romberg(double&, double&, int&)"
#define G_TRAPZD          "GIntegral::trapzd(double&, double&, int&, double)"
#define G_POLINT  "GIntegral::polint(double*, double*, int, double, double*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
namespace gammalib {

    // Gauss-Kronrod abscissae, common to the 10-, 21-, 43- and 87-point rule
    const double gkx1[5] = {
        0.973906528517171720077964012084452,
        0.865063366688984510732096688423493,
        0.679409568299024406234327365114874,
        0.433395394129247190799265943165784,
        0.148874338981631210884826001129720
    };

    // Gauss-Kronrod weights of the 10-point rule
    const double gkw10[5] = {
        0.066671344308688137593568809893332,
        0.149451349150580593145776339657697,
        0.219086362515982043995534934228163,
        0.269266719309996355091226921569469,
        0.295524224714752870173892994651338
    };

    // Gauss-Kronrod abscissae, common to the 21-, 43- and 87-point rule
    const double gkx2[5] = {
        0.995657163025808080735527280689003,
        0.930157491355708226001207180059508,
        0.780817726586416897063717578345042,
        0.562757134668604683339000099272694,
        0.294392862701460198131126603103866
    };

    // Gauss-Kronrod weights of the 21-point rule for abscissae gkx1
    const double gkw21a[5] = {
        0.032558162307964727478818972459390,
        0.075039674810919952767043140916190,
        0.109387158802297641899210590325805,
        0.134709217311473325928054001771707,
        0.147739104901338491374841515972068
    };

    // Gauss-Kronrod weights of the 21-point rule for abscissae gkx2
    const double gkw21b[6] = {
        0.011694638867371874278064396062192,
        0.054755896574351996031381300244580,
        0.093125454583697605535065465083366,
        0.123491976262065851077958109831074,
        0.142775938577060080797094273138717,
        0.149445554002916905664936468389821
    };

    // Gauss-Kronrod abscissae, common to the 43- and 87-point rule
    const double gkx3[11] = {
        0.999333360901932081394099323919911,
        0.987433402908088869795961478381209,
        0.954807934814266299257919200290473,
        0.900148695748328293625099494069092,
        0.825198314983114150847066732588520,
        0.732148388989304982612354848755461,
        0.622847970537725238641159120344323,
        0.499479574071056499952214885499755,
        0.364901661346580768043989548502644,
        0.222254919776601296498260928066212,
        0.074650617461383322043914435796506
    };

    // Gauss-Kronrod weights of the 43-point rule for abscissae gkx1, gkx3
    const double gkw43a[10] = {
        0.016296734289666564924281974617663,
        0.037522876120869501461613795898115,
        0.054694902058255442147212685465005,
        0.067355414609478086075553166302174,
        0.073870199632393953432140695251367,
        0.005768556059769796184184327908655,
        0.027371890593248842081276069289151,
        0.046560826910428830743339154433824,
        0.061744995201442564496240336030883,
        0.071387267268693397768559114425516
    };

    // Gauss-Kronrod weights of the 43-point formula for abscissae gkx3
    const double gkw43b[12] = {
        0.001844477640212414100389106552965,
        0.010798689585891651740465406741293,
        0.021895363867795428102523123075149,
        0.032597463975345689443882222526137,
        0.042163137935191811847627924327955,
        0.050741939600184577780189020092084,
        0.058379395542619248375475369330206,
        0.064746404951445885544689259517511,
        0.069566197912356484528633315038405,
        0.072824441471833208150939535192842,
        0.074507751014175118273571813842889,
        0.074722147517403005594425168280423
    };

    // Gauss-Kronrod abscissae, of the 87-point rule
    const double gkx4[22] = {
        0.999902977262729234490529830591582,
        0.997989895986678745427496322365960,
        0.992175497860687222808523352251425,
        0.981358163572712773571916941623894,
        0.965057623858384619128284110607926,
        0.943167613133670596816416634507426,
        0.915806414685507209591826430720050,
        0.883221657771316501372117548744163,
        0.845710748462415666605902011504855,
        0.803557658035230982788739474980964,
        0.757005730685495558328942793432020,
        0.706273209787321819824094274740840,
        0.651589466501177922534422205016736,
        0.593223374057961088875273770349144,
        0.531493605970831932285268948562671,
        0.466763623042022844871966781659270,
        0.399424847859218804732101665817923,
        0.329874877106188288265053371824597,
        0.258503559202161551802280975429025,
        0.185695396568346652015917141167606,
        0.111842213179907468172398359241362,
        0.037352123394619870814998165437704
    };

    // Gauss-Kronrod weights of the 87-point rule for abscissae gkx1, gkx2, gkx3
    const double gkw87a[21] = {
        0.008148377384149172900002878448190,
        0.018761438201562822243935059003794,
        0.027347451050052286161582829741283,
        0.033677707311637930046581056957588,
        0.036935099820427907614589586742499,
        0.002884872430211530501334156248695,
        0.013685946022712701888950035273128,
        0.023280413502888311123409291030404,
        0.030872497611713358675466394126442,
        0.035693633639418770719351355457044,
        0.000915283345202241360843392549948,
        0.005399280219300471367738743391053,
        0.010947679601118931134327826856808,
        0.016298731696787335262665703223280,
        0.021081568889203835112433060188190,
        0.025370969769253827243467999831710,
        0.029189697756475752501446154084920,
        0.032373202467202789685788194889595,
        0.034783098950365142750781997949596,
        0.036412220731351787562801163687577,
        0.037253875503047708539592001191226
    };

    // Gauss-Kronrod weights of the 87-point formula for abscissae gkx4
    const double gkw87b[23] = {
        0.000274145563762072350016527092881,
        0.001807124155057942948341311753254,
        0.004096869282759164864458070683480,
        0.006758290051847378699816577897424,
        0.009549957672201646536053581325377,
        0.012329447652244853694626639963780,
        0.015010447346388952376697286041943,
        0.017548967986243191099665352925900,
        0.019938037786440888202278192730714,
        0.022194935961012286796332102959499,
        0.024339147126000805470360647041454,
        0.026374505414839207241503786552615,
        0.028286910788771200659968002987960,
        0.030052581128092695322521110347341,
        0.031646751371439929404586051078883,
        0.033050413419978503290785944862689,
        0.034255099704226061787082821046821,
        0.035262412660156681033782717998428,
        0.036076989622888701185500318003895,
        0.036698604498456094498018047441094,
        0.037120549269832576114119958413599,
        0.037334228751935040321235449094698,
        0.037361073762679023410321241766599
    };

} // end gammalib namespace


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
 * @param[in] bounds Integration boundaries.
 * @param[in] order Integration order (default: 5)
 *
 * @exception GException::invalid_argument
 *            Integration order incompatible with number of iterations.
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
double GIntegral::romberg(std::vector<double> bounds, const int& order)
{
    // Sort integration boundaries in ascending order
    std::sort(bounds.begin(), bounds.end());

    // Initialise integral
    double value = 0.0;

    // Add integral of all intervals
    for (int i = 0; i < bounds.size()-1; ++i) {
        value += romb(bounds[i], bounds[i+1], order);
    }

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Perform Romberg integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 * @param[in] order Integration order (default: 5)
 *
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
double GIntegral::romb(const double& a, const double& b, const int& order)
{
    // Initialise result and status
    double result = 0.0;

    // Initialise integration status information
    m_isvalid    = true;
    m_calls      = 0;
    m_has_abserr = false;
    m_relerr     = false;
    
    // Continue only if integration range is valid
    if (b > a) {

        // Initialise variables
        bool   converged = false;
        double dss       = 0.0;

        // Determine (maximum) number of iterations
        int max_iter = (m_fix_iter > 0) ? m_fix_iter : m_max_iter;

        // Check whether maximum number of iterations is compliant with
        // order
        if (order > max_iter) {
            std::string msg = "Requested integration order "+
                              gammalib::str(order)+" is larger than the "
                              "maximum number of iterations "+
                              gammalib::str(max_iter)+". Either reduced the "
                              "integration order or increase the (maximum) "
                              "number of iterations.";
            throw GException::invalid_argument(G_ROMBERG, msg);
        }

        // Allocate temporal storage
        double* s = new double[max_iter+2];
        double* h = new double[max_iter+2];

        // Initialise step size
        h[1] = 1.0;
        s[0] = 0.0;

        // Iterative loop
        for (m_iter = 1; m_iter <= max_iter; ++m_iter) {

            // Integration using Trapezoid rule
            s[m_iter] = trapzd(a, b, m_iter, s[m_iter-1]);

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (is_notanumber(s[m_iter]) || is_infinite(s[m_iter])) {
                m_message = "*** ERROR: GIntegral::romberg"
                            "(a="+gammalib::str(a)+", b="+gammalib::str(b)+""
                            ", k="+gammalib::str(k)+"): NaN/Inf encountered"
                            " (s["+gammalib::str(m_iter)+"]="
                            ""+gammalib::str(s[m_iter])+")";
                std::cout << m_message << std::endl;
                m_isvalid = false;
            }
            #endif

            // Starting from iteration order on, use polynomial interpolation
            if (m_iter >= order) {

                // Compute result using polynom interpolation
                result = polint(&h[m_iter-order], &s[m_iter-order],
                                order, 0.0, &dss);

                // If a fixed number of iterations has been requested then
                // check whether we reached the final one; otherwise check
                // whether we reached the requested precision.
                if (((m_fix_iter > 0) && (m_iter == max_iter)) ||
                    (std::abs(dss) <= m_eps * std::abs(result))) {
                    converged    = true;
                    m_has_abserr = true;
                    m_abserr     = std::abs(dss);
                    if (std::abs(result) > 0) {
                        m_has_relerr = true;
                        m_relerr     = m_abserr / std::abs(result);
                    }
                    break;
                }

            } // endif: polynomial interpolation performed

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
                        gammalib::str(m_eps * std::abs(result))+
                        " after "+gammalib::str(m_iter)+
                        " iterations. Result "+
                        gammalib::str(result)+
                        " is inaccurate.";
            if (!m_silent) {
                std::string origin = "GIntegral::romberg("+
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
    m_isvalid    = true;
    m_calls      = 0;
    m_iter       = m_max_iter;
    m_has_abserr = false;
    m_relerr     = false;

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

    // Initialise absolute precision
    double epsilon = (std::abs(S) > 0) ? m_eps * S : m_eps;

    // Call recursive auxiliary function
    double value = adaptive_simpson_aux(a, b, epsilon, S, fa, fb, fc, m_max_iter);

    // Deduce the number of iterations from the iteration counter
    m_iter = m_max_iter - m_iter;

    // If result is not valid, set and output status message
    if (!m_isvalid) {
        m_message = "Integration uncertainty exceeds relative tolerance "
                    "of "+gammalib::str(m_eps)+" and absolute tolerance of "
                    ""+gammalib::str(epsilon)+" after "+gammalib::str(m_iter)+
                    " iterations and "+gammalib::str(m_calls)+" function "
                    "calls. Result "+gammalib::str(value)+" inaccurate.";
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
 * @brief Gauss-Kronrod integration
 *
 * @param[in] a Left integration boundary.
 * @param[in] b Right integration boundary.
 *
 ***************************************************************************/
double GIntegral::gauss_kronrod(const double& a, const double& b) const
{
    // Initialise integration status information
    m_isvalid    = true;
    m_iter       = 0;
    m_calls      = 0;
    m_has_abserr = false;
    m_relerr     = false;

    // Initialise integration result
    double result = 0.0;
    double error  = 0.0; 

    // Allocate some arrays
    double fv1[5];
    double fv2[5];
    double fv3[5];
    double fv4[5];
    double savfun[21];

    // Main code loop (so that we can exit using break)
    do {

        // Tolerance check
        if (m_eps < 1.12e-14) {
            m_isvalid = false;
            m_message = "Requested relative tolerance of "+gammalib::str(m_eps)+
                        " cannot be acheived. Please relax the integration "
                        "precision.";
            if (!m_silent) {
                std::string origin = "GIntegral::gauss_kronrod("+
                                     gammalib::str(a)+", "+
                                     gammalib::str(b)+")";
                gammalib::warning(origin, m_message);
            }
        }

        // Compute function at mid-point
        double h     = 0.5 * (b - a);
        double abs_h = std::abs(h);
        double c     = 0.5 * (b + a);
        double f_c   = m_kernel->eval(c);
        m_calls++;

        // Compute the integral using the 10- and 21-point formulae
        m_iter++;
        double res10  = 0;
        double res21  = gammalib::gkw21b[5] * f_c;
        double resabs = gammalib::gkw21b[5] * std::abs(f_c);
        for (int k = 0; k < 5; ++k) {
            double x     = h * gammalib::gkx1[k];
            double fval1 = m_kernel->eval(c+x);
            double fval2 = m_kernel->eval(c-x);
            double fval  = fval1 + fval2;
            m_calls     += 2;
            res10       += gammalib::gkw10[k]  * fval;
            res21       += gammalib::gkw21a[k] * fval;
            resabs      += gammalib::gkw21a[k] * (std::abs(fval1) + std::abs(fval2));
            savfun[k]    = fval;
            fv1[k]       = fval1;
            fv2[k]       = fval2;
        }
        for (int k = 0; k < 5; ++k) {
            double x = h * gammalib::gkx2[k];
            double fval1 = m_kernel->eval(c+x);
            double fval2 = m_kernel->eval(c-x);
            double fval  = fval1 + fval2;
            m_calls     += 2;
            res21       += gammalib::gkw21b[k] * fval;
            resabs      += gammalib::gkw21b[k] * (std::abs(fval1) + std::abs(fval2));
            savfun[k+5]  = fval;
            fv3[k]       = fval1;
            fv4[k]       = fval2;
        }
        resabs       *= abs_h;
        double mean   = 0.5 * res21;
        double resasc = gammalib::gkw21b[5] * std::abs(f_c - mean);
        for (int k = 0; k < 5; ++k) {
            resasc += (gammalib::gkw21a[k] *
                       (std::abs(fv1[k] - mean) + std::abs(fv2[k] - mean)) +
                       gammalib::gkw21b[k] *
                       (std::abs(fv3[k] - mean) + std::abs(fv4[k] - mean)));
        }
        resasc *= abs_h ;
        result  = res21 * h;
        error   = rescale_error((res21 - res10) * h, resabs, resasc);

        // Test for convergence */
        //if (err < epsabs || err < epsrel * fabs (result))
        if (error < m_eps * std::abs(result)) {
            m_has_abserr = true;
            m_abserr     = error;
            if (std::abs(result) > 0) {
                m_has_relerr = true;
                m_relerr     = error / std::abs(result);
            }
            break;
        }

        // Compute the integral using the 43-point formula
        m_iter++;
        double res43 = gammalib::gkw43b[11] * f_c;
        for (int k = 0; k < 10; ++k) {
            res43 += savfun[k] * gammalib::gkw43a[k];
        }
        for (int k = 0; k < 11; ++k) {
            double x     = h * gammalib::gkx3[k];
            double fval  = (m_kernel->eval(c+x) +
                            m_kernel->eval(c-x));
            m_calls     += 2;
            res43       += fval * gammalib::gkw43b[k];
            savfun[k+10] = fval;
        }
        result = res43 * h;

        // Test for convergence */
        error = rescale_error((res43 - res21) * h, resabs, resasc);
        //if (err < epsabs || err < epsrel * fabs (result))
        if (error < m_eps * std::abs(result)) {
            m_has_abserr = true;
            m_abserr     = error;
            if (std::abs(result) > 0) {
                m_has_relerr = true;
                m_relerr     = error / std::abs(result);
            }
            break;
        }

        // Compute the integral using the 87-point formula
        m_iter++;
        double res87 = gammalib::gkw87b[22] * f_c;
        for (int k = 0; k < 21; ++k) {
            res87 += savfun[k] * gammalib::gkw87a[k];
        }
        for (int k = 0; k < 22; ++k) {
            double x = h * gammalib::gkx4[k];
            res87   += gammalib::gkw87b[k] *
                       (m_kernel->eval(c+x) +
                        m_kernel->eval(c-x));
            m_calls += 2;
        }
        result = res87 * h ;

        // Test for convergence */
        error = rescale_error ((res87 - res43) * h, resabs, resasc);
        //if (err < epsabs || err < epsrel * fabs (result))
        if (error < m_eps * std::abs(result)) {
            m_has_abserr = true;
            m_abserr     = error;
            if (std::abs(result) > 0) {
                m_has_relerr = true;
                m_relerr     = error / std::abs(result);
            }
            break;
        }

        // Failed to converge
        m_isvalid = false;
        m_message = "Integration uncertainty "+gammalib::str(error)+" exceeds "
                    "absolute tolerance of "+
                    gammalib::str(m_eps * std::abs(result))+" after "+
                    gammalib::str(m_iter)+" iterations and "+
                    gammalib::str(m_calls)+" function calls. Result "+
                    gammalib::str(result)+" inaccurate.";
        if (!m_silent) {
            std::string origin = "GIntegral::gauss_kronrod("+
                                 gammalib::str(a)+", "+
                                 gammalib::str(b)+")";
            gammalib::warning(origin, m_message);
        }
        
    } while (false); // end of main loop
    
    // Return result
    return result;
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
        if (m_has_abserr) {
            result.append("\n"+gammalib::parformat("Absolute error"));
            result.append(gammalib::str(m_abserr));
        }
        if (m_has_relerr) {
            result.append("\n"+gammalib::parformat("Relative error"));
            result.append(gammalib::str(m_relerr));
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
void GIntegral::init_members(void)
{
    // Initialise members
    m_kernel    = NULL;
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
    m_abserr     = 0.0;
    m_relerr     = 0.0;

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
 * @brief Auxiliary function for adaptive Simpson's method
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
    if (bottom < m_iter) {
        m_iter = bottom;
    }

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
    if (std::abs(S2 - S) <= 15.0 * eps) {
//    if (std::abs(S2 - S) <= 15.0 * m_eps * std::abs(S2)) {
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


//#define GSL_DBL_EPSILON        2.2204460492503131e-16
//#define GSL_DBL_MIN        2.2250738585072014e-308
/***********************************************************************//**
 * @brief Rescale errors for Gauss-Kronrod integration
 *
 * @param[in] err Error estimate.
 * @param[in] result_abs ???.
 * @param[in] result_asc ???.
 * @return Rescaled error estimate.
 ***************************************************************************/
double GIntegral::rescale_error(double err, const double& result_abs, const double& result_asc) const
{
    // Take absolute value of error
    err = std::abs(err);

    // ...
    if (result_asc != 0.0 && err != 0.0) {
        double scale = std::pow((200.0 * err / result_asc), 1.5);
        if (scale < 1.0) {
            err = result_asc * scale ;
        }
        else {
            err = result_asc ;
        }
    }
    if (result_abs > 2.2250738585072014e-308 / (50.0 * 2.2204460492503131e-16)) {
        double min_err = 50.0 * 2.2204460492503131e-16 * result_abs;
        if (min_err > err) {
            err = min_err ;
        }
    }

    // Return error
    return err;
}
