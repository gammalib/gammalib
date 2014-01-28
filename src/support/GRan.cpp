/***************************************************************************
 *                 GRan.cpp - Random number generator class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GRan.cpp
 * @brief Random number generator class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GRan.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GRan::GRan(void)
{
    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Seed constructor
 *
 * @param[in] seed Random number generator seed.
 ***************************************************************************/
GRan::GRan(unsigned long long int seed)
{ 
    // Initialise private
    init_members(seed);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ran Random number generator.
 ***************************************************************************/
GRan::GRan(const GRan& ran)
{ 
    // Initialise private
    init_members();

    // Copy members
    copy_members(ran);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GRan::~GRan(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] ran Random number generator.
 * @return Random number generator.
 ***************************************************************************/
GRan& GRan::operator=(const GRan& ran)
{ 
    // Execute only if object is not identical
    if (this != &ran) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(ran);

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
 * @brief Clear random number generator
 *
 * This method properly resets the object to an initial state using the
 * default seed.
 ***************************************************************************/
void GRan::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone random number generator
 *
 * @return Pointer to deep copy of random number generator.
 **************************************************************************/
GRan* GRan::clone(void) const
{
    return new GRan(*this);
}


/***********************************************************************//**
 * @brief Set seed of random number generator
 *
 * @param[in] seed Random number generator seed.
 *
 * This method properly resets the object to an initial state using the
 * specified seed.
 ***************************************************************************/
void GRan::seed(unsigned long long int seed)
{
    // Free class members
    free_members();

    // Initialise members
    init_members(seed);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return 32-bit random unsigned integer
 ***************************************************************************/
unsigned long int GRan::int32(void)
{
    // Return random value
    return ((unsigned long int)int64());
}


/***********************************************************************//**
 * @brief Return 64-bit random unsigned integer
 ***************************************************************************/
unsigned long long int GRan::int64(void)
{
    // Update u
    m_value1 = m_value1 * 2862933555777941757LL + 7046029254386353087LL;

    // Shuffle v
    m_value2 ^= m_value2 >> 17;
    m_value2 ^= m_value2 << 31;
    m_value2 ^= m_value2 >> 8;

    // Update w
    m_value3 = 4294957665U * (m_value3 & 0xffffffff) + (m_value3 >> 32);

    // Compute x
    unsigned long long int x = m_value1 ^ (m_value1 << 21);
    x ^= x >> 35;
    x ^= x << 4;
    
    // Return random value
    return ((x + m_value2) ^ m_value3);
}


/***********************************************************************//**
 * @brief Returns random double precision floating value in range 0 to 1
 ***************************************************************************/
double GRan::uniform(void)
{
    // Return random value
    return (5.42101086242752217e-20 * int64());
}


/***********************************************************************//**
 * @brief Returns exponential deviates
 *
 * @param[in] lambda Mean rate.
 *
 * Returns exponential deviates from the probability distribution
 * \f[p(x) = \lambda \exp( -\lambda x )\f]
 * where
 * \f$\lambda\f$ is the parameter of the distribution.
 * This method may be used to simulate the occurence time of an event, where
 * \f$\lambda\f$ is the mean event rate. Convsersely, \f$1/\lambda\f$ is the
 * mean waiting time between events.
 *
 * @todo Check that \f$\lambda>0\f$.
 ***************************************************************************/
double GRan::exp(const double& lambda)
{
    // Allocate argument
    double x;

    // Get non-zero uniform deviate
    do {
        x = uniform();
    } while (x == 0.0);

    // Return random value
    return (-std::log(x)/lambda);
}


/***********************************************************************//**
 * @brief Returns Poisson deviates
 *
 * @param[in] lambda Expectation value.
 *
 * Returns Poisson deviates for an expectation value.
 *
 * This method may be used to simulate the number of events in case that a
 * given mean number of events is expected.
 ***************************************************************************/
double GRan::poisson(const double& lambda)
{
    // Declare result
    double value;

    // Use direct method for small numbers ...
    if (lambda < 12.0) {
        if (lambda != m_old_lambda) {
            m_old_lambda = lambda;
            m_exp_lambda = std::exp(-lambda);
        }
        value      = -1.0;
        double tmp =  1.0;
        do {
            value += 1.0;
            tmp   *= uniform();
        } while (tmp > m_exp_lambda);
    } // endif: direct method used

    // ... otherwise use rejection method        
    else {
        double tmp;
        if (lambda != m_old_lambda) {
            m_old_lambda  = lambda;
            m_sqrt_lambda = std::sqrt(2.0*lambda);
            m_log_lambda  = std::log(lambda);
            m_exp_lambda  = lambda * m_log_lambda - gammalib::gammln(lambda+1.0);
        }
        do {
            double factor;
            do {
                factor = std::tan(gammalib::pi * uniform());
                value  = m_sqrt_lambda * factor + lambda;
            } while (value < 0.0);
            value = floor(value);
            tmp   = 0.9*(1.0+factor*factor) *
                    std::exp(value*m_log_lambda - gammalib::gammln(value+1.0)-m_exp_lambda);
        } while (uniform() > tmp);
    }
    
    // Return random deviate
    return value;
}


/***********************************************************************//**
 * @brief Returns Chi2 deviates for 2 degrees of freedom
 *
 * Returns exponential deviates from the probability distribution
 * \f[p(x) = \frac{1}{2\pi} x \exp( -\frac{1}{2} x^2 )\f]
 * This method can be used to simulate the random radial offset of a measured
 * source position from the true source position, assuming an azimuthally
 * symmetric 2D Gaussian probability distribution.
 ***************************************************************************/
double GRan::chisq2(void)
{
    // Allocate argument
    double x;

    // Get uniform deviate < 1
    do {
        x = uniform();
    } while (x == 1.0);

    // Return random value
    return (std::sqrt(-2.0*std::log(1.0-x)));
}


/***********************************************************************//**
 * @brief Print random number generator information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing random number generator information.
 ***************************************************************************/
std::string GRan::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GRan ===");

        // Append information
        result.append("\n"+gammalib::parformat("Seed")+gammalib::str(m_seed));
        result.append("\n"+gammalib::parformat("Value 1")+gammalib::str(m_value1));
        result.append("\n"+gammalib::parformat("Value 2")+gammalib::str(m_value2));
        result.append("\n"+gammalib::parformat("Value 3")+gammalib::str(m_value3));

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
void GRan::init_members(unsigned long long int seed)
{
    // Initialise members
    m_seed        = seed;
    m_value2      = 4101842887655102017LL;
    m_value3      = 1;
    m_old_lambda  = -1.0;
    m_sqrt_lambda = 0.0;
    m_log_lambda  = 0.0;
    m_exp_lambda  = 0.0;

    // Initialise values
    m_value1 = m_seed ^ m_value2;
    int64();
    m_value2 = m_value1;
    int64();
    m_value3 = m_value2;
    int64();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ran Random number generator.
 ***************************************************************************/
void GRan::copy_members(const GRan& ran)
{
    // Copy members
    m_seed        = ran.m_seed;
    m_value1      = ran.m_value1;
    m_value2      = ran.m_value2;
    m_value3      = ran.m_value3;
    m_old_lambda  = ran.m_old_lambda;
    m_sqrt_lambda = ran.m_sqrt_lambda;
    m_log_lambda  = ran.m_log_lambda;
    m_exp_lambda  = ran.m_exp_lambda;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GRan::free_members(void)
{
    // Return
    return;
}
