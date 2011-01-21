/***************************************************************************
 *                 GRan.cpp - Randon number generator class                *
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
 * @file GRan.cpp
 * @brief Randon number generator class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GRan.hpp"
#include "GTools.hpp"
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
 ***************************************************************************/
GRan& GRan::operator= (const GRan& ran)
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
 * @brief Clear instance
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
 * @brief Clone instance
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
 * specified seed..
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
 *
 * This method is inspired from Numerical Recipes Third Edition, page 343.
 ***************************************************************************/
unsigned long int GRan::int32(void)
{
    // Return random value
    return ((unsigned long int)int64());
}


/***********************************************************************//**
 * @brief Return 64-bit random unsigned integer
 *
 * This method is inspired from Numerical Recipes Third Edition, page 343.
 ***************************************************************************/
unsigned long long int GRan::int64(void)
{
    // Update u
    m_u = m_u * 2862933555777941757LL + 7046029254386353087LL;

    // Shuffle v
    m_v ^= m_v >> 17;
    m_v ^= m_v << 31;
    m_v ^= m_v >> 8;

    // Update w
    m_w = 4294957665U * (m_w & 0xffffffff) + (m_w >> 32);

    // Compute x
    unsigned long long int x = m_u ^ (m_u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    
    // Return random value
    return ((x + m_v) ^ m_w);
}


/***********************************************************************//**
 * @brief Returns random double precision floating value in range 0 to 1
 *
 * This method is inspired from Numerical Recipes Third Edition, page 343.
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
 * This method is inspired from Numerical Recipes Third Edition, page 362.
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
    return (-log(x)/lambda);
}


/***********************************************************************//**
 * @brief Returns Chi2 deviates for 2 degrees of freedom
 *
 * @param[in] lambda Mean rate.
 *
 * Returns exponential deviates from the probability distribution
 * \f[p(x) = \frac{1}{2\pi} x \exp( -\frac{1}{2} x^2 )\f]
 * This method can be used to simulate the radom radial offset of a measured
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
    return (sqrt(-2.0*log(1.0-x)));
}


/***********************************************************************//**
 * @brief Print class information
 ***************************************************************************/
std::string GRan::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GRan ===");
    result.append("\n"+parformat("Seed")+str(m_seed));
    result.append("\n"+parformat("u")+str(m_u));
    result.append("\n"+parformat("v")+str(m_v));
    result.append("\n"+parformat("w")+str(m_w));

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
    m_seed = seed;
    m_v    = 4101842887655102017LL;
    m_w    = 1;

    // Initialise u, v, and w
    m_u = m_seed ^ m_v;
    int64();
    m_v = m_u;
    int64();
    m_w = m_v;
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
    m_seed = ran.m_seed;
    m_u    = ran.m_u;
    m_v    = ran.m_v;
    m_w    = ran.m_w;
    
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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] ran Column separated values table.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GRan& ran)
{
     // Write random number generator in output stream
    os << ran.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] ran Column separated values table.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GRan& ran)
{
    // Write random number generator into logger
    log << ran.print();

    // Return logger
    return log;
}
