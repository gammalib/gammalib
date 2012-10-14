/***************************************************************************
 *              GOptimizerPars.cpp  -  Parameter container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
 * @file GOptimizerPars.cpp
 * @brief Optimizer parameter container class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GOptimizerPars.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                     "GOptimizerPars::par(int&)"

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
GOptimizerPars::GOptimizerPars(void)
{
    // Initialise members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(const GOptimizerPars& pars)
{ 
    // Initialise members
    init_members();

    // Copy members
    copy_members(pars);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GOptimizerPars::~GOptimizerPars()
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
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
GOptimizerPars& GOptimizerPars::operator= (const GOptimizerPars& pars)
{ 
    // Execute only if object is not identical
    if (this != &pars) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(pars);

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
 * @brief Clear object
 ***************************************************************************/
void GOptimizerPars::clear(void)
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
 ***************************************************************************/
GOptimizerPars* GOptimizerPars::clone(void) const
{
    return new GOptimizerPars(*this);
}


/***********************************************************************//**
 * @brief Returns number of free parameters
 ***************************************************************************/
int GOptimizerPars::nfree(void) const
{
    // Initialise number of free parameters
    int nfree = 0;
    
    // Collect all free parameters
    for (int i = 0; i < npars(); ++i) {
        if (m_pars[i]->isfree())
            nfree++;
    }
    
    // Return
    return nfree;
}


/***********************************************************************//**
 * @brief Returns reference to model parameter
 *
 * @param[in] index Parameter index [0,...,npars()-1]
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
GModelPar& GOptimizerPars::par(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= npars())
        throw GException::out_of_range(G_PAR, index, 0, npars()-1);
    #endif
    
    // Return parameter reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,npars()-1]
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 ***************************************************************************/
const GModelPar& GOptimizerPars::par(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= npars())
        throw GException::out_of_range(G_PAR, index, 0, npars()-1);
    #endif
    
    // Return parameter reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Print models
 ***************************************************************************/
std::string GOptimizerPars::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GOptimizerPars ===");
    result.append("\n"+parformat("Number of parameters")+str(npars()));
    for (int i = 0; i < npars(); ++i)
        result.append("\n"+m_pars[i]->print());

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
void GOptimizerPars::init_members(void)
{
    // Initialise members
    m_pars.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
void GOptimizerPars::copy_members(const GOptimizerPars& pars)
{
    // Copy attributes
    m_pars = pars.m_pars;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerPars::free_members(void)
{
    // Return
    return;
}
