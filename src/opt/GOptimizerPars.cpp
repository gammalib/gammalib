/***************************************************************************
 *              GOptimizerPars.cpp  -  Parameter container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GOptimizerPars.cpp
 * @brief GOptimizerPars container class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GOptimizerPars.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                "GOptimizerPars::par(int) const"
#define G_COPY_MEMBERS  "GOptimizerPars::copy_members(const GOptimizerPars&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GOptimizerPars constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pars Parameters from which the instance should be built.
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(const GOptimizerPars& pars)
{ 
    // Initialise private members for clean destruction
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
 =                         GOptimizerPars operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pars Parameters which should be assigned.
 ***************************************************************************/
GOptimizerPars& GOptimizerPars::operator= (const GOptimizerPars& pars)
{ 
    // Execute only if object is not identical
    if (this != &pars) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pars);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                       GOptimizerPars public methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Compute number of free parameters
 ***************************************************************************/
int GOptimizerPars::nfree(void) const
{
    // Initialise number of free parameters
    int nfree = 0;
    
    // Collect all free parameters
    for (int i = 0; i < m_npars; ++i) {
        if (m_par[i]->isfree())
            nfree++;
    }
    
    // Return
    return nfree;
}


/***********************************************************************//**
 * @brief Get pointer to model parameter
 *
 * @param[in] index Parameter index.
 ***************************************************************************/
GModelPar* GOptimizerPars::par(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_PAR, index, 0, m_npars-1);
    
    // Return parameter pointer
    return m_par[index];
}


/*==========================================================================
 =                                                                         =
 =                       GOptimizerPars private methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizerPars::init_members(void)
{
    // Initialise members
    m_npars = 0;
    m_par   = NULL;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pars GOptimizerPars members which should be copied.
 ***************************************************************************/
void GOptimizerPars::copy_members(const GOptimizerPars& pars)
{
    // Copy attributes
    m_npars = pars.m_npars;
    
    // If there are parameters then copy them
    if (m_npars > 0 && pars.m_par != NULL) {
    
        // Allocate parameters
        m_par = new GModelPar*[m_npars];
        if (m_par == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_npars);
        
        // Copy parameters
        for (int i = 0; i < m_npars; ++i)
            m_par[i] = pars.m_par[i];

    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerPars::free_members(void)
{
    // Free memory
    if (m_par != NULL) delete [] m_par;

    // Signal free pointers
    m_par = NULL;
  
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GOptimizerPars friends                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put parameters in output stream
 *
 * @param[in] os Output stream into which the parameters will be dumped
 * @param[in] pars Parameters to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GOptimizerPars& pars)
{
    // Allocate filler
    std::string filler = " ..............";
    
    // Put model in stream
    os << "=== GOptimizerPars ===" << std::endl;
    os << " Number of parameters ......: " << pars.npars();
    for (int i = 0; i < pars.npars(); ++i) {
        os << std::endl;
        if (i == 10) filler = " .............";
        if (i == 100) filler = " ............";
        if (i == 1000) filler = " ...........";
        os << "  Parameter " << i << filler << ": " << *(pars.par(i));
    }

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GOptimizerPars                 =
 =                                                                         =
 ==========================================================================*/
