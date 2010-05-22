/***************************************************************************
 *                        GEnergy.hpp - Energy class                       *
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
 * @file GEnergy.hpp
 * @brief Energy value class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GEnergy.hpp"

/* __ Constants __________________________________________________________ */

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
GEnergy::GEnergy(void)
{
    // Initialise private members for clean destruction
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] eng Object from which the instance should be built.
 ***************************************************************************/
GEnergy::GEnergy(const GEnergy& eng)
{ 
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(eng);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEnergy::~GEnergy(void)
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
 * @param[in] eng Object which should be assigned.
 ***************************************************************************/
GEnergy& GEnergy::operator= (const GEnergy& eng)
{ 
    // Execute only if object is not identical
    if (this != &eng) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(eng);

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
 * @brief Return energy in keV
 ***************************************************************************/
double GEnergy::keV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e+3;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return energy in MeV
 ***************************************************************************/
double GEnergy::MeV(void) const
{
    // Return energy
    return m_energy;
}


/***********************************************************************//**
 * @brief Return energy in GeV
 ***************************************************************************/
double GEnergy::GeV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e-3;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Return energy in TeV
 ***************************************************************************/
double GEnergy::TeV(void) const
{
    // Compute energy
    double energy = m_energy * 1.0e-6;
    
    // Return energy
    return energy; 
}


/***********************************************************************//**
 * @brief Set energy in keV
 *
 * @param[in] energy Energy in keV.
 ***************************************************************************/
void GEnergy::keV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e-3;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in MeV
 *
 * @param[in] energy Energy in MeV.
 ***************************************************************************/
void GEnergy::MeV(const double& eng)
{
    // Set energy
    m_energy = eng;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in GeV
 *
 * @param[in] energy Energy in GeV.
 ***************************************************************************/
void GEnergy::GeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+3;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in TeV
 *
 * @param[in] energy Energy in TeV.
 ***************************************************************************/
void GEnergy::TeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+6;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEnergy::init_members(void)
{
    // Initialise members
    m_energy = 0.0;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] eng Object from which members which should be copied.
 ***************************************************************************/
void GEnergy::copy_members(const GEnergy& eng)
{
    // Copy time
    m_energy = eng.m_energy;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEnergy::free_members(void)
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
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] eng Object to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEnergy& eng)
{
    // Put object in stream
    os << eng.MeV() << std::endl;

    // Return output stream
    return os;
}


