/***************************************************************************
 *                        GEnergy.cpp - Energy class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GEnergy.cpp
 * @brief Energy value class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cfloat>
#include <cmath>
#include "GEnergy.hpp"
#include "GTools.hpp"

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
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] eng Energy.
 ***************************************************************************/
GEnergy::GEnergy(const GEnergy& eng)
{ 
    // Initialise private members
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
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] eng Energy.
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
 * @brief Clear instance
 ***************************************************************************/
void GEnergy::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GEnergy* GEnergy::clone(void) const
{
    // Clone this image
    return new GEnergy(*this);
}


/***********************************************************************//**
 * @brief Return energy in erg
 ***************************************************************************/
double GEnergy::erg(void) const
{
    // Compute energy
    double energy = m_energy * MeV2erg;
    
    // Return energy
    return energy; 
}


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
 * @brief Return log10 of energy in keV
 *
 * Returns the log10 of the energy in keV.
 ***************************************************************************/
double GEnergy::log10keV(void) const
{
    // Return log10 energy
    return (log10MeV()+3.0); 
}


/***********************************************************************//**
 * @brief Return log10 of energy in MeV
 *
 * Returns the log10 of the energy in MeV. The result is stored internally
 * and not recomputed when the method is called again with the same energy
 * value. This speeds up computation. In case that the energy is not positive
 * the method returns DBL_MIN.
 ***************************************************************************/
double GEnergy::log10MeV(void) const
{
    // If required compute log10 of energy.
    if (!m_has_log10) {
        m_elog10    = (m_energy > 0.0) ? std::log10(m_energy) : DBL_MIN;
        m_has_log10 = true;
    }
    
    // Return log10 energy
    return m_elog10; 
}


/***********************************************************************//**
 * @brief Return log10 of energy in GeV
 *
 * Returns the log10 of the energy in GeV.
 ***************************************************************************/
double GEnergy::log10GeV(void) const
{
    // Return log10 energy
    return (log10MeV()-3.0); 
}


/***********************************************************************//**
 * @brief Return log10 of energy in TeV
 *
 * Returns the log10 of the energy in TeV.
 ***************************************************************************/
double GEnergy::log10TeV(void) const
{
    // Return log10 energy
    return (log10MeV()-6.0); 
}


/***********************************************************************//**
 * @brief Set energy in erg
 *
 * @param[in] eng Energy in erg.
 ***************************************************************************/
void GEnergy::erg(const double& eng)
{
    // Set energy
    m_energy = eng * erg2MeV;

    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in keV
 *
 * @param[in] eng Energy in keV.
 ***************************************************************************/
void GEnergy::keV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e-3;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in MeV
 *
 * @param[in] eng Energy in MeV.
 ***************************************************************************/
void GEnergy::MeV(const double& eng)
{
    // Set energy
    m_energy = eng;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in GeV
 *
 * @param[in] eng Energy in GeV.
 ***************************************************************************/
void GEnergy::GeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+3;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy in TeV
 *
 * @param[in] eng Energy in TeV.
 ***************************************************************************/
void GEnergy::TeV(const double& eng)
{
    // Set energy
    m_energy = eng * 1.0e+6;
    
    // Signals that log10 of energy is not valid
    m_has_log10 = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in keV
 *
 * @param[in] eng log10 of energy in keV.
 ***************************************************************************/
void GEnergy::log10keV(const double& eng)
{
    // Set energy
    log10MeV(eng-3.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in MeV
 *
 * @param[in] eng log10 of energy in MeV.
 ***************************************************************************/
void GEnergy::log10MeV(const double& eng)
{
    // Set energy
    m_elog10    = eng;
    m_energy    = std::pow(10.0, eng);
    m_has_log10 = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in GeV
 *
 * @param[in] eng log10 of energy in GeV.
 ***************************************************************************/
void GEnergy::log10GeV(const double& eng)
{
    // Set energy
    log10MeV(eng+3.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set log10 of energy in TeV
 *
 * @param[in] eng log10 of energy in TeV.
 ***************************************************************************/
void GEnergy::log10TeV(const double& eng)
{
    // Set energy
    log10MeV(eng+6.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy
 ***************************************************************************/
std::string GEnergy::print(void) const
{
    // Initialise result string
    std::string result;

    // Append energy
    if (GeV() > 1000.0)
        result.append(str(TeV())+" TeV");
    else if (MeV() > 1000.0)
        result.append(str(GeV())+" GeV");
    else if (keV() > 1000.0)
        result.append(str(MeV())+" MeV");
    else
        result.append(str(keV())+" keV");

    // Return
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
void GEnergy::init_members(void)
{
    // Initialise members
    m_energy    = 0.0;
    m_elog10    = DBL_MIN;
    m_has_log10 = false;
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] eng Energy.
 ***************************************************************************/
void GEnergy::copy_members(const GEnergy& eng)
{
    // Copy time
    m_energy    = eng.m_energy;
    m_elog10    = eng.m_elog10;
    m_has_log10 = eng.m_has_log10;
    
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
