/***************************************************************************
 *                        GPhoton.hpp - Photon class                       *
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
 * @file GPhoton.hpp
 * @brief Photon class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <cfloat>
//#include <cmath>
#include "GPhoton.hpp"
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
GPhoton::GPhoton(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] ph Photon.
 ***************************************************************************/
GPhoton::GPhoton(const GPhoton& ph)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(ph);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GPhoton::~GPhoton(void)
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
 * @param[in] ph Photon.
 ***************************************************************************/
GPhoton& GPhoton::operator= (const GPhoton& ph)
{ 
    // Execute only if object is not identical
    if (this != &ph) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(ph);

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
void GPhoton::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Print photon
 *
 * @todo Implement print() method for GSkyDir
 ***************************************************************************/
std::string GPhoton::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GPhoton ===");

    // Append arrival direction
    result.append("\n"+parformat("Arrival direction"));
    result.append("RA="+str(m_dir.ra_deg())+", DEC="+str(m_dir.dec_deg()));

    // Append energy and time
    result.append("\n"+parformat("Energy")+m_energy.print());
    result.append("\n"+parformat("Arrival time (MET)")+m_time.print());

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
void GPhoton::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_energy.clear();
    m_time.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] ph Photon.
 ***************************************************************************/
void GPhoton::copy_members(const GPhoton& ph)
{
    // Copy time
    m_dir    = ph.m_dir;
    m_energy = ph.m_energy;
    m_time   = ph.m_time;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GPhoton::free_members(void)
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
 * @param[in] ph Photon.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GPhoton& ph)
{
     // Write photon in output stream
    os << ph.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] ph Photon.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GPhoton& ph)
{
    // Write photon into logger
    log << ph.print();

    // Return logger
    return log;
}
