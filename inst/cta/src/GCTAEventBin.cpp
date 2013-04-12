/***************************************************************************
 *                 GCTAEventBin.cpp - CTA event bin class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAEventBin.cpp
 * @brief CTA event bin class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <cmath>
#include "GCTAEventBin.hpp"
#include "GCTAException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DIR                                           "GCTAEventBin::dir()"
#define G_ENERGY                                     "GCTAEventBin::energy()"
#define G_TIME                                         "GCTAEventBin::time()"
#define G_COUNTS_GET                                 "GCTAEventBin::counts()"
#define G_COUNTS_SET                          "GCTAEventBin::counts(double&)"
#define G_OMEGA                                       "GCTAEventBin::omega()"
#define G_EWIDTH                                     "GCTAEventBin::ewidth()"
#define G_ONTIME                                     "GCTAEventBin::ontime()"

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
GCTAEventBin::GCTAEventBin(void) : GEventBin()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
GCTAEventBin::GCTAEventBin(const GCTAEventBin& bin) : GEventBin(bin)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(bin);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAEventBin::~GCTAEventBin(void)
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
 * @param[in] bin Event bin.
 ***************************************************************************/
GCTAEventBin& GCTAEventBin::operator= (const GCTAEventBin& bin)
{
    // Execute only if object is not identical
    if (this != &bin) {

        // Copy base class members
        this->GEventBin::operator=(bin);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(bin);

    } // endif: object was not identical

    // Return this object
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
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GCTAEventBin::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventBin::free_members();
    this->GEvent::free_members();

    // Initialise members
    this->GEvent::init_members();
    this->GEventBin::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAEventBin* GCTAEventBin::clone(void) const
{
    return new GCTAEventBin(*this);
}


/***********************************************************************//**
 * @brief Return size of event bin
 *
 * The size of the event bin (units sr MeV s) is given by
 * \f[size = \Omega \times \Delta E \times \Delta T\f]
 * where
 * \f$\Omega\f$ is the size of the spatial bin in sr,
 * \f$\Delta E\f$ is the size of the energy bin in MeV, and
 * \f$\Delta T\f$ is the ontime of the observation in seconds. 
 ***************************************************************************/
double GCTAEventBin::size(void) const
{
    // Compute bin size
    double size = omega() * ewidth().MeV() * ontime();

    // Return bin size
    return size;
}


/***********************************************************************//**
 * @brief Return instrument direction of event bin
 *
 * @return Instrument direction of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid instrument direction pointer.
 *
 * Returns reference to the instrument direction of the event bin.
 ***************************************************************************/
const GCTAInstDir& GCTAEventBin::dir(void) const
{
    // Throw an exception if instrument direction pointer is not valid
    if (m_dir == NULL) {
        throw GCTAException::no_member(G_DIR,
                                       "Invalid instrument direction pointer.");
    }

    // Return instrument direction
    return *m_dir;
}


/***********************************************************************//**
 * @brief Return energy of event bin
 *
 * @return Energy of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid energy pointer.
 *
 * Returns reference to the energy of the event bin.
 ***************************************************************************/
const GEnergy& GCTAEventBin::energy(void) const
{
    // Throw an exception if energy pointer is not valid
    if (m_energy == NULL) {
        throw GCTAException::no_member(G_ENERGY,
                                       "Invalid energy pointer.");
    }

    // Return energy
    return *m_energy;
}


/***********************************************************************//**
 * @brief Return time of event bin
 *
 * @return Time of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid time pointer.
 *
 * Returns reference to the time of the event bin.
 ***************************************************************************/
const GTime& GCTAEventBin::time(void) const
{
    // Throw an exception if time pointer is not valid
    if (m_time == NULL) {
        throw GCTAException::no_member(G_TIME,
                                       "Invalid time pointer.");
    }

    // Return time
    return *m_time;
}


/***********************************************************************//**
 * @brief Return number of counts in event bin
 *
 * @return Number of counts in event bin
 *
 * @exception GCTAException::no_member
 *            Invalid counts pointer.
 *
 * Returns reference to the number of counts in the event bin.
 ***************************************************************************/
double GCTAEventBin::counts(void) const
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        throw GCTAException::no_member(G_COUNTS_GET,
                                       "Invalid counts pointer.");
    }

    // Return counts
    return *m_counts;
}


/***********************************************************************//**
 * @brief Set number of counts in event bin
 *
 * @param[in] counts Number of counts.
 *
 * @exception GCTAException::no_member
 *            Invalid counts pointer.
 *
 * Set the number of counts in the event bin.
 ***************************************************************************/
void GCTAEventBin::counts(const double& counts)
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        throw GCTAException::no_member(G_COUNTS_SET,
                                       "Invalid counts pointer.");
    }

    // Set number of counts in event bin
    *m_counts = counts;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return error in number of counts
 *
 * @return Error in number of counts in event bin
 *
 * Returns \f$\sqrt(counts+delta)\f$ as the uncertainty in the number of
 * counts in the bin. Adding delta avoids uncertainties of 0 which will
 * lead in the optimisation step to the exlusion of the corresponding bin.
 * In the actual implementation delta=1e-50.
 *
 * @todo The choice of delta has been made somewhat arbitrary, mainly
 * because the optimizer routines filter error^2 below 1e-100.
 ***************************************************************************/
double GCTAEventBin::error(void) const
{
    // Compute uncertainty
    double error = sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Return solid angle of event bin
 *
 * @return Solid angle of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid solid angle pointer.
 *
 * Returns reference to the solid angle of the event bin.
 ***************************************************************************/
const double& GCTAEventBin::omega(void) const
{
    // Throw an exception if solid angle pointer is not valid
    if (m_omega == NULL) {
        throw GCTAException::no_member(G_OMEGA,
                                       "Invalid solid angle pointer.");
    }

    // Return solid angle
    return *m_omega;
}


/***********************************************************************//**
 * @brief Return energy width of event bin
 *
 * @return Energy width of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid energy width pointer.
 *
 * Returns reference to the energy width of the event bin.
 ***************************************************************************/
const GEnergy& GCTAEventBin::ewidth(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ewidth == NULL) {
        throw GCTAException::no_member(G_EWIDTH,
                                       "Invalid energy width pointer.");
    }

    // Return energy width
    return *m_ewidth;
}


/***********************************************************************//**
 * @brief Return ontime of event bin
 *
 * @return Ontime of event bin
 *
 * @exception GCTAException::no_member
 *            Invalid ontime pointer.
 *
 * Returns reference to the ontime of the event bin.
 ***************************************************************************/
const double& GCTAEventBin::ontime(void) const
{
    // Throw an exception if ontime pointer is not valid
    if (m_ontime == NULL) {
        throw GCTAException::no_member(G_ONTIME,
                                       "Invalid ontime pointer.");
    }

    // Return ontime
    return *m_ontime;
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event information.
 ***************************************************************************/
std::string GCTAEventBin::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append number of counts
        result.append(gammalib::str(counts()));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAEventBin::init_members(void)
{
    // Initialise members
    m_energy = NULL;
    m_dir    = NULL;
    m_time   = NULL;
    m_counts = NULL;
    m_omega  = NULL;
    m_ewidth = NULL;
    m_ontime = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
void GCTAEventBin::copy_members(const GCTAEventBin& bin)
{
    // Copy members
    m_energy = bin.m_energy;
    m_dir    = bin.m_dir;
    m_time   = bin.m_time;
    m_counts = bin.m_counts;
    m_omega  = bin.m_omega;
    m_ewidth = bin.m_ewidth;
    m_ontime = bin.m_ontime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventBin::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
