/***************************************************************************
 *              GLATEventBin.cpp - Fermi/LAT event bin class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GLATEventBin.cpp
 * @brief Fermi/LAT event bin class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <cmath>
#include "GLATEventBin.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DIR                                           "GLATEventBin::dir()"
#define G_ENERGY                                     "GLATEventBin::energy()"
#define G_TIME                                         "GLATEventBin::time()"
#define G_COUNTS_GET                                 "GLATEventBin::counts()"
#define G_COUNTS_SET                          "GLATEventBin::counts(double&)"
#define G_SOLIDANGLE                             "GLATEventBin::solidangle()"
#define G_EWIDTH                                     "GLATEventBin::ewidth()"
#define G_ONTIME                                     "GLATEventBin::ontime()"

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
GLATEventBin::GLATEventBin(void) : GEventBin()
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
GLATEventBin::GLATEventBin(const GLATEventBin& bin) : GEventBin(bin)
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
GLATEventBin::~GLATEventBin(void)
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
 * @param[in] bin LAT event bin.
 * @return Lat event bin.
 ***************************************************************************/
GLATEventBin& GLATEventBin::operator=(const GLATEventBin& bin)
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
 * @brief Clear event bin
 ***************************************************************************/
void GLATEventBin::clear(void)
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
 * @brief Clone event bin
 *
 * @return Pointer to deep copy of Fermi/LAT event bin.
 ***************************************************************************/
GLATEventBin* GLATEventBin::clone(void) const
{
    return new GLATEventBin(*this);
}


/***********************************************************************//**
 * @brief Return size of event bin
 *
 * @return Size of event bin in units of sr MeV s.
 *
 * The size of the event bin (units sr MeV s) is given by
 * \f[size = \Omega \times \Delta E \times \Delta T\f]
 * where
 * \f$\Omega\f$ is the size of the spatial bin in sr,
 * \f$\Delta E\f$ is the size of the energy bin in MeV, and
 * \f$\Delta T\f$ is the ontime of the observation in seconds. 
 ***************************************************************************/
double GLATEventBin::size(void) const
{
    // Compute bin size
    double size = solidangle() * ewidth().MeV() * ontime();

    // Return bin size
    return size;
}


/***********************************************************************//**
 * @brief Return instrument direction of event bin
 *
 * @return Instrument direction of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid instrument direction pointer.
 *
 * Returns reference to the instrument direction of the event bin.
 ***************************************************************************/
const GLATInstDir& GLATEventBin::dir(void) const
{
    // Throw an exception if instrument direction pointer is not valid
    if (m_dir == NULL) {
        throw GLATException::no_member(G_DIR,
                                       "Invalid instrument direction pointer.");
    }

    // Return instrument direction
    return *m_dir;
}


/***********************************************************************//**
 * @brief Return energy of event bin
 *
 * @return Energy of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid energy pointer.
 *
 * Returns reference to the energy of the event bin.
 ***************************************************************************/
const GEnergy& GLATEventBin::energy(void) const
{
    // Throw an exception if energy pointer is not valid
    if (m_energy == NULL) {
        throw GLATException::no_member(G_ENERGY,
                                       "Invalid energy pointer.");
    }

    // Return energy
    return *m_energy;
}


/***********************************************************************//**
 * @brief Return time of event bin
 *
 * @return Time of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid time pointer.
 *
 * Returns reference to the time of the event bin.
 ***************************************************************************/
const GTime& GLATEventBin::time(void) const
{
    // Throw an exception if time pointer is not valid
    if (m_time == NULL) {
        throw GLATException::no_member(G_TIME,
                                       "Invalid time pointer.");
    }

    // Return time
    return *m_time;
}


/***********************************************************************//**
 * @brief Return number of counts in event bin
 *
 * @return Number of counts in event bin.
 *
 * @exception GLATException::no_member
 *            Invalid counts pointer.
 *
 * Returns reference to the number of counts in the event bin.
 ***************************************************************************/
double GLATEventBin::counts(void) const
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        throw GLATException::no_member(G_COUNTS_GET,
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
 * @exception GLATException::no_member
 *            Invalid counts pointer.
 *
 * Set the number of counts in the event bin.
 ***************************************************************************/
void GLATEventBin::counts(const double& counts)
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        throw GLATException::no_member(G_COUNTS_SET,
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
 * @return Error in number of counts in event bin.
 *
 * Returns \f$\sqrt(counts+delta)\f$ as the uncertainty in the number of
 * counts in the bin. Adding delta avoids uncertainties of 0 which will
 * lead in the optimisation step to the exlusion of the corresponding bin.
 * In the actual implementation delta=1e-50.
 *
 * @todo The choice of delta has been made somewhat arbitrary, mainly
 * because the optimizer routines filter error^2 below 1e-100.
 ***************************************************************************/
double GLATEventBin::error(void) const
{
    // Compute uncertainty
    double error = sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Return solid angle of event bin
 *
 * @return Solid angle of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid solid angle pointer.
 *
 * Returns reference to the solid angle of the event bin.
 ***************************************************************************/
const double& GLATEventBin::solidangle(void) const
{
    // Throw an exception if solid angle pointer is not valid
    if (m_solidangle == NULL) {
        throw GLATException::no_member(G_SOLIDANGLE,
                                       "Invalid solid angle pointer.");
    }

    // Return solid angle
    return *m_solidangle;
}


/***********************************************************************//**
 * @brief Return energy width of event bin
 *
 * @return Energy width of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid energy width pointer.
 *
 * Returns reference to the energy width of the event bin.
 ***************************************************************************/
const GEnergy& GLATEventBin::ewidth(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ewidth == NULL) {
        throw GLATException::no_member(G_EWIDTH,
                                       "Invalid energy width pointer.");
    }

    // Return energy width
    return *m_ewidth;
}


/***********************************************************************//**
 * @brief Return ontime of event bin
 *
 * @return Ontime of event bin.
 *
 * @exception GLATException::no_member
 *            Invalid ontime pointer.
 *
 * Returns reference to the ontime of the event bin.
 ***************************************************************************/
const double& GLATEventBin::ontime(void) const
{
    // Throw an exception if ontime pointer is not valid
    if (m_ontime == NULL) {
        throw GLATException::no_member(G_ONTIME,
                                       "Invalid ontime pointer.");
    }

    // Return ontime
    return *m_ontime;
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing number of counts in event bin.
 ***************************************************************************/
std::string GLATEventBin::print(const GChatter& chatter) const
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
void GLATEventBin::init_members(void)
{
    // Initialise members
    m_cube       = NULL;
    m_index      = -1;
    m_ipix       = -1;
    m_ieng       = -1;
    m_energy     = NULL;
    m_dir        = NULL;
    m_time       = NULL;
    m_counts     = NULL;
    m_solidangle = NULL;
    m_ewidth     = NULL;
    m_ontime     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
void GLATEventBin::copy_members(const GLATEventBin& bin)
{
    // Copy members
    m_cube       = bin.m_cube;
    m_index      = bin.m_index;
    m_ipix       = bin.m_ipix;
    m_ieng       = bin.m_ieng;
    m_energy     = bin.m_energy;
    m_dir        = bin.m_dir;
    m_time       = bin.m_time;
    m_counts     = bin.m_counts;
    m_solidangle = bin.m_solidangle;
    m_ewidth     = bin.m_ewidth;
    m_ontime     = bin.m_ontime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventBin::free_members(void)
{
    // Return
    return;
}
