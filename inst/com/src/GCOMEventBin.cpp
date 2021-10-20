/***************************************************************************
 *                GCOMEventBin.cpp - COMPTEL event bin class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMEventBin.cpp
 * @brief COMPTEL event bin class implementation.
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <cmath>
#include "GCOMEventBin.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DIR_GET                                       "GCOMEventBin::dir()"
#define G_ENERGY_GET                                 "GCOMEventBin::energy()"
#define G_TIME_GET                                     "GCOMEventBin::time()"
#define G_POLARIZATION_GET                     "GCOMEventBin::polarization()"
#define G_COUNTS_GET                                 "GCOMEventBin::counts()"
#define G_SOLIDANGLE_GET                         "GCOMEventBin::solidangle()"
#define G_EWIDTH_GET                                 "GCOMEventBin::ewidth()"
#define G_ONTIME_GET                                 "GCOMEventBin::ontime()"
#define G_DIR_SET                           "GCOMEventBin::dir(GCOMInstDir&)"
#define G_ENERGY_SET                         "GCOMEventBin::energy(GEnergy&)"
#define G_TIME_SET                               "GCOMEventBin::time(GTime&)"
#define G_POLARIZATION_SET       "GCOMEventBin::polarization(GPolarization&)"
#define G_COUNTS_SET                          "GCOMEventBin::counts(double&)"
#define G_SOLIDANGLE_SET                  "GCOMEventBin::solidangle(double&)"
#define G_EWIDTH_SET                         "GCOMEventBin::ewidth(GEnergy&)"
#define G_ONTIME_SET                          "GCOMEventBin::ontime(double&)"

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
GCOMEventBin::GCOMEventBin(void) : GEventBin()
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
GCOMEventBin::GCOMEventBin(const GCOMEventBin& bin) : GEventBin(bin)
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
GCOMEventBin::~GCOMEventBin(void)
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
 * @return Event bin.
 ***************************************************************************/
GCOMEventBin& GCOMEventBin::operator= (const GCOMEventBin& bin)
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
void GCOMEventBin::clear(void)
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
 *
 * @return Pointer to deep copy of event bin.
 ***************************************************************************/
GCOMEventBin* GCOMEventBin::clone(void) const
{
    return new GCOMEventBin(*this);
}


/***********************************************************************//**
 * @brief Return size of event bin
 *
 * @return Size of event bin (sr MeV s)
 *
 * The size of the event bin (units: sr MeV s) is given by
 * \f[size = \Omega \times \Delta E \times \Delta T\f]
 * where
 * \f$\Omega\f$ is the size of the spatial bin in sr,
 * \f$\Delta E\f$ is the size of the energy bin in MeV, and
 * \f$\Delta T\f$ is the ontime of the observation in seconds. 
 ***************************************************************************/
double GCOMEventBin::size(void) const
{
    // Compute bin size
    double size = solidangle() * ewidth().MeV() * ontime();

    // Return bin size
    return size;
}


/***********************************************************************//**
 * @brief Return instrument direction of event bin
 *
 * @return Instrument direction of event bin
 *
 * @exception GException::invalid_value
 *            Invalid instrument direction pointer.
 *
 * Returns reference to the instrument direction of the event bin.
 ***************************************************************************/
const GCOMInstDir& GCOMEventBin::dir(void) const
{
    // Throw an exception if instrument direction pointer is not valid
    if (m_dir == NULL) {
        std::string msg = "No valid instrument direction found in event bin";
        throw GException::invalid_value(G_DIR_GET, msg);
    }

    // Return instrument direction
    return *m_dir;
}


/***********************************************************************//**
 * @brief Return energy of event bin
 *
 * @return Energy of event bin
 *
 * @exception GException::invalid_value
 *            Invalid energy pointer.
 *
 * Returns reference to the energy of the event bin.
 ***************************************************************************/
const GEnergy& GCOMEventBin::energy(void) const
{
    // Throw an exception if energy pointer is not valid
    if (m_energy == NULL) {
        std::string msg = "No valid energy found in event bin";
        throw GException::invalid_value(G_ENERGY_GET, msg);
    }

    // Return energy
    return *m_energy;
}


/***********************************************************************//**
 * @brief Return time of event bin
 *
 * @return Time of event bin
 *
 * @exception GException::invalid_value
 *            Invalid time pointer.
 *
 * Returns reference to the time of the event bin.
 ***************************************************************************/
const GTime& GCOMEventBin::time(void) const
{
    // Throw an exception if time pointer is not valid
    if (m_time == NULL) {
        std::string msg = "No valid time found in event bin";
        throw GException::invalid_value(G_TIME_GET, msg);
    }

    // Return time
    return *m_time;
}


/***********************************************************************//**
 * @brief Return polarization of event bin
 *
 * @return Polarization of event bin
 *
 * @exception GException::invalid_value
 *            Invalid polarization pointer.
 *
 * Returns reference to the polarization of the event bin.
 ***************************************************************************/
const GPolarization& GCOMEventBin::polarization(void) const
{
    // Throw an exception if polarization pointer is not valid
    if (m_polarization == NULL) {
        std::string msg = "No valid polarization found in event bin";
        throw GException::invalid_value(G_POLARIZATION_GET, msg);
    }

    // Return polarization
    return *m_polarization;
}


/***********************************************************************//**
 * @brief Return number of counts in event bin
 *
 * @return Number of counts in event bin
 *
 * @exception GCTAException::invalid_value
 *            Invalid counts pointer.
 *
 * Returns reference to the number of counts in the event bin.
 ***************************************************************************/
double GCOMEventBin::counts(void) const
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        std::string msg = "No valid counts found in event bin";
        throw GException::invalid_value(G_COUNTS_GET, msg);
    }

    // Return counts
    return *m_counts;
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
double GCOMEventBin::error(void) const
{
    // Compute uncertainty
    double error = std::sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Return solid angle of event bin
 *
 * @return Solid angle of event bin
 *
 * @exception GException::invalid_value
 *            Invalid solid angle pointer.
 *
 * Returns reference to the solid angle of the event bin.
 ***************************************************************************/
const double& GCOMEventBin::solidangle(void) const
{
    // Throw an exception if counts pointer is not valid
    if (m_solidangle == NULL) {
        std::string msg = "No valid solid angle found in event bin";
        throw GException::invalid_value(G_SOLIDANGLE_GET, msg);
    }

    // Return solid angle
    return *m_solidangle;
}


/***********************************************************************//**
 * @brief Return energy width of event bin
 *
 * @return Energy width of event bin
 *
 * @exception GException::invalid_value
 *            Invalid energy width pointer.
 *
 * Returns reference to the energy width of the event bin.
 ***************************************************************************/
const GEnergy& GCOMEventBin::ewidth(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ewidth == NULL) {
        std::string msg = "No valid energy width found in event bin";
        throw GException::invalid_value(G_EWIDTH_GET, msg);
    }

    // Return energy width
    return *m_ewidth;
}


/***********************************************************************//**
 * @brief Return ontime of event bin
 *
 * @return Ontime of event bin
 *
 * @exception GException::invalid_value
 *            Invalid ontime pointer.
 *
 * Returns reference to the ontime of the event bin.
 ***************************************************************************/
const double& GCOMEventBin::ontime(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ontime == NULL) {
        std::string msg = "No valid ontime found in event bin";
        throw GException::invalid_value(G_ONTIME_GET, msg);
    }

    // Return ontime
    return *m_ontime;
}


/***********************************************************************//**
 * @brief Set instrument direction of event bin
 *
 * @param[in] dir Instrument direction of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold instrument direction.
 *
 * Sets the instrument direction of the event bin.
 ***************************************************************************/
void GCOMEventBin::dir(const GCOMInstDir& dir)
{
    // Throw an exception if no memory has been allocated
    if (m_dir == NULL) {
        std::string msg = "No memory available to hold instrument direction.";
        throw GException::invalid_value(G_DIR_SET, msg);
    }

    // Set instrument direction
    *m_dir = dir;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy of event bin
 *
 * @param[in] energy Energy of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold energy.
 *
 * Sets the energy of the event bin.
 ***************************************************************************/
void GCOMEventBin::energy(const GEnergy& energy)
{
    // Throw an exception if no memory has been allocated
    if (m_energy == NULL) {
        std::string msg = "No memory available to hold energy.";
        throw GException::invalid_value(G_ENERGY_SET, msg);
    }

    // Set energy
    *m_energy = energy;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set time of event bin
 *
 * @param[in] time Time of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold time.
 *
 * Sets the time of the event bin.
 ***************************************************************************/
void GCOMEventBin::time(const GTime& time)
{
    // Throw an exception if no memory has been allocated
    if (m_time == NULL) {
        std::string msg = "No memory available to hold time.";
        throw GException::invalid_value(G_TIME_SET, msg);
    }

    // Set time
    *m_time = time;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set polarization of event bin
 *
 * @param[in] polarization Polarization of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold polarization.
 *
 * Sets the polarization of the event bin.
 ***************************************************************************/
void GCOMEventBin::polarization(const GPolarization& polarization)
{
    // Throw an exception if no memory has been allocated
    if (m_polarization == NULL) {
        std::string msg = "No memory available to hold polarization.";
        throw GException::invalid_value(G_POLARIZATION_SET, msg);
    }

    // Set polarization
    *m_polarization = polarization;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set number of counts in event bin
 *
 * @param[in] counts Number of counts.
 *
 * @exception GException::invalid_value
 *            No memory available to hold counts.
 *
 * Set the number of counts in the event bin.
 ***************************************************************************/
void GCOMEventBin::counts(const double& counts)
{
    // Throw an exception if counts pointer is not valid
    if (m_counts == NULL) {
        std::string msg = "No memory available to hold counts.";
        throw GException::invalid_value(G_COUNTS_SET, msg);
    }

    // Set number of counts in event bin
    *m_counts = counts;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set solid angle of event bin
 *
 * @param[in] solidangle Solid angle of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold solid angle.
 *
 * Sets the solid angle of the event bin.
 ***************************************************************************/
void GCOMEventBin::solidangle(const double& solidangle)
{
    // Throw an exception if no memory has been allocated
    if (m_solidangle == NULL) {
        std::string msg = "No memory available to hold solid angle.";
        throw GException::invalid_value(G_SOLIDANGLE_SET, msg);
    }

    // Set solid angle
    *m_solidangle = solidangle;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set energy width of event bin
 *
 * @param[in] ewidth Energy width of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold energy width.
 *
 * Sets the energy width of the event bin.
 ***************************************************************************/
void GCOMEventBin::ewidth(const GEnergy& ewidth)
{
    // Throw an exception if no memory has been allocated
    if (m_ewidth == NULL) {
        std::string msg = "No memory available to hold energy width.";
        throw GException::invalid_value(G_EWIDTH_SET, msg);
    }

    // Set energy width
    *m_ewidth = ewidth;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set ontime of event bin
 *
 * @param[in] ontime Ontime of event bin (sec).
 *
 * @exception GException::invalid_value
 *            No memory available to hold ontime.
 *
 * Sets the ontime of the event bin.
 ***************************************************************************/
void GCOMEventBin::ontime(const double& ontime)
{
    // Throw an exception if no memory has been allocated
    if (m_ontime == NULL) {
        std::string msg = "No memory available to hold ontime.";
        throw GException::invalid_value(G_ONTIME_SET, msg);
    }

    // Set solid angle
    *m_ontime = ontime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event information.
 ***************************************************************************/
std::string GCOMEventBin::print(const GChatter& chatter) const
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
 *
 * This method allocates memory for all event bin attributes and intialises
 * the attributes to well defined initial values.
 * 
 * The method assumes that on entry no memory is hold by the member pointers.
 ***************************************************************************/
void GCOMEventBin::init_members(void)
{
    // Allocate members
    m_alloc        = true;
    m_index        = -1;   // Signals that event bin does not correspond to cube
    m_dir          = new GCOMInstDir;
    m_time         = new GTime;
    m_polarization = new GPolarization;
    m_energy       = new GEnergy;
    m_ewidth       = new GEnergy;
    m_counts       = new double;
    m_solidangle   = new double;
    m_ontime       = new double;

    // Initialise members
    m_dir->clear();
    m_time->clear();
    m_polarization->clear();
    m_energy->clear();
    m_ewidth->clear();
    *m_counts     = 0.0;
    *m_solidangle = 0.0;
    *m_ontime     = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
void GCOMEventBin::copy_members(const GCOMEventBin& bin)
{
    // First de-allocate existing memory if needed
    free_members();

    // Copy members by cloning
    m_dir          = new GCOMInstDir(*bin.m_dir);
    m_time         = new GTime(*bin.m_time);
    m_polarization = new GPolarization(*bin.m_polarization);
    m_energy       = new GEnergy(*bin.m_energy);
    m_ewidth       = new GEnergy(*bin.m_ewidth);
    m_counts       = new double(*bin.m_counts);
    m_solidangle   = new double(*bin.m_solidangle);
    m_ontime       = new double(*bin.m_ontime);

    // Copy non-pointer members
    m_index = bin.m_index;

    // Signal memory allocation
    m_alloc = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 *
 * This method frees all memory of the class attributes and sets the member
 * pointers to NULL. This method should only be called if new memory is
 * allocated immediately afterwards (for example by cloning another event
 * bin), or upon destruction of the object.
 *
 * Note that some logic has been implemented that frees only memory that also
 * has indeed been allocated by the class. Thus if the class only serves as
 * container to hold memory pointer allocated by someone else (for example
 * the GCOMEventCube class), no memory is freed.
 ***************************************************************************/
void GCOMEventBin::free_members(void)
{
    // If memory was allocated then free members now
    if (m_alloc) {
        if (m_dir        != NULL) delete m_dir;
        if (m_time       != NULL) delete m_time;
        if (m_energy     != NULL) delete m_energy;
        if (m_ewidth     != NULL) delete m_ewidth;
        if (m_counts     != NULL) delete m_counts;
        if (m_solidangle != NULL) delete m_solidangle;
        if (m_ontime     != NULL) delete m_ontime;
    }

    // Signal member pointers as free
    m_dir        = NULL;
    m_time       = NULL;
    m_energy     = NULL;
    m_ewidth     = NULL;
    m_counts     = NULL;
    m_solidangle = NULL;
    m_ontime     = NULL;

    // Signal memory de-allocation
    m_alloc = false;

    // Return
    return;
}
