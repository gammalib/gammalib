/***************************************************************************
 *                 GCTAEventBin.cpp - CTA event bin class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
#define G_DIR_GET                                       "GCTAEventBin::dir()"
#define G_ENERGY                                     "GCTAEventBin::energy()"
#define G_TIME                                         "GCTAEventBin::time()"
#define G_COUNTS_GET                                 "GCTAEventBin::counts()"
#define G_SOLIDANGLE                             "GCTAEventBin::solidangle()"
#define G_EWIDTH                                     "GCTAEventBin::ewidth()"
#define G_ONTIME                                     "GCTAEventBin::ontime()"
#define G_WEIGHT                                     "GCTAEventBin::weight()"
#define G_DIR_SET                           "GCTAEventBin::dir(GCTAInstDir&)"
#define G_ENERGY_SET                         "GCTAEventBin::energy(GEnergy&)"
#define G_TIME_SET                               "GCTAEventBin::time(GTime&)"
#define G_COUNTS_SET                          "GCTAEventBin::counts(double&)"
#define G_SOLIDANGLE_SET                  "GCTAEventBin::solidangle(double&)"
#define G_EWIDTH_SET                         "GCTAEventBin::ewidth(GEnergy&)"
#define G_ONTIME_SET                          "GCTAEventBin::ontime(double&)"
#define G_WEIGHT_SET                          "GCTAEventBin::weight(double&)"

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
 * @return Event bin.
 ***************************************************************************/
GCTAEventBin& GCTAEventBin::operator=(const GCTAEventBin& bin)
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
 * @brief Clear eventbin
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
 * @brief Clone event bin
 *
 * @return Pointer to deep copy of event bin.
 ***************************************************************************/
GCTAEventBin* GCTAEventBin::clone(void) const
{
    return new GCTAEventBin(*this);
}


/***********************************************************************//**
 * @brief Return size of event bin
 *
 * The size of the event bin (units sr MeV s) is given by
 * \f[size = \Omega \times \Delta E \times \Delta T \times W\f]
 * where
 * \f$\Omega\f$ is the size of the spatial bin in sr,
 * \f$\Delta E\f$ is the size of the energy bin in MeV,
 * \f$\Delta T\f$ is the ontime of the observation in seconds, and
 * \f$W\f$ is the weight of the bin.
 ***************************************************************************/
double GCTAEventBin::size(void) const
{
    // Compute bin size
    double size = solidangle() * ewidth().MeV() * ontime() * weight();

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
const GCTAInstDir& GCTAEventBin::dir(void) const
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
const GEnergy& GCTAEventBin::energy(void) const
{
    // Throw an exception if energy pointer is not valid
    if (m_energy == NULL) {
        std::string msg = "No valid energy found in event bin";
        throw GException::invalid_value(G_ENERGY, msg);
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
const GTime& GCTAEventBin::time(void) const
{
    // Throw an exception if time pointer is not valid
    if (m_energy == NULL) {
        std::string msg = "No valid time found in event bin";
        throw GException::invalid_value(G_TIME, msg);
    }

    // Return time
    return *m_time;
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
double GCTAEventBin::counts(void) const
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
 * @exception GException::invalid_value
 *            Invalid solid angle pointer.
 *
 * Returns reference to the solid angle of the event bin.
 ***************************************************************************/
const double& GCTAEventBin::solidangle(void) const
{
    // Throw an exception if counts pointer is not valid
    if (m_solidangle == NULL) {
        std::string msg = "No valid solid angle found in event bin";
        throw GException::invalid_value(G_SOLIDANGLE, msg);
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
const GEnergy& GCTAEventBin::ewidth(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ewidth == NULL) {
        std::string msg = "No valid energy width found in event bin";
        throw GException::invalid_value(G_EWIDTH, msg);
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
const double& GCTAEventBin::ontime(void) const
{
    // Throw an exception if energy width pointer is not valid
    if (m_ontime == NULL) {
        std::string msg = "No valid ontime found in event bin";
        throw GException::invalid_value(G_ONTIME, msg);
    }

    // Return ontime
    return *m_ontime;
}


/***********************************************************************//**
 * @brief Return weight of event bin
 *
 * @return Weight of event bin
 *
 * @exception GException::invalid_value
 *            Invalid weight pointer.
 *
 * Returns reference to the weight of the event bin.
 ***************************************************************************/
const double& GCTAEventBin::weight(void) const
{
    // Throw an exception if weight pointer is not valid
    if (m_weight == NULL) {
        std::string msg = "No valid weight found in event bin";
        throw GException::invalid_value(G_WEIGHT, msg);
    }

    // Return weight
    return *m_weight;
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
void GCTAEventBin::dir(const GCTAInstDir& dir)
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
void GCTAEventBin::energy(const GEnergy& energy)
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
 *            No memory available to hold instrument direction.
 *
 * Sets the time of the event bin.
 ***************************************************************************/
void GCTAEventBin::time(const GTime& time)
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
 * @brief Set number of counts in event bin
 *
 * @param[in] counts Number of counts.
 *
 * @exception GException::invalid_value
 *            No memory available to hold counts.
 *
 * Set the number of counts in the event bin.
 ***************************************************************************/
void GCTAEventBin::counts(const double& counts)
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
void GCTAEventBin::solidangle(const double& solidangle)
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
void GCTAEventBin::ewidth(const GEnergy& ewidth)
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
void GCTAEventBin::ontime(const double& ontime)
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
 * @brief Set weight of event bin
 *
 * @param[in] weight Weight angle of event bin
 *
 * @exception GException::invalid_value
 *            No memory available to hold weight.
 *
 * Sets the weight of the event bin.
 ***************************************************************************/
void GCTAEventBin::weight(const double& weight)
{
    // Throw an exception if no memory has been allocated
    if (m_weight == NULL) {
        std::string msg = "No memory available to hold weight.";
        throw GException::invalid_value(G_WEIGHT_SET, msg);
    }

    // Set weight
    *m_weight = weight;

    // Return
    return;
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
 *
 * This method allocates memory for all event bin attributes and intialises
 * the attributes to well defined initial values.
 * 
 * The method assumes that on entry no memory is hold by the member pointers.
 ***************************************************************************/
void GCTAEventBin::init_members(void)
{
    // Initialise members
    m_alloc      = true;
    m_ipix       = -1;   //!< Not part of an event cube
    m_ieng       = -1;   //!< Not part of an event cube
    m_dir        = new GCTAInstDir;
    m_time       = new GTime;
    m_energy     = new GEnergy;
    m_ewidth     = new GEnergy;
    m_counts     = new double;
    m_solidangle = new double;
    m_ontime     = new double;
    m_weight     = new double;

    // Initialise members
    m_dir->clear();
    m_time->clear();
    m_energy->clear();
    m_ewidth->clear();
    *m_counts     = 0.0;
    *m_solidangle = 0.0;
    *m_ontime     = 0.0;
    *m_weight     = 0.0;

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
    // First de-allocate existing memory if needed
    free_members();

    // Copy members by cloning
    m_dir        = new GCTAInstDir(*bin.m_dir);
    m_time       = new GTime(*bin.m_time);
    m_energy     = new GEnergy(*bin.m_energy);
    m_ewidth     = new GEnergy(*bin.m_ewidth);
    m_counts     = new double(*bin.m_counts);
    m_solidangle = new double(*bin.m_solidangle);
    m_ontime     = new double(*bin.m_ontime);
    m_weight     = new double(*bin.m_weight);

    // Copy non-pointer members
    m_ipix = bin.m_ipix;
    m_ieng = bin.m_ieng;

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
 * the GCTAEventCube class), no memory is freed.
 ***************************************************************************/
void GCTAEventBin::free_members(void)
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
        if (m_weight     != NULL) delete m_weight;
    }

    // Signal member pointers as free
    m_dir        = NULL;
    m_time       = NULL;
    m_energy     = NULL;
    m_ewidth     = NULL;
    m_counts     = NULL;
    m_solidangle = NULL;
    m_ontime     = NULL;
    m_weight     = NULL;

    // Signal memory de-allocation
    m_alloc = false;

    // Return
    return;
}
