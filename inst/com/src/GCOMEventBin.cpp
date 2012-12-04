/***************************************************************************
 *               GCOMEventBin.cpp  -  COMPTEL event bin class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
    double size = omega() * ewidth().MeV() * ontime();

    // Return bin size
    return size;
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
    double error = sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @return String containing event information.
 ***************************************************************************/
std::string GCOMEventBin::print(void) const
{
    // Initialise result string
    std::string result;

    // Append number of counts
    result.append(str(counts()));

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
    m_alloc  = true;
    m_index  = -1;   // Signals that event bin does not correspond to cube
    m_counts = new double;
    m_dir    = new GCOMInstDir;
    m_omega  = new double;
    m_time   = new GTime;
    m_ontime = new double;
    m_energy = new GEnergy;
    m_ewidth = new GEnergy;

    // Initialise members
    *m_counts = 0.0;
    m_dir->clear();
    *m_omega = 0.0;
    m_time->clear();
    *m_ontime = 0.0;
    m_energy->clear();
    m_ewidth->clear();

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
    // Copy members by cloning
    m_alloc  = true;
    m_index  = bin.m_index;
    m_counts = new double(*bin.m_counts);
    m_dir    = new GCOMInstDir(*bin.m_dir);
    m_omega  = new double(*bin.m_omega);
    m_time   = new GTime(*bin.m_time);
    m_ontime = new double(*bin.m_ontime);
    m_energy = new GEnergy(*bin.m_energy);
    m_ewidth = new GEnergy(*bin.m_ewidth);

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
        if (m_counts != NULL) delete m_counts;
        if (m_dir    != NULL) delete m_dir;
        if (m_omega  != NULL) delete m_omega;
        if (m_time   != NULL) delete m_time;
        if (m_ontime != NULL) delete m_ontime;
        if (m_energy != NULL) delete m_energy;
        if (m_ewidth != NULL) delete m_ewidth;
    }

    // Signal member pointers as free
    m_alloc  = false;
    m_counts = NULL;
    m_dir    = NULL;
    m_omega  = NULL;
    m_time   = NULL;
    m_ontime = NULL;
    m_energy = NULL;
    m_ewidth = NULL;

    // Return
    return;
}
