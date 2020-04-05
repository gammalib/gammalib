/***************************************************************************
 *             GSPIEventBin.cpp - INTEGRAL/SPI event bin class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIEventBin.hpp
 * @brief INTEGRAL/SPI event bin class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <string>
#include "GTools.hpp"
#include "GSPIEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_MODEL                                   "GSPIEventBin::model(int&)"

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
 *
 * Creates an empty INTEGRAL/SPI event bin.
 ***************************************************************************/
GSPIEventBin::GSPIEventBin(void) : GEventBin()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bin INTEGRAL/SPI event bin.
 ***************************************************************************/
GSPIEventBin::GSPIEventBin(const GSPIEventBin& bin) : GEventBin(bin)
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
GSPIEventBin::~GSPIEventBin(void)
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
 * @param[in] bin INTEGRAL/SPI event bin.
 * @return INTEGRAL/SPI event bin.
 ***************************************************************************/
GSPIEventBin& GSPIEventBin::operator=(const GSPIEventBin& bin)
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
 * @brief Clear INTEGRAL/SPI event bin
 *
 * Clears INTEGRAL/SPI event bin by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GSPIEventBin::clear(void)
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
 * @return Pointer to deep copy of INTEGRAL/SPI event bin.
 ***************************************************************************/
GSPIEventBin* GSPIEventBin::clone(void) const
{
    return new GSPIEventBin(*this);
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
double GSPIEventBin::error(void) const
{
    // Compute uncertainty
    double error = std::sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Return model value
 *
 * @param[in] index Model index.
 * @return Model value.
 *
 * @exception GException::out_of_range
 *            Invalid model index
 ***************************************************************************/
const double& GSPIEventBin::model(const int& index) const
{
    // Optionally check if the model index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_models) {
        throw GException::out_of_range(G_MODEL, "Invalid model index",
                                       index, m_num_models);
    }
    #endif

    // Return reference to model
    return (m_models[index]);
}


/***********************************************************************//**
 * @brief Print event information
 *
 * @param[in] chatter Chattiness.
 * @return String containing event information.
 ***************************************************************************/
std::string GSPIEventBin::print(const GChatter& chatter) const
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
void GSPIEventBin::init_members(void)
{
    // Allocate members
    m_alloc      = true;
    m_index      = -1;  // Signals that event bin does not correspond to cube
    m_idir       = -1;  // Signals that event bin does not correspond to cube
    m_iebin      = -1;  // Signals that event bin does not correspond to cube
    m_num_models =  0;  // No models in event bin
    m_dir        = new GSPIInstDir;
    m_time       = new GTime;
    m_energy     = new GEnergy;
    m_counts     = new double;
    m_ontime     = new double;
    m_size       = new double;
    m_models     = NULL;

    // Initialise members
    m_dir->clear();
    m_time->clear();
    m_energy->clear();
    *m_counts = 0.0;
    *m_ontime = 0.0;
    *m_size   = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin INTEGRAL/SPI event bin.
 ***************************************************************************/
void GSPIEventBin::copy_members(const GSPIEventBin& bin)
{
    // First de-allocate existing memory if needed
    free_members();

    // Copy non-pointer members
    m_index      = bin.m_index;
    m_idir       = bin.m_idir;
    m_iebin      = bin.m_iebin;
    m_num_models = bin.m_num_models;

    // Copy members by cloning
    m_dir    = new GSPIInstDir(*bin.m_dir);
    m_time   = new GTime(*bin.m_time);
    m_energy = new GEnergy(*bin.m_energy);
    m_counts = new double(*bin.m_counts);
    m_ontime = new double(*bin.m_ontime);
    m_size   = new double(*bin.m_size);
    
    // Copy models
    if (m_num_models > 0) {
        m_models = new double[m_num_models];
        for (int i = 0; i < m_num_models; ++i) {
            m_models[i] = bin.m_models[i];
        }
    }

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
 * the GSPIEventCube class), no memory is freed.
 ***************************************************************************/
void GSPIEventBin::free_members(void)
{
    // If memory was allocated then free members now
    if (m_alloc) {
        if (m_dir    != NULL) delete m_dir;
        if (m_time   != NULL) delete m_time;
        if (m_energy != NULL) delete m_energy;
        if (m_counts != NULL) delete m_counts;
        if (m_ontime != NULL) delete m_ontime;
        if (m_size   != NULL) delete m_size;
        if (m_models != NULL) delete [] m_models;
    }

    // Signal member pointers as free
    m_dir    = NULL;
    m_time   = NULL;
    m_energy = NULL;
    m_counts = NULL;
    m_ontime = NULL;
    m_size   = NULL;
    m_models = NULL;

    // Signal memory de-allocation
    m_alloc = false;

    // Return
    return;
}
