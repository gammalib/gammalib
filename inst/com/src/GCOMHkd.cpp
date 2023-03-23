/***************************************************************************
 *          GCOMHkd.cpp - COMPTEL Housekeeping Data container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMHkd.cpp
 * @brief COMPTEL Housekeeping Data container class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GCOMTools.hpp"
#include "GCOMHkd.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_REMOVE                                      "GCOMHkd::remove(int&)"

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
GCOMHkd::GCOMHkd(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Name of Housekeeping parameter
 ***************************************************************************/
GCOMHkd::GCOMHkd(const std::string& name)
{
    // Initialise class members
    init_members();

    // Set name
    this->name(name);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] hkd Housekeeping Data container.
 ***************************************************************************/
GCOMHkd::GCOMHkd(const GCOMHkd& hkd)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(hkd);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMHkd::~GCOMHkd(void)
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
 * @param[in] hkd Housekeeping Data container.
 * @return Housekeeping Data container.
 ***************************************************************************/
GCOMHkd& GCOMHkd::operator=(const GCOMHkd& hkd)
{
    // Execute only if object is not identical
    if (this != &hkd) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(hkd);

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
 * @brief Clear Housekeeping Data container
 ***************************************************************************/
void GCOMHkd::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Housekeeping Data container
 *
 * @return Pointer to deep copy of Housekeeping Data container.
 ***************************************************************************/
GCOMHkd* GCOMHkd::clone(void) const
{
    return new GCOMHkd(*this);
}


/***********************************************************************//**
 * @brief Append Housekeeping Data to container
 *
 * @param[in] time Time stamp of Housekeeping Data.
 * @param[in] value Value of Housekeeping Data.
 *
 * Appends time stamp and corresponding value to Housekeeping Data container.
 ***************************************************************************/
void GCOMHkd::append(const GTime& time, const double& value)
{
    // Append time and value
    m_times.append(time);
    m_values.push_back(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove Housekeeping Data
 *
 * @param[in] index Housekeeping Data index [0,...,size()[.
 *
 * Removes Housekeeping Data at @p index from the container. All data after
 * the specified @p index are moved forward by one position.
 ***************************************************************************/
void GCOMHkd::remove(const int& index)
{
    #if defined(G_RANGE_CHECK)
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Housekeeping Data index",
                                       index, size());
    }
    #endif

    // Remove time and erase value
    m_times.remove(index);
    m_values.erase(m_values.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extend Housekeeping Data
 *
 * @param[in] hkd Housekeeping Data container.
 *
 * Append Housekeeping Data container to existing container.
 ***************************************************************************/
void GCOMHkd::extend(const GCOMHkd& hkd)
{
    // Do nothing if Housekeeping Data container is empty
    if (!hkd.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = hkd.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_times.append(hkd.m_times[i]);
            m_values.push_back(hkd.m_values[i]);
        }

    } // endif: Housekeeping Data container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Housekeeping Data container
 *
 * @param[in] chatter Chattiness.
 * @return String containing Housekeeping Data container information.
 ***************************************************************************/
std::string GCOMHkd::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMHkd ===");

        // Append information
        result.append("\n"+gammalib::parformat("Parameter"));
        result.append(m_name);
        result.append("\n"+gammalib::parformat("Records"));
        result.append(gammalib::str(size()));

        // If there are records then print time range
        if (size() > 0) {

            // Get start and stop time
            GTime tstart = m_times[0];
            GTime tstop  = m_times[size()-1];

            // Print time intervals
            result.append("\n"+gammalib::parformat("TJD interval"));
            result.append(gammalib::str(gammalib::com_tjd(tstart))+":");
            result.append(gammalib::str(gammalib::com_tics(tstart))+" - ");
            result.append(gammalib::str(gammalib::com_tjd(tstop))+":");
            result.append(gammalib::str(gammalib::com_tics(tstop)));
            result.append("\n"+gammalib::parformat("MJD interval"));
            result.append(gammalib::str(tstart.mjd())+" - ");
            result.append(gammalib::str(tstop.mjd())+" days");
            result.append("\n"+gammalib::parformat("UTC interval"));
            result.append(tstart.utc()+" - ");
            result.append(tstop.utc());

        }

    } // endif: chatter was not silent

    // Return result
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
void GCOMHkd::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_times.clear();
    m_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] hkd Housekeeping Data container.
 ***************************************************************************/
void GCOMHkd::copy_members(const GCOMHkd& hkd)
{
    // Copy members
    m_name   = hkd.m_name;
    m_times  = hkd.m_times;
    m_values = hkd.m_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMHkd::free_members(void)
{
    // Return
    return;
}
