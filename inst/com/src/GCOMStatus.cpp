/***************************************************************************
 *             GCOMStatus.cpp - COMPTEL instrument status class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMStatus.cpp
 * @brief COMPTEL instrument status class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GCaldb.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GCOMStatus.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_D1STATUS                         "GCOMStatus::d1status(int&, int&)"
#define G_D2STATUS                         "GCOMStatus::d2status(int&, int&)"
#define G_UPDATE_CACHE                       "GCOMStatus::update_cache(int&)"
#define G_LOAD_STATUS                             "GCOMStatus::load_status()"

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
GCOMStatus::GCOMStatus(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] status COMPTEL instrument status.
 ***************************************************************************/
GCOMStatus::GCOMStatus(const GCOMStatus& status)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(status);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMStatus::~GCOMStatus(void)
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
 * @param[in] status COMPTEL instrument status.
 * @return COMPTEL instrument status.
 ***************************************************************************/
GCOMStatus& GCOMStatus::operator=(const GCOMStatus& status)
{
    // Execute only if object is not identical
    if (this != &status) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(status);

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
 * @brief Clear COMPTEL instrument status
 ***************************************************************************/
void GCOMStatus::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL instrument status
 *
 * @return Pointer to deep copy of COMPTEL instrument status.
 ***************************************************************************/
GCOMStatus* GCOMStatus::clone(void) const
{
    return new GCOMStatus(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL instrument status database
 ***************************************************************************/
void GCOMStatus::load(void) const
{
    // Load database
    load_status();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return D1 module status
 *
 * @param[in] tjd TJD for status inquiry.
 * @param[in] module D1 module number (1-7).
 * @return D1 module status.
 ***************************************************************************/
int GCOMStatus::d1status(const int& tjd, const int& module) const
{
    // Load database if it's not yet available
    if (m_tjds.empty()) {
        load_status();
    }

    // Check that module number is in range
    if (module < 1 || module > 7) {
        std::string msg = "Invalid D1 module number ("+gammalib::str(module)+
                          "). Number needs to be comprised within [1,7].";
        throw GException::invalid_argument(G_D1STATUS, msg);
    }

    // Update cache
    update_cache(tjd);

    // Return module status
    return (m_last_d1status[module-1]);
}


/***********************************************************************//**
 * @brief Return D2 module status
 *
 * @param[in] tjd TJD for status inquiry.
 * @param[in] module D2 module number (1-14).
 * @return D2 module status.
 ***************************************************************************/
int GCOMStatus::d2status(const int& tjd, const int& module) const
{
    // Load database if it's not yet available
    if (m_tjds.empty()) {
        load_status();
    }

    // Check that module number is in range
    if (module < 1 || module > 14) {
        std::string msg = "Invalid D2 module number ("+gammalib::str(module)+
                          "). Number needs to be comprised within [1,14].";
        throw GException::invalid_argument(G_D2STATUS, msg);
    }

    // Update cache
    update_cache(tjd);

    // Return module status
    return (m_last_d2status[module-1]);
}


/***********************************************************************//**
 * @brief Print COMPTEL instrument status
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL instrument status information.
 ***************************************************************************/
std::string GCOMStatus::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMStatus ===");

        // Append information
        result.append("\n"+gammalib::parformat("Database"));
        if (m_tjds.empty()) {
            result.append("not loaded");
        }
        else {
            result.append("loaded");
            result.append("\n"+gammalib::parformat("Number of entries"));
            result.append(gammalib::str(m_tjds.size()));
            result.append("\n"+gammalib::parformat("TJD range"));
            result.append(gammalib::str(m_tjds[0]));
            result.append(" - ");
            result.append(gammalib::str(m_tjds[m_tjds.size()-1]));
            result.append(" days");
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
void GCOMStatus::init_members(void)
{
    // Initialise members
    m_tjds.clear();
    m_d1status.clear();
    m_d2status.clear();

    // Initialise cache
    m_last_tjd = 0;
    m_last_d1status.clear();
    m_last_d2status.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] status COMPTEL instrument status.
 ***************************************************************************/
void GCOMStatus::copy_members(const GCOMStatus& status)
{
    // Copy members
    m_tjds     = status.m_tjds;
    m_d1status = status.m_d1status;
    m_d2status = status.m_d2status;

    // Copy cache
    m_last_tjd      = status.m_last_tjd;
    m_last_d1status = status.m_last_d1status;
    m_last_d2status = status.m_last_d2status;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMStatus::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update module status cache
 *
 * @param[in] tjd TJD for status inquiry.
 ***************************************************************************/
void GCOMStatus::update_cache(const int& tjd) const
{
    // Update cache only if the TJD value has changed
    if (tjd != m_last_tjd) {

        // Update TJD value
        m_last_tjd = tjd;

        // Signal that cache update was successful
        bool success = false;

        // Find index
        for (int i = 0; i < m_tjds.size(); ++i) {
            if (m_tjds[i] == tjd) {
                m_last_d1status = m_d1status[i];
                m_last_d2status = m_d2status[i];
                success         = true;
                break;
            }
        }

        // Throw an exception if update was not successful
        if (!success) {
            std::string msg = "No COMPTEL status information found for TJD "+
                              gammalib::str(tjd)+" in module status database. "
                              "Please exclude that day from your analysis.";
            throw GException::invalid_value(G_UPDATE_CACHE, msg);
        }

    } // endif: TJD value has changed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load status information from database
 *
 * Loads the status information from the database.
 ***************************************************************************/
void GCOMStatus::load_status(void) const
{
    // Initialise status vectors
    m_tjds.clear();
    m_d1status.clear();
    m_d2status.clear();

    // Locate COMPTEL calibration database that hold the module status
    // information
    GCaldb caldb("cgro","comptel");

    // Get database filename
    GFilename filename = caldb.filename("","","MODULE_STATUS","","","DEFAULT");

    // If filename is empty then throw an exception
    if (filename.is_empty()) {
        std::string msg = "COMPTEL status information not found in COMPTEL "
                          "calibration database. Make sure that the COMPTEL "
                          "calibration database is in the path of the $CALDB "
                          "environment variable.";
        throw GException::invalid_value(G_LOAD_STATUS, msg);
    }

    // Open FITS file
    GFits fits(filename);

    // Get status table
    const GFitsTable& table = *fits.table(1);

    // Extract number of entries in table
    int num = table.nrows();

    // If there are entries then read them
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_tjd      = table["TJD"];
        const GFitsTableCol* ptr_d1status = table["D1STATUS"];
        const GFitsTableCol* ptr_d2status = table["D2STATUS"];

        // Reserve space
        m_tjds.reserve(num);
        m_d1status.reserve(num);
        m_d2status.reserve(num);
        
        // Copy data from table into vectors
        for (int i = 0; i < num; ++i) {

            // Copy TJD
            m_tjds.push_back(ptr_tjd->integer(i));

            // Copy D1 status
            std::vector<int> d1status;
            for (int k = 0; k < 7; ++k) {
                d1status.push_back(ptr_d1status->integer(i,k));
            }
            m_d1status.push_back(d1status);

            // Copy D2 status
            std::vector<int> d2status;
            for (int k = 0; k < 14; ++k) {
                d2status.push_back(ptr_d2status->integer(i,k));
            }
            m_d2status.push_back(d2status);

        } // endfor: looped over entries

    } // endif: there were entries

    // Close FITS file
    fits.close();

    // Return
    return;
}
