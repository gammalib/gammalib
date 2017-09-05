/***************************************************************************
 *             GCOMTim.cpp - COMPTEL Good Time Intervals class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMTim.cpp
 * @brief COMPTEL Good Time Intervals class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GCOMTim.hpp"
#include "GCOMSupport.hpp"

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
GCOMTim::GCOMTim(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] tim COMPTEL Good Time Intervals.
 ***************************************************************************/
GCOMTim::GCOMTim(const GCOMTim& tim)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(tim);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename TIM file name.
 * @param[in] usage Usage selection (blank: accept all usage strings).
 * @param[in] mode Mode selection (blank: accept all mode strings).
 *
 * Constructs COMPTEL Good Time Intervals from the information in a TIM
 * FITS file.
 ***************************************************************************/
GCOMTim::GCOMTim(const GFilename&   filename,
                 const std::string& usage,
                 const std::string& mode)
{
    // Initialise class members
    init_members();

    // Load TIM file
    load(filename, usage, mode);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMTim::~GCOMTim(void)
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
 * @param[in] tim COMPTEL Good Time Intervals.
 * @return COMPTEL Good Time Intervals.
 ***************************************************************************/
GCOMTim& GCOMTim::operator=(const GCOMTim& tim)
{
    // Execute only if object is not identical
    if (this != &tim) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(tim);

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
 * @brief Clear COMPTEL good time intervals
 ***************************************************************************/
void GCOMTim::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL good time intervals
 *
 * @return Pointer to deep copy of COMPTEL Good Time Intervals.
 ***************************************************************************/
GCOMTim* GCOMTim::clone(void) const
{
    return new GCOMTim(*this);
}


/***********************************************************************//**
 * @brief Load COMPTEL Good Time Intervals from FITS file
 *
 * @param[in] filename TIM file name.
 * @param[in] usage Usage selection (blank: accept all usage strings).
 * @param[in] mode Mode selection (blank: accept all mode strings).
 *
 * Load COMPTEL Good Time Intervals from the information in a TIM FITS file.
 ***************************************************************************/
void GCOMTim::load(const GFilename&   filename,
                   const std::string& usage,
                   const std::string& mode)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *fits.table(1);

    // Read TIM file
    read(hdu, usage, mode);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL Good Time Intervals from FITS table
 *
 * @param[in] table TIM FITS table.
 * @param[in] usage Usage selection (blank: accept all usage strings).
 * @param[in] mode Mode selection (blank: accept all mode strings).
 *
 * Reads COMPTEL Good Time Intervals from the information in a TIM FITS
 * table.
 ***************************************************************************/
void GCOMTim::read(const GFitsTable&  table,
                   const std::string& usage,
                   const std::string& mode)
{
    // Clear object
    clear();

    // Extract number of TIM rows
    int num = table.nrows();

    // If there is TIM information then load it
    if (num > 0) {

        // Get column pointers
        const GFitsTableCol* ptr_tjd_start = table["START_TJD"];
        const GFitsTableCol* ptr_tic_start = table["START_TIC"];
        const GFitsTableCol* ptr_tjd_end   = table["END_TJD"];
        const GFitsTableCol* ptr_tic_end   = table["END_TIC"];
        const GFitsTableCol* ptr_usage     = table["USAGE"];
        const GFitsTableCol* ptr_mode      = table["MODE"];
 
        // Convert data into GTI
        for (int i = 0; i < num; ++i) {

            // Skip if usage string does not match
            if (!usage.empty() && ptr_usage->string(i) != usage) {
                continue;
            }

            // Skip if mode string does not match
            if (!mode.empty() && ptr_mode->string(i) != mode) {
                continue;
            }

            // Convert times
            GTime tstart = com_time(ptr_tjd_start->integer(i),
                                    ptr_tic_start->integer(i));
            GTime tstop  = com_time(ptr_tjd_end->integer(i),
                                    ptr_tic_end->integer(i));

            // Append GTI
            m_gti.append(tstart, tstop);

        } // endfor: Looped over GTIs

    } // endif: there was TIM information

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Good Time Intervals
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Good Time Intervals information.
 ***************************************************************************/
std::string GCOMTim::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMTim ===");

        // Append GTI information
        result.append("\n"+gammalib::parformat("Number of intervals"));
        result.append(gammalib::str(m_gti.size()));
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(m_gti.ontime())+" sec");
        result.append("\n"+gammalib::parformat("Elapsed time"));
        result.append(gammalib::str(m_gti.telapse())+" sec");
        result.append("\n"+gammalib::parformat("MJD range"));
        result.append(gammalib::str(m_gti.tstart().mjd()));
        result.append(" - ");
        result.append(gammalib::str(m_gti.tstop().mjd()));
        result.append("\n"+gammalib::parformat("UTC range"));
        result.append(m_gti.tstart().utc());
        result.append(" - ");
        result.append(m_gti.tstop().utc());

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
void GCOMTim::init_members(void)
{
    // Initialise members
    m_gti.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] tim COMPTEL Good Time Intervals.
 ***************************************************************************/
void GCOMTim::copy_members(const GCOMTim& tim)
{
    // Copy members
    m_gti = tim.m_gti;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMTim::free_members(void)
{
    // Return
    return;
}
