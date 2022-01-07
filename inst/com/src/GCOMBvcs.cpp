/***************************************************************************
 *   GCOMBvcs.cpp - COMPTEL Solar System Barycentre Data container class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvcs.hpp
 * @brief COMPTEL Solar System Barycentre Data container class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMBvcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMBvcs::at(int&)"
#define G_INSERT                           "GCOMBvcs::insert(int&, GCOMBvc&)"
#define G_REMOVE                                     "GCOMBvcs::remove(int&)"
#define G_READ                                  "GCOMBvcs::read(GFitsTable&)"

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
 *
 * Constructs an empty COMPTEL Solar System Barycentre Data container
 ***************************************************************************/
GCOMBvcs::GCOMBvcs(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename COMPTEL Solar System Barycentre Data FITS file
 *
 * Constructs a COMPTEL Solar System Barycentre Data container from a BVC
 * FITS file.
 ***************************************************************************/
GCOMBvcs::GCOMBvcs(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load data
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bvcs COMPTEL Solar System Barycentre Data container.
 ***************************************************************************/
GCOMBvcs::GCOMBvcs(const GCOMBvcs& bvcs)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(bvcs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMBvcs::~GCOMBvcs(void)
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
 * @param[in] bvcs COMPTEL Solar System Barycentre Data container.
 * @return COMPTEL Solar System Barycentre Data container.
 ***************************************************************************/
GCOMBvcs& GCOMBvcs::operator=(const GCOMBvcs& bvcs)
{
    // Execute only if object is not identical
    if (this != &bvcs) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(bvcs);

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
 * @brief Clear COMPTEL Solar System Barycentre Data container
 ***************************************************************************/
void GCOMBvcs::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Solar System Barycentre Data container
 *
 * @return Pointer to deep copy of COMPTEL Solar System Barycentre Data
 *         container.
 ***************************************************************************/
GCOMBvcs* GCOMBvcs::clone(void) const
{
    return new GCOMBvcs(*this);
}


/***********************************************************************//**
 * @brief Return reference to Solar System Barycentre Data
 *
 * @param[in] index Solar System Barycentre Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Solar System Barycentre Data index is out of range.
 *
 * Returns a reference to the Solar System Barycentre Data with the specified
 * @p index.
 ***************************************************************************/
GCOMBvc& GCOMBvcs::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Solar System Barycentre Data index",
                                       index, size());
    }

    // Return reference
    return m_bvcs[index];
}


/***********************************************************************//**
 * @brief Return reference to Solar System Barycentre Data (const version)
 *
 * @param[in] index Solar System Barycentre Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Solar System Barycentre Data index is out of range.
 *
 * Returns a reference to the Solar System Barycentre Data with the specified
 * @p index.
 ***************************************************************************/
const GCOMBvc& GCOMBvcs::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Solar System Barycentre Data index",
                                       index, size());
    }

    // Return reference
    return m_bvcs[index];
}


/***********************************************************************//**
 * @brief Append Solar System Barycentre Data to container
 *
 * @param[in] bvc Solar System Barycentre Data.
 * @return Reference to appended Solar System Barycentre Data.
 *
 * Appends Solar System Barycentre Data to the container by making a deep
 * copy of the Solar System Barycentre Data.
 ***************************************************************************/
GCOMBvc& GCOMBvcs::append(const GCOMBvc& bvc)
{
    // Append bvc to list
    m_bvcs.push_back(bvc);

    // Return reference
    return m_bvcs[size()-1];
}


/***********************************************************************//**
 * @brief Insert Solar System Barycentre Data into container
 *
 * @param[in] index Solar System Barycentre Data index (0,...,size()-1).
 * @param[in] bvc Solar System Barycentre Data.
 *
 * @exception GException::out_of_range
 *            Solar System Barycentre Data index is out of range.
 *
 * Inserts Solar System Barycentre Data into the container before the Solar
 * System Barycentre Data with the specified @p index.
 ***************************************************************************/
GCOMBvc& GCOMBvcs::insert(const int& index, const GCOMBvc& bvc)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT,
                                           "Solar System Barycentre Data index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT,
                                           "Solar System Barycentre Data index",
                                           index, size());
        }
    }
    #endif

    // Inserts Solar System Barycentre Data
    m_bvcs.insert(m_bvcs.begin()+index, bvc);

    // Return reference
    return m_bvcs[index];
}


/***********************************************************************//**
 * @brief Remove Solar System Barycentre Data from container
 *
 * @param[in] index Solar System Barycentre Data index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Solar System Barycentre Data index is out of range.
 *
 * Remove Solar System Barycentre Data of specified @p index from container.
 ***************************************************************************/
void GCOMBvcs::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE,
                                       "Solar System Barycentre Data index",
                                       index, size());
    }
    #endif

    // Erase Solar System Barycentre Data from container
    m_bvcs.erase(m_bvcs.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Solar System Barycentre Data container
 *
 * @param[in] bvcs COMPTEL Solar System Barycentre Data container.
 *
 * Append COMPTEL Solar System Barycentre Data container to the container.
 ***************************************************************************/
void GCOMBvcs::extend(const GCOMBvcs& bvcs)
{
    // Do nothing if COMPTEL Solar System Barycentre Data container is empty
    if (!bvcs.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = bvcs.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_bvcs.push_back(bvcs[i]);
        }

    } // endif: COMPTEL Solar System Barycentre Data container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load COMPTEL Solar System Barycentre Data FITS file
 *
 * @param[in] filename COMPTEL BVC FITS file name.
 *
 * Loads an COMPTEL Solar System Barycentre FITS file in the container.
 ***************************************************************************/
void GCOMBvcs::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *fits.table(1);

    // Read Solar System Barycentre Data FITS table
    read(hdu);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL Solar System Barycentre Data FITS table
 *
 * @param[in] table COMPTEL BVC FITS table.
 *
 * Reads COMPTEL Solar System Barycentre Data FITS table into the container.
 ***************************************************************************/
void GCOMBvcs::read(const GFitsTable& table)
{
    // Clear
    clear();

    // Extract number of records in FITS file
    int num = table.nrows();

    // If there are records then load them
    if (num > 0) {

        // Reserve data
        m_bvcs.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_tjd    = table["TJD"];    // days
        const GFitsTableCol* ptr_tics   = table["TICS"];   // ticks
        const GFitsTableCol* ptr_ssbx   = table["SSB_X"];  // km
        const GFitsTableCol* ptr_ssby   = table["SSB_Y"];  // km
        const GFitsTableCol* ptr_ssbz   = table["SSB_Z"];  // km
        const GFitsTableCol* ptr_tdelta = table["TDELTA"]; // s

        // Copy data from columns into records
        for (int i = 0; i < num; ++i) {

            // Allocate BVC record
            GCOMBvc bvc;

            // Get time of record
            int tjd  = ptr_tjd->integer(i);
            int tics = ptr_tics->integer(i);

            // Store time information. The stop time is defined as the start
            // time plus 131071 tics, since the length of one superpacket is
            // 16.384 secs, i.e. 16.384 * 8000 = 131072 ticks
            bvc.tjd(tjd);
            bvc.tics(tics);
            bvc.tstart(gammalib::com_time(tjd, tics));
            bvc.tstop(gammalib::com_time(tjd, tics + 131071));

            // Set Solar System Barycentre vector in celestial system
            GVector ssb(ptr_ssbx->real(i), ptr_ssby->real(i), ptr_ssbz->real(i));
            bvc.ssb(ssb);

            // Set TDB-UTC time difference
            bvc.tdelta(ptr_tdelta->real(i));

            // Append record
            m_bvcs.push_back(bvc);

        } // endfor: looped over BVC records

    } // endif: there were records to load

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Solar System Barycentre Data container
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Solar System Barycentre Data container
 *         information.
 ***************************************************************************/
std::string GCOMBvcs::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMBvcs ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Superpackets"));
        result.append(gammalib::str(size()));
        if (size() > 0) {

            // Append time range
            result.append("\n"+gammalib::parformat("TJD range"));
            result.append(gammalib::str(m_bvcs[0].tjd()));
            result.append(":");
            result.append(gammalib::str(m_bvcs[0].tics()));
            result.append(" - ");
            result.append(gammalib::str(m_bvcs[size()-1].tjd()));
            result.append(":");
            result.append(gammalib::str(m_bvcs[size()-1].tics()));
            result.append("\n"+gammalib::parformat("MJD range"));
            result.append(gammalib::str(m_bvcs[0].tstart().mjd()));
            result.append(" - ");
            result.append(gammalib::str(m_bvcs[size()-1].tstop().mjd()));
            result.append(" days");
            result.append("\n"+gammalib::parformat("UTC range"));
            result.append(m_bvcs[0].tstart().utc());
            result.append(" - ");
            result.append(m_bvcs[size()-1].tstop().utc());

            // Append detailed information
            GChatter reduced_chatter = gammalib::reduce(chatter);
            if (reduced_chatter > SILENT) {

                // Append TJDs
                int tjd = 0;
                int num = 0;
                for (int i = 0; i < size(); ++i) {
                    if (m_bvcs[i].tjd() != tjd) {
                        if (num > 0) {
                            std::string key = "TJD "+gammalib::str(tjd);
                            result.append("\n"+gammalib::parformat(key));
                            result.append(gammalib::str(num)+" superpackets");
                        }
                        tjd = m_bvcs[i].tjd();
                        num = 1;
                    }
                    else {
                    num++;
                    }
                }
                std::string key = "TJD "+gammalib::str(tjd);
                result.append("\n"+gammalib::parformat(key));
                result.append(gammalib::str(num)+" superpackets");

            } // endif: detailed information requested

        } // endif: there were records

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
void GCOMBvcs::init_members(void)
{
    // Initialise members
    m_bvcs.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bvcs COMPTEL Solar System Barycentre Data container.
 ***************************************************************************/
void GCOMBvcs::copy_members(const GCOMBvcs& bvcs)
{
    // Copy members
    m_bvcs = bvcs.m_bvcs;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMBvcs::free_members(void)
{
    // Return
    return;
}
