/***************************************************************************
 *         GCOMOads.cpp - COMPTEL Orbit Aspect Data container class        *
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
 * @file GCOMOads.hpp
 * @brief COMPTEL Orbit Aspect Data container class implementation
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
#include "GCOMOads.hpp"
#include "GCOMSupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMOads::at(int&)"
#define G_INSERT                           "GCOMOads::insert(int&, GCOMOad&)"
#define G_REMOVE                                     "GCOMOads::remove(int&)"
#define G_READ                                  "GCOMOads::read(GFitsTable&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_GEORAD_WARNING

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty Orbit Aspect Data container
 ***************************************************************************/
GCOMOads::GCOMOads(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Filename constructor
 *
 * @param[in] filename Orbit Aspect Data FITS file
 *
 * Constructs an Orbit Aspect Data container from an OAD FITS file.
 ***************************************************************************/
GCOMOads::GCOMOads(const GFilename& filename)
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
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads::GCOMOads(const GCOMOads& oads)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(oads);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMOads::~GCOMOads(void)
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
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 * @return COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads& GCOMOads::operator=(const GCOMOads& oads)
{
    // Execute only if object is not identical
    if (this != &oads) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(oads);

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
 * @brief Clear COMPTEL Orbit Aspect Data container
 ***************************************************************************/
void GCOMOads::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Orbit Aspect Data container
 *
 * @return Pointer to deep copy of COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads* GCOMOads::clone(void) const
{
    return new GCOMOads(*this);
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
GCOMOad& GCOMOads::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Orbit Aspect Data index",
                                       index, size());
    }

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data (const version)
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
const GCOMOad& GCOMOads::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Orbit Aspect Data index",
                                       index, size());
    }

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Append Orbit Aspect Data to container
 *
 * @param[in] oad Orbit Aspect Data.
 * @return Reference to appended Orbit Aspect Data.
 *
 * Appends Orbit Aspect Data to the container by making a deep copy of the
 * Orbit Aspect Data.
 ***************************************************************************/
GCOMOad& GCOMOads::append(const GCOMOad& oad)
{
    // Append oad to list
    m_oads.push_back(oad);

    // Return reference
    return m_oads[size()-1];
}


/***********************************************************************//**
 * @brief Insert Orbit Aspect Data into container
 *
 * @param[in] index Orbit Aspect Data index (0,...,size()-1).
 * @param[in] oad Orbit Aspect Data.
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Inserts an @p Orbit Aspect Data into the container before the Orbit
 * Aspect Data with the specified @p index.
 ***************************************************************************/
GCOMOad& GCOMOads::insert(const int& index, const GCOMOad& oad)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Orbit Aspect Data index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Orbit Aspect Data index",
                                           index, size());
        }
    }
    #endif

    // Inserts Orbit Aspect Data
    m_oads.insert(m_oads.begin()+index, oad);

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Remove Orbit Aspect Data from container
 *
 * @param[in] index Orbit Aspect Data index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Remove Orbit Aspect Data of specified @p index from container.
 ***************************************************************************/
void GCOMOads::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Orbit Aspect Data index",
                                       index, size());
    }
    #endif

    // Erase Orbit Aspect Data from container
    m_oads.erase(m_oads.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Orbit Aspect Data container
 *
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 *
 * Append COMPTEL Orbit Aspect Data container to the container.
 ***************************************************************************/
void GCOMOads::extend(const GCOMOads& oads)
{
    // Do nothing if COMPTEL Orbit Aspect Data container is empty
    if (!oads.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = oads.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_oads.push_back(oads[i]);
        }

    } // endif: COMPTEL Orbit Aspect Data container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load COMPTEL Orbit Aspect Data FITS file
 *
 * @param[in] filename COMPTEL OAD FITS file name.
 *
 * Loads an COMPTEL Orbit Aspect Data FITS file in the container.
 ***************************************************************************/
void GCOMOads::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Get HDU (pointer is always valid)
    const GFitsTable& hdu = *fits.table(1);

    // Read Orbit Aspect Data FITS table
    read(hdu);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL Orbit Aspect Data FITS table
 *
 * @param[in] table COMPTEL OAD FITS table.
 *
 * Reads COMPTEL Orbit Aspect Data FITS table into the container.
 ***************************************************************************/
void GCOMOads::read(const GFitsTable& table)
{
    // Earth radius in km including atmosphere (see EVP routine EVSCDJ Rev.7)
    const double era = 6451.03;

    // Set Earth radius in km excluding atmosphere
    const double erd     = 6378.5;
    const double erd_min = erd + 300.0; // Lower limit on orbit
    const double erd_max = erd + 530.0; // Upper limit on orbit

    // Clear
    clear();

    // Extract number of records in FITS file
    int num = table.nrows();

    // If there are records then load them
    if (num > 0) {

        // Reserve data
        m_oads.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_tjd   = table["TJD"];  // days
        const GFitsTableCol* ptr_tics  = table["TICS"]; // ticks
        const GFitsTableCol* ptr_gcaz  = table["GCAZ"]; // rad
        const GFitsTableCol* ptr_gcel  = table["GCEL"]; // rad
        const GFitsTableCol* ptr_posx  = table["POSX"]; // km
        const GFitsTableCol* ptr_posy  = table["POSY"]; // km
        const GFitsTableCol* ptr_posz  = table["POSZ"]; // km
        const GFitsTableCol* ptr_zrasc = table["ZRASC"]; // rad
        const GFitsTableCol* ptr_zdecl = table["ZDECL"]; // rad
        const GFitsTableCol* ptr_xrasc = table["XRASC"]; // rad
        const GFitsTableCol* ptr_xdecl = table["XDECL"]; // rad

        // Initialise Earth radius angle
        double georad = 73.5;

        // Copy data from columns into records
        for (int i = 0; i < num; ++i) {

            // Allocate OAD record
            GCOMOad oad;

            // Get time of record
            int tjd  = ptr_tjd->integer(i);
            int tics = ptr_tics->integer(i);

            // Store time information. The stop time is defined as the start
            // time plus 131071 tics, since the length of one superpacket is
            // 16.384 secs, i.e. 16.384 * 8000 = 131072 ticks
            oad.tjd(tjd);
            oad.tics(tics);
            oad.tstart(com_time(tjd, tics));
            oad.tstop(com_time(tjd, tics + 131071));

            // Set geocentre azimuth and zenith angle in deg
            oad.gcaz(ptr_gcaz->real(i) * gammalib::rad2deg);
            oad.gcel(ptr_gcel->real(i) * gammalib::rad2deg);

            // Set telescope Z- and X-axes
            GSkyDir zaxis;
            GSkyDir xaxis;
            zaxis.radec(ptr_zrasc->real(i), ptr_zdecl->real(i));
            xaxis.radec(ptr_xrasc->real(i), ptr_xdecl->real(i));
            oad.zaxis(zaxis);
            oad.xaxis(xaxis);

            // Compute apparent radius of Earth
            double radius = std::sqrt(ptr_posx->real(i) * ptr_posx->real(i) +
                                      ptr_posy->real(i) * ptr_posy->real(i) +
                                      ptr_posz->real(i) * ptr_posz->real(i));
            if ((radius > erd_min) && (radius < erd_max)) {
                georad = std::asin(era/radius) * gammalib::rad2deg;
            }
            #if defined(G_GEORAD_WARNING)
            else {
                std::string msg = "Error in CGRO position. Distance from "
                                  "geocentre is "+gammalib::str(radius)+ " km "
                                  "while it should be in the interval ["+
                                  gammalib::str(erd_min)+","+
                                  gammalib::str(erd_max)+"] km. Use previous "
                                  "spacecraft altitude.";
                gammalib::warning(G_READ, msg);
            }
            #endif
            oad.georad(georad);

            // Append record
            m_oads.push_back(oad);

        } // endfor: looped over OAD records

    } // endif: there were records to load

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Orbit Aspect Data container
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Orbit Aspect Data container information.
 ***************************************************************************/
std::string GCOMOads::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMOads ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Number of records"));
        result.append(gammalib::str(size()));
        if (size() > 0) {
            result.append("\n"+gammalib::parformat("MJD range"));
            result.append(gammalib::str(m_oads[0].tstart().mjd()));
            result.append(" - ");
            result.append(gammalib::str(m_oads[size()-1].tstop().mjd()));
            result.append(" days");
            result.append("\n"+gammalib::parformat("UTC range"));
            result.append(m_oads[0].tstart().utc());
            result.append(" - ");
            result.append(m_oads[size()-1].tstop().utc());
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
void GCOMOads::init_members(void)
{
    // Initialise members
    m_oads.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
void GCOMOads::copy_members(const GCOMOads& oads)
{
    // Copy members
    m_oads = oads.m_oads;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMOads::free_members(void)
{
    // Return
    return;
}
