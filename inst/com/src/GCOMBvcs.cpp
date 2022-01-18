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
#include "GCOMOad.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMBvcs::at(int&)"
#define G_INSERT                           "GCOMBvcs::insert(int&, GCOMBvc&)"
#define G_REMOVE                                     "GCOMBvcs::remove(int&)"
#define G_READ                                  "GCOMBvcs::read(GFitsTable&)"
#define G_TDELTA                         "GCOMBvcs::tdelta(GSkyDir&, GTime&)"

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

            // Store time information, both as native COMPTEL time and as
            // GTime object
            bvc.tjd(tjd);
            bvc.tics(tics);
            bvc.time(gammalib::com_time(tjd, tics));

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
 * @brief Find Solar System Barycentre Data for Orbit Aspect Data
 *
 * @param[in] oad Orbit Aspect Data.
 * @return Pointer to Solar System Barycentre Data (NULL if no data were found)
 *
 * Finds Solar System Barycentre Data that correspond to the specified Orbit
 * Aspect Data. The method returns the Solar System Barycentre Data with the
 * same TJD as the Orbit Aspect Data and the smallest difference in the
 * number of tics.
 *
 * If this smallest difference is larger than 131072, which is the length
 * of one superpacket, the method returns a NULL pointer.
 *
 * The method also returns a NULL pointer in case that no matching Solar
 * System Barycentre Data was found.
 *
 * The client needs to verify the validity of the Solar System Barycentre
 * Data pointer. The client must not deallocate the associated memory.
 ***************************************************************************/
const GCOMBvc* GCOMBvcs::find(const GCOMOad& oad) const
{
    // Initialise Solar System Barycentre Data
    const GCOMBvc* bvc = NULL;

    // Compute TJD and tics at the centre of the OAD superpacket
    int oad_tjd  = oad.tjd();
    int oad_tics = oad.tics() + 65536;

    // Loop over all records and among those that have the same TJD find
    // the one with the closest number of tics
    int size                = this->size();
    int max_tics_difference = 700000000;
    int ibest               = -1;
    for (int i = 0; i < size; ++i) {
        if (m_bvcs[i].tjd() == oad_tjd) {
            int tics_difference = std::abs(m_bvcs[i].tics() - oad_tics);
            if (tics_difference < max_tics_difference) {
                ibest               = i;
                max_tics_difference = tics_difference;
            }
        }
    }

    // If matching Solar System Barycentre Data were found then check
    // whether the tics difference is acceptable. A maximum difference
    // of 131072 tics, corresponding to the duration of one superpacket
    // of 16.384 sec, is considered as acceptable.
    if (ibest >= 0) {
        if (max_tics_difference < 131072) {
            bvc = &(m_bvcs[ibest]);
        }
    }

    // Return Solar System Barycentre Data
    return bvc;
}


/***********************************************************************//**
 * @brief Return time difference between photon arrival time at CGRO and
 *        the Solar System Barycentre (SSB)
 *
 * @param[in] dir Source position.
 * @param[in] time CGRO photon arrival time.
 * @return Time difference between photon arrival times (s)
 *
 * @exception GException::invalid_value
 *            Not enough SSB vectors loaded for computation
 *
 * Returns the time difference between photon arrival time at CGRO and the
 * Solar System Barycentre (SSB). The arrival time at the SSB is computed
 * by adding the time difference to the photon arrival time as measured by
 * COMPTEL
 *
 * \f[
 *    T_{\rm SSB} = T_{\rm CGRO} + \Delta T
 * \f]
 *
 * The routine implements the algorithm PUL-AL-004 and is inspried from the
 * COMPASS code evpbin02.pulssb.f.
 *
 * It computes
 *
 * \f[
 *    \Delta T = \Delta T_{\rm travel} - \Delta T_{\rm rel} + \Delta T_{\rm BVC}
 * \f]
 *
 * where
 *
 * \f[
 *    \Delta T_{\rm travel} = \left( \vec{SSB} \cdot \vec{n} \right) \times 10^{-6}
 * \f]
 *
 * is the light travel time in seconds between CGRO and the SSB, with
 * \f$\vec{SSB}\f$ being the vector going from the SSB to CGRO, and
 * \f$\vec{n}\f$ is the normalised vector of the source position, provided
 * by the GSkyDir::celvector() method,
 *
 * \f[
 *    \Delta T_{\rm rel} =
 *    -2 R \log \left( 1 + \frac{\Delta T_{\rm travel}}{|\vec{SSB}| * 10^{-6}} \right)
 * \f]
 *
 * is the relativistic delay due to the Sun in seconds, with
 * \f$R=0.49254909 \times 10^{-5}\f$ s, and \f$\Delta T_{\rm BVC}\f$ is the
 * difference in seconds due to the time unit conversion.
 *
 * The values of \f$\vec{SSB}\f$ and \f$\Delta T_{\rm BVC}\f$ are linearly
 * interpolated from the tabulated values based on the specified @p time.
 ***************************************************************************/
double GCOMBvcs::tdelta(const GSkyDir& dir, const GTime& time) const
{
    // Set Schwartzschild radius of sun in lightseconds
    const double t_sun = 4.92549089483e-6;

    // Initialise SSB vector and time difference
    GVector ssb;
    double  tdelta = 0.0;

    // Get number of elements in container
    int size = this->size();

    // Throw an exception if no vectors were loaded
    if (size < 2) {
        std::string msg = gammalib::str(size)+" SSB vectors loaded but at "
                          "least two vectors are needed to compute the time "
                          "delay. Please load SSB vectors from a BVC file "
                          "before calling this method.";
        throw GException::invalid_value(G_TDELTA, msg);
    }

    // If time is before first BVC element then extrapolate using the first
    // two BVC elements
    if (time < m_bvcs[0].time()) {
        int    ilow = 0;
        int    iup  = 1;
        double wlow = (m_bvcs[iup].time() - time) / (m_bvcs[iup].time() - m_bvcs[ilow].time());
        double wup  = 1.0 - wlow;
        ssb    = wlow * m_bvcs[ilow].ssb()    + wup * m_bvcs[iup].ssb();
        tdelta = wlow * m_bvcs[ilow].tdelta() + wup * m_bvcs[iup].tdelta();
    }

    // ... else if time is after last BVC element then extrapolate using the last
    // two BVC elements
    else if (time > m_bvcs[size-1].time()) {
        int    ilow = size-2;
        int    iup  = size-1;
        double wlow = (m_bvcs[iup].time() - time) / (m_bvcs[iup].time() - m_bvcs[ilow].time());
        double wup  = 1.0 - wlow;
        ssb    = wlow * m_bvcs[ilow].ssb()    + wup * m_bvcs[iup].ssb();
        tdelta = wlow * m_bvcs[ilow].tdelta() + wup * m_bvcs[iup].tdelta();
    }

    // ... otherwise time is comprised in the elements hence we search for
    // the bracketing elements and perform a linear interpolation of the
    // vector and time delay
    else {

        // Compute TJD and tics from CGRO photon arrival time
        int tjd  = gammalib::com_tjd(time);
        int tics = gammalib::com_tics(time);

        // Loop over all records and among those that have the same TJD find
        // the one with the closest number of tics
        int max_tics_difference = 700000000;
        int ibest               = -1;
        for (int i = 0; i < size; ++i) {
            if (m_bvcs[i].tjd() == tjd) {
                int tics_difference = std::abs(m_bvcs[i].tics() - tics);
                if (tics_difference < max_tics_difference) {
                    ibest               = i;
                    max_tics_difference = tics_difference;
                }
            }
        }

        // Find bracketing records and interpolating weighting factors
        int    ilow =  -1;
        int    iup  =  -1;
        double wlow = 0.0;
        double wup  = 0.0;
        if (m_bvcs[ibest].tics() <= tics) {
            ilow = ibest;
            if (ibest < size-1) {
                iup  = ibest + 1;
                wlow = (m_bvcs[iup].time() - time) / (m_bvcs[iup].time() - m_bvcs[ilow].time());
                wup  = 1.0 - wlow;
            }
            else {
                iup  = ibest;
                wlow = 1.0;
            }
        }
        else {
            iup = ibest;
            if (ibest > 0) {
                ilow = ibest - 1;
                wlow = (m_bvcs[iup].time() - time) / (m_bvcs[iup].time() - m_bvcs[ilow].time());
                wup  = 1.0 - wlow;
            }
            else {
                ilow = ibest;
                wup  = 1.0;
            }
        }

        // Interpolate vector and time difference
        ssb    = wlow * m_bvcs[ilow].ssb()    + wup * m_bvcs[iup].ssb();
        tdelta = wlow * m_bvcs[ilow].tdelta() + wup * m_bvcs[iup].tdelta();

    } // endelse: performed linear interpolation

    // Get celestial vector
    GVector n = dir.celvector();

    // Compute the light travel time from the satellite to SSB along the
    // pulsar direction
    double travt = ssb * n * 1.0e-6;

    // Compute the relativistic delay due to the Sun
    double r     = norm(ssb) * 1.0e-6;
    double relat = -2.0 * t_sun * std::log(1.0 + (travt/r));

    // Compute the time difference at SSB
    tdelta = travt - relat + tdelta;

    // Return
    return tdelta;
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
            result.append(gammalib::str(m_bvcs[0].time().mjd()));
            result.append(" - ");
            result.append(gammalib::str(m_bvcs[size()-1].time().mjd()));
            result.append(" days");
            result.append("\n"+gammalib::parformat("UTC range"));
            result.append(m_bvcs[0].time().utc());
            result.append(" - ");
            result.append(m_bvcs[size()-1].time().utc());

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
