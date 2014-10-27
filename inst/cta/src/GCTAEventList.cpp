/***************************************************************************
 *            GCTAEventList.cpp - CTA event atom container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAEventList.cpp
 * @brief CTA event atom container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAEventList.hpp"
#include "GCTAException.hpp"
#include "GCTASupport.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableStringCol.hpp"
#include "GTime.hpp"
#include "GTimeReference.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                          "GCTAEventList::operator[](int&)"
#define G_ROI                                     "GCTAEventList::roi(GRoi&)"
#define G_READ_DS_EBOUNDS         "GCTAEventList::read_ds_ebounds(GFitsHDU*)"

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
GCTAEventList::GCTAEventList(void) : GEventList()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] list Event list.
 ***************************************************************************/
GCTAEventList::GCTAEventList(const GCTAEventList& list) : GEventList(list)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(list);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAEventList::~GCTAEventList(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] list Event list.
 * @return Event list.
 ***************************************************************************/
GCTAEventList& GCTAEventList::operator=(const GCTAEventList& list)
{
    // Execute only if object is not identical
    if (this != &list) {

        // Copy base class members
        this->GEventList::operator=(list);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(list);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Event atom access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to an event atom.
 ***************************************************************************/
GCTAEventAtom* GCTAEventList::operator[](const int& index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
    }
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/***********************************************************************//**
 * @brief Event atom access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to an event atom.
 ***************************************************************************/
const GCTAEventAtom* GCTAEventList::operator[](const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
    }
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear event list
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTAEventList::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventList::free_members();
    this->GEvents::free_members();

    // Initialise members
    this->GEvents::init_members();
    this->GEventList::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone event list
 *
 * @return Pointer to deep copy of event list.
 ***************************************************************************/
GCTAEventList* GCTAEventList::clone(void) const
{
    // Clone event list
    return new GCTAEventList(*this);
}


/***********************************************************************//**
 * @brief Load events from event FITS file.
 *
 * @param[in] filename Name of FITS file from which events are loaded.
 *
 * Load CTA events from the EVENTS extension EVENTS, read the region of
 * interest and the energy boundaries from the data selection keywords,
 * and read the Good Time Intervals from the GTI extension.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GCTAEventList::load(const std::string& filename)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Read event list
    read(file);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save CTA events into FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file (default=false).
 *
 * Write the CTA event list into FITS file.
 ***************************************************************************/
void GCTAEventList::save(const std::string& filename,
                         const bool& clobber) const
{
    // Create empty FITS file
    GFits fits;

    // Write event list
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS file.
 *
 * @param[in] fits FITS file.
 *
 * This method reads the CTA event list from a FITS file. The extension name
 * for the events is expected to be "EVENTS". The header of the HDU is
 * scanned for data selection keywords that define the energy boundaries and
 * the region of interest. So far, no scan is performed for time intervals.
 *
 * If present, Good Time Intervals will be read from an extension names 
 * "GTI". If no "GTI" extension is present, a single Good Time Interval will
 * be assumed based on the TSTART and TSTOP keywords.
 *
 * The method clears the object before reading, thus any information residing
 * in the event list prior to reading will be lost.
 *
 * @todo Ultimately, any events file should have a GTI extension, hence the
 *       extraction of GTIs from TSTART and TSTOP should not be necessary.
 ***************************************************************************/
void GCTAEventList::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get event list HDU
    const GFitsTable& events = *fits.table("EVENTS");

    // If we have a GTI extension, then read Good Time Intervals from that
    // extension
    if (fits.contains("GTI")) {
        const GFitsTable& gti = *fits.table("GTI");
        m_gti.read(gti);
    }

    // ... otherwise build GTI from TSTART and TSTOP
    else {

        // Read start and stop time
        double tstart = events.real("TSTART");
        double tstop  = events.real("TSTOP");

        // Create time reference from header information
        GTimeReference timeref(events);

        // Set start and stop time
        GTime start(tstart);
        GTime stop(tstop);

        // Append start and stop time as single time interval to GTI
        m_gti.append(start, stop);

        // Set GTI time reference
        m_gti.reference(timeref);

    } // endelse: GTI built from TSTART and TSTOP

    // Load event data
    read_events(events);

    // Read region of interest from data selection keyword
    m_roi = gammalib::read_ds_roi(events);

    // Read energy boundaries from data selection keyword
    m_ebounds = gammalib::read_ds_ebounds(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA events into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Write the CTA event list into FITS file.
 *
 * @todo The TELMASK column is allocated with a dummy length of 100.
 * @todo Implement agreed column format
 ***************************************************************************/
void GCTAEventList::write(GFits& file) const
{
    // Allocate FITS binary table HDU
    GFitsBinTable* events = new GFitsBinTable;

    // Write events
    write_events(*events);

    // Write data selection keywords
    write_ds_keys(*events);

    // Append event table to FITS file
    file.append(*events);

    // Free binary table
    delete events;

    // Append GTI to FITS file
    gti().write(file);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set ROI
 *
 * @param[in] roi ROI.
 *
 * @exception GCTAException::bad_roi_type
 *            ROI is not of type GCTARoi.
 ***************************************************************************/
void GCTAEventList::roi(const GRoi& roi)
{
    // Cast ROI dynamically
    const GCTARoi* ptr = dynamic_cast<const GCTARoi*>(&roi);

    // Throw exception if ROI is not of correct type
    if (ptr == NULL) {
        throw GCTAException::bad_roi_type(G_ROI);
    }

    // Set ROI
    m_roi = *ptr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event list information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing event list information.
 ***************************************************************************/
std::string GCTAEventList::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAEventList ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(size()));

        // Append GTI interval
        result.append("\n"+gammalib::parformat("Time interval"));
        if (gti().size() > 0) {
            result.append(gammalib::str(tstart().mjd()));
            result.append(" - ");
            result.append(gammalib::str(tstop().mjd())+" days");
        }
        else {
            result.append("not defined");
        }

        // Append energy intervals
        if (gammalib::reduce(chatter) > SILENT) {
            if (ebounds().size() > 0) {
                result.append("\n"+ebounds().print(gammalib::reduce(chatter)));
            }
            else {
                result.append("\n"+gammalib::parformat("Energy intervals") +
                              "not defined");
            }
        }

        // Append ROI
        if (gammalib::reduce(chatter) > SILENT) {
            if (roi().radius() > 0) {
                result.append("\n"+roi().print(gammalib::reduce(chatter)));
            }
            else {
                result.append("\n"+gammalib::parformat("Region of interest") +
                              "not defined");
            }
        }

        // EXPLICIT: Append IRF cache
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < m_irf_names.size(); ++i) {
                result.append("\n"+gammalib::parformat("IRF cache " +
                              gammalib::str(i)));
                result.append(m_irf_names[i]+" = ");
                int num   = 0;
                for (int k = 0; k < size(); ++k) {
                    if ((m_irf_values[i])[k] != -1.0) {
                        num++;
                    }
                }
                result.append(gammalib::str(num)+" values");
            }
        } // endif: chatter was explicit

    } // endif: chatter was not silent

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Append event atom to event list
 *
 * @param[in] event Event.
 *
 * Appends an event atom to the event list.
 ***************************************************************************/
void GCTAEventList::append(const GCTAEventAtom& event)
{
    // Append event
    m_events.push_back(event);

    // Set event index
    int index = m_events.size()-1;
    m_events[index].m_index = index;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAEventList::init_members(void)
{
    // Initialise members
    m_roi.clear();
    m_events.clear();

    // Initialise cache
    m_irf_names.clear();
    m_irf_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] list Event list.
 ***************************************************************************/
void GCTAEventList::copy_members(const GCTAEventList& list)
{
    // Copy members
    m_roi    = list.m_roi;
    m_events = list.m_events;

    // Copy cache
    m_irf_names  = list.m_irf_names;
    m_irf_values = list.m_irf_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventList::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS table
 *
 * @param[in] table FITS table.
 *
 * This method reads the CTA event list from a FITS table HDU into memory.
 * Depending on the columns existing in the file, it either selects v0 or
 * v1 of the event list reader.
 ***************************************************************************/
void GCTAEventList::read_events(const GFitsTable& table)
{
    // Clear existing events
    m_events.clear();

    // Extract number of events in FITS file
    int num = table.integer("NAXIS2");

    // Continue only if there are events
    if (num > 0) {

        // Read events for v1
        if (table.contains("SHWIDTH") && table.contains("SHLENGTH")) {
            read_events_v1(table);
        }

        // ... otherwise read events for v0
        else {
            read_events_v0(table);
        }

        // Read (optional) Hillas parameters
        read_events_hillas(table);

    } // endif: there were events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS table (version 0)
 *
 * @param[in] table FITS table.
 *
 * This method reads the CTA event list from a FITS table HDU into memory.
 * It is a minimal event reader that is compliant with the initial data
 * format distributed by Karl Kosack. Information that is not present in
 * that format is set to 0.
 ***************************************************************************/
void GCTAEventList::read_events_v0(const GFitsTable& table)
{
    // Clear existing events
    m_events.clear();

    // Extract number of events in FITS file
    int num = table.nrows();

    // If there are events then load them
    if (num > 0) {

        // Reserve data
        m_events.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_eid         = table["EVENT_ID"];
        const GFitsTableCol* ptr_time        = table["TIME"];
        const GFitsTableCol* ptr_multip      = table["MULTIP"];
        const GFitsTableCol* ptr_ra          = table["RA"];
        const GFitsTableCol* ptr_dec         = table["DEC"];
        const GFitsTableCol* ptr_dir_err     = table["DIR_ERR"];
        const GFitsTableCol* ptr_detx        = table["DETX"];
        const GFitsTableCol* ptr_dety        = table["DETY"];
        const GFitsTableCol* ptr_alt         = table["ALT"];
        const GFitsTableCol* ptr_az          = table["AZ"];
        const GFitsTableCol* ptr_corex       = table["COREX"];
        const GFitsTableCol* ptr_corey       = table["COREY"];
        const GFitsTableCol* ptr_core_err    = table["CORE_ERR"];
        const GFitsTableCol* ptr_xmax        = table["XMAX"];
        const GFitsTableCol* ptr_xmax_err    = table["XMAX_ERR"];
        const GFitsTableCol* ptr_energy      = table["ENERGY"];
        const GFitsTableCol* ptr_energy_err  = table["ENERGY_ERR"];

        // Copy data from columns into GCTAEventAtom objects
        GCTAEventAtom event;
        for (int i = 0; i < num; ++i) {
            event.m_index     = i;
            event.m_time.set(ptr_time->real(i), m_gti.reference());
            event.m_dir.dir().radec_deg(ptr_ra->real(i), ptr_dec->real(i));
            event.m_dir.detx(ptr_detx->real(i)*gammalib::deg2rad);
            event.m_dir.dety(ptr_dety->real(i)*gammalib::deg2rad);
            event.m_energy.TeV(ptr_energy->real(i));
            event.m_event_id    = ptr_eid->integer(i);
            event.m_obs_id      = 0;
            event.m_multip      = ptr_multip->integer(i);
            event.m_telmask     = 0;
            event.m_dir_err     = ptr_dir_err->real(i);
            event.m_alt         = ptr_alt->real(i);
            event.m_az          = ptr_az->real(i);
            event.m_corex       = ptr_corex->real(i);
            event.m_corey       = ptr_corey->real(i);
            event.m_core_err    = ptr_core_err->real(i);
            event.m_xmax        = ptr_xmax->real(i);
            event.m_xmax_err    = ptr_xmax_err->real(i);
            event.m_shwidth     = 0.0;
            event.m_shlength    = 0.0;
            event.m_energy_err  = ptr_energy_err->real(i);
            m_events.push_back(event);
        }

    } // endif: there were events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS table (version 1)
 *
 * @param[in] table FITS table.
 *
 * This method reads the CTA event list from a FITS table HDU into memory.
 *
 * @todo Implement agreed column format
 ***************************************************************************/
void GCTAEventList::read_events_v1(const GFitsTable& table)
{
    // Clear existing events
    m_events.clear();

    // Extract number of events in FITS file
    int num = table.nrows();

    // If there are events then load them
    if (num > 0) {

        // Reserve data
        m_events.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_eid         = table["EVENT_ID"];
        const GFitsTableCol* ptr_oid         = table["OBS_ID"];
        const GFitsTableCol* ptr_time        = table["TIME"];
        const GFitsTableCol* ptr_multip      = table["MULTIP"];
        const GFitsTableCol* ptr_ra          = table["RA"];
        const GFitsTableCol* ptr_dec         = table["DEC"];
        const GFitsTableCol* ptr_dir_err     = table["DIR_ERR"];
        const GFitsTableCol* ptr_detx        = table["DETX"];
        const GFitsTableCol* ptr_dety        = table["DETY"];
        const GFitsTableCol* ptr_alt         = table["ALT"];
        const GFitsTableCol* ptr_az          = table["AZ"];
        const GFitsTableCol* ptr_corex       = table["COREX"];
        const GFitsTableCol* ptr_corey       = table["COREY"];
        const GFitsTableCol* ptr_core_err    = table["CORE_ERR"];
        const GFitsTableCol* ptr_xmax        = table["XMAX"];
        const GFitsTableCol* ptr_xmax_err    = table["XMAX_ERR"];
        const GFitsTableCol* ptr_shw         = table["SHWIDTH"];
        const GFitsTableCol* ptr_shl         = table["SHLENGTH"];
        const GFitsTableCol* ptr_energy      = table["ENERGY"];
        const GFitsTableCol* ptr_energy_err  = table["ENERGY_ERR"];

        // Copy data from columns into GCTAEventAtom objects
        GCTAEventAtom event;
        for (int i = 0; i < num; ++i) {
            event.m_index     = i;
            event.m_time.set(ptr_time->real(i), m_gti.reference());
            event.m_dir.dir().radec_deg(ptr_ra->real(i), ptr_dec->real(i));
            event.m_dir.detx(ptr_detx->real(i)*gammalib::deg2rad);
            event.m_dir.dety(ptr_dety->real(i)*gammalib::deg2rad);
            event.m_energy.TeV(ptr_energy->real(i));
            event.m_event_id    = ptr_eid->integer(i);
            event.m_obs_id      = ptr_oid->integer(i);
            event.m_multip      = ptr_multip->integer(i);
            event.m_telmask     = 0;
            event.m_dir_err     = ptr_dir_err->real(i);
            event.m_alt         = ptr_alt->real(i);
            event.m_az          = ptr_az->real(i);
            event.m_corex       = ptr_corex->real(i);
            event.m_corey       = ptr_corey->real(i);
            event.m_core_err    = ptr_core_err->real(i);
            event.m_xmax        = ptr_xmax->real(i);
            event.m_xmax_err    = ptr_xmax_err->real(i);
            event.m_shwidth     = ptr_shw->real(i);
            event.m_shlength    = ptr_shl->real(i);
            event.m_energy_err  = ptr_energy_err->real(i);
            m_events.push_back(event);
        }

    } // endif: there were events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Hillas information for CTA events from FITS table
 *
 * @param[in] table FITS table.
 *
 * This method reads the Hillas reconstruction information for CTA events
 * from an EVENTS file. It searches for the columns HIL_MSW, HIL_MSW_ERR,
 * HIL_MSL, and HIL_MSL_ERR in the FITS table and extracts the relevant
 * columns from the FITS file. If a column is not found, no action is
 * performed.
 *
 * @todo Verify consistency of event list size
 ***************************************************************************/
void GCTAEventList::read_events_hillas(const GFitsTable& table)
{
    // Extract number of events in FITS file
    int num = table.nrows();

    //TODO: Make sure that dimension is consistent with existing event
    //      list

    // Continue only if there are events
    if (num > 0) {

        // HIL_MSW
        if (table.contains("HIL_MSW")) {
            const GFitsTableCol* ptr = table["HIL_MSW"];
            for (int i = 0; i < num; ++i) {
                m_events[i].m_hil_msw = ptr->real(i);
            }
        }

        // HIL_MSW_ERR
        if (table.contains("HIL_MSW_ERR")) {
            const GFitsTableCol* ptr = table["HIL_MSW_ERR"];
            for (int i = 0; i < num; ++i) {
                m_events[i].m_hil_msw_err = ptr->real(i);
            }
        }

        // HIL_MSL
        if (table.contains("HIL_MSL")) {
            const GFitsTableCol* ptr = table["HIL_MSL"];
            for (int i = 0; i < num; ++i) {
                m_events[i].m_hil_msl = ptr->real(i);
            }
        }

        // HIL_MSL_ERR
        if (table.contains("HIL_MSL_ERR")) {
            const GFitsTableCol* ptr = table["HIL_MSL_ERR"];
            for (int i = 0; i < num; ++i) {
                m_events[i].m_hil_msl_err = ptr->real(i);
            }
        }

    } // endif: there were events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA events into FITS table
 *
 * @param[in] hdu FITS table HDU.
 *
 * Write the CTA event list into FITS table.
 *
 * @todo The TELMASK column is allocated with a dummy length of 100.
 * @todo Implement agreed column format
 ***************************************************************************/
void GCTAEventList::write_events(GFitsBinTable& hdu) const
{
    // Set extension name
    hdu.extname("EVENTS");

    // If there are events then write them now
    if (size() > 0) {

        // Allocate columns
        GFitsTableULongCol  col_eid         = GFitsTableULongCol("EVENT_ID", size());
        GFitsTableULongCol  col_oid         = GFitsTableULongCol("OBS_ID", size());
        GFitsTableDoubleCol col_time        = GFitsTableDoubleCol("TIME", size());
        GFitsTableDoubleCol col_live        = GFitsTableDoubleCol("TLIVE", size());
        GFitsTableShortCol  col_multip      = GFitsTableShortCol("MULTIP", size());
        GFitsTableBitCol    col_telmask     = GFitsTableBitCol("TELMASK", size(), 100);
        GFitsTableFloatCol  col_ra          = GFitsTableFloatCol("RA", size());
        GFitsTableFloatCol  col_dec         = GFitsTableFloatCol("DEC", size());
        GFitsTableFloatCol  col_direrr      = GFitsTableFloatCol("DIR_ERR", size());
        GFitsTableFloatCol  col_detx        = GFitsTableFloatCol("DETX", size());
        GFitsTableFloatCol  col_dety        = GFitsTableFloatCol("DETY", size());
        GFitsTableFloatCol  col_alt         = GFitsTableFloatCol("ALT", size());
        GFitsTableFloatCol  col_az          = GFitsTableFloatCol("AZ", size());
        GFitsTableFloatCol  col_corex       = GFitsTableFloatCol("COREX", size());
        GFitsTableFloatCol  col_corey       = GFitsTableFloatCol("COREY", size());
        GFitsTableFloatCol  col_core_err    = GFitsTableFloatCol("CORE_ERR", size());
        GFitsTableFloatCol  col_xmax        = GFitsTableFloatCol("XMAX", size());
        GFitsTableFloatCol  col_xmax_err    = GFitsTableFloatCol("XMAX_ERR", size());
        GFitsTableFloatCol  col_shw         = GFitsTableFloatCol("SHWIDTH", size());
        GFitsTableFloatCol  col_shl         = GFitsTableFloatCol("SHLENGTH", size());
        GFitsTableFloatCol  col_energy      = GFitsTableFloatCol("ENERGY", size());
        GFitsTableFloatCol  col_energy_err  = GFitsTableFloatCol("ENERGY_ERR", size());
        GFitsTableFloatCol  col_hil_msw     = GFitsTableFloatCol("HIL_MSW", size());
        GFitsTableFloatCol  col_hil_msw_err = GFitsTableFloatCol("HIL_MSW_ERR", size());
        GFitsTableFloatCol  col_hil_msl     = GFitsTableFloatCol("HIL_MSL", size());
        GFitsTableFloatCol  col_hil_msl_err = GFitsTableFloatCol("HIL_MSL_ERR", size());

        // Fill columns
        for (int i = 0; i < size(); ++i) {
            col_eid(i)         = m_events[i].m_event_id;
            col_oid(i)         = m_events[i].m_obs_id;
            col_time(i)        = m_events[i].time().convert(m_gti.reference());
            col_live(i)        = 0.0;
            col_multip(i)      = 0;
            //col_telmask
            col_ra(i)          = m_events[i].dir().dir().ra_deg();
            col_dec(i)         = m_events[i].dir().dir().dec_deg();
            col_direrr(i)      = m_events[i].m_dir_err;
            col_detx(i)        = m_events[i].dir().detx() * gammalib::rad2deg;
            col_dety(i)        = m_events[i].dir().dety() * gammalib::rad2deg;
            col_alt(i)         = m_events[i].m_alt;
            col_az(i)          = m_events[i].m_az;
            col_corex(i)       = m_events[i].m_corex;
            col_corey(i)       = m_events[i].m_corey;
            col_core_err(i)    = m_events[i].m_core_err;
            col_xmax(i)        = m_events[i].m_xmax;
            col_xmax_err(i)    = m_events[i].m_xmax_err;
            col_shw(i)         = m_events[i].m_shwidth;
            col_shl(i)         = m_events[i].m_shlength;
            col_energy(i)      = m_events[i].energy().TeV();
            col_energy_err(i)  = m_events[i].m_energy_err;
            col_hil_msw(i)     = m_events[i].m_hil_msw;
            col_hil_msw_err(i) = m_events[i].m_hil_msw_err;
            col_hil_msl(i)     = m_events[i].m_hil_msl;
            col_hil_msl_err(i) = m_events[i].m_hil_msl_err;
        } // endfor: looped over rows

        // Append columns to table
        hdu.append(col_eid);
        hdu.append(col_oid);
        hdu.append(col_time);
        hdu.append(col_live);
        hdu.append(col_multip);
        hdu.append(col_telmask);
        hdu.append(col_ra);
        hdu.append(col_dec);
        hdu.append(col_direrr);
        hdu.append(col_detx);
        hdu.append(col_dety);
        hdu.append(col_alt);
        hdu.append(col_az);
        hdu.append(col_corex);
        hdu.append(col_corey);
        hdu.append(col_core_err);
        hdu.append(col_xmax);
        hdu.append(col_xmax_err);
        hdu.append(col_shw);
        hdu.append(col_shl);
        hdu.append(col_energy);
        hdu.append(col_energy_err);
        hdu.append(col_hil_msw);
        hdu.append(col_hil_msw_err);
        hdu.append(col_hil_msl);
        hdu.append(col_hil_msl_err);

    } // endif: there were events to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write data selection keywords into FITS HDU
 *
 * @param[in] hdu FITS HDU.
 *
 * This method does nothing if the HDU pointer is NULL.
 *
 * @todo This is a very dumb data selection keyword writing routine that does
 *       not take into account any existing keywords. We definitely want a
 *       more secure logic that checks for existing keywords and possible
 *       conflicts. But for the moment, this code does the job.
 ***************************************************************************/
void GCTAEventList::write_ds_keys(GFitsHDU& hdu) const
{
    // Set ROI parameters
    double ra  = roi().centre().dir().ra_deg();
    double dec = roi().centre().dir().dec_deg();
    double rad = roi().radius();

    // Set energy range parameters
    double e_min = emin().TeV();
    double e_max = emax().TeV();

    // Set cone selection string
    std::string dsval2 = "CIRCLE(" +
                         gammalib::str(ra) + "," +
                         gammalib::str(dec) + "," +
                         gammalib::str(rad) + ")";

    // Set energy selection string
    std::string dsval3 = gammalib::str(e_min) + ":" +
                         gammalib::str(e_max);

    // Add time selection keywords
    hdu.card("DSTYP1", "TIME",  "Data selection type");
    hdu.card("DSUNI1", "s",     "Data selection unit");
    hdu.card("DSVAL1", "TABLE", "Data selection value");
    hdu.card("DSREF1", ":GTI",  "Data selection reference");

    // Add acceptance cone selection
    hdu.card("DSTYP2", "POS(RA,DEC)", "Data selection type");
    hdu.card("DSUNI2", "deg",         "Data selection unit");
    hdu.card("DSVAL2", dsval2,        "Data selection value");

    // Add energy range selection
    hdu.card("DSTYP3", "ENERGY", "Data selection type");
    hdu.card("DSUNI3", "TeV",    "Data selection unit");
    hdu.card("DSVAL3", dsval3,   "Data selection value");

    // Set number of data selection keys
    hdu.card("NDSKEYS", 3,  "Number of data selections");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialize IRF cache for a given model
 *
 * @param[in] name Model name.
 * @return Cache index (-1 if invalid).
 ***************************************************************************/
int GCTAEventList::irf_cache_init(const std::string& name) const
{
    // Initialise cache index
    int index = irf_cache_index(name);

    // Continue only if model does not yet exist
    if (index == -1) {

        // Add model name and vector to cache. The vector is initialized
        // to -1, which signals that no cache values exist
        m_irf_names.push_back(name);
        m_irf_values.push_back(std::vector<double>(size(), -1.0));

        // Set index
        index = m_irf_names.size()-1;

    } // endif: model cache did not yet exist

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Determines the cache index for a given model name
 *
 * @param[in] name Model name.
 * @return Cache index (-1 if model has not been found).
 ***************************************************************************/
int GCTAEventList::irf_cache_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Continue only if there are models in cache
    if (!m_irf_names.empty()) {

         // Search for model name
         for (int i = 0; i < m_irf_names.size(); ++i) {
             if (m_irf_names[i] == name) {
                 index = i;
                 break;
             }
         }

    } // endif: there were models in cache

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Get cache IRF value
 *
 * @param[in] name Model name.
 * @param[in] index Event index [0,...,size()-1].
 * @return IRF value (-1 if no cache value found).
 ***************************************************************************/
double GCTAEventList::irf_cache(const std::string& name, const int& index) const
{
    // Initialise IRF value to invalid value
    double irf = -1.0;

    // Get cache index. Continue only if index is valid
    int icache = irf_cache_index(name);
    if (icache != -1) {
        irf = (m_irf_values[icache])[index];
    }

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Set cache IRF value
 *
 * @param[in] name Model name.
 * @param[in] index Event index [0,...,size()-1].
 * @param[in] irf IRF value.
 ***************************************************************************/
void GCTAEventList::irf_cache(const std::string& name, const int& index,
                              const double& irf) const
{
    // Initialize cache index. Continue only if index is valid
    int icache = irf_cache_init(name);
    if (icache != -1) {
        (m_irf_values[icache])[index] = irf;
    }

    // Return
    return;
}
