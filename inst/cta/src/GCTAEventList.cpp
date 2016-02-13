/***************************************************************************
 *                GCTAEventList.cpp - CTA event list class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @brief CTA event list class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GTime.hpp"
#include "GTimeReference.hpp"
#include "GCTAEventList.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                          "GCTAEventList::operator[](int&)"
#define G_ROI                                     "GCTAEventList::roi(GRoi&)"
#define G_FETCH                                      "GCTAEventList::fetch()"

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
 * Constructs empty event list.
 ***************************************************************************/
GCTAEventList::GCTAEventList(void) : GEventList()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File name constructor
 *
 * @param[in] filename Counts cube filename.
 *
 * Constructs event list by loading the events from a FITS file.
 ***************************************************************************/
GCTAEventList::GCTAEventList(const GFilename& filename) : GEventList()
{
    // Initialise members
    init_members();

    // Load event list
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] list Event list.
 *
 * Constructs event list by coping from another event list.
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
 *
 * Destructs event list.
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
 *
 * Assigns events from another event list.
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
 * @return Pointer to event atom.
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to an event atom.
 ***************************************************************************/
GCTAEventAtom* GCTAEventList::operator[](const int& index)
{
    // Make sure that the events are online
    fetch();

    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, "Event index",
                                       index, size());
    }
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/***********************************************************************//**
 * @brief Event atom access operator
 *
 * @param[in] index Event index [0,...,size()-1].
 * @return Pointer to event atom.
 *
 * @exception GException::out_of_range
 *            Event index outside valid range.
 *
 * Returns pointer to an event atom.
 ***************************************************************************/
const GCTAEventAtom* GCTAEventList::operator[](const int& index) const
{
    // Make sure that the events are online
    fetch();

    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OPERATOR, "Event index",
                                       index, size());
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
 * Clear event list.
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
 *
 * Returns a pointer to a deep copy of the event list.
 ***************************************************************************/
GCTAEventList* GCTAEventList::clone(void) const
{
    // Clone event list
    return new GCTAEventList(*this);
}


/***********************************************************************//**
 * @brief Load events from FITS file.
 *
 * @param[in] filename FITS file name.
 *
 * Loads the event list from a FITS file. See the read() method for details.
 ***************************************************************************/
void GCTAEventList::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename.filename());

    // Read event list
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save events into FITS file.
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing FITS file (default: false).
 *
 * Saves the event list into a FITS file. See the write() method for details.
 ***************************************************************************/
void GCTAEventList::save(const GFilename& filename,
                         const bool&      clobber) const
{
    // Open or create FITS file
    GFits fits(filename.filename(), true);

    // Write event list
    write(fits);

    // Save FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read events from FITS file.
 *
 * @param[in] fits FITS file.
 *
 * Reads the event list from a FITS file.
 *
 * The events will be read by default from the extension "EVENTS" unless
 * an extension name is explicitly specified in the FITS file name. The
 * FITS header of the extension will be scanned for data sub-space keywords
 * to extract the energy boundaries, the region of interest, as well as the
 * extension name for the Good Time Intervals. If no extension name for
 * Good Time Intervals is found it is expected that the Good Time Intervals
 * reside in the "GTI" extension.
 *
 * If a Good Time Intervals is found in the same FITS file, the Good Time
 * Intervals will be loaded. Otherwise, a single Good Time Interval will
 * be assumed based on the "TSTART" and "TSTOP" keywords found in the header
 * of the event list.
 *
 * The method clears the event list before reading, thus any events that
 * existed before in the event list will be lost. Note that the method will
 * actually not read the events but only the event metadata. The events are
 * only read from the FITS file when they are actually need. This reduces
 * the memory requirements. Events can be manually loaded into memory using
 * the fetch() method. Events can be unloaded using the dispose() method, but
 * the user has to take care of saving the events to disk in case he wants
 * to keep them for further use.
 ***************************************************************************/
void GCTAEventList::read(const GFits& fits)
{
    // Clear object
    clear();

    // Store the FITS file name
    m_filename = fits.filename();

    // Initialise events extension name
    std::string extname = fits.filename().extname("EVENTS");

    // Get event list HDU
    const GFitsTable& events = *fits.table(extname);

    // Determine number of events from event list HDU
    m_num_events = events.nrows();

    // Signal presence of "DETX" and "DETY" columns
    if (events.contains("DETX") && events.contains("DETY")) {
        m_has_detxy = true;
    }
    else {
        m_has_detxy = false;
    }

    // Signal presence of "PHASE" column
    if (events.contains("PHASE")) {
        m_has_phase = true;
    }
    else {
        m_has_phase = false;
    }

    // If the file name contains an expression then load the events now
    // (we apparently can't do this later properly, not fully sure why this
    // is the case)
    if (m_filename.has_expression()) {
        read_events(events);
    }

    // Read GTI extension name from data sub-space keyword
    std::string gti_extname = gammalib::read_ds_gti_extname(events);

    // If no GTI extension name was found then
    if (gti_extname.empty()) {
        gti_extname = "GTI";
    }

    // If GTI extension is present in FITS file then read Good Time Intervals
    // from that extension
    if (fits.contains(gti_extname)) {
        const GFitsTable& gti = *fits.table(gti_extname);
        m_gti_extname         = gti_extname;
        m_gti.read(gti);
    }

    // ... otherwise build GTI from TSTART and TSTOP
    else {

        // Read start and stop time
        double tstart = events.real("TSTART");
        double tstop  = events.real("TSTOP");

        // Create time reference from header information
        GTimeReference timeref(events);

        // Set GTI time reference
        m_gti.reference(timeref);

        // Set start and stop time
        GTime start(tstart, m_gti.reference());
        GTime stop(tstop, m_gti.reference());

        // Append start and stop time as single time interval to GTI
        m_gti.append(start, stop);


    } // endelse: GTI built from TSTART and TSTOP

    // Read region of interest from data sub-space keyword
    m_roi = gammalib::read_ds_roi(events);

    // Read energy boundaries from data sub-space keyword
    m_ebounds = gammalib::read_ds_ebounds(events);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA events and Good Time Intervals into FITS file
 *
 * @param[in] fits FITS file.
 *
 * Writes the CTA event list and the Good Time intervals into a FITS file.
 *
 * The events will be written by default into the extension "EVENTS" unless
 * an extension name is explicitly specified in the FITS file name. The
 * method also writes the data sub-space keywords in the FITS header of the
 * events table.
 *
 * In addition, the method will also append a table containing the Good Time
 * Intervals of the events to the FITS file. The extension name for the Good
 * Time Intervals is either taken from the m_gti_extname member, or if empty,
 * is set to "GTI".
 ***************************************************************************/
void GCTAEventList::write(GFits& fits) const
{
    // Set event extension name
    std::string evtname = fits.filename().extname("EVENTS");

    // Set GTI extension name
    std::string gtiname = (m_gti_extname.empty()) ? "GTI" : m_gti_extname;

    // Write events and GTIs
    write(fits, evtname, gtiname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA events and Good Time Intervals into FITS file
 *
 * @param[in] fits FITS file.
 * @param[in] evtname Event FITS extension name.
 * @param[in] gtiname Good Time Intervals FITS extension name.
 *
 * Writes the CTA event list and the Good Time intervals into a FITS file.
 *
 * The events will be written into the extension @p evtname while the Good
 * Time Intervals will be written into the extension @p gtiname.
 ***************************************************************************/
void GCTAEventList::write(GFits& fits, const std::string& evtname,
                                       const std::string& gtiname) const
{
    // Allocate empty FITS binary table
    GFitsBinTable table;

    // Write events into binary table
    write_events(table);

    // Set FITS extension name
    table.extname(evtname);

    // Write data selection keywords
    write_ds_keys(table, gtiname);

    // Append event table to FITS file
    fits.append(table);

    // Write GTI extension
    gti().write(fits, gtiname);

    // Store FITS file name
    m_filename = fits.filename();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set ROI
 *
 * @param[in] roi ROI.
 *
 * @exception GException::invalid_argument
 *            Invalid region of interest class specified.
 *
 * Sets the region of interest of the event list.
 ***************************************************************************/
void GCTAEventList::roi(const GRoi& roi)
{
    // Cast ROI dynamically
    const GCTARoi* ptr = dynamic_cast<const GCTARoi*>(&roi);

    // Throw exception if ROI is not of correct type
    if (ptr == NULL) {
        std::string cls = std::string(typeid(&roi).name());
        std::string msg = "Invalid region of interest type \""+cls+"\" "
                          "provided on input. Please specify a \"GCTARoi\" "
                          "object as argument.";
        throw GException::invalid_argument(G_ROI, msg);
    }

    // Set ROI
    m_roi = *ptr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append event to event list
 *
 * @param[in] event Event.
 *
 * Appends an event to the end of the event list.
 ***************************************************************************/
void GCTAEventList::append(const GCTAEventAtom& event)
{
    // Make sure that the events are online
    fetch();

    // Append event
    m_events.push_back(event);

    // Append an element to all additional columns
    for (int i = 0; i < m_columns.size(); ++i) {
        m_columns[i]->insert(m_columns[i]->length(),1);
    }

    // Set number of events
    m_num_events = m_events.size();

    // Set event index
    int index               = m_num_events - 1;
    m_events[index].m_index = index;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove events from event list
 *
 * @param[in] index Index from which on events should be removed.
 * @param[in] number Number of event to remove (default: 1).
 *
 * Removes events from the event list. This method does nothing if @p index
 * points beyond the event list. The method does also gently accept
 * @p number arguments where @p index + @p number reach beyond the event
 * list. In that case, all events from event @p index on will be removed.
 ***************************************************************************/
void GCTAEventList::remove(const int& index, const int& number)
{
    // Make sure that the events are online
    fetch();

    // Continue only if index is valid
    if (index < size()) {

        // Determine number of elements to remove
        int n_remove = (index + number > size()) ? size() - index : number;

        // Remove events
        m_events.erase(m_events.begin() + index,
                       m_events.begin() + index + n_remove);

        // Remove elements from additional columns
        for (int i = 0; i < m_columns.size(); ++i) {
            m_columns[i]->remove(index, n_remove);
        }

        // Set number of events
        m_num_events = m_events.size();

    } // endif: index was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch events
 *
 * @exception GException::file_error
 *            File not found.
 *            Unable to load events.
 * @exception GException::invalid_value
 *            No file name has been specified.
 *            Fetched number of events mismatches the expectation.
 *
 * Fetches the events by reading them from a FITS file. This method does
 * nothing if the events are already loaded, if there is nothing to fetch,
 * or if the m_filename member is empty.
 *
 * The method is thread save. The method furthermore checks whether the
 * loaded number of events corresponds to the expectation. If also checks
 * whether the file from which events should be loaded actually exists.
 ***************************************************************************/
void GCTAEventList::fetch(void) const
{
    // Continue only if events are not loaded but there are events to fetch
    if (m_events.empty() && size() > 0) {

        // Continue only if the file name is not empty
        if (!m_filename.empty()) {

            // Throw an exception if the file does not exist
            if (!gammalib::file_exists_gzip(m_filename.filename())) {
                std::string msg = "File \""+m_filename.filename()+"\" not "
                                  "found. Cannot fetch events. Maybe the "
                                  "file has been deleted in the meantime.";
                GException::file_error(G_FETCH, msg);
            }

            // Initialise exception flag
            bool has_exception = false;

            // Load events. Catch any exception. Put the code into a critical
            // zone as it might be called from within a parallelized thread.
            #pragma omp critical
            {
            try {

                // Open FITS file
                GFits fits(m_filename.fullname());

                // Initialise events extension name
                std::string extname = fits.filename().extname("EVENTS");

                // Get event list HDU
                const GFitsTable& events = *fits.table(extname);

                // Load event data
                read_events(events);

                // Close FITS file
                fits.close();


            }
            catch (...) {
                has_exception = true;
            }
            }

            // Throw an exception if an exception has occured
            if (has_exception) {
                std::string msg = "Unable to load events from file \""+
                                  m_filename.filename()+"\"."; 
                throw GException::file_error(G_FETCH, msg);
            }

            // Throw an exception if the number of events does not correspond
            // to the expectation
            if (m_events.size() != size()) {
                std::string msg = "Loaded "+gammalib::str(m_events.size())+
                                  " events from FITS file while "+
                                  gammalib::str(size())+" events were "
                                  "expected. This should never happen!"; 
                throw GException::invalid_value(G_FETCH, msg);
            }

        
        } // endif: filename was not empty

        // Throw an exception if the FITS file name is not known
        else {
            std::string msg = "Unable to fetch "+gammalib::str(size())+
                              " events since the name of the FITS file from "
                              "which to fetch the events is not known.";
            throw GException::invalid_value(G_FETCH, msg);
        }

    } // endif: there were no events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Dispose events
 *
 * Disposes the events to save memory. The method free's all memory
 * associated to the events. The method should be used with care as it does
 * not check whether the events can actually be recovered using the fetch()
 * method. Eventually, the events should be saved before so that fetch()
 * can recover them.
 ***************************************************************************/
void GCTAEventList::dispose(void) const
{
    // Free additional columns
    for (int i = 0; i < m_columns.size(); ++i) {
        delete m_columns[i];
        m_columns[i] = NULL;
    }

    // Clear events and additional columns
    m_events.clear();
    m_columns.clear();

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

        // Append number of events and disposal information
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(size()));
        if (m_events.empty()) {
            if (!m_filename.empty()) {
                result.append(" (disposed in \"");
                result.append(m_filename.filename());
                result.append("\")");
            }
            else {
                result.append(" (disposed without possibility to recover)");
            }
        }
        else {
            result.append(" (loaded)");
        }

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
        else {
            result.append("\n"+gammalib::parformat("Energy interval"));
            if (ebounds().size() > 0) {
                result.append(gammalib::str(emin().TeV()));
                result.append(" - ");
                result.append(gammalib::str(emax().TeV())+" TeV");
            }
            else {
                result.append("not defined");
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
        else {
            result.append("\n"+gammalib::parformat("Region of interest"));
            if (roi().radius() > 0) {
                result.append(roi().centre().print());
                result.append(" Radius=");
                result.append(gammalib::str(roi().radius())+" deg");
            }
            else {
                result.append("not defined");
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
    m_columns.clear();
    m_filename.clear();
    m_num_events  = 0;
    m_gti_extname = "GTI"; //!< Default GTI extension name
    m_has_phase   = false;
    m_has_detxy   = true;

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
    m_roi         = list.m_roi;
    m_events      = list.m_events;
    m_filename    = list.m_filename;
    m_num_events  = list.m_num_events;
    m_gti_extname = list.m_gti_extname;
    m_has_phase   = list.m_has_phase;
    m_has_detxy   = list.m_has_detxy;

    // Copy cache
    m_irf_names  = list.m_irf_names;
    m_irf_values = list.m_irf_values;

    // Copy column pointers
    m_columns.clear();
    for (int i = 0; i < list.m_columns.size(); ++i) {
        m_columns.push_back((list.m_columns[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventList::free_members(void)
{
    // Free columns
    for (int i = 0; i < m_columns.size(); ++i) {
        delete m_columns[i];
        m_columns[i] = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS table
 *
 * @param[in] table FITS table.
 *
 * This method reads the CTA event list from a FITS table into memory.
 * The following columns are mandatory in the FITS table:
 *
 *      - EVENT_ID
 *      - TIME
 *      - RA
 *      - DEC
 *      - ENERGY
 *
 * Optionally, the
 *
 *      - DETX
 *      - DETY
 *      - PHASE
 *
 * columns are read. Any other columns present in the FITS table will also be
 * read, and written into the FITS table when the write_events() method is
 * called.
 ***************************************************************************/
void GCTAEventList::read_events(const GFitsTable& table) const
{
    // Dispose any existing events
    dispose();

    // Extract number of events in FITS file
    int num = table.nrows();

    // If there are events then load them
    if (num > 0) {

        // Reserve data
        m_events.reserve(num);

        // Get column pointers
        const GFitsTableCol* ptr_eid    = table["EVENT_ID"];
        const GFitsTableCol* ptr_time   = table["TIME"];
        const GFitsTableCol* ptr_ra     = table["RA"];
        const GFitsTableCol* ptr_dec    = table["DEC"];
        const GFitsTableCol* ptr_energy = table["ENERGY"];

        // Optionally get pointers to DETX and DETY columns
        const GFitsTableCol* ptr_detx = (m_has_detxy) ? table["DETX"] : NULL;
        const GFitsTableCol* ptr_dety = (m_has_detxy) ? table["DETY"] : NULL;

        // Optionally get pointer to PHASE column
        const GFitsTableCol* ptr_phase = (m_has_phase) ? table["PHASE"] : NULL;

        // Copy data from columns into GCTAEventAtom objects
        for (int i = 0; i < num; ++i) {

            // Allocate event
            GCTAEventAtom event;

            // Set mandatory information
            event.m_time.set(ptr_time->real(i), m_gti.reference());
            event.m_dir.dir().radec_deg(ptr_ra->real(i), ptr_dec->real(i));
            event.m_energy.TeV(ptr_energy->real(i));
            event.m_event_id = ptr_eid->integer(i);
            event.m_index    = i;

            // If available, set detector coordinates in radians
            if (m_has_detxy) {
                event.m_dir.detx(ptr_detx->real(i)*gammalib::deg2rad);
                event.m_dir.dety(ptr_dety->real(i)*gammalib::deg2rad);
            }

            // If available, set pulse phase
            if (m_has_phase) {
                event.m_phase = ptr_phase->real(i);
            }

            // Append event
            m_events.push_back(event);

        } // endfor: looped over all events

        // Loop over table and find optional columns
        for (int i = 0; i < table.ncols(); ++i) {

            // Get column pointer
            const GFitsTableCol* col = table[i];

            // Get column name
            const std::string name = col->name();

            // If column was mandatory or optional then skip it ...
            if (name == "EVENT_ID" ||
                name == "TIME"     ||
                name == "RA"       ||
                name == "DEC"      ||
                name == "ENERGY"   ||
                name == "DETX"     ||
                name == "DETY"     ||
                name == "PHASE") {
                continue;
            }

            // ... otherwise keep a copy
            else {
                m_columns.push_back(col->clone());
            }

        } // endfor: loop over optional columns

    } // endif: there were events

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write CTA events into FITS table
 *
 * @param[in] hdu FITS binary table.
 *
 * Write the CTA event list into a FITS binary table.
 *
 * The following mandatory columns will be written into the FITS table:
 *
 *      - EVENT_ID
 *      - TIME
 *      - RA
 *      - DEC
 *      - ENERGY
 *
 * If available, also the
 *
 *      - DETX
 *      - DETY
 *      - PHASE
 *
 * columns are written. Any other columns that were read by the read_events()
 * method will be also written into the table.
 ***************************************************************************/
void GCTAEventList::write_events(GFitsBinTable& hdu) const
{
    // Make sure that the events are online
    fetch();

    // Set extension name
    hdu.extname("EVENTS");

    // If there are events then write them now
    if (size() > 0) {

        // Allocate mandatory columns
        GFitsTableULongCol  col_eid    = GFitsTableULongCol("EVENT_ID", size());
        GFitsTableDoubleCol col_time   = GFitsTableDoubleCol("TIME", size());
        GFitsTableFloatCol  col_ra     = GFitsTableFloatCol("RA", size());
        GFitsTableFloatCol  col_dec    = GFitsTableFloatCol("DEC", size());
        GFitsTableFloatCol  col_energy = GFitsTableFloatCol("ENERGY", size());

        // Set units of columns
        // (see http://fits.gsfc.nasa.gov/standard30/fits_standard30aa.pdf)
        col_time.unit("s");
        col_ra.unit("deg");
        col_dec.unit("deg");
        col_energy.unit("TeV");

        // Fill mandatory columns
        for (int i = 0; i < size(); ++i) {
            col_eid(i)    = m_events[i].m_event_id;
            col_time(i)   = m_events[i].time().convert(m_gti.reference());
            col_ra(i)     = m_events[i].dir().dir().ra_deg();
            col_dec(i)    = m_events[i].dir().dir().dec_deg();
            col_energy(i) = m_events[i].energy().TeV();
        }

        // Append mandatory columns to table
        hdu.append(col_eid);
        hdu.append(col_time);
        hdu.append(col_ra);
        hdu.append(col_dec);
        hdu.append(col_energy);

        // If available, add detector coordinates in degrees
        if (m_has_detxy) {

            // Allocate columns
            GFitsTableFloatCol col_detx = GFitsTableFloatCol("DETX", size());
            GFitsTableFloatCol col_dety = GFitsTableFloatCol("DETY", size());

            // Set units of columns
            // (see http://fits.gsfc.nasa.gov/standard30/fits_standard30aa.pdf)
            col_detx.unit("deg");
            col_dety.unit("deg");

            // Fill columns
            for (int i = 0; i < size(); ++i) {
                col_detx(i) = m_events[i].dir().detx() * gammalib::rad2deg;
                col_dety(i) = m_events[i].dir().dety() * gammalib::rad2deg;
            }

            // Append columns to table
            hdu.append(col_detx);
            hdu.append(col_dety);

        }

        // If available, add event phase
        if (m_has_phase) {

            // Allocate columns
            GFitsTableFloatCol col_phase = GFitsTableFloatCol("PHASE", size());

            // Fill columns
            for (int i = 0; i < size(); ++i) {
                col_phase(i) = m_events[i].m_phase;
            }

            // Append columns to table
            hdu.append(col_phase);

        }

        // Append other columns to table
        for (int i = 0; i < m_columns.size(); ++i) {
            hdu.append(*m_columns[i]);
        }

    } // endif: there were events to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write data sub-space keywords into FITS HDU
 *
 * @param[in] hdu FITS HDU.
 * @param[in] gtiname Good Time Interval FITS extension.
 *
 * Writes the data sub-space keywords for an event list into the FITS HDU.
 * The following keywords will be written:
 *
 *      DSTYP1 = "TIME"                     / Data sub-space type
 *      DSUNI1 = "s"                        / Data sub-space unit
 *      DSVAL1 = "TABLE"                    / Data sub-space value
 *      DSREF1 = ":[extname]"               / Data sub-space reference
 *      DSTYP2 = "ENERGY"                   / Data sub-space type
 *      DSUNI2 = "TeV"                      / Data sub-space unit
 *      DSVAL2 = "[emin]:[emax]"            / Data sub-space value
 *      DSTYP3 = "POS(RA,DEC)"              / Data sub-space type
 *      DSUNI3 = "deg"                      / Data sub-space unit
 *      DSVAL3 = "CIRCLE([ra],[dec],[rad])" / Data sub-space value
 *
 * where
 *
 *      [extname] is the GTI extension name @p gtiname
 *      [emin] is the minimum event energy in TeV
 *      [emax] is the maximum event energy in TeV
 *      [ra] is the Right Ascension of the Region of Interest centre in degrees
 *      [dec] is the Declination of the Region of Interest centre in degrees
 *      [rad] is the radius of the Region of Interest in degrees
 ***************************************************************************/
void GCTAEventList::write_ds_keys(GFitsHDU& hdu, const std::string& gtiname) const
{
    // Set ROI parameters
    double ra  = roi().centre().dir().ra_deg();
    double dec = roi().centre().dir().dec_deg();
    double rad = roi().radius();

    // Set energy range parameters
    double e_min = emin().TeV();
    double e_max = emax().TeV();

    // Set energy selection string
    std::string dsval2 = gammalib::str(e_min) + ":" +
                         gammalib::str(e_max);

    // Set Good Time Intervals extension name
    std::string dsref1 = ":"+gtiname;

    // Add time selection keywords
    hdu.card("DSTYP1", "TIME",  "Data sub-space type");
    hdu.card("DSUNI1", "s",     "Data sub-space unit");
    hdu.card("DSVAL1", "TABLE", "Data sub-space value");
    hdu.card("DSREF1", dsref1,  "Data sub-space reference");

    // Add energy range selection
    hdu.card("DSTYP2", "ENERGY", "Data sub-space type");
    hdu.card("DSUNI2", "TeV",    "Data sub-space unit");
    hdu.card("DSVAL2", dsval2,   "Data sub-space value");

    // Initialise number of NDSKEYS
    int ndskeys = 2;

    // Add acceptance cone only if RoI information is valid
    if (m_roi.is_valid()) {

        // Set cone selection string
        std::string dsval3 = "CIRCLE(" +
                             gammalib::str(ra) + "," +
                             gammalib::str(dec) + "," +
                             gammalib::str(rad) + ")";

        // Write DS keywords
        hdu.card("DSTYP3", "POS(RA,DEC)", "Data sub-space type");
        hdu.card("DSUNI3", "deg",         "Data sub-space unit");
        hdu.card("DSVAL3", dsval3,        "Data sub-space value");
        ndskeys++;
        
    } // endif: RoI was valid

    // Set number of data selection keys
    hdu.card("NDSKEYS", ndskeys,  "Number of data sub-space keys");

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
