/***************************************************************************
 *                GCTAEventList.cpp  -  CTA event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAEventList.cpp
 * @brief GCTAEventList class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAEventList.hpp"
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableStringCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_POINTER                               "GCTAEventList::pointer(int)"

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
    // Initialise class members for clean destruction
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
    // Initialise class members for clean destruction
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
 ***************************************************************************/
GCTAEventList& GCTAEventList::operator= (const GCTAEventList& list)
{
    // Execute only if object is not identical
    if (this != &list) {

        // Copy base class members
        this->GEventList::operator=(list);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(list);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
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
 * @brief Clone object
***************************************************************************/
GCTAEventList* GCTAEventList::clone(void) const
{
    return new GCTAEventList(*this);
}


/***********************************************************************//**
 * @brief Load events from event FITS file.
 *
 * @param[in] filename Name of FITS file from which events are loaded.
 *
 * Load events from CTA event file (extension EVENTS). The method clears the
 * object before loading, thus any events residing in the object before
 * loading will be lost.
 ***************************************************************************/
void GCTAEventList::load(const std::string& filename)
{
    // Clear object
    clear();

    // Open FITS file
    GFits file(filename);

    // Get HDU
    GFitsTable* hdu = file.table("EVENTS");

    // Load columns
    read(hdu);

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
void GCTAEventList::save(const std::string& filename, bool clobber) const
{
    // Open/create FITS file
    GFits fits(filename);

    // Write CTA events into FITS file
    write(&fits);
    
    // Save FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read CTA events from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * This method read the CTA event list from a FITS table HDU into memory.
 *
 * @todo Implement agreed column format
 ***************************************************************************/
void GCTAEventList::read(GFitsTable* hdu)
{
    // Clear existing events
    m_events.clear();

    // Allocate space for keyword name
    char keyword[10];

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract number of events in FITS file
        int num = hdu->integer("NAXIS2");

        // If there are events then load them
        if (num > 0) {

            // Reserve data
            m_events.reserve(num);

            // Get column pointers
            GFitsTableULongCol*  ptr_eid         = (GFitsTableULongCol*)hdu->column("EVENT_ID");
            GFitsTableDoubleCol* ptr_time        = (GFitsTableDoubleCol*)hdu->column("TIME");
            GFitsTableShortCol*  ptr_multip      = (GFitsTableShortCol*)hdu->column("MULTIP");
            GFitsTableFloatCol*  ptr_ra          = (GFitsTableFloatCol*)hdu->column("RA");
            GFitsTableFloatCol*  ptr_dec         = (GFitsTableFloatCol*)hdu->column("DEC");
            GFitsTableFloatCol*  ptr_dir_err     = (GFitsTableFloatCol*)hdu->column("DIR_ERR");
            GFitsTableFloatCol*  ptr_detx        = (GFitsTableFloatCol*)hdu->column("DETX");
            GFitsTableFloatCol*  ptr_dety        = (GFitsTableFloatCol*)hdu->column("DETY");
            GFitsTableFloatCol*  ptr_alt_pnt     = (GFitsTableFloatCol*)hdu->column("ALT_PNT");
            GFitsTableFloatCol*  ptr_az_pnt      = (GFitsTableFloatCol*)hdu->column("AZ_PNT");
            GFitsTableFloatCol*  ptr_alt         = (GFitsTableFloatCol*)hdu->column("ALT");
            GFitsTableFloatCol*  ptr_az          = (GFitsTableFloatCol*)hdu->column("AZ");
            GFitsTableFloatCol*  ptr_corex       = (GFitsTableFloatCol*)hdu->column("COREX");
            GFitsTableFloatCol*  ptr_corey       = (GFitsTableFloatCol*)hdu->column("COREY");
            GFitsTableFloatCol*  ptr_core_err    = (GFitsTableFloatCol*)hdu->column("CORE_ERR");
            GFitsTableFloatCol*  ptr_xmax        = (GFitsTableFloatCol*)hdu->column("XMAX");
            GFitsTableFloatCol*  ptr_xmax_err    = (GFitsTableFloatCol*)hdu->column("XMAX_ERR");
            GFitsTableFloatCol*  ptr_energy      = (GFitsTableFloatCol*)hdu->column("ENERGY");
            GFitsTableFloatCol*  ptr_energy_err  = (GFitsTableFloatCol*)hdu->column("ENERGY_ERR");
            GFitsTableFloatCol*  ptr_hil_msw     = (GFitsTableFloatCol*)hdu->column("HIL_MSW");
            GFitsTableFloatCol*  ptr_hil_msw_err = (GFitsTableFloatCol*)hdu->column("HIL_MSW_ERR");
            GFitsTableFloatCol*  ptr_hil_msl     = (GFitsTableFloatCol*)hdu->column("HIL_MSL");
            GFitsTableFloatCol*  ptr_hil_msl_err = (GFitsTableFloatCol*)hdu->column("HIL_MSL_ERR");

            // Copy data from columns into GCTAEventAtom objects
            GCTAEventAtom event;
            for (int i = 0; i < num; ++i) {
                event.m_time.met((*ptr_time)(i));
                event.m_dir.radec_deg((*ptr_ra)(i), (*ptr_dec)(i));
                event.m_energy.TeV((*ptr_energy)(i));
                event.m_event_id    = (*ptr_eid)(i);
                event.m_multip      = (*ptr_multip)(i);
                event.m_dir_err     = (*ptr_dir_err)(i);
                event.m_detx        = (*ptr_detx)(i);
                event.m_dety        = (*ptr_dety)(i);
                event.m_alt_pnt     = (*ptr_alt_pnt)(i);
                event.m_az_pnt      = (*ptr_az_pnt)(i);
                event.m_alt         = (*ptr_alt)(i);
                event.m_az          = (*ptr_az)(i);
                event.m_corex       = (*ptr_corex)(i);
                event.m_corey       = (*ptr_corey)(i);
                event.m_core_err    = (*ptr_core_err)(i);
                event.m_xmax        = (*ptr_xmax)(i);
                event.m_xmax_err    = (*ptr_xmax_err)(i);
                event.m_energy_err  = (*ptr_energy_err)(i);
                event.m_hil_msw     = (*ptr_hil_msw)(i);
                event.m_hil_msw_err = (*ptr_hil_msw_err)(i);
                event.m_hil_msl     = (*ptr_hil_msl)(i);
                event.m_hil_msl_err = (*ptr_hil_msl_err)(i);
                m_events.push_back(event);
            }

        } // endif: there were events

    } // endif: HDU was valid

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
void GCTAEventList::write(GFits* file) const
{
    // Allocate FITS binary table HDU
    GFitsBinTable* hdu = new GFitsBinTable;

    // Set extension name
    hdu->extname("EVENTS");

    // Write header
    write_header(hdu);

    // If there are events then write them now
    if (size() > 0) {

        // Allocate columns
        GFitsTableULongCol  col_eid         = GFitsTableULongCol("EVENT_ID", size());
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
        GFitsTableFloatCol  col_alt_pnt     = GFitsTableFloatCol("ALT_PNT", size());
        GFitsTableFloatCol  col_az_pnt      = GFitsTableFloatCol("AZ_PNT", size());
        GFitsTableFloatCol  col_corex       = GFitsTableFloatCol("COREX", size());
        GFitsTableFloatCol  col_corey       = GFitsTableFloatCol("COREY", size());
        GFitsTableFloatCol  col_core_err    = GFitsTableFloatCol("CORE_ERR", size());
        GFitsTableFloatCol  col_xmax        = GFitsTableFloatCol("XMAX", size());
        GFitsTableFloatCol  col_xmax_err    = GFitsTableFloatCol("XMAX_ERR", size());
        GFitsTableFloatCol  col_energy      = GFitsTableFloatCol("ENERGY", size());
        GFitsTableFloatCol  col_energy_err  = GFitsTableFloatCol("ENERGY_ERR", size());
        GFitsTableFloatCol  col_hil_msw     = GFitsTableFloatCol("HIL_MSW", size());
        GFitsTableFloatCol  col_hil_msw_err = GFitsTableFloatCol("HIL_MSW_ERR", size());
        GFitsTableFloatCol  col_hil_msl     = GFitsTableFloatCol("HIL_MSL", size());
        GFitsTableFloatCol  col_hil_msl_err = GFitsTableFloatCol("HIL_MSL_ERR", size());

        // Fill columns
        for (int i = 0; i < size(); ++i) {
            col_eid(i)         = i;
            col_time(i)        = m_events[i].time().met();
            col_live(i)        = 0.0;
            col_multip(i)      = 0;
            //col_telmask
            col_ra(i)          = m_events[i].dir().ra_deg();
            col_dec(i)         = m_events[i].dir().dec_deg();
            col_direrr(i)      = 0.0;
            col_detx(i)        = 0.0;
            col_dety(i)        = 0.0;
            col_alt(i)         = 0.0;
            col_az(i)          = 0.0;
            col_alt_pnt(i)     = 0.0;
            col_az_pnt(i)      = 0.0;
            col_corex(i)       = 0.0;
            col_corey(i)       = 0.0;
            col_core_err(i)    = 0.0;
            col_xmax(i)        = 0.0;
            col_xmax_err(i)    = 0.0;
            col_energy(i)      = m_events[i].energy().TeV();
            col_energy_err(i)  = 0.0;
            col_hil_msw(i)     = 0.0;
            col_hil_msw_err(i) = 0.0;
            col_hil_msl(i)     = 0.0;
            col_hil_msl_err(i) = 0.0;
        } // endfor: looped over rows

        // Append columns to table
        hdu->append_column(col_eid);
        hdu->append_column(col_time);
        hdu->append_column(col_live);
        hdu->append_column(col_multip);
        hdu->append_column(col_telmask);
        hdu->append_column(col_ra);
        hdu->append_column(col_dec);
        hdu->append_column(col_direrr);
        hdu->append_column(col_detx);
        hdu->append_column(col_dety);
        hdu->append_column(col_alt);
        hdu->append_column(col_az);
        hdu->append_column(col_alt_pnt);
        hdu->append_column(col_az_pnt);
        hdu->append_column(col_corex);
        hdu->append_column(col_corey);
        hdu->append_column(col_core_err);
        hdu->append_column(col_xmax);
        hdu->append_column(col_xmax_err);
        hdu->append_column(col_energy);
        hdu->append_column(col_energy_err);
        hdu->append_column(col_hil_msw);
        hdu->append_column(col_hil_msw_err);
        hdu->append_column(col_hil_msl);
        hdu->append_column(col_hil_msl_err);
        
    } // endif: there were events to write

    // Append HDU to FITS file. The FITS class will later handle the
    // proper deallocation of the HDU
    file->append(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to an event
 *
 * @param[in] index Event index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Event index not in valid range.
 *
 * Returns a pointer to an event atom.
 ***************************************************************************/
GCTAEventAtom* GCTAEventList::pointer(int index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_POINTER, index, 0, size()-1);
    #endif

    // Return pointer
    return (&(m_events[index]));
}


/***********************************************************************//**
 * @brief Print event list information
 ***************************************************************************/
std::string GCTAEventList::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAEventList ===\n");
    result.append(parformat("Number of events")+str(size()));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Append event to event list
 *
 * @param[in] event Event.
 *
 * Appends an event atom to the event list.
 ***************************************************************************/
void GCTAEventList::append(const GCTAEventAtom& event)
{
    // Append event
    m_events.push_back(event);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserves space for events
 *
 * @param[in] number Number of events.
 *
 * Reserves space for number events in the event list.
 ***************************************************************************/
void GCTAEventList::reserve(const int& number)
{
    // Reserve space
    m_events.reserve(number);

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
    m_events.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] list GCTAEventList members which should be copied.
 ***************************************************************************/
void GCTAEventList::copy_members(const GCTAEventList& list)
{
    // Copy member
    m_events = list.m_events;

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
 * @brief Write CTA events file header
 *
 * @param[in] hdu FITS events table HDU.
 ***************************************************************************/
void GCTAEventList::write_header(GFitsBinTable* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Set observation information
        hdu->card("CREATOR",  "GammaLib", "Program which created the file");
        hdu->card("TELESCOP", "CTA",      "Telescope");
        hdu->card("EXTNAME",  "EVENTS",   "Extension name");
        hdu->card("OBS_ID",   0,          "Observation identifier");
        hdu->card("DATE_OBS", "string",   "Observation start date");
        hdu->card("TIME_OBS", "string",   "Observation start time");
        hdu->card("DATE_END", "string",   "Observation end date");
        hdu->card("TIME_END", "string",   "Observation end time");

        // Set observation time information
        hdu->card("TSTART",   0.0,     "[s] Mission time of start of observation");
        hdu->card("TSTOP",    0.0,     "[s] Mission time of end of observation");
        hdu->card("MJDREFI",  51910,   "[days] Integer part of mission time reference MJD");
        hdu->card("MJDREFF",  7.428703703703703e-14, "[days] Fractional part of mission time reference MJD");
        hdu->card("TIMEUNIT", "s",     "Time unit");
        hdu->card("TIMESYS",  "TT",    "Time system");
        hdu->card("TIMEREF",  "LOCAL", "Time reference");
        hdu->card("TELAPSE",  0.0,     "[s] Mission elapsed time");
        hdu->card("ONTIME",   0.0,     "[s] Total good time including deadtime");
        hdu->card("LIVETIME", 0.0,     "[s] Total livetime");
        hdu->card("DEADC",    0.0,     "Deadtime fraction");
        hdu->card("TIMEDEL",  1.0,     "Time resolution");

        // Set pointing information
        hdu->card("OBJECT",   "string", "Observed object");
        hdu->card("RA_OBJ",   0.0,      "[deg] Target Right Ascension");
        hdu->card("DEC_OBJ",  0.0,      "[deg] Target Declination");
        hdu->card("RA_PNT",   0.0,      "[deg] Pointing Right Ascension");
        hdu->card("DEC_PNT",  0.0,      "[deg] Pointing Declination");
        hdu->card("ALT_PNT",  0.0,      "[deg] Average altitude of pointing");
        hdu->card("AZ_PNT",   0.0,      "[deg] Average azimuth of pointing");
        hdu->card("RADECSYS", "FK5",    "Coordinate system");
        hdu->card("EQUINOX",  2000.0,   "Epoch");
        hdu->card("CONV_DEP", 0.0,      "Convergence depth of telescopes");
        hdu->card("CONV_RA",  0.0,      "[deg] Convergence Right Ascension");
        hdu->card("CONV_DEC", 0.0,      "[deg] Convergence Declination");
        hdu->card("OBSERVER", "string", "Observer");

        // Telescope information
        hdu->card("N_TELS",   100,      "Number of telescopes in event list");
        hdu->card("TELLIST",  "string", "Telescope IDs");
        hdu->card("GEOLAT",   0.0,      "[deg] Geographic latitude of array centre");
        hdu->card("GEOLON",   0.0,      "[deg] Geographic longitude of array centre");
        hdu->card("ALTITUDE", 0.0,      "[km] Altitude of array centre");

        // Other information
        hdu->card("EUNIT",    "TeV",    "Energy unit");
        hdu->card("EVTVER",   "draft1", "Event list version number");
        
    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
