/***************************************************************************
 *           GCTAEventList.cpp  -  CTA event atom container class          *
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
 * @brief CTA event atom container class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAEventList.hpp"
#include "GCTAException.hpp"
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableStringCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                          "GCTAEventList::operator[](int&)"
#define G_ROI                                     "GCTAEventList::roi(GRoi&)"
#define G_READ_DS_EBOUNDS         "GCTAEventList::read_ds_ebounds(GFitsHDU*)"
#define G_READ_DS_ROI                 "GCTAEventList::read_ds_roi(GFitsHDU*)"

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
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
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
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_OPERATOR, index, 0, size()-1);
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
void GCTAEventList::save(const std::string& filename, bool clobber) const
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
 * @param[in] file FITS file.
 *
 * This method read the CTA event list from a FITS file.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GCTAEventList::read(const GFits& file)
{
    // Clear object
    clear();

    // Get HDUs
    GFitsTable* events = file.table("EVENTS");
    GFitsTable* gti    = file.table("GTI");

    // Load event data
    read_events(events);

    // Read region of interest from data selection keyword
    read_ds_roi(events);

    // Read energy boundaries from data selection keyword
    read_ds_ebounds(events);

    // Read Good Time Intervals
    m_gti.read(gti);

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
    write_events(events);

    // Write data selection keywords
    write_ds_keys(events);

    // Append event table to FITS file
    file.append(*events);

    // Free binary table
    delete events;

    // Append GTI to FITS file
    gti().write(&file);

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
    if (ptr == NULL)
        throw GCTAException::bad_roi_type(G_ROI);

    // Set ROI
    m_roi = *ptr;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print event list information
 ***************************************************************************/
std::string GCTAEventList::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAEventList ===");
    result.append("\n"+parformat("Number of events")+str(size()));
    result.append("\n"+parformat("Time interval"));
    if (gti().size() > 0)
        result.append(str(tstart().met())+" - "+str(tstop().met())+" sec");
    else
        result.append("not defined");
    result.append("\n"+parformat("Energy range"));
    if (ebounds().size() > 0)
        result.append(emin().print()+" - "+emax().print());
    else
        result.append("not defined");
    if (roi().radius() > 0)
        result.append("\n"+roi().print());
    else
        result.append("\n"+parformat("Region of interest")+"not defined");

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
    m_roi.clear();
    m_events.clear();

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
 * @param[in] table FITS table pointer.
 *
 * This method reads the CTA event list from a FITS table HDU into memory.
 *
 * @todo Implement agreed column format
 * @todo Declare GFitsTable::column() as const, or at least, add a const
 *       version. Even better, add column access operators instead of the
 *       column() method as GFitsTable is a container for columns, and
 *       each container should have an access operator.
 ***************************************************************************/
void GCTAEventList::read_events(const GFitsTable* table)
{
    // Clear existing events
    m_events.clear();

    // Allocate space for keyword name
    char keyword[10];

    // Continue only if HDU is valid
    if (table != NULL) {

        // Circumvent const correctness. We need this because the column()
        // method is not declared const. This should be corrected.
        GFitsTable* hdu = (GFitsTable*)table;

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
 * @brief Read energy boundary data selection keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GCTAException::no_ebds
 *            Invalid energy data selection encountered
 *
 * Reads the energy boundary data selection keywords by searching for a 
 * DSTYPx keyword named "ENERGY". The data selection information is expected
 * to be in the format "200:50000", where the 2 arguments are the minimum
 * and maximum energy. The energy unit is given by the keyword DSUNIx, which
 * supports keV, MeV, GeV and TeV (case independent). No detailed syntax
 * checking is performed.
 *
 * This method clears the energy boundaries, irrespectively of whether energy
 * boundary information has been found in the HDU.
 ***************************************************************************/
void GCTAEventList::read_ds_ebounds(const GFitsHDU* hdu)
{
    // Reset energy boundaries
    m_ebounds.clear();

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get number of data selection keywords
        int ndskeys = 0;
        try {
            ndskeys = hdu->integer("NDSKEYS");
        }
        catch (GException::fits_key_not_found) {
            ;
        }

        // Loop over all data selection keys
        for (int i = 1; i <= ndskeys; ++i) {
            std::string type_key  = "DSTYP"+str(i);
            std::string unit_key  = "DSUNI"+str(i);
            std::string value_key = "DSVAL"+str(i);
            try {
                if (hdu->string(type_key) == "ENERGY") {
                    std::string unit                 = toupper(hdu->string(unit_key));
                    std::string value                = hdu->string(value_key);
                    std::vector<std::string> ebounds = split(value, ":");
                    if (ebounds.size() == 2) {
                        double  emin = todouble(ebounds[0]);
                        double  emax = todouble(ebounds[1]);
                        GEnergy e_min;
                        GEnergy e_max;
                        if (unit == "KEV") {
                            e_min.keV(emin);
                            e_max.keV(emax);
                        }
                        else if (unit == "MEV") {
                            e_min.MeV(emin);
                            e_max.MeV(emax);
                        }
                        else if (unit == "GEV") {
                            e_min.GeV(emin);
                            e_max.GeV(emax);
                        }
                        else if (unit == "TEV") {
                            e_min.TeV(emin);
                            e_max.TeV(emax);
                        }
                        else {
                            throw GCTAException::no_ebds(G_READ_DS_EBOUNDS,
                                  "Invalid energy unit \""+unit+
                                  "\" encountered in data selection key \""+
                                  unit_key+"\"");
                        }
                        m_ebounds.append(e_min, e_max);
                    }
                    else {
                        throw GCTAException::no_ebds(G_READ_DS_EBOUNDS,
                              "Invalid energy value \""+value+
                              "\" encountered in data selection key \""+
                              value_key+"\"");
                    }
                } // endif: ENERGY type found
            }
            catch (GException::fits_key_not_found) {
                ;
            }
        } // endfor: looped over data selection keys

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read ROI data selection keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GException::no_roi
 *            Invalid ROI data selection encountered
 *
 * Reads the ROI data selection keywords by searching for a DSTYPx keyword
 * named "POS(RA,DEC)". The data selection information is expected to be
 * in the format "CIRCLE(267.0208,-24.78,4.5)", where the 3 arguments are
 * Right Ascension, Declination and radius in units of degrees. No detailed
 * syntax checking is performed.
 *
 * This method clears the ROI, irrespectively of whether ROI information has
 * been found in the HDU.
 ***************************************************************************/
void GCTAEventList::read_ds_roi(const GFitsHDU* hdu)
{
    // Reset ROI
    m_roi.clear();

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get number of data selection keywords
        int ndskeys = 0;
        try {
            ndskeys = hdu->integer("NDSKEYS");
        }
        catch (GException::fits_key_not_found) {
            ;
        }

        // Loop over all data selection keys
        for (int i = 1; i <= ndskeys; ++i) {
            std::string type_key  = "DSTYP"+str(i);
            std::string unit_key  = "DSUNI"+str(i);
            std::string value_key = "DSVAL"+str(i);
            try {
                if (hdu->string(type_key) == "POS(RA,DEC)") {
                    std::string unit              = toupper(hdu->string(unit_key));
                    std::string value             = hdu->string(value_key);
                    value                         = strip_chars(value, "CIRCLE(");
                    value                         = strip_chars(value, ")");
                    std::vector<std::string> args = split(value, ",");
                    if (args.size() == 3) {
                        double  ra  = todouble(args[0]);
                        double  dec = todouble(args[1]);
                        double  rad = todouble(args[2]);
                        GCTAInstDir dir;
                        dir.radec_deg(ra, dec);
                        m_roi.centre(dir);
                        m_roi.radius(rad);
                    }
                    else {
                        throw GException::no_roi(G_READ_DS_ROI,
                              "Invalid acceptance cone value \""+value+
                              "\" encountered in data selection key \""+
                              value_key+"\"");
                    }
                } // endif: POS(RA,DEC) type found
            }
            catch (GException::fits_key_not_found) {
                ;
            }
        } // endfor: looped over data selection keys

    } // endif: HDU was valid

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
void GCTAEventList::write_events(GFitsBinTable* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Set extension name
        hdu->extname("EVENTS");

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

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write data selection keywords into FITS HDU
 *
 * @param[in] hdu FITS HDU pointer.
 *
 * This method does nothing if the HDU pointer is NULL.
 *
 * @todo This is a very dumb data selection keyword writing routine that does
 *       not take into account any existing keywords. We definitely want a
 *       more secure logic that checks for existing keywords and possible
 *       conflicts. But for this prototype software, this code does the job.
 ***************************************************************************/
void GCTAEventList::write_ds_keys(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Set ROI parameters
        double ra  = roi().centre().ra_deg();
        double dec = roi().centre().dec_deg();
        double rad = roi().radius();

        // Set energy range parameters
        double e_min = emin().TeV();
        double e_max = emax().TeV();

        // Set cone selection string
        std::string dsval2 = "CIRCLE("+str(ra)+","+str(dec)+","+str(rad)+")";

        // Set energy selection string
        std::string dsval3 = str(e_min)+":"+str(e_max);

        // Add time selection keywords
        hdu->card("DSTYP1", "TIME",  "Data selection type");
        hdu->card("DSUNI1", "s",     "Data selection unit");
        hdu->card("DSVAL1", "TABLE", "Data selection value");
        hdu->card("DSREF1", ":GTI",  "Data selection reference");

        // Add acceptance cone selection
        hdu->card("DSTYP2", "POS(RA,DEC)", "Data selection type");
        hdu->card("DSUNI2", "deg",         "Data selection unit");
        hdu->card("DSVAL2", dsval2,        "Data selection value");
        
        // Add energy range selection
        hdu->card("DSTYP3", "ENERGY", "Data selection type");
        hdu->card("DSUNI3", "TeV",    "Data selection unit");
        hdu->card("DSVAL3", dsval3,   "Data selection value");

        // Set number of data selection keys
        hdu->card("NDSKEYS", 3,  "Number of data selections");

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
