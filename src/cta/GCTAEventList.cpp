/***************************************************************************
 *                GCTAEventList.cpp  -  CTA event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
#include <iostream>
#include "GException.hpp"
#include "GCTAEventList.hpp"
#include "GFitsTableFltCol.hpp"
#include "GFitsTableDblCol.hpp"
#include "GFitsTableLngCol.hpp"
#include "GFitsTableShtCol.hpp"
#include "GFitsTableStrCol.hpp"

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
 * @brief Constructor
 *
 * Creates empty instance of GCTAEventList.
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
 * @param[in] list Event list from which the instance should be built.
 *
 * Creates instance of GCTAEventList by copying information from another
 * instance.
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
 * @param[in] list Event list to be assigned.
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
 * @brief Load events from event FITS file.
 *
 * @param[in] filename Name of FITS file from which events are loaded.
 *
 * Load events from CTA event file (extension EVENTS).
 ***************************************************************************/
void GCTAEventList::load(const std::string& filename)
{
    // Declarse FITS file
    GFits file;

    // Open FITS file
    file.open(filename);

    // Get HDU
    GFitsHDU* hdu = file.hdu("EVENTS");

    // Load columns
    load_events(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get pointer to element
 *
 * @param[in] index Event index for which pointer will be returned.
 *
 * A valid pointer is only returned if index is in the valid range. Otherwise
 * a NULL pointer is returned.
 ***************************************************************************/
GCTAEventAtom* GCTAEventList::pointer(int index)
{
    // Assign pointer if in range, NULL otherwise
    GCTAEventAtom* ptr = (m_events != NULL && index >=0 && index < m_num)
                       ? &(((GCTAEventAtom*)m_events)[index]) : NULL;

    // Return pointer
    return ptr;
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
    // Initialise base class members
    m_num    = 0;
    m_events = NULL;

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
    // Copy attributes
    m_num = list.m_num;

    // If there are events then copy them
    if (m_num > 0 && list.m_events != NULL) {

        // Allocate memory for events
        m_events = new GCTAEventAtom[m_num];

        // Copy events using the correct casts
        GCTAEventAtom* dst = (GCTAEventAtom*)m_events;
        GCTAEventAtom* src = (GCTAEventAtom*)list.m_events;
        for (int i = 0; i < m_num; ++i)
            *dst++ = *src++;

    } // endif: there were events to copy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEventList::free_members(void)
{
    // Free memory.
    if (m_events != NULL) delete [] m_events;

    // Signal free pointers
    m_events = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GCTAEventList* GCTAEventList::clone(void) const
{
    return new GCTAEventList(*this);
}


/***********************************************************************//**
 * @brief Load LAT events from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which events are loaded.
 *
 * Note that this method does not handle memory deallocation.
 ***************************************************************************/
void GCTAEventList::load_events(GFitsHDU* hdu)
{
    // Allocate space for keyword name
    char keyword[10];

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract number of events in FITS file
        m_num = hdu->card("NAXIS2")->integer();

        // If there are events then load them
        if (m_num > 0) {

            // Allocate data
            m_events = new GCTAEventAtom[m_num];

            // Get column pointers
            GFitsTableLngCol* ptr_eid         = (GFitsTableLngCol*)hdu->column("EVENT_ID");
            GFitsTableDblCol* ptr_time        = (GFitsTableDblCol*)hdu->column("TIME");
            GFitsTableShtCol* ptr_multip      = (GFitsTableShtCol*)hdu->column("MULTIP");
            GFitsTableFltCol* ptr_ra          = (GFitsTableFltCol*)hdu->column("RA");
            GFitsTableFltCol* ptr_dec         = (GFitsTableFltCol*)hdu->column("DEC");
            GFitsTableFltCol* ptr_dir_err     = (GFitsTableFltCol*)hdu->column("DIR_ERR");
            GFitsTableFltCol* ptr_detx        = (GFitsTableFltCol*)hdu->column("DETX");
            GFitsTableFltCol* ptr_dety        = (GFitsTableFltCol*)hdu->column("DETY");
            GFitsTableFltCol* ptr_alt_pnt     = (GFitsTableFltCol*)hdu->column("ALT_PNT");
            GFitsTableFltCol* ptr_az_pnt      = (GFitsTableFltCol*)hdu->column("AZ_PNT");
            GFitsTableFltCol* ptr_alt         = (GFitsTableFltCol*)hdu->column("ALT");
            GFitsTableFltCol* ptr_az          = (GFitsTableFltCol*)hdu->column("AZ");
            GFitsTableFltCol* ptr_corex       = (GFitsTableFltCol*)hdu->column("COREX");
            GFitsTableFltCol* ptr_corey       = (GFitsTableFltCol*)hdu->column("COREY");
            GFitsTableFltCol* ptr_core_err    = (GFitsTableFltCol*)hdu->column("CORE_ERR");
            GFitsTableFltCol* ptr_xmax        = (GFitsTableFltCol*)hdu->column("XMAX");
            GFitsTableFltCol* ptr_xmax_err    = (GFitsTableFltCol*)hdu->column("XMAX_ERR");
            GFitsTableFltCol* ptr_energy      = (GFitsTableFltCol*)hdu->column("ENERGY");
            GFitsTableFltCol* ptr_energy_err  = (GFitsTableFltCol*)hdu->column("ENERGY_ERR");
            GFitsTableFltCol* ptr_hil_msw     = (GFitsTableFltCol*)hdu->column("HIL_MSW");
            GFitsTableFltCol* ptr_hil_msw_err = (GFitsTableFltCol*)hdu->column("HIL_MSW_ERR");
            GFitsTableFltCol* ptr_hil_msl     = (GFitsTableFltCol*)hdu->column("HIL_MSL");
            GFitsTableFltCol* ptr_hil_msl_err = (GFitsTableFltCol*)hdu->column("HIL_MSL_ERR");

            // Copy data from columns into GCTAEventAtom objects
            GCTAEventAtom* ptr = (GCTAEventAtom*)m_events;
            for (int i = 0; i < m_num; ++i) {
                ptr[i].m_time.met((*ptr_time)(i));
                ptr[i].m_dir.radec_deg((*ptr_ra)(i), (*ptr_dec)(i));
                ptr[i].m_energy.TeV((*ptr_energy)(i));
                ptr[i].m_event_id    = (*ptr_eid)(i);
                ptr[i].m_multip      = (*ptr_multip)(i);
                ptr[i].m_dir_err     = (*ptr_dir_err)(i);
                ptr[i].m_detx        = (*ptr_detx)(i);
                ptr[i].m_dety        = (*ptr_dety)(i);
                ptr[i].m_alt_pnt     = (*ptr_alt_pnt)(i);
                ptr[i].m_az_pnt      = (*ptr_az_pnt)(i);
                ptr[i].m_alt         = (*ptr_alt)(i);
                ptr[i].m_az          = (*ptr_az)(i);
                ptr[i].m_corex       = (*ptr_corex)(i);
                ptr[i].m_corey       = (*ptr_corey)(i);
                ptr[i].m_core_err    = (*ptr_core_err)(i);
                ptr[i].m_xmax        = (*ptr_xmax)(i);
                ptr[i].m_xmax_err    = (*ptr_xmax_err)(i);
                ptr[i].m_energy_err  = (*ptr_energy_err)(i);
                ptr[i].m_hil_msw     = (*ptr_hil_msw)(i);
                ptr[i].m_hil_msw_err = (*ptr_hil_msw_err)(i);
                ptr[i].m_hil_msl     = (*ptr_hil_msl)(i);
                ptr[i].m_hil_msl_err = (*ptr_hil_msl_err)(i);
            }

        }
        else
            m_num = 0;

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put LAT event list in output stream
 *
 * @param[in] os Output stream into which the event list will be dumped
 * @param[in] list Event list to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTAEventList& list)
{
    // Put LAT event list in output stream
    os << "=== GCTAEventList ===" << std::endl;
    os << " Number of events in list ..: " << list.number();

    // Return output stream
    return os;
}
