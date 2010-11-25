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
#include "GCTAException.hpp"
#include "GCTAEventList.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponse.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
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
    load_events(hdu);

    // Close FITS file
    file.close();

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
 * This method returns a pointer on an event atom.
 ***************************************************************************/
GCTAEventAtom* GCTAEventList::pointer(int index)
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num)
        throw GException::out_of_range(G_POINTER, index, 0, m_num-1);
    #endif

    // Return pointer
    return ((GCTAEventAtom*)m_events + index);
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
    //m_obs    = NULL;

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
 * @brief Load CTA events from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS table.
 *
 * This method loads the CTA event list from a FITS file into memory. Only
 * quantities that are relevant for science analysis will be loaded. This
 * method can therefore not be used to modify the contents of a CTA events
 * file. For this purpose specific code should be written that accesses the
 * data directly through the GFits class.
 ***************************************************************************/
void GCTAEventList::load_events(GFitsTable* hdu)
{
    // Allocate space for keyword name
    char keyword[10];

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract number of events in FITS file
        m_num = hdu->integer("NAXIS2");

        // If there are events then load them
        if (m_num > 0) {

            // Allocate data
            m_events = new GCTAEventAtom[m_num];

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
            GFitsTableDoubleCol* ptr_energy      = (GFitsTableDoubleCol*)hdu->column("ENERGY");
            GFitsTableFloatCol*  ptr_energy_err  = (GFitsTableFloatCol*)hdu->column("ENERGY_ERR");
            GFitsTableFloatCol*  ptr_hil_msw     = (GFitsTableFloatCol*)hdu->column("HIL_MSW");
            GFitsTableFloatCol*  ptr_hil_msw_err = (GFitsTableFloatCol*)hdu->column("HIL_MSW_ERR");
            GFitsTableFloatCol*  ptr_hil_msl     = (GFitsTableFloatCol*)hdu->column("HIL_MSL");
            GFitsTableFloatCol*  ptr_hil_msl_err = (GFitsTableFloatCol*)hdu->column("HIL_MSL_ERR");

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
