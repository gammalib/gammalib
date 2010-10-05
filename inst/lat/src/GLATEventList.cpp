/***************************************************************************
 *                GLATEventList.cpp  -  LAT event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventList.cpp
 * @brief GLATEventList class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GLATEventList.hpp"
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
 ***************************************************************************/
GLATEventList::GLATEventList(void) : GEventList()
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
 ***************************************************************************/
GLATEventList::GLATEventList(const GLATEventList& list) : GEventList(list)
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
GLATEventList::~GLATEventList(void)
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
 * @param[in] list Event list to be assigned.
 ***************************************************************************/
GLATEventList& GLATEventList::operator= (const GLATEventList& list)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load LAT events from FT1 file.
 *
 * @param[in] ft1name FT1 FITS filename from which events are loaded.
 ***************************************************************************/
void GLATEventList::load(const std::string& ft1name)
{
    // Allocate FITS file
    GFits file;

    // Open FT1 FITS file
    file.open(ft1name);

    // Get HDU
    GFitsHDU* hdu = file.hdu("EVENTS");

    // Load columns
    load_ft1(hdu);

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
 *
 * @todo Needs assignment of pointing and response function pointers.
 * @todo Should we really return a NULL pointer in case that the list or
 *       the index is not valid? Should we not better throw an exception?
 ***************************************************************************/
GLATEventAtom* GLATEventList::pointer(int index)
{
    // Preset pointer with NULL
    GLATEventAtom* ptr = NULL;

    // Set pointer if index is in range
    if (m_events != NULL && index >=0 && index < m_num) {

        // Point to the requested event atom
        ptr = &(((GLATEventAtom*)m_events)[index]);

        // Set instrument pointing
        ptr->m_pnt = NULL; // DUMMY

        // Set instrument response function
        ptr->m_rsp = NULL; // DUMMY

    } // endif: valid index

    // Return pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATEventList::init_members(void)
{
    // Initialise members
    m_num    = 0;
    m_events = NULL;

    // Initialise diffuse models
    m_num_difrsp   = 0;
    m_difrsp_label = NULL;

    // Initialise data selection keys
    m_ds_num       = 0;
    m_ds_type      = NULL;
    m_ds_unit      = NULL;
    m_ds_value     = NULL;
    m_ds_reference = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] list GLATEventList members which should be copied.
 ***************************************************************************/
void GLATEventList::copy_members(const GLATEventList& list)
{
    // Copy attributes
    m_num        = list.m_num;
    m_num_difrsp = list.m_num_difrsp;
    m_ds_num     = list.m_ds_num;

    // If there are events then copy them
    if (m_num > 0 && list.m_events != NULL) {

        // Allocate memory for events
        m_events = new GLATEventAtom[m_num];

        // Copy events using the correct casts
        GLATEventAtom* dst = (GLATEventAtom*)m_events;
        GLATEventAtom* src = (GLATEventAtom*)list.m_events;
        for (int i = 0; i < m_num; ++i)
            *dst++ = *src++;

    }

    // If there are diffuse response model labels then copy them
    if (m_num_difrsp > 0) {

        // Allocate memory for keys
        m_difrsp_label = new std::string[m_num_difrsp];

        // Copy keys
        for (int i = 0; i < m_num_difrsp; ++i)
            m_difrsp_label[i] = list.m_difrsp_label[i];

    } // endif: there were diffuse response model labels

    // If there are data selection keys then copy them
    if (m_ds_num > 0) {

        // Allocate memory for keys
        m_ds_type      = new std::string[m_ds_num];
        m_ds_unit      = new std::string[m_ds_num];
        m_ds_value     = new std::string[m_ds_num];
        m_ds_reference = new std::string[m_ds_num];

        // Copy keys
        for (int i = 0; i < m_ds_num; ++i) {
            m_ds_type[i]      = list.m_ds_type[i];
            m_ds_unit[i]      = list.m_ds_unit[i];
            m_ds_value[i]     = list.m_ds_value[i];
            m_ds_reference[i] = list.m_ds_reference[i];
        }

    } // endif: there were data selection keys to copy

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventList::free_members(void)
{
    // Free memory.
    if (m_events       != NULL) delete [] m_events;
    if (m_difrsp_label != NULL) delete [] m_difrsp_label;
    if (m_ds_type      != NULL) delete [] m_ds_type;
    if (m_ds_unit      != NULL) delete [] m_ds_unit;
    if (m_ds_value     != NULL) delete [] m_ds_value;
    if (m_ds_reference != NULL) delete [] m_ds_reference;

    // Signal free pointers
    m_events       = NULL;
    m_difrsp_label = NULL;
    m_ds_type      = NULL;
    m_ds_unit      = NULL;
    m_ds_value     = NULL;
    m_ds_reference = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GLATEventList* GLATEventList::clone(void) const
{
    return new GLATEventList(*this);
}


/***********************************************************************//**
 * @brief Load LAT events from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which events are loaded.
 ***************************************************************************/
void GLATEventList::load_ft1(GFitsHDU* hdu)
{
    // Free and initialise base class members
    this->GEvents::free_members();
    this->GEvents::init_members();

    // Free and initialise base class members
    this->GEventList::free_members();
    this->GEventList::init_members();

    // Free and initialise class members
    free_members();
    init_members();

    // Load event data
    load_events(hdu);

    // Load data selection keywords
    load_ds_keys(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load LAT events from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which events are loaded.
 *
 * Note that this method does not handle memory deallocation.
 ***************************************************************************/
void GLATEventList::load_events(GFitsHDU* hdu)
{
    // Allocate space for keyword name
    char keyword[10];

    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Extract number of events in FT1 file
        m_num = hdu->card("NAXIS2")->integer();

        // If there are events then load them
        if (m_num > 0) {

            // Allocate data
            m_events = new GLATEventAtom[m_num];

            // Get column pointers
            GFitsTableDblCol* ptr_time    = (GFitsTableDblCol*)hdu->column("TIME");
            GFitsTableFltCol* ptr_energy  = (GFitsTableFltCol*)hdu->column("ENERGY");
            GFitsTableFltCol* ptr_ra      = (GFitsTableFltCol*)hdu->column("RA");
            GFitsTableFltCol* ptr_dec     = (GFitsTableFltCol*)hdu->column("DEC");
            GFitsTableFltCol* ptr_theta   = (GFitsTableFltCol*)hdu->column("THETA");
            GFitsTableFltCol* ptr_phi     = (GFitsTableFltCol*)hdu->column("PHI");
            GFitsTableFltCol* ptr_zenith  = (GFitsTableFltCol*)hdu->column("ZENITH_ANGLE");
            GFitsTableFltCol* ptr_azimuth = (GFitsTableFltCol*)hdu->column("EARTH_AZIMUTH_ANGLE");
            GFitsTableLngCol* ptr_eid     = (GFitsTableLngCol*)hdu->column("EVENT_ID");
            GFitsTableLngCol* ptr_rid     = (GFitsTableLngCol*)hdu->column("RUN_ID");
            GFitsTableShtCol* ptr_recon   = (GFitsTableShtCol*)hdu->column("RECON_VERSION");
            GFitsTableShtCol* ptr_calib   = (GFitsTableShtCol*)hdu->column("CALIB_VERSION");
            GFitsTableShtCol* ptr_class   = (GFitsTableShtCol*)hdu->column("EVENT_CLASS");
            GFitsTableShtCol* ptr_conv    = (GFitsTableShtCol*)hdu->column("CONVERSION_TYPE");
            GFitsTableDblCol* ptr_ltime   = (GFitsTableDblCol*)hdu->column("LIVETIME");

            // Copy data from columns into GLATEventAtom objects
            GLATEventAtom* ptr = (GLATEventAtom*)m_events;
            for (int i = 0; i < m_num; ++i) {
                ptr[i].m_time.met((*ptr_time)(i));
                ptr[i].m_energy.MeV((*ptr_energy)(i));
                ptr[i].m_dir.radec_deg((*ptr_ra)(i), (*ptr_dec)(i));
                ptr[i].m_theta               = (*ptr_theta)(i);
                ptr[i].m_phi                 = (*ptr_phi)(i);
                ptr[i].m_zenith_angle        = (*ptr_zenith)(i);
                ptr[i].m_earth_azimuth_angle = (*ptr_azimuth)(i);
                ptr[i].m_event_id            = (*ptr_eid)(i);
                ptr[i].m_run_id              = (*ptr_rid)(i);
                ptr[i].m_recon_version       = (*ptr_recon)(i);
                ptr[i].m_calib_version[0]    = (*ptr_calib)(i,0);
                ptr[i].m_calib_version[1]    = (*ptr_calib)(i,1);
                ptr[i].m_calib_version[2]    = (*ptr_calib)(i,2);
                ptr[i].m_event_class         = (*ptr_class)(i);
                ptr[i].m_conversion_type     = (*ptr_conv)(i);
                ptr[i].m_livetime            = (*ptr_ltime)(i);
            }

            // Extract number of diffuse response labels
            m_num_difrsp = hdu->card("NDIFRSP")->integer();

            // Allocate diffuse response components
            if (m_num_difrsp > 0) {

                // Allocate labels
                m_difrsp_label = new std::string[m_num_difrsp];

                // Allocate components
                GLATEventAtom* ptr = (GLATEventAtom*)m_events;
                for (int i = 0; i < m_num; ++i)
                    ptr[i].m_difrsp = new double[m_num_difrsp];

                // Load diffuse columns
                for (int k = 0; k < m_num_difrsp; ++k) {

                    // Set keyword name
                    sprintf(keyword, "DIFRSP%d", k);

                    // Get DIFRSP label
                    try {
                        m_difrsp_label[k] = hdu->card(std::string(keyword))->string();
                    }
                    catch (GException::fits_key_not_found &e) {
                        m_difrsp_label[k] = "NONE";
                    }

                    // Get column pointer
                    GFitsTableFltCol* ptr_dif = (GFitsTableFltCol*)hdu->column(std::string(keyword));

                    // Copy data from columns into GLATEventAtom objects
                    GLATEventAtom* ptr = (GLATEventAtom*)m_events;
                    for (int i = 0; i < m_num; ++i)
                        ptr[i].m_difrsp[k] = (*ptr_dif)(i);

                } // endfor: looped over diffuse columns

            } // endif: allocated memory for diffuse components
            else
                m_num_difrsp = 0;

        }
        else
            m_num = 0;

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Data Selection keywords from FITS HDU.
 *
 * @param[in] hdu Pointer to FITS HDU from which DS keywords are loaded.
 *
 * Note that this method does not handle memory deallocation.
 ***************************************************************************/
void GLATEventList::load_ds_keys(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get number of data selection keys
        m_ds_num = hdu->card("NDSKEYS")->integer();

        // Get data selection keys
        if (m_ds_num > 0) {

            // Allocate data selection key information
            m_ds_type      = new std::string[m_ds_num];
            m_ds_unit      = new std::string[m_ds_num];
            m_ds_value     = new std::string[m_ds_num];
            m_ds_reference = new std::string[m_ds_num];

            // Allocate space for the keyword
            char keyword[10];

            // Get columns
            for (int i = 0; i < m_ds_num; ++i) {

                // Get DSTYPnn
                try {
                    sprintf(keyword, "DSTYP%d", i+1);
                    m_ds_type[i] = hdu->card(std::string(keyword))->string();
                }
                catch (GException::fits_key_not_found &e) {
                    m_ds_type[i] = "";
                }

                // Get DSUNInn
                try {
                    sprintf(keyword, "DSUNI%d", i+1);
                    m_ds_unit[i] = hdu->card(std::string(keyword))->string();
                }
                catch (GException::fits_key_not_found &e) {
                    m_ds_unit[i] = "";
                }

                // Get DSVALnn
                try {
                    sprintf(keyword, "DSVAL%d", i+1);
                    m_ds_value[i] = hdu->card(std::string(keyword))->string();
                }
                catch (GException::fits_key_not_found &e) {
                    m_ds_value[i] = "";
                }

                // Get DSREFnn
                try {
                    sprintf(keyword, "DSREF%d", i+1);
                    m_ds_reference[i] = hdu->card(std::string(keyword))->string();
                }
                catch (GException::fits_key_not_found &e) {
                    m_ds_reference[i] = "";
                }

            } // endfor: looped over data selection keys

        } // endif: there were data selection keys
        else
            m_ds_num = 0;

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GLATEventList friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put LAT event list in output stream
 *
 * @param[in] os Output stream into which the event list will be dumped
 * @param[in] list Event list to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATEventList& list)
{
    // Put LAT event list in output stream
    os << "=== GLATEventList ===" << std::endl;
    os << " Number of DS keywords .....: " << list.m_ds_num << std::endl;
    for (int i = 0; i < list.m_ds_num; ++i) {
        os << "  Data selection ...........: Type=" << list.m_ds_type[i];
        os << " Value=" << list.m_ds_value[i];
        if (list.m_ds_unit[i] != "")
            os << " Unit=" << list.m_ds_unit[i];
        if (list.m_ds_reference[i] != "")
            os << " Reference=" << list.m_ds_reference[i];
        os << std::endl;
    }
    os << " Number of diffuse models ..: " << list.m_num_difrsp << std::endl;
    for (int i = 0; i < list.m_num_difrsp; ++i)
        os << "  Diffuse response component: " << list.m_difrsp_label[i] << std::endl;
    os << " Number of events in list ..: " << list.number();

    // Return output stream
    return os;
}
