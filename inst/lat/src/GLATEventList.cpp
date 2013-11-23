/***************************************************************************
 *               GLATEventList.cpp - Fermi/LAT event list class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GLATEventList.cpp
 * @brief Fermi/LAT event list class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::sprintf
#include "GException.hpp"
#include "GTools.hpp"
#include "GLATEventList.hpp"
#include "GLATException.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableShortCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR                          "GLATEventList::operator[](int&)"
#define G_ROI                                     "GLATEventList::roi(GRoi&)"

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
 * @param[in] list LAT event list.
 ***************************************************************************/
GLATEventList::GLATEventList(const GLATEventList& list) : GEventList(list)
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
 * @param[in] list LAT event list.
 ***************************************************************************/
GLATEventList& GLATEventList::operator= (const GLATEventList& list)
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
GLATEventAtom* GLATEventList::operator[](const int& index)
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
const GLATEventAtom* GLATEventList::operator[](const int& index) const
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATEventList::clear(void)
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
GLATEventList* GLATEventList::clone(void) const
{
    return new GLATEventList(*this);
}


/***********************************************************************//**
 * @brief Load LAT events from FITS file
 *
 * @param[in] filename FITS filename.
 *
 * This method loads LAT events from a FT1 file.
 ***************************************************************************/
void GLATEventList::load(const std::string& filename)
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
 * @brief Save LAT events
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file (default=false).
 *
 * This method saves LAT events into a FT1 file.
 *
 * @todo To be implemented.
 ***************************************************************************/
void GLATEventList::save(const std::string& filename, bool clobber) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read LAT events from FITS file.
 *
 * @param[in] file FITS file.
 *
 * This method read the LAT event list from a FITS file.
 *
 * The method clears the object before loading, thus any events residing in
 * the object before loading will be lost.
 ***************************************************************************/
void GLATEventList::read(const GFits& file)
{
    // Clear object
    clear();

    // Get HDU
    GFitsTable* hdu = file.table("EVENTS");

    // Continue only if event HDU is valid
    if (hdu != NULL) {

        // Read event data
        read_events(*hdu);

        // Read data selection keywords
        read_ds_keys(*hdu);
    
    } // endif: HDU was valid

    // If we have a GTI extension, then read Good Time Intervals from that
    // extension
    if (file.hashdu("GTI")) {
        GFitsTable* gti = file.table("GTI");
        m_gti.read(gti);
    }

    // ... otherwise build GTI from TSTART and TSTOP
    else if (hdu != NULL) {

        // Read start and stop time
        double tstart = hdu->real("TSTART");
        double tstop  = hdu->real("TSTOP");

        // Create time reference from header information
        GTimeReference timeref(hdu);
        
        // Set start and stop time
        GTime start(tstart);
        GTime stop(tstop);

        // Append start and stop time as single time interval to GTI
        m_gti.append(start, stop);

        // Set GTI time reference
        m_gti.reference(timeref);

    } // endelse: GTI built from TSTART and TSTOP

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write LAT event list into FITS file.
 *
 * @param[in] file FITS file.
 *
 * Write the LAT event list into FITS file.
 *
 * @todo To be implemented.
 ***************************************************************************/
void GLATEventList::write(GFits& file) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set region of interest
 *
 * @param[in] roi Region of interest.
 *
 * @exception GLATException::bad_roi_type
 *            ROI is not of type GLATRoi.
 ***************************************************************************/
void GLATEventList::roi(const GRoi& roi)
{
    // Cast ROI dynamically
    const GLATRoi* ptr = dynamic_cast<const GLATRoi*>(&roi);

    // Throw exception if ROI is not of correct type
    if (ptr == NULL)
        throw GLATException::bad_roi_type(G_ROI);

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
std::string GLATEventList::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATEventList ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of events") +
                      gammalib::str(number()));

        // Append DS keywords
        result.append("\n"+gammalib::parformat("Number of DS keywords"));
        result.append(gammalib::str(m_ds_type.size()));
        for (int i = 0; i < m_ds_type.size(); ++i) {
            result.append("\n"+gammalib::parformat(" Data selection"));
            result.append("Type="+m_ds_type[i]);
            if (m_ds_unit[i].length() > 0) {
                result.append(" Unit="+m_ds_unit[i]);
            }
            if (m_ds_reference[i].length() > 0) {
                result.append(" Reference="+m_ds_reference[i]);
            }
        }

        // Append diffuse keywords
        result.append("\n"+gammalib::parformat("Number of diffuse models"));
        result.append(gammalib::str(m_difrsp_label.size()));
        for (int i = 0; i < m_difrsp_label.size(); ++i) {
            result.append("\n"+gammalib::parformat(" Diffuse component"));
            result.append(m_difrsp_label[i]);
        }

    } // endif: chatter was not silent

    // Return result
    return result;
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
    m_roi.clear();
    m_events.clear();
    m_difrsp_label.clear();
    m_ds_type.clear();
    m_ds_unit.clear();
    m_ds_value.clear();
    m_ds_reference.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] list LAT event list.
 ***************************************************************************/
void GLATEventList::copy_members(const GLATEventList& list)
{
    // Copy members
    m_roi          = list.m_roi;
    m_events       = list.m_events;
    m_difrsp_label = list.m_difrsp_label;
    m_ds_type      = list.m_ds_type;
    m_ds_unit      = list.m_ds_unit;
    m_ds_value     = list.m_ds_value;
    m_ds_reference = list.m_ds_reference;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventList::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read LAT events from FITS table.
 *
 * @param[in] table Event table.
 *
 * Read the LAT events from the event table.
 ***************************************************************************/
void GLATEventList::read_events(const GFitsTable& table)
{
    // Clear existing events
    m_events.clear();

    // Allocate space for keyword name
    char keyword[10];

    // Extract number of events in FT1 file
    int num = table.integer("NAXIS2");

    // If there are events then load them
    if (num > 0) {

        // Reserve data
        m_events.reserve(num);

        // Get column pointers
        GFitsTableDoubleCol* ptr_time    = (GFitsTableDoubleCol*)table["TIME"];
        GFitsTableFloatCol*  ptr_energy  = (GFitsTableFloatCol*)table["ENERGY"];
        GFitsTableFloatCol*  ptr_ra      = (GFitsTableFloatCol*)table["RA"];
        GFitsTableFloatCol*  ptr_dec     = (GFitsTableFloatCol*)table["DEC"];
        GFitsTableFloatCol*  ptr_theta   = (GFitsTableFloatCol*)table["THETA"];
        GFitsTableFloatCol*  ptr_phi     = (GFitsTableFloatCol*)table["PHI"];
        GFitsTableFloatCol*  ptr_zenith  = (GFitsTableFloatCol*)table["ZENITH_ANGLE"];
        GFitsTableFloatCol*  ptr_azimuth = (GFitsTableFloatCol*)table["EARTH_AZIMUTH_ANGLE"];
        GFitsTableLongCol*   ptr_eid     = (GFitsTableLongCol*)table["EVENT_ID"];
        GFitsTableLongCol*   ptr_rid     = (GFitsTableLongCol*)table["RUN_ID"];
        GFitsTableShortCol*  ptr_recon   = (GFitsTableShortCol*)table["RECON_VERSION"];
        GFitsTableShortCol*  ptr_calib   = (GFitsTableShortCol*)table["CALIB_VERSION"];
        GFitsTableShortCol*  ptr_class   = (GFitsTableShortCol*)table["EVENT_CLASS"];
        GFitsTableShortCol*  ptr_conv    = (GFitsTableShortCol*)table["CONVERSION_TYPE"];
        GFitsTableDoubleCol* ptr_ltime   = (GFitsTableDoubleCol*)table["LIVETIME"];

        // Copy data from columns into GLATEventAtom objects
        GLATEventAtom event;
        for (int i = 0; i < num; ++i) {
            event.m_time.set((*ptr_time)(i), m_gti.reference());
            event.m_energy.MeV((*ptr_energy)(i));
            event.m_dir.dir().radec_deg((*ptr_ra)(i), (*ptr_dec)(i));
            event.m_theta               = (*ptr_theta)(i);
            event.m_phi                 = (*ptr_phi)(i);
            event.m_zenith_angle        = (*ptr_zenith)(i);
            event.m_earth_azimuth_angle = (*ptr_azimuth)(i);
            event.m_event_id            = (*ptr_eid)(i);
            event.m_run_id              = (*ptr_rid)(i);
            event.m_recon_version       = (*ptr_recon)(i);
            event.m_calib_version[0]    = (*ptr_calib)(i,0);
            event.m_calib_version[1]    = (*ptr_calib)(i,1);
            event.m_calib_version[2]    = (*ptr_calib)(i,2);
            event.m_event_class         = (*ptr_class)(i);
            event.m_conversion_type     = (*ptr_conv)(i);
            event.m_livetime            = (*ptr_ltime)(i);
            m_events.push_back(event);
        }

        // Extract number of diffuse response labels
        int num_difrsp = table.integer("NDIFRSP");

        // Allocate diffuse response components
        if (num_difrsp > 0) {

            // Reserve space
            m_difrsp_label.reserve(num_difrsp);

            // Allocate components
            for (int i = 0; i < num; ++i) {
                m_events[i].m_difrsp = new double[num_difrsp];
            }

            // Load diffuse columns
            for (int k = 0; k < num_difrsp; ++k) {

                // Set keyword name
                std::sprintf(keyword, "DIFRSP%d", k);

                // Get DIFRSP label
                try {
                    GFitsTable* ptr = const_cast<GFitsTable*>(&table); // Kluge to const
                    m_difrsp_label.push_back(ptr->card(std::string(keyword))->string());
                }
                catch (GException::fits_key_not_found &e) {
                    m_difrsp_label.push_back("NONE");
                }

                // Get column pointer
                GFitsTableFloatCol* ptr_dif = 
                    static_cast<GFitsTableFloatCol*>(const_cast<GFitsTableCol*>(table[std::string(keyword)]));

                // Copy data from columns into GLATEventAtom objects
                for (int i = 0; i < num; ++i) {
                    m_events[i].m_difrsp[k] = (*ptr_dif)(i);
                }

            } // endfor: looped over diffuse columns

        } // endif: diffuse components found

    } // endif: events found

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data selection keywords from FITS HDU.
 *
 * @param[in] hdu FITS HDU pointer.
 *
 * @todo Declared header card const in to GFitsHDU.
 * @todo Add check key method to GFitsHDU to avoid unneccesary try/catch
 *       blocks.
 ***************************************************************************/
void GLATEventList::read_ds_keys(const GFitsHDU& hdu)
{
    // Get number of data selection keys
    int ds_num = hdu.integer("NDSKEYS");

    // Get data selection keys
    if (ds_num > 0) {

        // Circumvent const correctness. We need this because the header()
        // card method is not declared const. This should be corrected.
        GFitsHDU* ptr = (GFitsHDU*)&hdu;

        // Reserve space
        m_ds_type.reserve(ds_num);
        m_ds_unit.reserve(ds_num);
        m_ds_value.reserve(ds_num);
        m_ds_reference.reserve(ds_num);

        // Allocate space for the keyword
        char keyword[10];

        // Get columns
        for (int i = 0; i < ds_num; ++i) {

            // Get DSTYPnn
            try {
                std::sprintf(keyword, "DSTYP%d", i+1);
                m_ds_type.push_back(ptr->card(std::string(keyword))->string());
            }
            catch (GException::fits_key_not_found &e) {
                m_ds_type.push_back("");
            }

            // Get DSUNInn
            try {
                std::sprintf(keyword, "DSUNI%d", i+1);
                m_ds_unit.push_back(ptr->card(std::string(keyword))->string());
            }
            catch (GException::fits_key_not_found &e) {
                m_ds_unit.push_back("");
            }

            // Get DSVALnn
            try {
                std::sprintf(keyword, "DSVAL%d", i+1);
                m_ds_value.push_back(ptr->card(std::string(keyword))->string());
            }
            catch (GException::fits_key_not_found &e) {
                m_ds_value.push_back("");
            }

            // Get DSREFnn
            try {
                std::sprintf(keyword, "DSREF%d", i+1);
                m_ds_reference.push_back(ptr->card(std::string(keyword))->string());
            }
            catch (GException::fits_key_not_found &e) {
                m_ds_reference.push_back("");
            }

        } // endfor: looped over data selection keys

    } // endif: there were data selection keys

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Friends                                   =
 =                                                                         =
 ==========================================================================*/
