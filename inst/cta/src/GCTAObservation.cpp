/***************************************************************************
 *               GCTAObservation.cpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GCTAObservation.cpp
 * @brief CTA observation class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAObservation      g_obs_cta_seed("CTA");
const GCTAObservation      g_obs_hess_seed("HESS");
const GCTAObservation      g_obs_magic_seed("MAGIC");
const GCTAObservation      g_obs_veritas_seed("VERITAS");
const GObservationRegistry g_obs_cta_registry(&g_obs_cta_seed);
const GObservationRegistry g_obs_hess_registry(&g_obs_hess_seed);
const GObservationRegistry g_obs_magic_registry(&g_obs_magic_seed);
const GObservationRegistry g_obs_veritas_registry(&g_obs_veritas_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GCTAObservation::response(GResponse&)"
#define G_READ                          "GCTAObservation::read(GXmlElement&)"
#define G_WRITE                        "GCTAObservation::write(GXmlElement&)"
#define G_READ_DS_EBOUNDS       "GCTAObservation::read_ds_ebounds(GFitsHDU*)"
#define G_READ_DS_ROI               "GCTAObservation::read_ds_roi(GFitsHDU*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates empty class instance.
 ***************************************************************************/
GCTAObservation::GCTAObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument constructor
 *
 * @param[in] instrument Instrument name.
 *
 * Creates empty CTA observation class for a given instrument. This enables
 * using the CTA specific interface for any other THE instrument. Note that
 * each other THE instruments needs a specific registry at the beginning
 * of the GCTAObservation.cpp file.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const std::string& instrument) : GObservation()
{
    // Initialise members
    init_members();

    // Set instrument
    m_instrument = instrument;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs CTA observation.
 *
 * Creates class instance by copying an existing CTA observation.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const GCTAObservation& obs) : GObservation(obs)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAObservation::~GCTAObservation(void)
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
 * @param[in] obs CTA observation.
 *
 * Assign CTA observation to this object.
 ***************************************************************************/
GCTAObservation& GCTAObservation::operator= (const GCTAObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAObservation::clear(void)
{
    // Free members
    free_members();
    this->GObservation::free_members();

    // Initialise members
    this->GObservation::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAObservation* GCTAObservation::clone(void) const
{
    return new GCTAObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GCTAException::bad_response_type
 *            Specified response in not of type GCTAResponse.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GCTAResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GCTAObservation::response(const GResponse& rsp)
{
    // Get pointer on CTA response
    const GCTAResponse* ctarsp = dynamic_cast<const GCTAResponse*>(&rsp);
    if (ctarsp == NULL)
        throw GCTAException::bad_response_type(G_RESPONSE);

    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = ctarsp->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA response function
 *
 * @param[in] irfname Name of CTA response function.
 * @param[in] caldb Optional name of calibration database.
 ***************************************************************************/
void GCTAObservation::response(const std::string& irfname, std::string caldb)
{
    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Allocate new CTA response function
    m_response = new GCTAResponse;

    // Set calibration database
    m_response->caldb(caldb);

    // Load instrument response function
    m_response->load(irfname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to CTA response function
 ***************************************************************************/
GCTAResponse* GCTAObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to CTA pointing direction
 *
 * Returns pointer to pointing direction.
 ***************************************************************************/
GCTAPointing* GCTAObservation::pointing(void) const
{
    // Return pointing pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Set CTA pointing direction
 *
 * @param[in] pointing Pointing.
 ***************************************************************************/
void GCTAObservation::pointing(const GCTAPointing& pointing)
{
    // Free any existing pointing
    if (m_pointing != NULL) delete m_pointing;
    m_pointing = NULL;

    // Clone pointing
    m_pointing = pointing.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::xml_invalid_parnum
 *            Invalid number of parameters found in XML element.
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Reads information for a CTA observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="EventList" file="..."/>
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *     </observation>
 *
 * for an unbinned observation and
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="CountsMap" file="..."/>
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *     </observation>
 *
 * for a binned observation.
 *
 * @todo Still supports old ARF, PSF and RMF parameter names.
 ***************************************************************************/
void GCTAObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Initialise response parameters
    std::string aeff;
    std::string psf;
    std::string rmf;
    double      thetacut = 0.0;
    double      scale    = 1.0;

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || npars != 4) {
        throw GException::xml_invalid_parnum(G_READ, xml,
              "CTA observation requires exactly 4 parameters.");
    }

    // Extract parameters
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

        // Handle EventList
        if (par->attribute("name") == "EventList") {

            // Read eventlist file name
            std::string filename = par->attribute("file");
            
            // Load unbinned observation (sets also m_evenfile member)
            load_unbinned(filename);
            
            // Increment parameter counter
            npar[0]++;
        }

        // Handle CountsMap
        else if (par->attribute("name") == "CountsMap") {

            // Read countsmap file name
            std::string filename = par->attribute("file");
            
            // Load binned observation (sets also m_evenfile member)
            load_binned(filename);
            
            // Increment parameter counter
            npar[0]++;
        }

        // Handle effective area
        else if ((par->attribute("name") == "EffectiveArea") ||
                 (par->attribute("name") == "ARF")) {

            // Store filename
            aeff = par->attribute("file");

            // Optionally extract thetacut (0.0 if no thetacut)
            std::string s_thetacut = par->attribute("thetacut");
            if (s_thetacut.length() > 0) {
                thetacut = todouble(s_thetacut);
            }

            // Optionally extract scale factor (1.0 if no scale)
            std::string s_scale = par->attribute("scale");
            if (s_scale.length() > 0) {
                scale = todouble(s_scale);
            }

            // Increase number of parameters
            npar[1]++;
        }

        // Handle PSF
        else if ((par->attribute("name") == "PointSpreadFunction") ||
                 (par->attribute("name") == "PSF")) {
            psf = par->attribute("file");
            npar[2]++;
        }


        // Handle RMF
        else if ((par->attribute("name") == "EnergyDispersion") ||
                 (par->attribute("name") == "RMF")) {
            rmf = par->attribute("file");
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require \"EventList\" or \"CountsMap\" and \"EffectiveArea\""
              ", \"PointSpreadFunction\" and \"EnergyDispersion\" parameters.");
    }

    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Allocate new CTA response function
    m_response = new GCTAResponse;

    // Load PSF (need this before Aeff to handle thetacut)
    if (strip_whitespace(psf).length() > 0) {
        m_response->load_psf(psf);
    }

    // Load effective area
    if (strip_whitespace(aeff).length() > 0) {
        m_response->arf_thetacut(thetacut);
        m_response->arf_scale(scale);
        m_response->load_arf(aeff);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::no_list
 *            No valid CTA event list or event cube found.
 * @exception GException::xml_invalid_parnum
 *            Invalid number of parameters found in XML element.
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Writes information for a CTA observation into an XML element. The expected
 * format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="EventList" file="..."/>
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *     </observation>
 *
 * for an unbinned observation and
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="CountsMap" file="..."/>
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *     </observation>
 *
 * for a binned observation.
 *
 * @todo We should create a special exception that informs that there is
 *       neither a valid CTA event list nor a valid CTA counts map in this
 *       observations.
 ***************************************************************************/
void GCTAObservation::write(GXmlElement& xml) const
{
    // Determine if we deal with a binned or unbinned observation
    const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(m_events);
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(m_events);
    if (list == NULL && cube == NULL) {
        throw GException::no_list(G_WRITE);
    }

    // Set event list flag
    bool is_list = (list != NULL);

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        if (is_list) {
            xml.append(new GXmlElement("parameter name=\"EventList\""));
        }
        else {
            xml.append(new GXmlElement("parameter name=\"CountsMap\""));
        }
        xml.append(new GXmlElement("parameter name=\"EffectiveArea\""));
        xml.append(new GXmlElement("parameter name=\"PointSpreadFunction\""));
        xml.append(new GXmlElement("parameter name=\"EnergyDispersion\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::xml_invalid_parnum(G_WRITE, xml,
              "CTA observation requires exactly 4 parameters.");
    }

    // Set or update parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

        // Handle Eventlist
        if (par->attribute("name") == "EventList") {
            par->attribute("file", m_eventfile);
            npar[0]++;
        }

        // Handle Countsmap
        else if (par->attribute("name") == "CountsMap") {
            par->attribute("file", m_eventfile);
            npar[0]++;
        }

        // Handle ARF
        else if (par->attribute("name") == "EffectiveArea") {
            std::string aeff     = "";
            double      thetacut = 0.0;
            double      scale    = 1.0;
            if (m_response != NULL) {
                aeff     = m_response->arffile();
                thetacut = m_response->arf_thetacut();
                scale    = m_response->arf_scale();
            }
            par->attribute("file", aeff);
            if (thetacut > 0.0) {
                par->attribute("thetacut", str(thetacut));
            }
            if (scale != 1.0) {
                par->attribute("scale", str(scale));
            }
            npar[1]++;
        }

        // Handle PSF
        else if (par->attribute("name") == "PointSpreadFunction") {
            std::string psf = "";
            if (m_response != NULL) {
                psf = m_response->psffile();
            }
            par->attribute("file", psf);
            npar[2]++;
        }

        // Handle RMF
        else if (par->attribute("name") == "EnergyDispersion") {
            std::string rmf = "";
            if (m_response != NULL) {
                rmf = m_response->rmffile();
            }
            par->attribute("file", rmf);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::xml_invalid_parnames(G_WRITE, xml,
              "Require \"EventList\" or \"CountsMap\" and \"ARF\", \"RMF\""
              " \"PSF\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA observation information
 ***************************************************************************/
std::string GCTAObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAObservation ===");
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Identifier")+id());
    result.append("\n"+parformat("Instrument")+instrument());
    result.append("\n"+parformat("Statistics")+statistics());
    result.append("\n"+parformat("Ontime")+str(ontime()));
    result.append("\n"+parformat("Livetime")+str(livetime()));
    result.append("\n"+parformat("Deadtime correction")+str(m_deadc));

    // Append pointing
    if (m_pointing != NULL) {
        result.append("\n"+m_pointing->print());
    }
    else {
        result.append("\n"+parformat("CTA pointing")+"undefined");
    }

    // Append response
    if (m_response != NULL) {
        result.append("\n"+response()->print());
    }
    else {
        result.append("\n"+parformat("CTA response")+"undefined");
    }

    // Append events
    if (m_events != NULL) {
        result.append("\n"+m_events->print());
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] filename Event FITS file name.
 ***************************************************************************/
void GCTAObservation::load_unbinned(const std::string& filename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event list
    GCTAEventList* events = new GCTAEventList;

    // Assign event list as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(filename);

    // Read event list
    events->read(file);

    // Read observation attributes from EVENTS extension
    GFitsHDU* hdu = file.hdu("EVENTS");
    read_attributes(hdu);

    // Close FITS file
    file.close();

    // Store event filename
    m_eventfile = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] filename Counts map FITS file name.
 ***************************************************************************/
void GCTAObservation::load_binned(const std::string& filename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event cube
    GCTAEventCube* events = new GCTAEventCube;

    // Assign event cube as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(filename);

    // Read event cube
    events->read(file);

    // Read observation attributes from primary extension
    GFitsHDU* hdu = file.hdu(0);
    read_attributes(hdu);

    // Close FITS file
    file.close();

    // Store event filename
    m_eventfile = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save CTA observation into FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file (default=false).
 ***************************************************************************/
void GCTAObservation::save(const std::string& filename, bool clobber) const
{
    // Create FITS file
    GFits fits;

    // Get pointers on event list
    GCTAEventList* list = dynamic_cast<GCTAEventList*>(m_events);
    GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(m_events);

    // Case A: Observation contains an event list
    if (list != NULL) {

        // Write event list into FITS file. This method also writes
        // the GTI as they are part of the event list.
        list->write(fits);

        // Write observation attributes into EVENTS header
        GFitsHDU* hdu = fits.hdu("EVENTS");
        write_attributes(hdu);

    } // endif: observation contained an event list

    // Case B: Observation contains an event cube
    else if (cube != NULL) {

        // Write events cube into FITS file. This method also writes
        // the energy boundaries and the GTI as they are also part
        // of the event cube.
        cube->write(fits);

        // Write observation attributes into primary header
        GFitsHDU* hdu = fits.hdu(0);
        write_attributes(hdu);

    } // endelse: observation contained an event cube

    // Save FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAObservation::init_members(void)
{
    // Initialise members
    m_instrument = "CTA";
    m_eventfile.clear();
    m_response   = NULL;
    m_pointing   = NULL;
    m_obs_id     = 0;
    m_ontime     = 0.0;
    m_livetime   = 0.0;
    m_deadc      = 0.0;
    m_ra_obj     = 0.0;
    m_dec_obj    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs CTA observation.
 ***************************************************************************/
void GCTAObservation::copy_members(const GCTAObservation& obs)
{
    // Clone members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Copy members
    m_instrument = obs.m_instrument;
    m_eventfile  = obs.m_eventfile;
    m_obs_id     = obs.m_obs_id;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;
    m_ra_obj     = obs.m_ra_obj;
    m_dec_obj    = obs.m_dec_obj;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAObservation::free_members(void)
{
    // Free memory
    if (m_response != NULL) delete m_response;
    if (m_pointing != NULL) delete m_pointing;

    // Mark memory as free
    m_response = NULL;
    m_pointing = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Reads CTA observation attributes from HDU. Mandatory attributes are
 *
 * RA_PNT   - Right Ascension of pointing
 * DEC_PNT  - Declination of pointing
 * ONTIME   - Exposure time
 * LIVETIME - Livetime
 *
 * and optional attributes are
 *
 * OBJECT   - Name of observed object
 * DEADC    - Deadtime correction
 * RA_OBJ   - Right Ascension of observed object,
 * DEC_OBJ  - Declination of observed object,
 * OBS_ID   - Observation identifier
 * ALT_PNT  - Altitude of pointing above horizon
 * AZ_PNT   - Azimuth of pointing
 *
 * Based on RA_PNT and DEC_PNT, the CTA pointing direction is set. Note that
 * DEADC is computed using DEADC=LIVETIME/ONTIME
 *
 * Nothing is done if the HDU pointer is NULL.
 *
 * @todo The actual reader is a minimal reader to accomodate as many
 *       different datasets as possible. Once the CTA data format is fixed
 *       the reader should have more mandatory attributes.
 ***************************************************************************/
void GCTAObservation::read_attributes(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read mandatory attributes
        double ra_pnt  = hdu->real("RA_PNT");
        double dec_pnt = hdu->real("DEC_PNT");
        m_ontime   = (hdu->hascard("ONTIME"))   ? hdu->real("ONTIME") : 0.0;
        m_livetime = (hdu->hascard("LIVETIME")) ? hdu->real("LIVETIME") : 0.0;

        // Read optional attributes
        m_name     = (hdu->hascard("OBJECT"))   ? hdu->string("OBJECT") : "unknown";
        m_deadc    = (hdu->hascard("DEADC"))    ? hdu->real("DEADC") : 0.0;
        m_ra_obj   = (hdu->hascard("RA_OBJ"))   ? hdu->real("RA_OBJ") : 0.0;
        m_dec_obj  = (hdu->hascard("DEC_OBJ"))  ? hdu->real("DEC_OBJ") : 0.0;
        m_obs_id   = (hdu->hascard("OBS_ID"))   ? hdu->integer("OBS_ID") : 0;
        //double alt = (hdu->hascard("ALT_PNT"))  ? hdu->real("ALT_PNT") : 0.0;
        //double az  = (hdu->hascard("AZ_PNT"))   ? hdu->real("AZ_PNT") : 0.0;

        // Kluge: compute DEADC from livetime and ontime instead of using the
        // keyword value as the original event lists had this values badly
        // assigned
        if (m_ontime > 0) {
            m_deadc = m_livetime / m_ontime;
        }
        else {
            m_deadc = 0.0;
        }

        // Set pointing information
        GSkyDir pnt;
        pnt.radec_deg(ra_pnt, dec_pnt);
        if (m_pointing != NULL) delete m_pointing;
        m_pointing = new GCTAPointing;
        m_pointing->dir(pnt);

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Nothing is done if the HDU pointer is NULL.
 ***************************************************************************/
void GCTAObservation::write_attributes(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Compute some attributes
        double ra_pnt  = (m_pointing != NULL) ? m_pointing->dir().ra_deg() : 0.0;
        double dec_pnt = (m_pointing != NULL) ? m_pointing->dir().dec_deg() : 0.0;
        double tstart  = events()->tstart().met();
        double tstop   = events()->tstop().met();
        double telapse = events()->gti().telapse();
        double ontime  = events()->gti().ontime();
        double deadc   = (ontime > 0.0) ? livetime() / ontime : 0.0;

        // Set observation information
        hdu->card("CREATOR",  "GammaLib",   "Program which created the file");
        hdu->card("TELESCOP", instrument(), "Telescope");
        hdu->card("OBS_ID",   obs_id(),     "Observation identifier");
        hdu->card("DATE_OBS", "string",     "Observation start date");
        hdu->card("TIME_OBS", "string",     "Observation start time");
        hdu->card("DATE_END", "string",     "Observation end date");
        hdu->card("TIME_END", "string",     "Observation end time");

        // Set observation time information
        hdu->card("TSTART",   tstart, "[s] Mission time of start of observation");
        hdu->card("TSTOP",    tstop, "[s] Mission time of end of observation");
        hdu->card("MJDREFI",  51910, "[days] Integer part of mission time reference MJD");
        hdu->card("MJDREFF",  7.428703703703703e-14, "[days] Fractional part of mission time reference MJD");
        hdu->card("TIMEUNIT", "s", "Time unit");
        hdu->card("TIMESYS",  "TT", "Time system");
        hdu->card("TIMEREF",  "LOCAL", "Time reference");
        hdu->card("TELAPSE",  telapse, "[s] Mission elapsed time");
        hdu->card("ONTIME",   ontime, "[s] Total good time including deadtime");
        hdu->card("LIVETIME", livetime(), "[s] Total livetime");
        hdu->card("DEADC",    deadc, "Deadtime correction factor");
        hdu->card("TIMEDEL",  1.0, "Time resolution");

        // Set pointing information
        hdu->card("OBJECT",   name(),    "Observed object");
        hdu->card("RA_OBJ",   ra_obj(),  "[deg] Target Right Ascension");
        hdu->card("DEC_OBJ",  dec_obj(), "[deg] Target Declination");
        hdu->card("RA_PNT",   ra_pnt,    "[deg] Pointing Right Ascension");
        hdu->card("DEC_PNT",  dec_pnt,   "[deg] Pointing Declination");
        hdu->card("ALT_PNT",  0.0,       "[deg] Average altitude of pointing");
        hdu->card("AZ_PNT",   0.0,       "[deg] Average azimuth of pointing");
        hdu->card("RADECSYS", "FK5",     "Coordinate system");
        hdu->card("EQUINOX",  2000.0,    "Epoch");
        hdu->card("CONV_DEP", 0.0,       "Convergence depth of telescopes");
        hdu->card("CONV_RA",  0.0,       "[deg] Convergence Right Ascension");
        hdu->card("CONV_DEC", 0.0,       "[deg] Convergence Declination");
        hdu->card("OBSERVER", "string",  "Observer");

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
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Temporally integrate spatially & spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * Implement the temporal integration as a simple multiplication by the
 * ontime. This assumes that the source is non-variable during the
 * observation and that the CTA pointing is stable.
 *
 * Note that we use here the ontime as the GResponse::irf method returns
 * already deadtime corrected response values.
 ***************************************************************************/
double GCTAObservation::npred_temp(const GModel& model) const
{
    // Initialise result
    double result = 0.0;

    // Determine ontime
    double ontime = events()->gti().ontime();

    // Integrate only if ontime is positive
    if (ontime > 0.0) {

        // Integration is a simple multiplication by the time
        result = npred_spec(model, events()->gti().tstart()) * ontime;

    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
