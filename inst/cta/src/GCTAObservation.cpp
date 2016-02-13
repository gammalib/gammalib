/***************************************************************************
 *                GCTAObservation.cpp - CTA Observation class              *
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
 * @file GCTAObservation.cpp
 * @brief CTA observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GObservationRegistry.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GFilename.hpp"
#include "GGti.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GCTASupport.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAResponseCube.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"

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
#define G_RESPONSE_SET                "GCTAObservation::response(GResponse&)"
#define G_RESPONSE_GET                          "GCTAObservation::response()"
#define G_ROI                                        "GCTAObservation::roi()"
#define G_GTI                                        "GCTAObservation::gti()"
#define G_EBOUNDS                                "GCTAObservation::ebounds()"
#define G_READ                          "GCTAObservation::read(GXmlElement&)"
#define G_WRITE                        "GCTAObservation::write(GXmlElement&)"
#define G_LOAD           "GCTAObservation::load(std::string&, std::string&, "\
                                                "std::string&, std::string&)"

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
 * Constructs an empty CTA observation.
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
 * Constructs an empty CTA observation for a given instrument. This enables
 * using the CTA specific interface for any other VHE instrument. Note that
 * each other VHE instruments needs a specific registry at the beginning
 * of the GCTAObservation.cpp file. So far the following instruments are
 * supported: CTA, HESS, VERITAS, MAGIC.
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
 * @brief Stacked cube analysis constructor
 *
 * @param[in] cntcube Counts cube file name.
 * @param[in] expcube Exposure cube file name.
 * @param[in] psfcube Psf cube file name.
 * @param[in] bkgcube Backgorund cube file name.
 *
 * Constructs a CTA observation from a counts cube, an exposure cube, a Psf
 * cube and a background cube.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const std::string& cntcube,
                                 const std::string& expcube,
                                 const std::string& psfcube,
                                 const std::string& bkgcube) : GObservation()
{
    // Initialise members
    init_members();

    // Load data
    load(cntcube, expcube, psfcube, bkgcube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs CTA observation.
 *
 * Constructs a CTA observation by copying an existing CTA observation.
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
 *
 * Destructs CTA observation.
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
GCTAObservation& GCTAObservation::operator=(const GCTAObservation& obs)
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
 * @brief Clear CTA observation
 *
 * Clear CTA observation.
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
 * @brief Clone CTA observation
 *
 * @return Pointer to deep copy of CTA observation.
 *
 * Returns a pointer to a deep copy of a CTA observation.
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
 * @exception GException::invalid_argument
 *            Invalid response class specified.
 *
 * Sets the response function for the observation.
 ***************************************************************************/
void GCTAObservation::response(const GResponse& rsp)
{
    // Cast response dynamically
    const GCTAResponse* ptr = dynamic_cast<const GCTAResponse*>(&rsp);

    // Throw exception if response is not of correct type
    if (ptr == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Invalid response type \""+cls+"\" provided on "
                          "input. Please specify a \"GCTAResponse\" "
                          "object as argument.";
        throw GException::invalid_argument(G_RESPONSE_SET, msg);
    }

    // Free response
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = ptr->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer to CTA response function
 *
 * @return Pointer to CTA response function.
 *
 * @exception GException::invalid_value
 *            No valid response found in CTA observation.
 *
 * Returns a pointer to the CTA response function. An exception is thrown if
 * the pointer is not valid, hence the user does not need to verify the
 * validity of the pointer.
 ***************************************************************************/
const GCTAResponse* GCTAObservation::response(void) const
{
    // Throw an exception if the response pointer is not valid
    if (m_response == NULL) {
        std::string msg = "No valid response function found in CTA "
                          "observation.\n";
        throw GException::invalid_value(G_RESPONSE_GET, msg);
    }

    // Return pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Set CTA response function
 *
 * @param[in] rspname Response name.
 * @param[in] caldb Calibration database.
 *
 * Sets the CTA response function by specifying a response name and a
 * calibration database. This method also loads the response function so that
 * it is available for data analysis.
 ***************************************************************************/
void GCTAObservation::response(const std::string& rspname, const GCaldb& caldb)
{
    // Free response
    if (m_response != NULL) delete m_response;
    m_response = NULL;

    // Allocate fresh response function
    GCTAResponseIrf* rsp = new GCTAResponseIrf;

    // Set calibration database
    rsp->caldb(caldb);

    // Load instrument response function
    rsp->load(rspname);

    // Store pointer
    m_response = rsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set CTA response function
 *
 * @param[in] expcube Exposure cube.
 * @param[in] psfcube Psf cube.
 * @param[in] bkgcube Background cube.
 *
 * Sets the CTA response function fur cube analysis by specifying the
 * exposure cube, the Psf cube and the background cube. The method also
 * copies over the ontime, the livetime and the deadtime correction factor
 * from the exposure cube.
 ***************************************************************************/
void GCTAObservation::response(const GCTACubeExposure&   expcube,
                               const GCTACubePsf&        psfcube,
                               const GCTACubeBackground& bkgcube)
{
    // Free response
    if (m_response != NULL) delete m_response;
    m_response = NULL;

    // Allocate fresh response function
    GCTAResponseCube* rsp = new GCTAResponseCube(expcube, psfcube, bkgcube);

    // Store pointer
    m_response = rsp;

    // Copy over time information from exposure cube
    ontime(expcube.ontime());
    livetime(expcube.livetime());
    deadc(expcube.deadc());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get Region of Interest
 *
 * @return Region of Interest.
 *
 * @exception GException::invalid_value
 *            Observation does not contain events.
 *            Observation does not contain an event list.
 *
 * Extracts the Region of Interest from the event list. An exception is
 * thrown if the observation does not contain an event list.
 ***************************************************************************/
GCTARoi GCTAObservation::roi(void) const
{
    // Throw an exception if no events exist
    if (m_events == NULL) {
        std::string msg = "Region of Interest is not defined since the "
                          "observation does not contain events.";
        throw GException::invalid_value(G_ROI, msg);
    }

    // Get pointer to event list. Throw an exception if no event list is found
    const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(m_events);
    if (list == NULL) {
        std::string msg = "Region of Interest is not defined since the "
                          "observation does not contain an event list.";
        throw GException::invalid_value(G_ROI, msg);
    }

    // Return ROI
    return (list->roi());
}


/***********************************************************************//**
 * @brief Get Good Time Intervals
 *
 * @return Good Time Intervals.
 *
 * @exception GException::invalid_value
 *            Observation does not contain events.
 *
 * Extracts the Good Time Intervals from the events. An exception is thrown
 * if the observation does not contain events.
 ***************************************************************************/
GGti GCTAObservation::gti(void) const
{
    // Throw an exception if no events exist
    if (m_events == NULL) {
        std::string msg = "Good Time Intervals are not defined since the "
                          "observation does not contain events.";
        throw GException::invalid_value(G_GTI, msg);
    }

    // Return GTI
    return (m_events->gti());
}


/***********************************************************************//**
 * @brief Get energy boundaries
 *
 * @return Energy boundaries.
 *
 * @exception GException::invalid_value
 *            Observation does not contain events.
 *
 * Extract the energy boundaries from the events. An exception is thrown if
 * the observation does not contain events.
 ***************************************************************************/
GEbounds GCTAObservation::ebounds(void) const
{
    // Throw an exception if no events exist
    if (m_events == NULL) {
        std::string msg = "Energy boundaries are not defined since the "
                          "observation does not contain events.";
        throw GException::invalid_value(G_EBOUNDS, msg);
    }

    // Return energy boundaries
    return (m_events->ebounds());
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
 * Reads information for a CTA observation from an XML element. This method
 * handles to variants: a first where an event list of counts cube is
 * given and a second where the observation definition information is
 * provided by parameter tags.
 *
 * The XML format for an event list is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="EventList" file="..."/>
 *       ...
 *     </observation>
 *
 * and for a counts cube is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="CountsCube" file="..."/>
 *       ...
 *     </observation>
 *
 * The second variant without event information has the following format
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="Pointing" ra="..." dec="..."/>
 *       <parameter name="Deadtime" deadc="..."/>
 *       ...
 *     </observation>
 *
 * In addition, calibration information can be specified using the format
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="Calibration" database="..." response="..."/>
 *     </observation>
 *
 * If even more control is required over the response, individual file names
 * can be specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *       <parameter name="Background"          file="..."/>
 *     </observation>
 *
 * In case that a @a CountsCube is handled, the stacked response can also be
 * provided in the format
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *        <parameter name="BkgCube"      file="..."/>
 *      </observation>
 *   
 * Optional user energy thresholds can be specified by adding the @a emin
 * and @a emax attributes to the @p observation tag:
 *
 *     <observation name="..." id="..." instrument="..." emin="..." emax="...">
 *
 * The method does no load the events into memory but stores the file name
 * of the event file. The events are only loaded when required. This reduces
 * the memory needs for an CTA observation object and allows for loading
 * of event information upon need.
 ***************************************************************************/
void GCTAObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Read in user defined energy boundaries of this observation
    if (xml.attribute("emin") != "") {
        m_lo_user_thres = gammalib::todouble(xml.attribute("emin"));
    }
    if (xml.attribute("emax") != "") {
        m_hi_user_thres = gammalib::todouble(xml.attribute("emax"));
    }

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has at least 1 parameter
    if (xml.elements() < 1 || npars < 1) {
        throw GException::xml_invalid_parnum(G_READ, xml,
              "CTA observation requires at least 1 parameter.");
    }

    // First try to extract observation parameters for an event file
    int n_evfile = 0;
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle EventList or CountsCube
        if ((par->attribute("name") == "EventList") ||
            (par->attribute("name") == "CountsCube")) {

            // Read eventlist file name
            std::string filename = par->attribute("file");

            // Load events
            load(filename);

            // Increment parameter counter
            n_evfile++;
        }

        // Read Background filename (needed by GCTAModelCubeBackground)
        else if (par->attribute("name") == "Background") {

            // Read background file name
            m_bgdfile = par->attribute("file");

        }

    } // endfor: looped over observation parameters

    // Analyse parameters
    bool has_evfile = (n_evfile > 0);

    // If we have an event file then verify that all required parameters
    // were found
    if (has_evfile) {
        gammalib::xml_check_par(G_READ, "EventList\" or \"CountsCube", n_evfile);
    }

    // ... otherwise extract information from observation definition parameters
    else {

        // Initialise event list definition
        GEbounds ebounds;
        GGti     gti;
        GCTARoi  roi;
        bool     has_ebounds = false;
        bool     has_gti     = false;
        bool     has_roi     = false;

        // Extract parameters (do nothing if parameter does not exist)
        if (gammalib::xml_has_par(xml, "Pointing")) {
            m_pointing.read(xml);
        }
        if (gammalib::xml_has_par(xml, "EnergyBoundaries")) {
            ebounds.read(xml);
            has_ebounds = true;
        }
        if (gammalib::xml_has_par(xml, "GoodTimeIntervals")) {
            gti.read(xml);
            has_gti = true;
        }
        if (gammalib::xml_has_par(xml, "RegionOfInterest")) {
            roi.read(xml);
            has_roi = true;
        }
        if (gammalib::xml_has_par(xml, "Deadtime")) {
            const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Deadtime");
            m_deadc                = gammalib::todouble(par->attribute("deadc"));
        }
        else {
            m_deadc = 1.0;
        }

        // If we have at least one event list definition then create event
        // list and attach it to the observation
        if (has_ebounds || has_gti || has_roi) {

            // Create an empty event list
            GCTAEventList events;

            // Handle energy boundaries
            if (has_ebounds) {
                events.ebounds(ebounds);
            }

            // Handle GTIs
            if (has_gti) {
                events.gti(gti);
                ontime(gti.ontime());
                livetime(gti.ontime() * m_deadc);
            }

            // Handle ROIs
            if (has_roi) {
                events.roi(roi);
            }

            // Attach event list
            this->events(events);

        } // endif: handled event list information

    } // endelse: extracted information from observation definition

    // Determine response type as function of the information that is
    // provided in the XML file. The response type is determined from the
    // parameter names. An exception is thrown if the response type cannot
    // be unambigously determined
    int response_type = 0;
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Check for response type 1 (GCTAResponseIrf)
        if ((par->attribute("name") == "EffectiveArea")       ||
            (par->attribute("name") == "ARF")                 ||
            (par->attribute("name") == "PointSpreadFunction") ||
            (par->attribute("name") == "PSF")                 ||
            (par->attribute("name") == "EnergyDispersion")    ||
            (par->attribute("name") == "RMF")                 ||
            (par->attribute("name") == "Background")          ||
            (par->attribute("name") == "Calibration")) {
            if (response_type == 2) {
                throw GException::xml_invalid_parnames(G_READ, xml,
                      "Incompatible parameter names encountered in the "
                      "response definition of a CTA observation.\n");
            }
            response_type = 1;
        }
        
        // Check for response type 2 (GCTAResponseCube)
        else if ((par->attribute("name") == "ExposureCube") ||
                 (par->attribute("name") == "PsfCube")      ||
                 (par->attribute("name") == "BkgCube")) {
            if (response_type == 1) {
                throw GException::xml_invalid_parnames(G_READ, xml,
                      "Incompatible parameter names encountered in the "
                      "response definition of a CTA observation.\n");
            }
            response_type = 2;
        }

    } // endfor: looped over all parameters
    
    // Allocate response
    switch (response_type) {
        case 1:
            m_response = new GCTAResponseIrf;
            break;
        case 2:
            m_response = new GCTAResponseCube;
            break;
        default:
            break;
    }

    // Extract response information
    if (m_response != NULL) {
        m_response->read(xml);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            No valid events found in observation.
 *
 * Writes information for a CTA observation into an XML element. This method
 * handles to variants: a first where an event list of counts cube is
 * given and a second where the observation definition information is
 * provided by parameter tags. Note that in both bases, a valid event
 * type needs to be set (either @a EventList or @a CountsCube).
 *
 * The XML format for an event list is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="EventList" file="..."/>
 *       ...
 *     </observation>
 *
 * and for a counts cube is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="CountsCube" file="..."/>
 *       ...
 *     </observation>
 *
 * The second variant without event information has the following format
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="Pointing" ra="..." dec="..."/>
 *       <parameter name="Deadtime" deadc="..."/>
 *       ...
 *     </observation>
 *
 * In addition, calibration information can be specified using the format
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="Calibration" database="..." response="..."/>
 *     </observation>
 *
 * If even more control is required over the response, individual file names
 * can be specified using
 *
 *     <observation name="..." id="..." instrument="...">
 *       ...
 *       <parameter name="EffectiveArea"       file="..."/>
 *       <parameter name="PointSpreadFunction" file="..."/>
 *       <parameter name="EnergyDispersion"    file="..."/>
 *       <parameter name="Background"          file="..."/>
 *     </observation>
 *
 * In case that a @a CountsCube is handled, the stacked response can also be
 * provided in the format
 *
 *      <observation name="..." id="..." instrument="...">
 *        ...
 *        <parameter name="ExposureCube" file="..."/>
 *        <parameter name="PsfCube"      file="..."/>
 *        <parameter name="BkgCube"      file="..."/>
 *      </observation>
 *   
 * If user energy thresholds are defined (i.e. threshold values are >0) the
 * additional @a emin and @a emax attributes will be written to the
 * @p observation tag:
 *
 *     <observation name="..." id="..." instrument="..." emin="..." emax="...">
 *
 ***************************************************************************/
void GCTAObservation::write(GXmlElement& xml) const
{
    // Throw an exception if eventtype() is neither "EventList" nor
    std::string evttype = eventtype();
    if ((evttype != "EventList") && (evttype != "CountsCube")) {
        std::string msg;
        if (evttype.empty()) {
            msg = "The observation does not contain any events, hence "
                  "it cannot be written into an XML file.";
        }
        else {
            msg = "The observation contains an unknown event type \""+
                  evttype+"\". The event type needs to be either "
                  "\"EventList\" or \"CountsCube\".";
        }
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Add user energy threshold attributes (if set)
    if (m_lo_user_thres > 0.0) {
        xml.attribute("emin", gammalib::str(m_lo_user_thres));
    }
    if (m_hi_user_thres > 0.0) {
        xml.attribute("emax", gammalib::str(m_hi_user_thres));
    }

    // If there is an event filename then write the event information to the
    // XML file ...
    if (!m_eventfile.empty()) {

        // Write event file name
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, evttype);
        par->attribute("file", m_eventfile);

    }

    // ... otherwise write the observation definition information
    else {

        // Write pointing, energy bounds, GTIs and ROI
        m_pointing.write(xml);
        events()->ebounds().write(xml);
        events()->gti().write(xml);

        // Write ROI
        const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(events());
        if (list != NULL) {
            list->roi().write(xml);
        }

        // Write deadtime correction factor
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "Deadtime");
        par->attribute("deadc",  gammalib::str(m_deadc));

    }

    // Write response information
    if (m_response != NULL) {
        m_response->write(xml);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read data from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads event data from a FITS file and sets the observation attributes.
 *
 * The method automatically switches between an event list and a counts cube,
 * depending of the information provided in the FITS file. If an extension
 * name is specified, the method checks whether the extension exists and
 * loads the extension as event list. Otherwise, it checks whether the file
 * contains an "EVENTS" extension and loads the extension as event list.
 * If none of the above are satistified, the method loads a counts cube.
 ***************************************************************************/
void GCTAObservation::read(const GFits& fits)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Initialise file name from FITS file
    GFilename fname(fits.filename());

    // Get extension name
    std::string extname = fname.extname("EVENTS");

    // If FITS file contains an EVENTS extension we have an unbinned
    // observation ...
    if (fits.contains(extname)) {

        // Allocate event list
        GCTAEventList* events = new GCTAEventList;

        // Assign event list as the observation's event container
        m_events = events;

        // Read event list
        events->read(fits);

        // Read observation attributes from EVENTS extension
        const GFitsHDU& hdu = *fits.at(extname);
        read_attributes(hdu);

    }

    // ... otherwise we have a binned observation
    else {

        // Allocate event cube
        GCTAEventCube* events = new GCTAEventCube;

        // Assign event cube as the observation's event container
        m_events = events;

        // Read event cube
        events->read(fits);

        // Read observation attributes from primary extension
        const GFitsHDU& hdu = *fits.at(0);
        read_attributes(hdu);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into FITS file.
 *
 * @param[in] fits FITS file.
 * @param[in] evtname Events FITS extension name (default: "EVENTS").
 * @param[in] gtiname Good Time Intervals FITS extension name (default: "GTI").
 *
 * Writes the observation into a FITS file.
 *
 * If the observation contains an event list, the event list and Good Time
 * Intervals will be written to the FITS file. The FITS extension name for
 * the event list and the Good Time Intervals can be optionally specified
 * thru the @p evtname and @p gtiname arguments.
 *
 * If the observation contains an event cube, the event cube will be written
 * into the primary extension of the FITS file.
 *
 * This method does nothing if no events are in the CTA observation.
 ***************************************************************************/
void GCTAObservation::write(GFits& fits, const std::string& evtname,
                                         const std::string& gtiname) const
{
    // Get pointers on event list
    const GCTAEventList* list = dynamic_cast<const GCTAEventList*>(events());
    const GCTAEventCube* cube = dynamic_cast<const GCTAEventCube*>(events());

    // Case A: Observation contains an event list
    if (list != NULL) {

        // Write event list and Good Time Intervals into FITS file
        list->write(fits, evtname, gtiname);

        // Get reference to events extension
        GFitsHDU& hdu = *fits.at(evtname);

        // Write observation attributes to events extension
        write_attributes(hdu);

    } // endif: observation contained an event list

    // Case B: Observation contains an event cube
    else if (cube != NULL) {

        // Write events cube into FITS file. This method also writes
        // the energy boundaries and the GTI as they are also part
        // of the event cube.
        cube->write(fits);

        // Write observation attributes into primary header
        GFitsHDU& hdu = *fits.at(0);
        write_attributes(hdu);

    } // endelse: observation contained an event cube

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load event list or counts cube from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads either an event list or a counts cube from a FITS file.
 ***************************************************************************/
void GCTAObservation::load(const std::string& filename)
{
    // Store event filename
    m_eventfile = filename;

    // Open FITS file
    GFits fits(filename);

    // Read data
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load stacked cube from FITS files
 *
 * @param[in] cntcube Counts cube file name.
 * @param[in] expcube Exposure cube file name.
 * @param[in] psfcube Psf cube file name.
 * @param[in] bkgcube Background cube file name.
 *
 * Loads a counts map, an exposure cube, a Psf cube and a background cube
 * for stacked cube analysis.
 ***************************************************************************/
void GCTAObservation::load(const std::string& cntcube,
                           const std::string& expcube,
                           const std::string& psfcube,
                           const std::string& bkgcube) {

    // Load counts cube FITS file
    load(cntcube);

    // Check whether we have an event cube
    if (dynamic_cast<const GCTAEventCube*>(events()) == NULL) {
        std::string msg = "Specified file \""+cntcube+"\" is not a CTA "
                          "counts cube. Please provide a counts cube.";
        throw GException::invalid_argument(G_LOAD, msg);
    }

    // Load exposure cube
    GCTACubeExposure exposure(expcube);

    // Load Psf cube
    GCTACubePsf psf(psfcube);

    // Load background cube
    GCTACubeBackground background(bkgcube);

    // Attach exposure cube, Psf cube and background cube as response
    response(exposure, psf, background);
 
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save CTA observation into FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file? (default: false)
 *
 * Saves the CTA observation into a FITS file.
 *
 * If the CTA observation contains an event list, the method will write an
 * events and a Good Time Intervals extension to the FITS file. The names
 * for both extension can be optionally specified in the filename using
 * the format
 *
 *      filename[EVENTS;GTI]
 *
 * where the string in the squared bracket defines the extension names. The
 * part before the semi-colon is the events extension name and the part after
 * the semi-colon is the Good Time Intervals extension name.
 *
 * If the CTA observation contains an event cube, the method will write the
 * cube into the primary image, followed by binary tables containing the
 * energy boundaries and the Good Time Intervals. The extension names of
 * these binary tables are "EBOUNDS" and "GTI", and cannot be modified.
 ***************************************************************************/
void GCTAObservation::save(const std::string& filename, const bool& clobber) const
{
    // Initialise filename
    GFilename fname(filename);

    // Set default events and Good Time Intervals extension name and extract
    // a possible overwrite from the extension name argument of the filename.
    // The specific format that is implemented is [events;gti], where the
    // part before the semi-colon is the events extension name and the part
    // after the semi-colon is the Good Time Intervals extension name.
    std::string evtname = "EVENTS";
    std::string gtiname = "GTI";
    if (fname.has_extname()) {
        std::vector<std::string> extnames = gammalib::split(fname.extname(), ";");
        if (extnames.size() > 0) {
            evtname = gammalib::strip_whitespace(extnames[0]);
        }
        if (extnames.size() > 1) {
            gtiname = gammalib::strip_whitespace(extnames[1]);
        }
    }

    // Create FITS file
    GFits fits;

    // Write data into FITS file
    write(fits, evtname, gtiname);

    // Save FITS file
    fits.saveto(fname.filename(), clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return event type string
 *
 * @return Event type string.
 *
 * Returns "EventList" if the observation contains an event list,
 * "CountsCube" if it contains a counts cube, "Events" if it container an
 * unknown type of events (which should never occur), and an empty string if
 * no events have been allocated.
 ***************************************************************************/
std::string GCTAObservation::eventtype(void) const
{
    // Initialise empty event type string
    std::string eventtype = "";

    // Continue only if events are allocated
    if (m_events != NULL) {

        // Case A: we have a list
        GCTAEventList* list = dynamic_cast<GCTAEventList*>(m_events);
        if (list != NULL) {
            eventtype = "EventList";
        }

        // Case B: we have a cube
        else {
            GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(m_events);
            if (cube != NULL) {
                eventtype = "CountsCube";
            }

            // Case C: we don't know what we have
            else {
                eventtype = "Events";
            }
        }

    } // endif: events were allocated

    // Return event type
    return eventtype;
}


/***********************************************************************//**
 * @brief Dispose events
 *
 * Disposes the events to save memory. The method only applies to event
 * lists. If does nothing in case that the observation does not contain an
 * event list.
 ***************************************************************************/
void GCTAObservation::dispose_events(void)
{
    // Get pointer to event list
    GCTAEventList* list = dynamic_cast<GCTAEventList*>(m_events);

    // Dispose event list if pointer is valid
    if (list != NULL) {
        list->dispose();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print CTA observation information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing observation information.
 ***************************************************************************/
std::string GCTAObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAObservation ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Identifier")+id());
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Event file")+eventfile());
        result.append("\n"+gammalib::parformat("Event type")+eventtype());
        result.append("\n"+gammalib::parformat("Statistics")+statistics());
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" s");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" s");
        result.append("\n"+gammalib::parformat("Deadtime correction"));
        result.append(gammalib::str(m_deadc));

        // Append user energy threshold information
        result.append("\n"+gammalib::parformat("User energy range"));
        if (m_lo_user_thres > 0.0 && m_hi_user_thres) {
            result.append(gammalib::str(m_lo_user_thres));
            result.append(" - ");
            result.append(gammalib::str(m_hi_user_thres));
            result.append(" TeV");
        }
        else if (m_lo_user_thres > 0.0) {
            result.append("> ");
            result.append(gammalib::str(m_lo_user_thres));
            result.append(" TeV");
        }
        else if (m_hi_user_thres > 0.0) {
            result.append("< ");
            result.append(gammalib::str(m_hi_user_thres));
            result.append(" TeV");
        }
        else {
            result.append("undefined");
        }

        // Append detailed information
        GChatter reduced_chatter = gammalib::reduce(chatter);
        if (reduced_chatter > SILENT) {

            // Append pointing
            result.append("\n"+pointing().print(reduced_chatter));

            // Append response
            if (m_response != NULL) {
                result.append("\n"+m_response->print(reduced_chatter));
            }
            else {
                result.append("\n"+gammalib::parformat("Response function"));
                result.append("undefined");
            }

            // Append events
            if (m_events != NULL) {
                result.append("\n"+m_events->print(reduced_chatter));
            }
        
        } // endif: appended detailed information

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
void GCTAObservation::init_members(void)
{
    // Initialise members
    m_instrument = "CTA";
    m_object.clear();
    m_eventfile.clear();
    m_bgdfile.clear();
    m_response = NULL;
    m_pointing.clear();
    m_obs_id        = 0;
    m_ontime        = 0.0;
    m_livetime      = 0.0;
    m_deadc         = 1.0;
    m_ra_obj        = 0.0;
    m_dec_obj       = 0.0;
    m_lo_user_thres = 0.0;
    m_hi_user_thres = 0.0;
    m_n_tels        = 0;

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
    // Copy members
    m_instrument    = obs.m_instrument;
    m_object        = obs.m_object;
    m_eventfile     = obs.m_eventfile;
    m_bgdfile       = obs.m_bgdfile;
    m_pointing      = obs.m_pointing;
    m_obs_id        = obs.m_obs_id;
    m_ontime        = obs.m_ontime;
    m_livetime      = obs.m_livetime;
    m_deadc         = obs.m_deadc;
    m_ra_obj        = obs.m_ra_obj;
    m_dec_obj       = obs.m_dec_obj;
    m_lo_user_thres = obs.m_lo_user_thres;
    m_hi_user_thres = obs.m_hi_user_thres;
    m_n_tels        = obs.m_n_tels;

    // Clone members
    m_response = (obs.m_response != NULL) ? obs.m_response->clone() : NULL;

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

    // Initialise pointers
    m_response = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation attributes
 *
 * @param[in] hdu FITS HDU.
 *
 * Reads CTA observation attributes from HDU. Mandatory attributes are
 *
 *     RA_PNT   - Right Ascension of pointing
 *     DEC_PNT  - Declination of pointing
 *     ONTIME   - Exposure time
 *     LIVETIME - Livetime
 *
 * and optional attributes are
 *
 *     OBJECT   - Name of observed object
 *     DEADC    - Deadtime correction
 *     RA_OBJ   - Right Ascension of observed object,
 *     DEC_OBJ  - Declination of observed object,
 *     OBS_ID   - Observation identifier
 *     ALT_PNT  - Altitude of pointing above horizon
 *     AZ_PNT   - Azimuth of pointing
 *     TELESCOP - Telescope name
 *     N_TELS   - Number of telescopes
 *
 * Based on RA_PNT and DEC_PNT, the CTA pointing direction is set. Note that
 * DEADC is computed using DEADC=LIVETIME/ONTIME
 *
 * @todo The actual reader is a minimal reader to accomodate as many
 *       different datasets as possible. Once the CTA data format is fixed
 *       the reader should have more mandatory attributes.
 ***************************************************************************/
void GCTAObservation::read_attributes(const GFitsHDU& hdu)
{
    // Read mandatory attributes
    double ra_pnt  = hdu.real("RA_PNT");
    double dec_pnt = hdu.real("DEC_PNT");
    m_ontime   = (hdu.has_card("ONTIME"))   ? hdu.real("ONTIME") : 0.0;
    m_livetime = (hdu.has_card("LIVETIME")) ? hdu.real("LIVETIME") : 0.0;

    // Read optional attributes
    m_object     = (hdu.has_card("OBJECT"))   ? hdu.string("OBJECT") : "unknown";
    m_deadc      = (hdu.has_card("DEADC"))    ? hdu.real("DEADC") : 0.0;
    m_ra_obj     = (hdu.has_card("RA_OBJ"))   ? hdu.real("RA_OBJ") : 0.0;
    m_dec_obj    = (hdu.has_card("DEC_OBJ"))  ? hdu.real("DEC_OBJ") : 0.0;
    m_obs_id     = (hdu.has_card("OBS_ID"))   ? hdu.integer("OBS_ID") : 0;
    double alt   = (hdu.has_card("ALT_PNT"))  ? hdu.real("ALT_PNT") : 0.0;
    double az    = (hdu.has_card("AZ_PNT"))   ? hdu.real("AZ_PNT") : 0.0;
    m_instrument = (hdu.has_card("TELESCOP")) ? hdu.string("TELESCOP") : "CTA";
    m_n_tels     = (hdu.has_card("N_TELS"))   ? hdu.integer("N_TELS") : 0;

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
    m_pointing.dir(pnt);
    m_pointing.zenith(90.0-alt);
    m_pointing.azimuth(az);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation attributes
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GCTAObservation::write_attributes(GFitsHDU& hdu) const
{
    // Get time reference
    GTimeReference timeref = events()->gti().reference();

    // Compute some attributes
    double      ra_pnt   = m_pointing.dir().ra_deg();
    double      dec_pnt  = m_pointing.dir().dec_deg();
    double      alt      = 90.0 - m_pointing.zenith();
    double      az       = m_pointing.azimuth();
    double      tstart   = events()->tstart().convert(timeref);
    double      tstop    = events()->tstop().convert(timeref);
    double      telapse  = events()->gti().telapse();
    double      ontime   = events()->gti().ontime();
    double      deadc    = (ontime > 0.0 && livetime() > 0.0) ? 
                           livetime() / ontime : m_deadc;
    std::string utc_obs  = events()->tstart().utc();
    std::string utc_end  = events()->tstop().utc();
    std::string date_obs = utc_obs.substr(0, 10);
    std::string time_obs = utc_obs.substr(11, 8);
    std::string date_end = utc_end.substr(0, 10);
    std::string time_end = utc_end.substr(11, 8);

    // Set observation information
    hdu.card("CREATOR",  "GammaLib",   "Program which created the file");
    hdu.card("TELESCOP", instrument(), "Telescope");
    hdu.card("OBS_ID",   obs_id(),     "Observation identifier");
    hdu.card("DATE_OBS", date_obs,     "Observation start date");
    hdu.card("TIME_OBS", time_obs,     "Observation start time");
    hdu.card("DATE_END", date_end,     "Observation end date");
    hdu.card("TIME_END", time_end,     "Observation end time");

    // Set observation time information
    hdu.card("TSTART",   tstart, "[s] Mission time of start of observation");
    hdu.card("TSTOP",    tstop, "[s] Mission time of end of observation");
    timeref.write(hdu);
    hdu.card("TELAPSE",  telapse, "[s] Mission elapsed time");
    hdu.card("ONTIME",   ontime, "[s] Total good time including deadtime");
    hdu.card("LIVETIME", livetime(), "[s] Total livetime");
    hdu.card("DEADC",    deadc, "Deadtime correction factor");
    hdu.card("TIMEDEL",  1.0, "Time resolution");

    // Set pointing information
    hdu.card("OBJECT",   object(),  "Observed object");
    hdu.card("RA_OBJ",   ra_obj(),  "[deg] Target Right Ascension");
    hdu.card("DEC_OBJ",  dec_obj(), "[deg] Target Declination");
    hdu.card("RA_PNT",   ra_pnt,    "[deg] Pointing Right Ascension");
    hdu.card("DEC_PNT",  dec_pnt,   "[deg] Pointing Declination");
    hdu.card("ALT_PNT",  alt,       "[deg] Average altitude of pointing");
    hdu.card("AZ_PNT",   az,        "[deg] Average azimuth of pointing");
    hdu.card("RADECSYS", "FK5",     "Coordinate system");
    hdu.card("EQUINOX",  2000.0,    "Epoch");
    hdu.card("CONV_DEP", 0.0,       "Convergence depth of telescopes");
    hdu.card("CONV_RA",  0.0,       "[deg] Convergence Right Ascension");
    hdu.card("CONV_DEC", 0.0,       "[deg] Convergence Declination");
    hdu.card("OBSERVER", "string",  "Observer");

    // Telescope information
    hdu.card("N_TELS",   n_tels(), "Number of telescopes in event list");
    hdu.card("TELLIST",  "string", "Telescope IDs");
    hdu.card("GEOLAT",   0.0,      "[deg] Geographic latitude of array centre");
    hdu.card("GEOLON",   0.0,      "[deg] Geographic longitude of array centre");
    hdu.card("ALTITUDE", 0.0,      "[km] Altitude of array centre");

    // Other information
    hdu.card("EUNIT",    "TeV",    "Energy unit");
    hdu.card("EVTVER",   "draft1", "Event list version number");

    // Return
    return;
}
