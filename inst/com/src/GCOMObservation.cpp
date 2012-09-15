/***************************************************************************
 *            GCOMObservation.cpp  -  COMPTEL Observation class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMObservation.cpp
 * @brief COMPTEL observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventCube.hpp"

/* __ Globals ____________________________________________________________ */
const GCOMObservation      g_obs_com_seed;
const GObservationRegistry g_obs_com_registry(&g_obs_com_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GCOMObservation::response(GResponse&)"
#define G_READ                          "GCOMObservation::read(GXmlElement&)"
#define G_WRITE                        "GCOMObservation::write(GXmlElement&)"

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
GCOMObservation::GCOMObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Instrument constructor
 *
 * @param[in] drename Event cube name.
 * @param[in] drbname Background cube name.
 * @param[in] drgname Geometry cube name.
 * @param[in] drxname Exposure map name.
 *
 * Creates COMPTEL observation by loading the relevant FITS files.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const std::string& drename,
                                 const std::string& drbname,
                                 const std::string& drgname,
                                 const std::string& drxname) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(drename, drbname, drgname, drxname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs COMPTEL observation.
 *
 * Creates class instance by copying an existing COMPTEL observation.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GCOMObservation& obs) : GObservation(obs)
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
GCOMObservation::~GCOMObservation(void)
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
 * @param[in] obs COMPTEL observation.
 *
 * Assign COMPTEL observation to this object.
 ***************************************************************************/
GCOMObservation& GCOMObservation::operator= (const GCOMObservation& obs)
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
void GCOMObservation::clear(void)
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
GCOMObservation* GCOMObservation::clone(void) const
{
    return new GCOMObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GException::rsp_invalid_type
 *            Specified response is not of type GCOMResponse.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GCOMResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GCOMObservation::response(const GResponse& rsp)
{
    // Get pointer on COM response
    const GCOMResponse* comrsp = dynamic_cast<const GCOMResponse*>(&rsp);
    if (comrsp == NULL) {
        throw GException::rsp_invalid_type(G_RESPONSE,
              typeid(&rsp).name(), "Expected type \"GCOMResponse\".");
    }

    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = comrsp->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set COM response function
 *
 * @param[in] iaqname Name of COMPTEL IAQ response function.
 * @param[in] caldb Optional name of calibration database.
 ***************************************************************************/
void GCOMObservation::response(const std::string& iaqname,
                               const std::string& caldb)
{
    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Allocate new COM response function
    m_response = new GCOMResponse;

    // Set calibration database
    m_response->caldb(caldb);

    // Load instrument response function
    m_response->load(iaqname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to COMPTEL response function
 ***************************************************************************/
GCOMResponse* GCOMObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to COMPTEL pointing direction
 *
 * Returns pointer to pointing direction.
 ***************************************************************************/
GCOMPointing* GCOMObservation::pointing(void) const
{
    // Return pointing pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Set COMPTEL pointing direction
 *
 * @param[in] pointing Pointing.
 ***************************************************************************/
/*
void GCOMObservation::pointing(const GCOMPointing& pointing)
{
    // Free any existing pointing
    if (m_pointing != NULL) delete m_pointing;
    m_pointing = NULL;

    // Clone pointing
    m_pointing = pointing.clone();

    // Return
    return;
}
*/


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
 * Reads information for a COMPTEL observation from an XML element.
 * The expected format of the XML element is
 *
 *  <observation name="..." id="..." instrument="...">
 *    <parameter name="DRE" file="..."/>
 *    <parameter name="DRB" file="..."/>
 *    <parameter name="DRG" file="..."/>
 *    <parameter name="DRX" file="..."/>
 *  </observation>
 *
 * for a binned observation.
 ***************************************************************************/
void GCOMObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Initialise FITS file names
    std::string drename;
    std::string drbname;
    std::string drgname;
    std::string drxname;

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || npars != 4) {
        throw GException::xml_invalid_parnum(G_READ, xml,
              "COMPTEL observation requires exactly 4 parameters.");
    }

    // Extract parameters
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

        // Handle EventList
        if (par->attribute("name") == "DRE") {
            drename = par->attribute("file");
            npar[0]++;
        }

        // Handle DRB
        else if (par->attribute("name") == "DRB") {
            drbname = par->attribute("file");
            npar[1]++;
        }

        // Handle DRG
        else if (par->attribute("name") == "DRG") {
            drgname = par->attribute("file");
            npar[2]++;
        }

        // Handle PSF
        else if (par->attribute("name") == "DRX") {
            drxname = par->attribute("file");
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require \"DRE\", \"DRB\", \"DRG\" and \"DRX\""
              " parameters.");
    }

    // Load observation
    load(drename, drbname, drgname, drxname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::xml_invalid_parnum
 *            Invalid number of parameters found in XML element.
 * @exception GException::xml_invalid_parnames
 *            Invalid parameter names found in XML element.
 *
 * Writes information for a COMPTEL observation into an XML element. The
 * expected format of the XML element is
 *
 *  <observation name="..." id="..." instrument="...">
 *    <parameter name="DRE" file="..."/>
 *    <parameter name="DRB" file="..."/>
 *    <parameter name="DRG" file="..."/>
 *    <parameter name="DRX" file="..."/>
 *  </observation>
 *
 * for a binned observation.
 ***************************************************************************/
void GCOMObservation::write(GXmlElement& xml) const
{
    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"DRE\""));
        xml.append(new GXmlElement("parameter name=\"DRB\""));
        xml.append(new GXmlElement("parameter name=\"DRG\""));
        xml.append(new GXmlElement("parameter name=\"DRX\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::xml_invalid_parnum(G_WRITE, xml,
              "COMPTEL observation requires exactly 4 parameters.");
    }

    // Set or update parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = static_cast<GXmlElement*>(xml.element("parameter", i));

        // Handle DRE
        if (par->attribute("name") == "DRE") {
            par->attribute("file", m_drename);
            npar[0]++;
        }

        // Handle DRB
        else if (par->attribute("name") == "DRB") {
            par->attribute("file", m_drbname);
            npar[1]++;
        }

        // Handle DRG
        else if (par->attribute("name") == "DRG") {
            par->attribute("file", m_drgname);
            npar[2]++;
        }

        // Handle DRX
        else if (par->attribute("name") == "DRX") {
            par->attribute("file", m_drxname);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::xml_invalid_parnames(G_WRITE, xml,
              "Require \"DRE\", \"DRB\", \"DRG\" and \"DRX\""
              " parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] drename Event cube name.
 * @param[in] drbname Background cube name.
 * @param[in] drgname Geometry cube name.
 * @param[in] drxname Exposure map name.
 ***************************************************************************/
void GCOMObservation::load(const std::string& drename,
                           const std::string& drbname,
                           const std::string& drgname,
                           const std::string& drxname)
{
    // Load DRE
    load_dre(drename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COM observation information
 ***************************************************************************/
std::string GCOMObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCOMObservation ===");
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
        result.append("\n"+parformat("Pointing")+"undefined");
    }

    // Append response
    if (m_response != NULL) {
        result.append("\n"+response()->print());
    }
    else {
        result.append("\n"+parformat("Response")+"undefined");
    }

    // Append events
    if (m_events != NULL) {
        result.append("\n"+m_events->print());
    }
    else {
        result.append("\n"+parformat("Events")+"undefined");
    }

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
void GCOMObservation::init_members(void)
{
    // Initialise members
    m_instrument = "COM";
    m_drename.clear();
    m_drbname.clear();
    m_drgname.clear();
    m_drxname.clear();
    m_pointing   = NULL;
    m_response   = NULL;
    m_obs_id     = 0;
    m_ontime     = 0.0;
    m_livetime   = 0.0;
    m_deadc      = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs COMPTEL observation.
 ***************************************************************************/
void GCOMObservation::copy_members(const GCOMObservation& obs)
{
    // Clone members. Note that the events are cloned by the base class.
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Copy members
    m_instrument = obs.m_instrument;
    m_drename    = obs.m_drename;
    m_drbname    = obs.m_drbname;
    m_drgname    = obs.m_drgname;
    m_drxname    = obs.m_drxname;
    m_obs_id     = obs.m_obs_id;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMObservation::free_members(void)
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
 * @brief Load event cube data from DRE file
 *
 * @param[in] drename DRE filename.
 ***************************************************************************/
void GCOMObservation::load_dre(const std::string& drename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event cube
    m_events = new GCOMEventCube;

    // Open FITS file
    GFits file(drename);

    // Read event cube
    m_events->read(file);

    // Read observation attributes from primary extension
    GFitsHDU* hdu = file.hdu(0);
    read_attributes(hdu);

    // Close FITS file
    file.close();

    // Store event filename
    m_drename = drename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Reads COM observation attributes from HDU. Mandatory attributes are
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
 * Based on RA_PNT and DEC_PNT, the COM pointing direction is set. Note that
 * DEADC is computed using DEADC=LIVETIME/ONTIME
 *
 * Nothing is done if the HDU pointer is NULL.
 *
 * @todo The actual reader is a minimal reader to accomodate as many
 *       different datasets as possible. Once the COM data format is fixed
 *       the reader should have more mandatory attributes.
 ***************************************************************************/
void GCOMObservation::read_attributes(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read observation information
        m_name   = (hdu->hascard("OBJECT"))   ? hdu->string("OBJECT") : "unknown";
        m_obs_id = (hdu->hascard("OBS_ID"))   ? hdu->real("OBS_ID") : 0;

        // Read ontime, livetime, deadtime correction
        m_ontime   = (hdu->hascard("ONTIME"))   ? hdu->real("ONTIME") : 0.0;
        m_livetime = (hdu->hascard("LIVETIME")) ? hdu->real("LIVETIME") : 0.0;
        m_deadc    = (hdu->hascard("DEADC"))    ? hdu->real("DEADC") : 0.0;

        // Set pointing information
        GSkyDir pnt;
        double  ra_scz  = hdu->real("RA_SCZ");
        double  dec_scz = hdu->real("DEC_SCZ");
        pnt.radec_deg(ra_scz, dec_scz);
        if (m_pointing != NULL) delete m_pointing;
        m_pointing = new GCOMPointing;
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
void GCOMObservation::write_attributes(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Compute some attributes
        double ra_scz  = (m_pointing != NULL) ? m_pointing->dir().ra_deg() : 0.0;
        double dec_scz = (m_pointing != NULL) ? m_pointing->dir().dec_deg() : 0.0;
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
        hdu->card("RA_OBJ",   0.0,  "[deg] Target Right Ascension");
        hdu->card("DEC_OBJ",  0.0, "[deg] Target Declination");
        hdu->card("RA_SCZ",   ra_scz,    "[deg] Pointing Right Ascension");
        hdu->card("DEC_SCZ",  dec_scz,   "[deg] Pointing Declination");
        hdu->card("RADECSYS", "FK5",     "Coordinate system");
        hdu->card("EQUINOX",  2000.0,    "Epoch");
        hdu->card("OBSERVER", "string",  "Observer");

        // Other information
        hdu->card("EUNIT",    "MeV",    "Energy unit");

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
