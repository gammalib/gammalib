/***************************************************************************
 *             GLATObservation.cpp - Fermi/LAT observation class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GLATObservation.cpp
 * @brief Fermi/LAT observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GObservationRegistry.hpp"
#include "GLATObservation.hpp"
#include "GLATEventList.hpp"
#include "GLATEventCube.hpp"
#include "GLATRoi.hpp"
#include "GLATException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GEnergy.hpp"

/* __ Globals ____________________________________________________________ */
const GLATObservation      g_obs_lat_seed;
const GObservationRegistry g_obs_lat_registry(&g_obs_lat_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GLATObservation::response(GResponse&)"
#define G_READ                          "GLATObservation::read(GXmlElement&)"
#define G_WRITE                        "GLATObservation::write(GXmlElement&)"

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
GLATObservation::GLATObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs LAT observation.
 ***************************************************************************/
GLATObservation::GLATObservation(const GLATObservation& obs) : GObservation(obs)
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
GLATObservation::~GLATObservation(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] obs Fermi/LAT observation.
 * @return Fermi/LAT observation.
 ***************************************************************************/
GLATObservation& GLATObservation::operator=(const GLATObservation& obs)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear Fermi/LAT observation
 ***************************************************************************/
void GLATObservation::clear(void)
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
 * @brief Clone Fermi/LAT observation
 *
 * @return Pointer to deep copy of Fermi/LAT observation.
 ***************************************************************************/
GLATObservation* GLATObservation::clone(void) const
{
    return new GLATObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GLATException::bad_response_type
 *            Specified response in not of type GLATResponse.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GLATResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GLATObservation::response(const GResponse& rsp)
{
    // Get pointer on LAT response
    const GLATResponse* latrsp = dynamic_cast<const GLATResponse*>(&rsp);
    if (latrsp == NULL) {
        throw GLATException::bad_response_type(G_RESPONSE);
    }

    // Copy response function
    m_response = *latrsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] irfname Name of instrument response function.
 * @param[in] caldb Optional path to calibration database.
 *
 * Set the LAT response function using the IRF name and the path to the
 * calibration database. The IFR name has to be one of
 * name
 * name::front
 * name::back
 * where name is the response name (e.g. P6_v3). Note that the name is case
 * sensitive for the moment.
 ***************************************************************************/
void GLATObservation::response(const std::string& irfname,
                               const std::string& caldb)
{
    // Clear LAT response function
    m_response.clear();

    // Set calibration database
    m_response.caldb(caldb);

    // Load instrument response function
    m_response.load(irfname);

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
 * Reads information for a LAT observation from an XML element. The expected
 * format of the XML element is
 *
 *     <observation name="..." id="..." instrument="LAT">
 *       <parameter name="FT1" file="..."/>
 *       <parameter name="FT2" file="..."/>
 *       <parameter name="LiveTimeCube" file="..."/>
 *       <parameter name="IRF" file="..."/>
 *     </observation>
 *
 * for an unbinned observation and
 *
 *     <observation name="..." id="..." instrument="LAT">
 *       <parameter name="CountsMap" file="..."/>
 *       <parameter name="ExposureMap" file="..."/>
 *       <parameter name="LiveTimeCube" file="..."/>
 *       <parameter name="IRF" value="..."/>
 *     </observation>
 *
 * for a binned observation.
 ***************************************************************************/
void GLATObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Initialise attributes
    std::string ft1file = "";
    std::string ft2file = "";
    std::string ltfile  = "";
    std::string cntfile = "";
    std::string expfile = "";
    std::string irfname = "";

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || npars != 4) {
        throw GException::xml_invalid_parnum(G_READ, xml,
              "LAT observation requires exactly 4 parameters.");
    }

    // Extract parameters
    int npar1[] = {0, 0, 0, 0};
    int npar2[] = {0, 0, 0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle Unbinned format
        if (par->attribute("name") == "FT1") {
            ft1file = par->attribute("file");
            npar1[0]++;
        }
        else if (par->attribute("name") == "FT2") {
            ft2file = par->attribute("file");
            npar1[1]++;
        }

        // Handle Binned format
        else if (par->attribute("name") == "CountsMap") {
            cntfile = par->attribute("file");
            npar2[0]++;
        }
        else if (par->attribute("name") == "ExposureMap") {
            expfile = par->attribute("file");
            npar2[1]++;
        }

        // Handle common parameters
        else if (par->attribute("name") == "LiveTimeCube") {
            ltfile = par->attribute("file");
            npar1[2]++;
            npar2[2]++;
        }
        else if (par->attribute("name") == "IRF") {
            irfname = par->attribute("value");
            npar1[3]++;
            npar2[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    bool unbin_ok = (npar1[0] == 1 && npar1[1] == 1 && npar1[2] == 1 && npar1[3] == 1);
    bool bin_ok   = (npar2[0] == 1 && npar2[1] == 1 && npar2[2] == 1 && npar2[3] == 1);
    if (!bin_ok && !unbin_ok) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require either \"FT1\", \"FT2\", \"LiveTimeCube\", and \"IRF\""
              " or \"CountsMap\", \"ExposureMap\", \"LiveTimeCube\", and"
              " \"IRF\" parameters.");
    }

    // Load data
    if (unbin_ok) {
        load_unbinned(ft1file, ft2file, ltfile);
    }
    else {
        load_binned(cntfile, expfile, ltfile);
    }
    
    // Set response function
    response(irfname);

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
 * Writes information for a LAT observation into an XML element. The expected
 * format of the XML element is
 *
 *     <observation name="..." id="..." instrument="LAT">
 *       <parameter name="FT1" file="..."/>
 *       <parameter name="FT2" file="..."/>
 *       <parameter name="LiveTimeCube" file="..."/>
 *       <parameter name="IRF" file="..."/>
 *     </observation>
 *
 * for an unbinned observation and
 *
 *     <observation name="..." id="..." instrument="LAT">
 *       <parameter name="CountsMap" file="..."/>
 *       <parameter name="ExposureMap" file="..."/>
 *       <parameter name="LiveTimeCube" file="..."/>
 *       <parameter name="IRF" value="..."/>
 *     </observation>
 *
 * for a binned observation.
 *
 * @todo We should create a special exception that informs that there is
 *       neither a valid LAT event list nor a valid LAT counts map in this
 *       observations.
 ***************************************************************************/
void GLATObservation::write(GXmlElement& xml) const
{
    // Determine if we deal with a binned or unbinned observation
    const GLATEventList* list = dynamic_cast<const GLATEventList*>(m_events);
    const GLATEventCube* cube = dynamic_cast<const GLATEventCube*>(m_events);
    if (list == NULL && cube == NULL) {
        throw GException::no_list(G_WRITE);
    }

    // Set event list flag
    bool is_list = (list != NULL);

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        if (is_list) {
            xml.append(GXmlElement("parameter name=\"FT1\""));
            xml.append(GXmlElement("parameter name=\"FT2\""));
            xml.append(GXmlElement("parameter name=\"LiveTimeCube\""));
            xml.append(GXmlElement("parameter name=\"IRF\""));
        }
        else {
            xml.append(GXmlElement("parameter name=\"CountsMap\""));
            xml.append(GXmlElement("parameter name=\"ExposureMap\""));
            xml.append(GXmlElement("parameter name=\"LiveTimeCube\""));
            xml.append(GXmlElement("parameter name=\"IRF\""));
        }
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::xml_invalid_parnum(G_WRITE, xml,
              "LAT observation requires exactly 4 parameters.");
    }

    // Set or update parameter attributes
    int npar[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle FT1
        if (par->attribute("name") == "FT1") {
            par->attribute("file", m_ft1file);
            npar[0]++;
        }

        // Handle CountsMap
        else if (par->attribute("name") == "CountsMap") {
            par->attribute("file", m_cntfile);
            npar[0]++;
        }

        // Handle FT2
        else if (par->attribute("name") == "FT2") {
            par->attribute("file", m_ft2file);
            npar[1]++;
        }

        // Handle ExposureMap
        else if (par->attribute("name") == "ExposureMap") {
            par->attribute("file", m_expfile);
            npar[1]++;
        }

        // Handle LiveTimeCube
        else if (par->attribute("name") == "LiveTimeCube") {
            par->attribute("file", m_ltfile);
            npar[2]++;
        }

        // Handle IRF
        else if (par->attribute("name") == "IRF") {
            std::string irfname = "";
            irfname = m_response.rspname();
            par->attribute("value", irfname);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require either \"FT1\", \"FT2\", \"LiveTimeCube\", and \"IRF\""
              " or \"CountsMap\", \"ExposureMap\", \"LiveTimeCube\", and"
              " \"IRF\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print LAT observation information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing LAT observation information.
 ***************************************************************************/
std::string GLATObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATObservation ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Identifier")+id());
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Statistics")+statistics());
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" s");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" s");

        // Append detailed information
        GChatter reduced_chatter = gammalib::reduce(chatter);
        if (reduced_chatter > SILENT) {

            // Append response
            result.append("\n"+m_response.print(reduced_chatter));

            // Append livetime cube
            if (m_ltcube != NULL) {
                result.append("\n"+m_ltcube->print(reduced_chatter));
            }
            else {
                result.append("\n"+gammalib::parformat("LAT livetime cube"));
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


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] ft1name FT1 FITS filename.
 * @param[in] ft2name FT2 FITS filename.
 * @param[in] ltcube_name Livetime cube FITS filename
 *
 * @todo So far nothing is done with the ft2 file and the ltcube file.
 *       Loading of the relevant information needs to be implemented.
 ***************************************************************************/
void GLATObservation::load_unbinned(const std::string& ft1name,
                                    const std::string& ft2name,
                                    const std::string& ltcube_name)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event list
    GLATEventList* events = new GLATEventList;

    // Assign event list as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(ft1name);

    // Read event list
    events->read(file);

    // Read observation attributes from EVENTS extension
    //GFitsHDU* hdu = file.hdu("EVENTS");
    //read_attributes(hdu);

    // Close FITS file
    file.close();

    // Optionally allocate and load livetime cube
    if (ltcube_name.length() > 0) {
        m_ltcube = new GLATLtCube;
        m_ltcube->load(ltcube_name);
    }

    // Store filenames
    m_ft1file = ft1name;
    m_ft2file = ft2name;
    m_ltfile  = ltcube_name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] cntmap_name Counts map or Source map FITS filename
 * @param[in] expmap_name Binned explosure map FITS filename
 * @param[in] ltcube_name Livetime cube FITS filename
 *
 * @todo So far nothing is done with the expmap file.
 *       Approriate loading needs to be implemented.
 ***************************************************************************/
void GLATObservation::load_binned(const std::string& cntmap_name,
                                  const std::string& expmap_name,
                                  const std::string& ltcube_name)
{
    // Delete old events and livetime cube.  We do not call clear() here
    // since we want to preserve any existing response function.
    if (m_events != NULL) delete m_events;
    if (m_ltcube != NULL) delete m_ltcube;
    m_events = NULL;
    m_ltcube = NULL;

    // Allocate event cube
    GLATEventCube* events = new GLATEventCube;

    // Assign event cube as the observation's event container
    m_events = events;

    // Load event list
    events->load(cntmap_name);

    // Optionally allocate and load livetime cube
    if (ltcube_name.length() > 0) {
        m_ltcube = new GLATLtCube;
        m_ltcube->load(ltcube_name);
    }

    // Store filenames
    m_cntfile = cntmap_name;
    m_expfile = expmap_name;
    m_ltfile  = ltcube_name;

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
void GLATObservation::init_members(void)
{
    // Initialise members
    m_ft1file.clear();
    m_ft2file.clear();
    m_ltfile.clear();
    m_cntfile.clear();
    m_expfile.clear();
    m_response.clear();
    m_ltcube = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs LAT observation.
 ***************************************************************************/
void GLATObservation::copy_members(const GLATObservation& obs)
{
    // Copy members
    m_ft1file  = obs.m_ft1file;
    m_ft2file  = obs.m_ft2file;
    m_ltfile   = obs.m_ltfile;
    m_cntfile  = obs.m_cntfile;
    m_expfile  = obs.m_expfile;
    m_response = obs.m_response;
    
    // Clone members
    if (obs.m_ltcube != NULL) m_ltcube = obs.m_ltcube->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATObservation::free_members(void)
{
    // Free memory
    if (m_ltcube  != NULL) delete m_ltcube;

    // Mark memory as free
    m_ltcube = NULL;

    // Return
    return;
}
