/***************************************************************************
 *         GMWLObservation.cpp - Multi-wavelength observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GMWLObservation.cpp
 * @brief Multi-wavelength observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GTools.hpp"
#include "GException.hpp"
#include "GMWLObservation.hpp"
#include "GMWLSpectrum.hpp"
#include "GMWLException.hpp"

/* __ Globals ____________________________________________________________ */
const GMWLObservation      g_obs_mwl_seed;
const GObservationRegistry g_obs_mwl_registry(&g_obs_mwl_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GMWLObservation::response(GResponse&)"
#define G_READ                          "GMWLObservation::read(GXmlElement&)"
#define G_WRITE                        "GMWLObservation::write(GXmlElement&)"

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
 * Creates instance of an undefined observation.
 ***************************************************************************/
GMWLObservation::GMWLObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 * @param[in] extno FITS file extension number.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename,
                                 const int& extno) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename, extno);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename File name.
 * @param[in] extname FITS file extension name.
 *
 * Creates instance from file.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const std::string& filename,
                                 const std::string& extname) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename, extname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Observation.
 *
 * Creates instance by copying an existing observation.
 ***************************************************************************/
GMWLObservation::GMWLObservation(const GMWLObservation& obs) : GObservation(obs)
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
 * Destroy instance.
 ***************************************************************************/
GMWLObservation::~GMWLObservation(void)
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
 * @param[in] obs Observation.
 *
 * Copies observation into the instance.
 ***************************************************************************/
GMWLObservation& GMWLObservation::operator= (const GMWLObservation& obs)
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
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GMWLObservation::clear(void)
{
    // Free class members (base and derived classes, derived class first)
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
GMWLObservation* GMWLObservation::clone(void) const
{
    return new GMWLObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GMWLException::bad_response_type
 *            Specified response in not of type GMWLResponse.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GMWLResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GMWLObservation::response(const GResponse& rsp)
{
    // Get pointer on MWL response
    const GMWLResponse* mwlrsp = dynamic_cast<const GMWLResponse*>(&rsp);
    if (mwlrsp == NULL) {
        throw GMWLException::bad_response_type(G_RESPONSE);
    }

    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = mwlrsp->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to response function
 ***************************************************************************/
GMWLResponse* GMWLObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to pointing
 ***************************************************************************/
GMWLPointing* GMWLObservation::pointing(void) const
{
    // Return pointing pointer
    return m_pointing;
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
 * Reads information for a multi-wavelength observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="MWL">
 *       <parameter name="Instrument" value="..."/>
 *       <parameter name="Data" file="..." [extno="..."] [extname="..."]/>
 *     </observation>
 *
 * The extno and extname attributes of the data parameter are optional, and
 * can be used to indicate the extension number or name from which the
 * multi-wavelength data should be loaded. If both are given, the extension
 * number will take precedence and the extension name is ignored.
 ***************************************************************************/
void GMWLObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Determine number of parameter nodes in XML element
    int npars = xml.elements("parameter");

    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || npars != 2) {
        throw GException::xml_invalid_parnum(G_READ, xml,
              "MWL observation requires exactly 2 parameters.");
    }

    // Extract parameters
    int npar[2] = {0, 0};
    for (int i = 0; i < npars; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle Instrument name
        if (par->attribute("name") == "Instrument") {

            // Read instrument name
            m_instrument = par->attribute("value");
            
            // Increment parameter counter
            npar[0]++;
        }

        // Handle data filename
        else if (par->attribute("name") == "Data") {

            // Read filename
            std::string filename = par->attribute("file");
            
            // Read (optional) extension number and name
            std::string extno   = par->attribute("extno");
            std::string extname = par->attribute("extname");
            
            // Load file (this also stores the filename, extno and
            // extname)
            if (gammalib::strip_whitespace(extno).length() > 0) {
                load(filename, gammalib::toint(extno));
            }
            else if (gammalib::strip_whitespace(extname).length() > 0) {
                load(filename, extname);
            }
            else {
                load(filename);
            }
            
            // Increment parameter counter
            npar[1]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1) {
        throw GException::xml_invalid_parnames(G_READ, xml,
              "Require \"Instrument\" and \"Data\" parameters.");
    }

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
 * Writes information for a multi-wavelength observation into an XML element.
 * The format of the XML element is
 *
 *     <observation name="..." id="..." instrument="MWL">
 *       <parameter name="Instrument" value="..."/>
 *       <parameter name="Data" file="..." [extno="..."] [extname="..."]/>
 *     </observation>
 ***************************************************************************/
void GMWLObservation::write(GXmlElement& xml) const
{
    // If XML element has 0 nodes then append 2 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Instrument\""));
        xml.append(GXmlElement("parameter name=\"Data\""));
    }

    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2) {
        throw GException::xml_invalid_parnum(G_WRITE, xml,
              "MWL observation requires exactly 2 parameters.");
    }

    // Set or update parameter attributes
    int npar[] = {0, 0};
    for (int i = 0; i < 2; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle Instrument
        if (par->attribute("name") == "Instrument") {
            npar[0]++;
            par->attribute("value", m_instrument);
        }

        // Handle Data
        else if (par->attribute("name") == "Data") {
            npar[1]++;
            par->attribute("file", m_filename);
            if (gammalib::strip_whitespace(m_extno).length() > 0) {
                par->attribute("extno", m_extno);
            }
            if (gammalib::strip_whitespace(m_extname).length() > 0) {
                par->attribute("extname", m_extname);
            }
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1)
        throw GException::xml_invalid_parnames(G_WRITE, xml,
              "Require \"Instrument\" and \"Data\" parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 ***************************************************************************/
void GMWLObservation::load(const std::string& filename)
{
    // Clear observation
    clear();

    // Allocate spectrum
    GMWLSpectrum* spec = new GMWLSpectrum;
    m_events = spec;

    // Load spectrum
    spec->load(filename);

    // Set attributes
    name("Multi-wavelength observation");
    id(filename);
    m_filename   = filename;
    m_instrument = spec->instrument();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 * @param[in] extno FITS extension number.
 ***************************************************************************/
void GMWLObservation::load(const std::string& filename,
                           const int&         extno)
{
    // Clear observation
    clear();

    // Allocate spectrum
    GMWLSpectrum* spec = new GMWLSpectrum;
    m_events = spec;

    // Load spectrum
    spec->load(filename, extno);

    // Set attributes
    name("Multi-wavelength observation");
    id(filename+"["+gammalib::str(extno)+"]");
    m_filename   = filename;
    m_extno      = gammalib::str(extno);
    m_instrument = spec->instrument();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 * @param[in] extname FITS extension name.
 ***************************************************************************/
void GMWLObservation::load(const std::string& filename,
                           const std::string& extname)
{
    // Clear observation
    clear();

    // Allocate spectrum
    GMWLSpectrum* spec = new GMWLSpectrum;
    m_events = spec;

    // Load spectrum
    spec->load(filename, extname);

    // Set attributes
    name("Multi-wavelength observation");
    id(filename+"["+extname+"]");
    m_filename   = filename;
    m_extname    = extname;
    m_instrument = spec->instrument();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print multi-wavelength information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing multi-wavelength information.
 ***************************************************************************/
std::string GMWLObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GMWLObservation ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Identifier")+id());
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Statistics")+statistics());

        // EXPLICIT: Append events
        if (chatter >=  EXPLICIT) {
            if (m_events != NULL) {
                result.append("\n"+m_events->print(chatter));
            }
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
 *
 * The instrument name is set here to "MWL" so that the registry has an
 * instrument type with that name. This may be later overwritten by a
 * specific instrument.
 ***************************************************************************/
void GMWLObservation::init_members(void)
{
    // Initialise members
    m_instrument = "MWL";
    m_filename.clear();
    m_extno.clear();
    m_extname.clear();
    m_response   = new GMWLResponse;
    m_pointing   = new GMWLPointing;

    // Overwrite base class statistics
    m_statistics = "Gaussian";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 ***************************************************************************/
void GMWLObservation::copy_members(const GMWLObservation& obs)
{
    // Copy members
    m_instrument = obs.m_instrument;
    m_filename   = obs.m_filename;
    m_extno      = obs.m_extno;
    m_extname    = obs.m_extname;
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMWLObservation::free_members(void)
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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
