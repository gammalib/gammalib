/***************************************************************************
 *         GMWLObservation.cpp - Multi-wavelength observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
#include <typeinfo>
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GObservationRegistry.hpp"
#include "GException.hpp"
#include "GMWLObservation.hpp"
#include "GMWLSpectrum.hpp"

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
GMWLObservation::GMWLObservation(const GFilename& filename) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(filename);

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
 * @exception GException::invalid_argument
 *            Response @p rsp in not a MWL response.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GMWLResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GMWLObservation::response(const GResponse& rsp)
{
    // Get pointer on MWL response
    const GMWLResponse* mwlrsp = dynamic_cast<const GMWLResponse*>(&rsp);

    // If pointer is not valid then throw an exception
    if (mwlrsp == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Invalid response type \""+cls+"\" specified. "
                          "Please specify a \"GMWLResponse\" instance as "
                          "argument.";
        throw GException::invalid_argument(G_RESPONSE, msg);
    }

    // Copy response function
    m_response = *mwlrsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for a multi-wavelength observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="MWL">
 *       <parameter name="Instrument" value="..."/>
 *       <parameter name="Data" file="..."/>
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

    // Get parameters
    m_instrument = gammalib::xml_get_attr(G_READ, xml, "Instrument", "value");
    std::string filename = gammalib::xml_get_attr(G_READ, xml, "Data", "file");

    // Expand file names
    filename = gammalib::xml_file_expand(xml, filename);

    // Load file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a multi-wavelength observation into an XML element.
 * The format of the XML element is
 *
 *     <observation name="..." id="..." instrument="MWL">
 *       <parameter name="Instrument" value="..."/>
 *       <parameter name="Data" file="..."/>
 *     </observation>
 ***************************************************************************/
void GMWLObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // Set Instrument parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Instrument");
    par->attribute("value", gammalib::xml_file_reduce(xml, m_instrument));

    // Set Data parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "Data");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load observation
 *
 * @param[in] filename File name.
 ***************************************************************************/
void GMWLObservation::load(const GFilename& filename)
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
        result.append("\n"+gammalib::parformat("Statistic")+statistic());

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
    m_response.clear();

    // Overwrite base class statistic
    m_statistic = "Gaussian";

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
    m_response   = obs.m_response;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMWLObservation::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
