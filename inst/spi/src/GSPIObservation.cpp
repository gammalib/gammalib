/***************************************************************************
 *           GSPIObservation.cpp - INTEGRAL/SPI observation class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIObservation.cpp
 * @brief INTEGRAL/SPI observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo>
#include "GException.hpp"
#include "GObservationRegistry.hpp"
#include "GSPIObservation.hpp"
#include "GSPITools.hpp"
#include "GSPIEventCube.hpp"

/* __ Globals ____________________________________________________________ */
const GSPIObservation      g_obs_spi_seed;
const GObservationRegistry g_obs_spi_registry(&g_obs_spi_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GSPIObservation::response(GResponse&)"
#define G_READ                          "GSPIObservation::read(GXmlElement&)"
#define G_READ_FITS                           "GSPIObservation::read(GFits&)"
#define G_WRITE                        "GSPIObservation::write(GXmlElement&)"

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
 * Creates an empty INTEGRAL/SPI observation.
 ***************************************************************************/
GSPIObservation::GSPIObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs an INTEGRAL/SPI observation from the information that is found
 * in the XML element.
 ***************************************************************************/
GSPIObservation::GSPIObservation(const GXmlElement& xml) : GObservation()
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation Group constructor
 *
 * @param[in] filename Observation Group FITS file name.
 *
 * Constructs an INTEGRAL/SPI observation from an Observation Group.
 ***************************************************************************/
GSPIObservation::GSPIObservation(const GFilename& filename) : GObservation()
{
    // Initialise members
    init_members();

    // Load data
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs INTEGRAL/SPI observation.
 *
 * Creates INTEGRAL/SPI observation by copying from an existing INTEGRAL/SPI
 * observation.
 ***************************************************************************/
GSPIObservation::GSPIObservation(const GSPIObservation& obs) : GObservation(obs)
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
GSPIObservation::~GSPIObservation(void)
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
 * @param[in] obs INTEGRAL/SPI observation.
 * @return INTEGRAL/SPI observation.
 *
 * Assign INTEGRAL/SPI observation to this object. The assignment performs
 * a deep copy of all information, hence the original object from which the
 * assignment was made can be destroyed after this operation without any loss
 * of information.
 ***************************************************************************/
GSPIObservation& GSPIObservation::operator=(const GSPIObservation& obs)
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
 * @brief Clear INTEGRAL/SPI observation
 *
 * Clears INTEGRAL/SPI observation by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GSPIObservation::clear(void)
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
 * @brief Clone INTEGRAL/SPI observation
 *
 * @return Pointer to deep copy of INTEGRAL/SPI observation.
 ***************************************************************************/
GSPIObservation* GSPIObservation::clone(void) const
{
    return new GSPIObservation(*this);
}


/***********************************************************************//**
 * @brief Set INTEGRAL/SPI response function
 *
 * @param[in] rsp INTEGRAL/SPI response function.
 *
 * @exception GException::invalid_argument
 *            Specified response is not a INTEGRAL/SPI response.
 *
 * Sets the response function for the observation.
 ***************************************************************************/
void GSPIObservation::response(const GResponse& rsp)
{
    // Get pointer on INTEGRAL/SPI response
    const GSPIResponse* spirsp = dynamic_cast<const GSPIResponse*>(&rsp);
    if (spirsp == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Response of type \""+cls+"\" is "
                          "not a INTEGRAL/SPI response. Please specify a "
                          "INTEGRAL/SPI response as argument.";
        throw GException::invalid_argument(G_RESPONSE, msg);
    }

    // Clone response function
    m_response = *spirsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read INTEGRAL/SPI observation from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for a INTEGRAL/SPI observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="Crab" id="0044" instrument="SPI">
 *       <parameter name="ObservationGroup" file="og_spi.fits"/>
 *     </observation>
 *
 ***************************************************************************/
void GSPIObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Get parameters
    std::string filename = gammalib::xml_get_attr(G_READ, xml, "ObservationGroup", "file");

    // Expand file names
    filename = gammalib::xml_file_expand(xml, filename);

    // Load observation
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write INTEGRAL/SPI observation into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a INTEGRAL/SPI observation into an XML element. The
 * format of the XML element is
 *
 *     <observation name="Crab" id="0044" instrument="SPI">
 *       <parameter name="ObservationGroup" file="og_spi.fits"/>
 *     </observation>
 *
 ***************************************************************************/
void GSPIObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // Set ObservationGroup parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "ObservationGroup");
    par->attribute("file", gammalib::xml_file_reduce(xml, m_filename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Observation Group from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads Observation Group from a FITS file.
 ***************************************************************************/
void GSPIObservation::read(const GFits& fits)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;

    // Allocate new event cube
    GSPIEventCube* events = new GSPIEventCube;

    // Assign event cube as the observation's event container
    m_events = events;

    // Read in Observation Group
    events->read(fits);

    // Store ontime and livetime and compute deadtime correction
    m_ontime   = events->ontime();
    m_livetime = events->livetime();
    m_deadc    = (m_ontime > 0.0) ? m_livetime/m_ontime : 1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load Observation Group
 *
 * @param[in] filename Observation Group FITS file name.
 *
 * Loads data from an Observation Group FITS file into an INTEGRAL/SPI
 * observation.
 ***************************************************************************/
void GSPIObservation::load(const GFilename& filename)
{
    #pragma omp critical(GSPIObservation_load)
    {
        // Store event filename
        m_filename = filename;

        // Open FITS file
        GFits fits(filename);

        // Read data
        read(fits);

        // Close FITS file
        fits.close();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print observation information
 *
 * @param[in] chatter Chattiness.
 * @return String containing INTEGRAL/SPI observation information.
 ***************************************************************************/
std::string GSPIObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSPIObservation ===");

        // Append standard information for observation
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Identifier")+id());
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Statistic")+statistic());
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" sec");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" sec");
        result.append("\n"+gammalib::parformat("Deadtime correction"));
        result.append(gammalib::str(m_deadc));

        // Append detailed information
        GChatter reduced_chatter = gammalib::reduce(chatter);
        if (reduced_chatter > SILENT) {

            // Append response
            result.append("\n"+m_response.print(reduced_chatter));

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
void GSPIObservation::init_members(void)
{
    // Initialise members
    m_instrument = "SPI";
    m_filename.clear();
    m_response.clear();
    m_ontime   = 0.0;
    m_livetime = 0.0;
    m_deadc    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs INTEGRAL/SPI observation.
 ***************************************************************************/
void GSPIObservation::copy_members(const GSPIObservation& obs)
{
    // Copy members
    m_instrument = obs.m_instrument;
    m_filename   = obs.m_filename;
    m_response   = obs.m_response;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIObservation::free_members(void)
{
    // Return
    return;
}
