/***************************************************************************
 *           GXXXObservation.cpp - [INSTRUMENT] observation class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXObservation.cpp
 * @brief [INSTRUMENT] observation class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo> 
#include "GException.hpp"
#include "GObservationRegistry.hpp"
#include "GXXXObservation.hpp"

/* __ Globals ____________________________________________________________ */
const GXXXObservation      g_obs_xxx_seed;
const GObservationRegistry g_obs_xxx_registry(&g_obs_xxx_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GXXXObservation::response(GResponse&)"
#define G_READ                          "GXXXObservation::read(GXmlElement&)"
#define G_WRITE                        "GXXXObservation::write(GXmlElement&)"

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
 * Creates an empty [INSTRUMENT] observation.
 ***************************************************************************/
GXXXObservation::GXXXObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs [INSTRUMENT] observation.
 *
 * Creates [INSTRUMENT] observation by copying from an existing [INSTRUMENT]
 * observation.
 ***************************************************************************/
GXXXObservation::GXXXObservation(const GXXXObservation& obs) : GObservation(obs)
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
GXXXObservation::~GXXXObservation(void)
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
 * @param[in] obs [INSTRUMENT] observation.
 * @return [INSTRUMENT] observation.
 *
 * Assign [INSTRUMENT] observation to this object. The assignment performs
 * a deep copy of all information, hence the original object from which the
 * assignment was made can be destroyed after this operation without any loss
 * of information.
 ***************************************************************************/
GXXXObservation& GXXXObservation::operator=(const GXXXObservation& obs)
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
 * @brief Clear [INSTRUMENT] observation
 *
 * Clears [INSTRUMENT] observation by resetting all class members to an
 * initial state. Any information that was present before will be lost.
 ***************************************************************************/
void GXXXObservation::clear(void)
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
 * @brief Clone [INSTRUMENT] observation
 *
 * @return Pointer to deep copy of [INSTRUMENT] observation.
 ***************************************************************************/
GXXXObservation* GXXXObservation::clone(void) const
{
    return new GXXXObservation(*this);
}


/***********************************************************************//**
 * @brief Set [INSTRUMENT] response function
 *
 * @param[in] rsp [INSTRUMENT] response function.
 *
 * @exception GException::invalid_argument
 *            Specified response is not a [INSTRUMENT] response.
 *
 * Sets the response function for the observation.
 ***************************************************************************/
void GXXXObservation::response(const GResponse& rsp)
{
    // Get pointer on [INSTRUMENT] response
    const GXXXResponse* xxxrsp = dynamic_cast<const GXXXResponse*>(&rsp);
    if (xxxrsp == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Response of type \""+cls+"\" is "
                          "not a [INSTRUMENT] response. Please specify a "
                          "[INSTRUMENT] response as argument.";
        throw GException::invalid_argument(G_RESPONSE, msg);
    }

    // Clone response function
    m_response = *xxxrsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read [INSTRUMENT] observation from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for a [INSTRUMENT] observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="Crab" id="00001" instrument="XXX">
 *       ...
 *       @todo Define XML format for a [INSTRUMENT] observation
 *       ...
 *     </observation>
 *
 ***************************************************************************/
void GXXXObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Get parameters
    // TODO: add relevant parameters. An example is given below
    //std::string drename = gammalib::xml_get_attr(G_READ, xml, "DRE", "file");

    // Expand file names
    // TODO: the following example shows how to expand environment variables
    //drename = gammalib::xml_file_expand(xml, drename);

    // Load observation
    // TODO: you may load the observation here. An example below.
    //load(drename, drbname, drgname, drxname);

    // Load IAQ
    // TODO: you may load the response here.  An example below.
    //response(GCaldb("cgro", "comptel"), iaqname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write [INSTRUMENT] observation into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a [INSTRUMENT] observation into an XML element. The
 * format of the XML element is
 *
 *     <observation name="Crab" id="00001" instrument="XXX">
 *       ...
 *       @todo Define XML format for a [INSTRUMENT] observation
 *       ...
 *     </observation>
 *
 ***************************************************************************/
void GXXXObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // TODO: Write out all parameters. An example is shown below.
    //par = gammalib::xml_need_par(G_WRITE, xml, "DRE");
    //par->attribute("file", gammalib::xml_file_reduce(xml, m_drename));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print observation information
 *
 * @param[in] chatter Chattiness.
 * @return String containing [INSTRUMENT] observation information.
 ***************************************************************************/
std::string GXXXObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GXXXObservation ===");

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

        // Append additional information
        // TODO: Add code to append any additional information that might
        // be relevant.

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
void GXXXObservation::init_members(void)
{
    // Initialise members
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
 * @param[in] obs [INSTRUMENT] observation.
 ***************************************************************************/
void GXXXObservation::copy_members(const GXXXObservation& obs)
{
    // Copy members
    m_response = obs.m_response;
    m_ontime   = obs.m_ontime;
    m_livetime = obs.m_livetime;
    m_deadc    = obs.m_deadc;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXObservation::free_members(void)
{
    // Return
    return;
}
