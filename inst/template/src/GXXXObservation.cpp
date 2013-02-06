/***************************************************************************
 *              GXXXObservation.cpp  -  XXX Observation class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @brief XXX observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GObservationRegistry.hpp"
#include "GXXXException.hpp"
#include "GXXXObservation.hpp"
#include "GXXXEventCube.hpp"

/* __ Globals ____________________________________________________________ */
const GXXXObservation      g_obs_spi_seed;
const GObservationRegistry g_obs_spi_registry(&g_obs_spi_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GXXXObservation::response(GResponse&)"

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
 * @param[in] obs XXX observation.
 *
 * Creates class instance by copying an existing XXX observation.
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
 * @param[in] obs XXX observation.
 * @return XXX observation.
 *
 * Assign XXX observation to this object.
 ***************************************************************************/
GXXXObservation& GXXXObservation::operator= (const GXXXObservation& obs)
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
 * @brief Clone instance
 *
 * @return Pointer to deep copy of XXX observation.
 ***************************************************************************/
GXXXObservation* GXXXObservation::clone(void) const
{
    return new GXXXObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GException::rsp_invalid_type
 *            Specified response is not of type GXXXResponse.
 *
 * Sets the response function for the observation. The argument has to be of
 * type GXXXResponse, otherwise an exception is thrown.
 ***************************************************************************/
void GXXXObservation::response(const GResponse& rsp)
{
    // Get pointer on XXX response
    const GXXXResponse* spirsp = dynamic_cast<const GXXXResponse*>(&rsp);
    if (spirsp == NULL) {
        throw GException::rsp_invalid_type(G_RESPONSE,
              typeid(&rsp).name(), "Expected type \"GXXXResponse\".");
    }

    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Clone response function
    m_response = spirsp->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to response function
 *
 * @return Pointer to response function.
 ***************************************************************************/
GXXXResponse* GXXXObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to pointing direction
 *
 * @return Pointer to pointing direction.
 ***************************************************************************/
GXXXPointing* GXXXObservation::pointing(void) const
{
    // Return pointing pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Read observation from XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GXXXObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // TODO: Implement method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GXXXObservation::write(GXmlElement& xml) const
{
    // TODO: Implement method

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print observation information
 *
 * @return String containing observation information.
 ***************************************************************************/
std::string GXXXObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GXXXObservation ===");

    // Append observation information
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Identifier")+id());
    result.append("\n"+parformat("Instrument")+instrument());
    result.append("\n"+parformat("Statistics")+statistics());
    result.append("\n"+parformat("Ontime")+str(ontime())+" sec");
    result.append("\n"+parformat("Livetime")+str(livetime())+" sec");
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
void GXXXObservation::init_members(void)
{
    // Initialise members
    m_instrument = "XXX";
    m_ontime     = 0.0;
    m_livetime   = 0.0;
    m_deadc      = 0.0;
    m_pointing   = NULL;
    m_response   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs XXX observation.
 ***************************************************************************/
void GXXXObservation::copy_members(const GXXXObservation& obs)
{
    // Clone members. Note that the events are cloned by the base class.
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Copy members
    m_instrument = obs.m_instrument;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXXXObservation::free_members(void)
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
