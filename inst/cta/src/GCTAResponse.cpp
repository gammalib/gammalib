/***************************************************************************
 *                   GCTAResponse.cpp - CTA Response class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse.cpp
 * @brief CTA response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GCTAResponse.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs void CTA response.
 ***************************************************************************/
GCTAResponse::GCTAResponse(void) : GResponse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp CTA response.
 *
 * Constructs CTA response by making a deep copy of an existing object.
 **************************************************************************/
GCTAResponse::GCTAResponse(const GCTAResponse& rsp) : GResponse(rsp)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of CTA response object.
 ***************************************************************************/
GCTAResponse::~GCTAResponse(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] rsp CTA response.
 * @return CTA response.
 *
 * Assigns CTA response object to another CTA response object. The assignment
 * performs a deep copy of all information, hence the original object from
 * which the assignment has been performed can be destroyed after this
 * operation without any loss of information.
 ***************************************************************************/
GCTAResponse& GCTAResponse::operator=(const GCTAResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(rsp);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAResponse::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
void GCTAResponse::copy_members(const GCTAResponse& rsp)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Retrieve CTA observation from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs is not a CTA observations.
 *
 * Dynamically casts generic observation into a CTA observation. If the
 * generic observation is not a CTA observation, an exception is thrown.
 ***************************************************************************/
const GCTAObservation& GCTAResponse::retrieve_obs(const std::string& origin,
                                                  const GObservation& obs) const
{
    // Get pointer on CTA observation
    const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>(&obs);

    // If pointer is not valid then throw an exception
    if (cta == NULL) {
        std::string msg = "Specified observation is not a CTA observation.\n"
                          "Please specify a CTA observation when calling"
                          " this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return *cta;
}


/***********************************************************************//**
 * @brief Retrieve CTA pointing from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * Extract CTA pointing from a CTA observation.
 ***************************************************************************/
const GCTAPointing& GCTAResponse::retrieve_pnt(const std::string& origin,
                                               const GObservation& obs) const
{
    // Retrieve CTA observation and pointing
    const GCTAObservation& cta = retrieve_obs(origin, obs);

    // Return CTA pointing
    return (cta.pointing());
}


/***********************************************************************//**
 * @brief Retrieve CTA ROI from generic observation
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] obs Generic observation.
 *
 * @exception GException::invalid_argument
 *            Observation @p obs does not contain a CTA event list.
 *
 * Extract CTA Region of Interest from a CTA observation.
 ***************************************************************************/
const GCTARoi& GCTAResponse::retrieve_roi(const std::string& origin,
                                          const GObservation& obs) const
{
    // Retrieve CTA observation
    const GCTAObservation& cta = retrieve_obs(origin, obs);

    // Get pointer on CTA events list
    const GCTAEventList* events = dynamic_cast<const GCTAEventList*>(cta.events());

    // If pointer is not valid then throw an exception
    if (events == NULL) {
        std::string msg = "Specified observation does not contain a CTA event"
                          " list.\n"
                          "Please specify a CTA observation containing a CTA"
                          " event list when calling this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return CTA ROI
    return (events->roi());
}


/***********************************************************************//**
 * @brief Retrieve CTA instrument direction from generic event
 *
 * @param[in] origin Method asking for pointer retrieval.
 * @param[in] event Generic event.
 *
 * @exception GException::invalid_argument
 *            @p event does not contain a CTA instrument direction.
 *
 * Extract CTA Instrument Direction from an event.
 ***************************************************************************/
const GCTAInstDir& GCTAResponse::retrieve_dir(const std::string& origin,
                                              const GEvent&      event) const
{
    // Get pointer on CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&(event.dir()));

    // If pointer is not valid then throw an exception
    if (dir == NULL) {
        std::string msg = "Specified event is not a CTA event.\n"
                          "Please specify a CTA event when calling this method.";
        throw GException::invalid_argument(origin, msg);
    }

    // Return reference
    return *dir;
}
