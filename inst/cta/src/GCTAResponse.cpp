/***************************************************************************
 *                   GCTAResponse.cpp - CTA Response class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
#include "GEvent.hpp"
#include "GSource.hpp"
#include "GObservations.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialComposite.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponse.hpp"

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

/***********************************************************************//**
 * @brief Return instrument response
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response.
 *
 * Returns the instrument response for a given event, source and observation.
 ***************************************************************************/
double GCTAResponse::irf(const GEvent&       event,
                         const GSource&      source,
                         const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Set IRF value attributes
    std::string    name  = obs.id() + "::" + source.name();
    const GEnergy& ereco = event.energy();
    const GEnergy& etrue = source.energy();

    // Signal if spatial model has free parameters
    bool has_free_pars = source.model()->has_free_pars();

    // If the spatial model component has free parameters, or the response
    // cache should not be used, or the cache does not contain the requested
    // IRF value then compute the IRF value for the spatial model.
    if (has_free_pars    ||
        !m_use_irf_cache ||
        !m_irf_cache.contains(name, ereco, etrue, &irf)) {

        // Compute IRF for spatial model
        switch (source.model()->code()) {
            case GMODEL_SPATIAL_POINT_SOURCE:
                irf = irf_ptsrc(event, source, obs);
                break;
            case GMODEL_SPATIAL_RADIAL:
                irf = irf_radial(event, source, obs);
                break;
            case GMODEL_SPATIAL_ELLIPTICAL:
                irf = irf_elliptical(event, source, obs);
                break;
            case GMODEL_SPATIAL_DIFFUSE:
                irf = irf_diffuse(event, source, obs);
                break;
            case GMODEL_SPATIAL_COMPOSITE:
                irf = irf_composite(event, source, obs);
                break;
            default:
                break;
        }

    } // endif: computed spatial model

    // If the spatial model has no free parameters and the response cache
    // should be used then put the IRF value in the response cache.
    if (!has_free_pars && m_use_irf_cache) {
        m_irf_cache.set(name, ereco, etrue, irf);
    }

    // Return IRF value
    return irf;
}


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
    // Initialize members
    m_use_irf_cache = true;
    m_irf_cache.clear();

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
    // Copy members
    m_use_irf_cache = rsp.m_use_irf_cache;
    m_irf_cache     = rsp.m_irf_cache;

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
 * @brief Return instrument response to composite source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to composite source.
 *
 * Returns the instrument response to a specified composite source.
 ***************************************************************************/
double GCTAResponse::irf_composite(const GEvent&       event,
                                   const GSource&      source,
                                   const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to composite model
    const GModelSpatialComposite* model =
        dynamic_cast<const GModelSpatialComposite*>(source.model());

    // Loop over model components
    for (int i = 0; i < model->components(); ++i) {

        // Get pointer to spatial component
        GModelSpatial* spat = const_cast<GModelSpatial*>(model->component(i));

        // Create new GSource object
        GSource src(source.name(), spat, source.energy(), source.time());

        // Compute irf value
        irf += this->irf(event, src, obs) * model->scale(i);

    }

    // Divide by number of model components
    double sum = model->sum_of_scales();
    if (sum > 0.0) {
        irf /= sum;
    }

    // Return IRF value
    return irf;
}
