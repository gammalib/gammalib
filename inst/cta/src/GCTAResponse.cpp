/***************************************************************************
 *                   GCTAResponse.cpp - CTA Response class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
#include "GCTAEventCube.hpp"               // Kludge
#include "GCTAEventBin.hpp"                // Kludge
#include "GCTAResponse.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF           "GCTAResponse::irf(GEvent&, GSource&, GObservation&)"

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
 * @brief Return instrument response integrated over the spatial model
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to a spatial model.
 *
 * Returns the instrument response for a given event, source and observation
 * integrated over the spatial model component. The method computes
 *
 * \f[
 *    {\tt irf}(p', E', t') = \int_p M_{\rm S}(p | E, t) \,
 *                                   R(p', E', t' | p, E, t) \, d\,p
 * \f]
 *
 * where
 * * \f$M_{\rm S}(p | E, t)\f$ is the spatial model component,
 * * \f$R(p', E', t' | p, E, t)\f$ is the Instrument Response Function (IRF),
 * * \f$p'\f$ is the measured instrument direction,
 * * \f$E'\f$ is the measured or reconstructed energy,
 * * \f$t'\f$ is the measured arrival time,
 * * \f$p\f$ is the true photon arrival direction,
 * * \f$E\f$ is the true photon energy, and
 * * \f$t\f$ is the true trigger time.
 *
 * The integration is done over all relevant true sky directions \f$p\f$.
 *
 * Depending on the type of the source model the method branches to the
 * following methods to perform the actual computations
 *
 *      irf_ptsrc() - for the handling of a point source
 *      irf_radial() - for radial models
 *      irf_elliptical() - for elliptical models
 *      irf_diffuse() - for diffuse models
 *      irf_composite() - for composite models
 *
 * The method implements a caching mechanism for spatial models that have all
 * parameters fixed. For those models the instrument response for a given
 * event and observation is only computed once and then stored in an internal
 * cache from which it is fetched back in case that the method is called
 * again for the same event and observation.
 ***************************************************************************/
double GCTAResponse::irf_spatial(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Set IRF value attributes
    std::string     name  = obs.id() + "::" + source.name();
    const GInstDir& dir   = event.dir();
    const GEnergy&  ereco = event.energy();
    const GEnergy&  etrue = source.energy();

    // Signal if spatial model has free parameters
    bool has_free_pars = source.model()->has_free_pars();

    // Kludge: if the response cache should be used, the event is a bin, the
    // model has no free parameters, and there is no energy dispersion then
    // use the special cache for binned or stacked analysis
    bool use_kludge = false;
    if (m_use_irf_cache && event.is_bin() && !has_free_pars && !use_edisp()) {

        // Get reference to CTA event bin
        const GCTAEventBin& bin = static_cast<const GCTAEventBin&>(event);

        // If pixel and energy indices are valid then proceed with the kludge
        if (bin.ipix() >= 0 && bin.ieng() >=0) {

            // Signal that kludge will be used
            use_kludge = true;

            // Build unique cache name
            std::string name  = obs.id() + "::" + source.name();

            // Get index in cache, returns -1 if name is not found in cache
            int index = -1;
            for (int i = 0; i < m_cache_names.size(); ++i) {
                if (m_cache_names[i] == name) {
                    index = i;
                    break;
                }
            }

            // If index was not found then allocate a new cache map
            if (index == -1) {

                // Get pointer to event cube
                const GCTAEventCube* cube =
                      static_cast<const GCTAEventCube*>(obs.events());

                // Allocate cache
                GNdarray cache(cube->nx()*cube->ny(), cube->ebins());

                // Initialise all cache values with -1 (not set)
                for (int i = 0; i < cache.size(); ++i) {
                    cache(i) = -1.0;
                }

                // Insert cache
                m_cache_names.push_back(name);
                m_cache_values.push_back(cache);

                // Set index
                index = m_cache_names.size()-1;

            } // endif: allocated new cache

            // Get cache value
            irf = m_cache_values[index](bin.ipix(), bin.ieng());

            // If cache value is not valid then compute IRF
            if (irf < 0.0) {

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

                // Set cache value
                m_cache_values[index](bin.ipix(), bin.ieng()) = irf;

            } // endif: computed IRF

        } // endif: event bin was part of cube

    } // endif: kludge

    // If kludge was not used then use the release 1.7 response cache
    if (!use_kludge) {

        // If the spatial model component has free parameters, or the response
        // cache should not be used, or the cache does not contain the requested
        // IRF value then compute the IRF value for the spatial model.
        if (has_free_pars    ||
            !m_use_irf_cache ||
            !m_irf_cache.contains(name, dir, ereco, etrue, &irf)) {

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
            m_irf_cache.set(name, dir, ereco, etrue, irf);
        }

    } // endelse: used release 1.7 response cache

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
    // Kludge: Initialise cube response cache
    m_cache_names.clear();
    m_cache_values.clear();

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
    // Kludge: Copy cube response cache
    m_cache_names  = rsp.m_cache_names;
    m_cache_values = rsp.m_cache_values;

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
