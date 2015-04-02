/***************************************************************************
 *                GResponse.cpp - Abstract response base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2015 by Juergen Knoedlseder                         *
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
 * @file GResponse.hpp
 * @brief Abstract response base class interface definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <unistd.h>           // access() function
#include "GTools.hpp"
#include "GMath.hpp"
#include "GException.hpp"
#include "GResponse.hpp"
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GSource.hpp"        // will become obsolete
#include "GEbounds.hpp"       // will become obsolete
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF_RADIAL               "GResponse::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL       "GResponse::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE             "GResponse::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_EBOUNDS_SRC                      "GResponse::ebounds_src(GEnergy&)"

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
GResponse::GResponse(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response.
 ***************************************************************************/
GResponse::GResponse(const GResponse& rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GResponse::~GResponse(void)
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
 * @param[in] rsp Response.
 * @return Response.
 ***************************************************************************/
GResponse& GResponse::operator=(const GResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @brief Convolve sky model with the instrument response
 *
 * @param[in] model Sky model.
 * @param[in] event Event.
 * @param[in] obs Observation.
 * @param[in] grad Should model gradients be computed? (default: true)
 * @return Event probability.
 *
 * Computes the event probability
 *
 * \f[
 *    P(p',E',t') = \int \int \int
 *                  S(p,E,t) \times R(p',E',t'|p,E,t) \, dp \, dE \, dt
 * \f]
 *
 * without taking into account any energy or time dispersion. If energy or
 * time dispersion should be considered the method needs to be reimplemented
 * on the level of the instrument specific response class.
 ***************************************************************************/
double GResponse::convolve(const GModelSky&    model,
                           const GEvent&       event,
                           const GObservation& obs,
                           const bool&         grad) const
{
    // Initialise result
    double prob = 0.0;

    // Continue only if the model has a spatial component
    if (model.spatial() != NULL) {

        // Get source energy and time (no dispersion)
        GEnergy srcEng  = event.energy();
        GTime   srcTime = event.time();

        // Set source
        GSource source(model.name(), model.spatial(), srcEng, srcTime);
        
        // Get IRF value. This method returns the spatial component of the
        // source model.
        double irf = this->irf(event, source, obs);

        // If required, apply instrument specific model scaling
        if (model.has_scales()) {
            irf *= model.scale(obs.instrument()).value();
        }

        // Case A: evaluate gradients
        if (grad) {

            // Evaluate source model
            double spec = (model.spectral() != NULL) ? model.spectral()->eval_gradients(srcEng, srcTime) : 1.0;
            double temp = (model.temporal() != NULL) ? model.temporal()->eval_gradients(srcTime) : 1.0;

            // Set probability
            prob = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
                std::cout << "*** ERROR: GResponse::convolve:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (prob=" << prob;
                std::cout << ", spec=" << spec;
                std::cout << ", temp=" << temp;
                std::cout << ", irf=" << irf;
                std::cout << ")" << std::endl;
            }
            #endif

            // Multiply factors to spectral gradients
            if (model.spectral() != NULL) {
                double fact = temp * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < model.spectral()->size(); ++i) {
                        (*model.spectral())[i].factor_gradient((*model.spectral())[i].factor_gradient() * fact);
                    }
                }
            }

            // Multiply factors to temporal gradients
            if (model.temporal() != NULL) {
                double fact = spec * irf;
                if (fact != 1.0) {
                    for (int i = 0; i < model.temporal()->size(); ++i) {
                        (*model.temporal())[i].factor_gradient((*model.temporal())[i].factor_gradient() * fact);
                    }
                }
            }

        } // endif: gradient evaluation has been requested

        // Case B: evaluate no gradients
        else {

            // Evaluate source model
            double spec = (model.spectral() != NULL) ? model.spectral()->eval(srcEng, srcTime) : 1.0;
            double temp = (model.temporal() != NULL) ? model.temporal()->eval(srcTime) : 1.0;

            // Set probability
            prob = spec * temp * irf;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
                std::cout << "*** ERROR: GResponse::convolve:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (prob=" << prob;
                std::cout << ", spec=" << spec;
                std::cout << ", temp=" << temp;
                std::cout << ", irf=" << irf;
                std::cout << ")" << std::endl;
            }
            #endif

        }

    } // endif: Gamma-ray source model had a spatial component

    // Return event probability
    return prob;
}


/*==========================================================================
 =                                                                         =
 =                          Obsolete public methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * Returns the instrument response function for a given event, source and
 * observation.
 ***************************************************************************/
double GResponse::irf(const GEvent&       event,
                      const GSource&      source,
                      const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Select IRF depending on the spatial model type
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
        default:
            break;
    }

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of point source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Value of instrument response function for a point source.
 *
 * This method returns the value of the instrument response function for a
 * point source. The method assumes that source.model() is of type
 * GModelSpatialPointSource.
 ***************************************************************************/
double GResponse::irf_ptsrc(const GEvent&       event,
                            const GSource&      source,
                            const GObservation& obs) const
{
    // Get point source spatial model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Set Photon
    GPhoton photon(src->dir(), source.energy(), source.time());
    
    // Compute IRF
    double irf = this->irf(event, photon, obs);

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Return IRF value for radial source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_radial(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_RADIAL,
          "IRF computation not implemented for radial models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return IRF value for elliptical source model
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_elliptical(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_ELLIPTICAL,
          "IRF computation not implemented for elliptical models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return value of diffuse source instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GResponse::irf_diffuse(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_IRF_DIFFUSE,
          "IRF computation not implemented for diffuse models.");

    // Return IRF
    return 0.0;
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented.
 ***************************************************************************/
GEbounds GResponse::ebounds_src(const GEnergy& obsEnergy) const
{
    // Feature not yet implemented
    throw GException::feature_not_implemented(G_EBOUNDS_SRC,
          "Energy boundary computation not implemented.");

    // Allocate dummy energy boundaries
    GEbounds ebounds;

    // Return energy boundaries
    return (ebounds);
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GResponse::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response.
 ***************************************************************************/
void GResponse::copy_members(const GResponse& rsp)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GResponse::free_members(void)
{
    // Return
    return;
}
