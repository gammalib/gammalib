/***************************************************************************
 *                GResponse.cpp - Abstract response base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
#include "GIntegral.hpp"
#include "GResponse.hpp"
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GSource.hpp"
#include "GEbounds.hpp"
#include "GObservation.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialComposite.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF_RADIAL               "GResponse::irf_radial(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_ELLIPTICAL       "GResponse::irf_elliptical(GEvent&, GSource&,"\
                                                            " GObservation&)"
#define G_IRF_DIFFUSE             "GResponse::irf_diffuse(GEvent&, GSource&,"\
                                                            " GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_CONVOLVE_EDISP    //!< Debug convolve for energy dispersion


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
 * without taking into account any time dispersion. Energy dispersion is
 * correctly handled by this method. If time dispersion is indeed needed,
 * an instrument specific method needs to be provided.
 ***************************************************************************/
double GResponse::convolve(const GModelSky&    model,
                           const GEvent&       event,
                           const GObservation& obs,
                           const bool&         grad) const
{
    // Set number of iterations for Romberg integration.
    static const int iter = 6;

    // Initialise result
    double prob = 0.0;

    // Continue only if the model has a spatial component
    if (model.spatial() != NULL) {

        // Get source time (no dispersion)
        GTime srcTime = event.time();

        // Case A: Integration
        if (use_edisp()) {

            // Retrieve true energy boundaries
            GEbounds ebounds = this->ebounds(event.energy());

            // Loop over all boundaries
            for (int i = 0; i < ebounds.size(); ++i) {

                // Get true energy boundaries in MeV
                double etrue_min = ebounds.emin(i).MeV();
                double etrue_max = ebounds.emax(i).MeV();

                // Continue only if valid
                if (etrue_max > etrue_min) {

                    // Setup integration function
                    edisp_kern integrand(this, &obs, &model, &event, srcTime, grad);
                    GIntegral  integral(&integrand);

                    // Set number of iterations
                    integral.fixed_iter(iter);

                    // Do Romberg integration
                    prob += integral.romberg(etrue_min, etrue_max);

                } // endif: interval was valid

            } // endfor: looped over intervals

        }

        // Case B: No integration (assume no energy dispersion)
        else {

            // Get source energy (no dispersion)
            GEnergy srcEng  = event.energy();

            // Evaluate probability
            prob = eval_prob(model, event, srcEng, srcTime, obs, grad);

        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
            std::cout << "*** ERROR: GResponse::convolve:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (prob=" << prob;
            std::cout << ", event=" << event;
            std::cout << ", srcTime=" << srcTime;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: spatial component valid

    // Return probability
    return prob;
}


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
double GResponse::irf(const GEvent&       event,
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

    // Return IRF value
    return irf;
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
    // Initialize members
    m_use_irf_cache  = false;   //!< Switched off by default
    m_use_nroi_cache = false;   //!< Switched off by default
    m_irf_cache.clear();
    m_nroi_cache.clear();

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
    // Copy members
    m_use_irf_cache  = rsp.m_use_irf_cache;
    m_use_nroi_cache = rsp.m_use_nroi_cache;
    m_irf_cache      = rsp.m_irf_cache;
    m_nroi_cache     = rsp.m_nroi_cache;

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


/***********************************************************************//**
 * @brief Return instrument response to point source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to point source.
 *
 * Returns the instrument response to a point source.
 ***************************************************************************/
double GResponse::irf_ptsrc(const GEvent&       event,
                            const GSource&      source,
                            const GObservation& obs) const
{
    // Get pointer to point source model
    const GModelSpatialPointSource* src =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Setup photon
    GPhoton photon(src->dir(), source.energy(), source.time());

    // Get IRF
    double irf = this->irf(event, photon, obs);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to radial source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to radial source.
 *
 * Returns the instrument response to a radial source.
 *
 * @todo Implement method
 ***************************************************************************/
double GResponse::irf_radial(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Throw exception
    std::string msg = "Response computation not yet implemented for spatial "
                      "model type \""+source.model()->type()+"\".";
    throw GException::feature_not_implemented(G_IRF_RADIAL, msg);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to elliptical source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to elliptical source.
 *
 * Returns the instrument response to a elliptical source.
 *
 * @todo Implement method
 ***************************************************************************/
double GResponse::irf_elliptical(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Throw exception
    std::string msg = "Response computation not yet implemented for spatial "
                      "model type \""+source.model()->type()+"\".";
    throw GException::feature_not_implemented(G_IRF_ELLIPTICAL, msg);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response to diffuse source
 *
 * @param[in] event Observed event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response to diffuse source.
 *
 * Returns the instrument response to a diffuse source.
 *
 * @todo Implement method
 ***************************************************************************/
double GResponse::irf_diffuse(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Throw exception
    std::string msg = "Response computation not yet implemented for spatial "
                      "model type \""+source.model()->type()+"\".";
    throw GException::feature_not_implemented(G_IRF_DIFFUSE, msg);

    // Return IRF value
    return irf;
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
double GResponse::irf_composite(const GEvent&       event,
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


/***********************************************************************//**
 * @brief Convolve sky model with the instrument response
 *
 * @param[in] model Sky model.
 * @param[in] event Event.
 * @param[in] srcEng Source energy.
 * @param[in] srcTime Source time.
 * @param[in] obs Observation.
 * @param[in] grad Should model gradients be computed? (default: true)
 * @return Event probability.
 *
 * Computes the event probability
 *
 * \f[
 *    P(p',E',t'|E,t) = \int S(p,E,t) \times R(p',E',t'|p,E,t) \, dp
 * \f]
 *
 * for a given true energy \f$E\f$ and time \f$t\f$.
 ***************************************************************************/
double GResponse::eval_prob(const GModelSky&    model,
                            const GEvent&       event,
                            const GEnergy&      srcEng,
                            const GTime&        srcTime,
                            const GObservation& obs,
                            const bool&         grad) const
{
    // Initialise result
    double prob = 0.0;

    // Continue only if the model has a spatial component
    if (model.spatial() != NULL) {

        // Set source
        GSource source(model.name(), model.spatial(), srcEng, srcTime);

        // Compute IRF value
        double irf = this->irf(event, source, obs);

        // Continue only if IRF value is positive
        if (irf > 0.0) {

            // Optionally get scaling
            double scale = (model.has_scales())
                           ? model.scale(obs.instrument()).value() : 1.0;

            // Evaluate spectral and temporal components
            double spec = (model.spectral() != NULL)
                          ? model.spectral()->eval(srcEng, srcTime, grad) : 1.0;
            double temp = (model.temporal() != NULL)
                          ? model.temporal()->eval(srcTime, grad) : 1.0;

            // Compute probability
            prob = spec * temp * irf * scale;

            // Optionally compute partial derivatives
            if (grad) {

                // Multiply factors to spectral gradients
                if (model.spectral() != NULL) {
                    double fact = temp * irf * scale;
                    if (fact != 1.0) {
                        for (int i = 0; i < model.spectral()->size(); ++i) {
                            (*model.spectral())[i].factor_gradient((*model.spectral())[i].factor_gradient() * fact);
                        }
                    }
                }

                // Multiply factors to temporal gradients
                if (model.temporal() != NULL) {
                    double fact = spec * irf * scale;
                    if (fact != 1.0) {
                        for (int i = 0; i < model.temporal()->size(); ++i) {
                            (*model.temporal())[i].factor_gradient((*model.temporal())[i].factor_gradient() * fact);
                        }
                    }
                }

                // Optionally compute scale gradient for instrument
                if (model.has_scales()) {
                    for (int i = 0; i < model.scales(); ++i) {
                        double g_scale = 0.0;
                        if (model.scale(i).name() == obs.instrument()) {
                            if (model.scale(i).is_free()) {
                                g_scale = model.scale(i).scale() * spec * temp * irf;
                            }
                        }
                        model.scale(i).factor_gradient(g_scale);
                    }
                }

            } // endif: computed partial derivatives

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(prob) || gammalib::is_infinite(prob)) {
                std::cout << "*** ERROR: GResponse::eval_prob:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (prob=" << prob;
                std::cout << ", spec=" << spec;
                std::cout << ", temp=" << temp;
                std::cout << ", irf=" << irf;
                std::cout << ")" << std::endl;
            }
            #endif

        } // endif: IRF value was positive

        // ... otherwise if gradient computation is requested then set the
        // spectral and temporal gradients to zero
        else if (grad) {

            // Reset spectral gradients
            if (model.spectral() != NULL) {
                for (int i = 0; i < model.spectral()->size(); ++i) {
                    (*model.spectral())[i].factor_gradient(0.0);
                }
            }

            // Reset temporal gradients
            if (model.temporal() != NULL) {
                for (int i = 0; i < model.temporal()->size(); ++i) {
                    (*model.temporal())[i].factor_gradient(0.0);
                }
            }

            // Optionally reset scale gradients
            if (model.has_scales()) {
                for (int i = 0; i < model.scales(); ++i) {
                    model.scale(i).factor_gradient(0.0);
                }
            }

        } // endelse: IRF value was not positive

    } // endif: Gamma-ray source model had a spatial component

    // Return event probability
    return prob;
}


/***********************************************************************//**
 * @brief Integration kernel for GResponse::edisp_kern() class
 *
 * @param[in] etrue True photon energy in MeV.
 *
 * This method implements the integration kernel needed for the
 * GResponse::edisp_kern() class.
 ***************************************************************************/
double GResponse::edisp_kern::eval(const double& etrue)
{
    // Set true energy
    GEnergy srcEng;
    srcEng.MeV(etrue);

    // Get function value
    double value = m_parent->eval_prob(*m_model, *m_event,
                                       srcEng, m_srcTime, *m_obs, m_grad);

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GResponse::edisp_kern::eval";
        std::cout << "(etrue=" << etrue << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}
