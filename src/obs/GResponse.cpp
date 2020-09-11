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
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GIntegral.hpp"
#include "GIntegrals.hpp"
#include "GResponse.hpp"
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GSource.hpp"
#include "GEbounds.hpp"
#include "GObservation.hpp"
#include "GSkyMap.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GModelSpatialDiffuse.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpatialDiffuseCube.hpp"
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
 * @param[in] grad Should model gradients be computed?
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

            // Initialise Ndarray array
            GNdarray array;

            // Loop over all boundaries
            for (int i = 0; i < ebounds.size(); ++i) {

                // Get true energy boundaries in MeV
                double etrue_min = ebounds.emin(i).MeV();
                double etrue_max = ebounds.emax(i).MeV();

                // Continue only if valid
                if (etrue_max > etrue_min) {

                    // Setup integration function
                    edisp_kern integrand(this, &obs, &model, &event, srcTime, grad);
                    GIntegrals integral(&integrand);

                    // Set number of iterations
                    integral.fixed_iter(iter);

                    // Do Romberg integration
                    if (array.size() == 0) {
                        array = integral.romberg(etrue_min, etrue_max);
                    }
                    else {
                        array += integral.romberg(etrue_min, etrue_max);
                    }

                } // endif: interval was valid

            } // endfor: looped over intervals

            // Initialise array index
            int index = 0;

            // Get probability
            prob = array(index++);

            // Set gradients
            if (grad) {

                // Set spectral gradients
                if (model.spectral() != NULL) {
                    for (int i = 0; i < model.spectral()->size(); ++i) {
                        GModelPar& par = (*(model.spectral()))[i];
                        if (par.is_free() && par.has_grad()) {
                            par.factor_gradient(array(index++));
                        }
                    }
                }

                // Set temporal gradients
                if (model.temporal() != NULL) {
                    for (int i = 0; i < model.temporal()->size(); ++i) {
                        GModelPar& par = (*(model.temporal()))[i];
                        if (par.is_free() && par.has_grad()) {
                            par.factor_gradient(array(index++));
                        }
                    }
                }

                // Optionally set scale gradient for instrument
                if (model.has_scales()) {
                    for (int i = 0; i < model.scales(); ++i) {
                        const GModelPar& par = model.scale(i);
                        if (par.name() == obs.instrument()) {
                            if (par.is_free() && par.has_grad()) {
                                par.factor_gradient(array(index++));
                            }
                        }
                    }
                }

            } // endif: set gradients

        } // endif: Case A

        // Case B: No integration (assume no energy dispersion)
        else {

            // Get source energy (no dispersion)
            GEnergy srcEng  = event.energy();

            // Evaluate probability
            prob = eval_prob(model, event, srcEng, srcTime, obs, grad);

        } // endelse: Case B

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
 * * irf_ptsrc() - for the handling of a point source
 * * irf_radial() - for radial models
 * * irf_elliptical() - for elliptical models
 * * irf_diffuse() - for diffuse models
 * * irf_composite() - for composite models
 *
 * The method implements a caching mechanism for spatial models that have all
 * parameters fixed. For those models the instrument response for a given
 * event and observation is only computed once and then stored in an internal
 * cache from which it is fetched back in case that the method is called
 * again for the same event and observation.
 ***************************************************************************/
double GResponse::irf_spatial(const GEvent&       event,
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


/***********************************************************************//**
 * @brief Remove response cache for model
 *
 * @param[in] name Model name.
 *
 * Remove response cache for model @p name from response cache.
 ***************************************************************************/
void GResponse::remove_response_cache(const std::string& name)
{
    // Remove model from response caches
    m_irf_cache.remove(name);
    m_nroi_cache.remove(name);

    // Return
    return;
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
    m_use_irf_cache             = true;   //!< Switched on by default
    m_use_nroi_cache            = true;   //!< Switched on by default
    m_irf_radial_iter_theta     = 6;
    m_irf_radial_iter_phi       = 6;
    m_irf_elliptical_iter_theta = 6;
    m_irf_elliptical_iter_phi   = 6;
    m_irf_diffuse_resolution    = 0.1;    //!< Angular resolution (deg)

    // Clear cache
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
    m_use_irf_cache             = rsp.m_use_irf_cache;
    m_use_nroi_cache            = rsp.m_use_nroi_cache;
    m_irf_radial_iter_theta     = rsp.m_irf_radial_iter_theta;
    m_irf_radial_iter_phi       = rsp.m_irf_radial_iter_phi;
    m_irf_elliptical_iter_theta = rsp.m_irf_elliptical_iter_theta;
    m_irf_elliptical_iter_phi   = rsp.m_irf_elliptical_iter_phi;
    m_irf_diffuse_resolution    = rsp.m_irf_diffuse_resolution;
    m_irf_cache                 = rsp.m_irf_cache;
    m_nroi_cache                = rsp.m_nroi_cache;

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
    const GModelSpatialPointSource* model =
          static_cast<const GModelSpatialPointSource*>(source.model());

    // Setup photon
    GPhoton photon(model->dir(), source.energy(), source.time());

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
 * Returns the instrument response to a radial source for a given event,
 * source and observation. The method computes
 *
 * \f[
 *    {\tt irf}(p', E', t') = \int_0^{\theta_{\rm max}}
 *                            M_{\rm S}(\theta | E, t) \, sin \, \theta
 *                            \int_0^{2\pi}
 *                            R(p', E', t' | \theta, \phi, E, t) \,
 *                            d\,\phi \, d\,\theta
 * \f]
 *
 * where
 * * \f$M_{\rm S}(\theta | E, t)\f$ is the radial model component,
 * * \f$R(p', E', t' | \theta, \phi, E, t)\f$ is the Instrument Response
 *   Function (IRF),
 * * \f$p'\f$ is the measured instrument direction,
 * * \f$E'\f$ is the measured or reconstructed energy,
 * * \f$t'\f$ is the measured arrival time,
 * * \f$\theta\f$ is the radial distance from the model centre,
 * * \f$\phi\f$ is the azimuthal angle around the model centre,
 * * \f$E\f$ is the true photon energy, and
 * * \f$t\f$ is the true trigger time.
 *
 * The azimuth integration is done over \f$2\pi\f$ using the
 * irf_radial_kern_phi::eval() method.
 * The radial integration is done out to a maximum angle
 * \f$\theta_{\rm max}\f$, given by GModelSpatialRadial::theta_max(), by the
 * irf_radial_kern_theta::eval() method.
 ***************************************************************************/
double GResponse::irf_radial(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to radial model
    const GModelSpatialRadial* model =
          static_cast<const GModelSpatialRadial*>(source.model());

    // Get source attributes
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Set radial integration range (radians)
    double theta_min = 0.0;
    double theta_max = model->theta_max() - 1.0e-12; // Kludge to stay inside
                                                     // the boundary of a
                                                     // sharp edged model

    // Integrate if interval is valid
    if (theta_max > theta_min) {

        // Allocate Euler and rotation matrices
        GMatrix ry;
        GMatrix rz;

        // Set up rotation matrix to rotate from native model coordinates to
        // celestial coordinates
        ry.eulery( model->dir().dec_deg() - 90.0);
        rz.eulerz(-model->dir().ra_deg());
        GMatrix rot = (ry * rz).transpose();

        // Setup integration kernel
        irf_radial_kern_theta integrand(this,
                                        &event,
                                        &obs,
                                        model,
                                        &rot,
                                        &srcEng,
                                        &srcTime,
                                        m_irf_radial_iter_phi);

        // Integrate over model's radial coordinate
        GIntegral integral(&integrand);
        integral.fixed_iter(m_irf_radial_iter_theta);

        // Integrate kernel
        irf = integral.romberg(theta_min, theta_max, m_irf_radial_iter_theta);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GResponse::irf_radial:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", theta_min=" << theta_min;
            std::cout << ", theta_max=" << theta_max << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: integration interval was valid

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
 * Returns the instrument response to an elliptical source for a given event,
 * source and observation. The method computes
 *
 * \f[
 *    {\tt irf}(p', E', t') = \int_0^{\theta_{\rm max}} \int_0^{2\pi}
 *                            M_{\rm S}(\theta, \phi | E, t) \,
 *                            R(p', E', t' | \theta, \phi, E, t) \,
 *                            sin \, \theta \,
 *                            d\,\phi \, d\,\theta
 * \f]
 *
 * where
 * * \f$M_{\rm S}(\theta, \phi | E, t)\f$ is the elliptical model component,
 * * \f$R(p', E', t' | \theta, \phi, E, t)\f$ is the Instrument Response
 *   Function (IRF),
 * * \f$p'\f$ is the measured instrument direction,
 * * \f$E'\f$ is the measured or reconstructed energy,
 * * \f$t'\f$ is the measured arrival time,
 * * \f$\theta\f$ is the radial distance from the model centre,
 * * \f$\phi\f$ is the azimuthal angle around the model centre,
 * * \f$E\f$ is the true photon energy, and
 * * \f$t\f$ is the true trigger time.
 *
 * The azimuth integration is done over \f$2\pi\f$ using the
 * irf_elliptical_kern_phi::eval() method.
 * The radial integration is done out to a maximum angle
 * \f$\theta_{\rm max}\f$, given by GModelSpatialElliptical::theta_max(), by
 * the irf_elliptical_kern_theta::eval() method.
 ***************************************************************************/
double GResponse::irf_elliptical(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get pointer to elliptical model
    const GModelSpatialElliptical* model =
          static_cast<const GModelSpatialElliptical*>(source.model());

    // Get source attributes
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Set radial integration range (radians)
    double theta_min = 0.0;
    double theta_max = model->theta_max() - 1.0e-6; // Kludge to stay inside
                                                    // the boundary of a
                                                    // sharp edged model

    // Integrate if interval is valid
    if (theta_max > theta_min) {

        // Allocate Euler and rotation matrices
        GMatrix ry;
        GMatrix rz;

        // Set up rotation matrix to rotate from native model coordinates to
        // celestial coordinates
        ry.eulery( model->dir().dec_deg() - 90.0);
        rz.eulerz(-model->dir().ra_deg());
        GMatrix rot = (ry * rz).transpose();

        // Setup integration kernel
        irf_elliptical_kern_theta integrand(this,
                                            &event,
                                            &obs,
                                            model,
                                            &rot,
                                            &srcEng,
                                            &srcTime,
                                            m_irf_elliptical_iter_phi);

        // Integrate over model's radial coordinate
        GIntegral integral(&integrand);
        integral.fixed_iter(m_irf_elliptical_iter_theta);

        // Integrate kernel
        irf = integral.romberg(theta_min, theta_max, m_irf_elliptical_iter_theta);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: GResponse::irf_elliptical:";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", theta_min=" << theta_min;
            std::cout << ", theta_max=" << theta_max << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: integration interval was valid

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
 * Returns the instrument response to a diffuse source for a given event,
 * source and observation. The method computes
 *
 * \f[
 *    {\tt irf}(p', E', t') = \sum_p M_{\rm S}(p | E, t) \, \Delta(p)
 *                                   R(p', E', t' | p, E, t)
 * \f]
 *
 * where
 * * \f$M_{\rm S}(p | E, t)\f$ is the diffuse model component,
 * * \f$\Delta(p)\f$ is the solid angle of the map pixel
 * * \f$R(p', E', t' | p, E, t)\f$ is the Instrument Response Function (IRF),
 * * \f$p'\f$ is the measured instrument direction,
 * * \f$E'\f$ is the measured or reconstructed energy,
 * * \f$t'\f$ is the measured arrival time,
 * * \f$p\f$ is the true photon arrival direction,
 * * \f$E\f$ is the true photon energy, and
 * * \f$t\f$ is the true trigger time.
 *
 * The method simply adds up all sky map pixels multiplied with the IRF for
 * all diffuse model pixels. For models that do not contain sky map pixels
 * (such as the GModelSpatialDiffuseConst model) the method allocates an
 * all-sky sky map with an angular resolution that is specified by the
 * m_irf_diffuse_resolution member.
 ***************************************************************************/
double GResponse::irf_diffuse(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Get source attributes
    const GEnergy& srcEng  = source.energy();
    const GTime&   srcTime = source.time();

    // Initialise map
    //GSkyMap skymap;

    // If model is a diffuse map model then extract sky map
    const GModelSpatialDiffuseMap* map =
          dynamic_cast<const GModelSpatialDiffuseMap*>(source.model());
    if (map != NULL) {

        // Get reference to sky map
        const GSkyMap& skymap = map->map();

        // Loop over all sky map pixels
        for (int i = 0; i < skymap.npix(); ++i) {

            // Get map value
            double value = skymap(i);

            // Go to next pixel if map value is not positive
            if (value <= 0.0) {
                continue;
            }

            // Get sky direction
            GSkyDir srcDir = skymap.inx2dir(i);

            // Setup photon
            GPhoton photon(srcDir, srcEng, srcTime);

            // Add IRF multiplied by flux in map pixel. Since the map value
            // is per steradian we have to multiply with the solid angle of
            // the map pixel
            irf += this->irf(event, photon, obs) * value * skymap.solidangle(i);

        } // endfor: looped over all sky map pixels

    } // endif: model was diffuse map

    // ... otherwise if model is a cube map model then extract sky map
    else {
        const GModelSpatialDiffuseCube* cube =
              dynamic_cast<const GModelSpatialDiffuseCube*>(source.model());
        if (cube != NULL) {

            // Get reference to sky cube
            const GSkyMap& skymap = cube->cube();

            // Loop over all sky map pixels
            for (int i = 0; i < skymap.npix(); ++i) {

                // Get sky direction
                GSkyDir srcDir = skymap.inx2dir(i);

                // Setup photon
                GPhoton photon(srcDir, srcEng, srcTime);

                // Get map value
                double value = cube->eval(photon);

                // Go to next pixel if map value is not positive
                if (value <= 0.0) {
                    continue;
                }

                // Add IRF multiplied by flux in map pixel. Since the map
                // value is per steradian we have to multiply with the solid
                // angle of the map pixel
                irf += this->irf(event, photon, obs) * value * skymap.solidangle(i);

            } // endfor: looped over all sky map pixels

        } // endelse: model was diffuse cube

        // ... otherwise allocate all sky map with specified resolution
        else {

            // Allocate sky map
            int     nx     = int(360.0/m_irf_diffuse_resolution + 0.5);
            int     ny     = int(180.0/m_irf_diffuse_resolution + 0.5);
            double  dx     = 360.0 / double(nx);
            double  dy     = 180.0 / double(ny);
            GSkyMap skymap = GSkyMap("CAR", "CEL", 0.0, 0.0, -dx, dy, nx, ny);

            // Get pointer to diffuse model
            const GModelSpatialDiffuse* model =
                  static_cast<const GModelSpatialDiffuse*>(source.model());

            // Loop over all sky map pixels
            for (int i = 0; i < skymap.npix(); ++i) {

                // Get sky direction
                GSkyDir srcDir = skymap.inx2dir(i);

                // Setup photon
                GPhoton photon(srcDir, srcEng, srcTime);

                // Get model value
                double value = model->eval(photon);

                // Go to next pixel if map value is not positive
                if (value <= 0.0) {
                    continue;
                }

                // Add IRF multiplied by flux in map pixel. Since the map
                // value is per steradian we have to multiply with the solid
                // angle of the map pixel
                irf += this->irf(event, photon, obs) * value * skymap.solidangle(i);

            } // endfor: looped over sky map pixels

        } // endelse: model was not a cube

    } // endelse: model was not a map

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
        irf += irf_spatial(event, src, obs) * model->scale(i);

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
        double irf = irf_spatial(event, source, obs);

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
 * @brief Constructor for energy dispersion integration kernel class
 *
 * @param[in] parent Pointer to response.
 * @param[in] obs Pointer to observation.
 * @param[in] model Pointer to sky model.
 * @param[in] event Pointer to event.
 * @param[in] srcTime True time.
 * @param[in] grad Evaluate gradients?
 *
 * This method constructs the integration kernel needed for the energy
 * dispersion computation.
 ***************************************************************************/
GResponse::edisp_kern::edisp_kern(const GResponse*    parent,
                                  const GObservation* obs,
                                  const GModelSky*    model,
                                  const GEvent*       event,
                                  const GTime&        srcTime,
                                  const bool&         grad)
{
    // Set members
    m_parent  = parent;
    m_obs     = obs;
    m_model   = model;
    m_event   = event;
    m_srcTime = srcTime;
    m_grad    = grad;

    // If gradients are requested then put pointers to all relevant parameter
    // in the parameter pointer vector
    m_pars.clear();
    if (m_grad) {
        if (m_model->spectral() != NULL) {
            for (int i = 0; i < m_model->spectral()->size(); ++i) {
                GModelPar& par = (*(m_model->spectral()))[i];
                if (par.is_free() && par.has_grad()) {
                    m_pars.push_back(&par);
                }
            }
        }
        if (m_model->temporal() != NULL) {
            for (int i = 0; i < m_model->temporal()->size(); ++i) {
                GModelPar& par = (*(m_model->temporal()))[i];
                if (par.is_free() && par.has_grad()) {
                    m_pars.push_back(&par);
                }
            }
        }
        if (m_model->has_scales()) {
            for (int i = 0; i < m_model->scales(); ++i) {
                GModelPar& par = const_cast<GModelPar&>(m_model->scale(i));
                if (par.name() == m_obs->instrument()) {
                    if (par.is_free() && par.has_grad()) {
                        m_pars.push_back(&par);
                    }
                }
            }
        }
    }

    // Allocate Ndarray to be returned as return type
    m_array = GNdarray(m_pars.size()+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return energy dispersion integration kernel
 *
 * @return Ndarray.
 *
 * This method returns an Ndarray with the last values that were computed
 * by the eval() method. In case that eval() was not called before, the
 * method returns an Ndarray with the right dimension and with all values
 * set to zero.
 ***************************************************************************/
const GNdarray& GResponse::edisp_kern::array(void) const
{
    // Return Ndarray
    return m_array;
}


/***********************************************************************//**
 * @brief Evaluate energy dispersion integration kernel
 *
 * @param[in] etrue True photon energy in MeV.
 *
 * This method implements the integration kernel needed for the
 * GResponse::edisp_kern() class.
 ***************************************************************************/
GNdarray GResponse::edisp_kern::eval(const double& etrue)
{
    // Initialise Ndarray array index
    int index = 0;

    // Set true energy
    GEnergy srcEng;
    srcEng.MeV(etrue);

    // Get function value
    m_array(index++) = m_parent->eval_prob(*m_model, *m_event,
                                           srcEng, m_srcTime, *m_obs,
                                           m_grad);

    // If gradients are requested then extract them and put them into the
    // array
    if (m_grad) {
        for (int i = 0; i < m_pars.size(); ++i, ++index) {
            m_array(index) = m_pars[i]->factor_gradient();
        }
    }

    // Compile option: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(m_array(0)) || gammalib::is_infinite(m_array(0))) {
        std::cout << "*** ERROR: GResponse::edisp_kern::eval";
        std::cout << "(etrue=" << etrue << "): ";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << m_array(0);
        std::cout << ")" << std::endl;
    }
    #endif

    // Return array
    return m_array;
}


/***********************************************************************//**
 * @brief Zenith angle integration kernel for radial model
 *
 * @param[in] theta Zenith angle (radians).
 * @return IRF value.
 *
 * Integrated the IRF multiplied by the model value for a given zenith angle
 * over all azimuth angles.
 ***************************************************************************/
double GResponse::irf_radial_kern_theta::eval(const double& theta)
{
    // Initialise result
    double irf = 0.0;

    // Continue only for positive zenith angles (otherwise the integral will
    // be zero)
    if (theta > 0.0) {

        // Evaluate sky model
        double model = m_model->eval(theta, *m_srcEng, *m_srcTime);

        // Continue only if model is positive
        if (model > 0.0) {

            // Set azimuthal integration range
            double phi_min = 0.0;
            double phi_max = gammalib::twopi;

            // Setup integration kernel
            irf_radial_kern_phi integrand(m_rsp,
                                          m_event,
                                          m_obs,
                                          m_rot,
                                          theta,
                                          m_srcEng,
                                          m_srcTime);

            // Setup integration
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter_phi);

            // Integrate over phi
            irf = integral.romberg(phi_min, phi_max, m_iter_phi) *
                  model * std::sin(theta);

        } // endif: model was positive

    } // endif: theta was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Azimuth angle integration kernel for radial model
 *
 * @param[in] phi Azimuth angle (radians).
 * @return IRF value.
 *
 * Computes the IRF at a given zenith and azimuth angle with respect to a
 * specified centre.
 ***************************************************************************/
double GResponse::irf_radial_kern_phi::eval(const double& phi)
{
    // Set up native coordinate vector
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate vector into celestial coordinates
    GVector vector = (*m_rot) * native;

    // Convert vector into sky direction
    GSkyDir dir(vector);

    // Setup photon
    GPhoton photon(dir, *m_srcEng, *m_srcTime);

    // Evaluate IRF for photon
    double irf = m_rsp->irf(*m_event, photon, *m_obs);

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Zenith angle integration kernel for elliptical model
 *
 * @param[in] theta Zenith angle (radians).
 * @return IRF value.
 *
 * Integrated the IRF multiplied by the model value for a given zenith angle
 * over all azimuth angles.
 ***************************************************************************/
double GResponse::irf_elliptical_kern_theta::eval(const double& theta)
{
    // Initialise result
    double irf = 0.0;

    // Continue only for positive zenith angles (otherwise the integral will
    // be zero)
    if (theta > 0.0) {

        // Set azimuthal integration range
        double phi_min = 0.0;
        double phi_max = gammalib::twopi;

        // Setup integration kernel
        irf_elliptical_kern_phi integrand(m_rsp,
                                          m_event,
                                          m_obs,
                                          m_model,
                                          m_rot,
                                          theta,
                                          m_srcEng,
                                          m_srcTime);

        // Setup integration
        GIntegral integral(&integrand);
        integral.fixed_iter(m_iter_phi);

        // Integrate over phi
        irf = integral.romberg(phi_min, phi_max, m_iter_phi) * std::sin(theta);

    } // endif: theta was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Azimuth angle integration kernel for elliptical model
 *
 * @param[in] phi Azimuth angle (radians).
 * @return IRF value.
 *
 * Computes the IRF at a given zenith and azimuth angle with respect to a
 * specified centre.
 ***************************************************************************/
double GResponse::irf_elliptical_kern_phi::eval(const double& phi)
{
    // Initialise result
    double irf = 0.0;

    // Evaluate sky model
    double model = m_model->eval(m_theta, phi, *m_srcEng, *m_srcTime);

    // Continue only if model is positive
    if (model > 0.0) {

        // Set up native coordinate vector
        double  cos_phi = std::cos(phi);
        double  sin_phi = std::sin(phi);
        GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

        // Rotate vector into celestial coordinates
        GVector vector = (*m_rot) * native;

        // Convert vector into sky direction
        GSkyDir dir(vector);

        // Setup photon
        GPhoton photon(dir, *m_srcEng, *m_srcTime);

        // Evaluate IRF for photon
        irf = m_rsp->irf(*m_event, photon, *m_obs) * model;

    } // endif: model was positive

    // Return
    return irf;
}
