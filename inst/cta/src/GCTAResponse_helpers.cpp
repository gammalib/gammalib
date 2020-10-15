/***************************************************************************
 *         GCTAResponse_helpers.cpp - CTA response helper classes          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse_helpers.cpp
 * @brief CTA response helper classes implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GIntegrals.hpp"
#include "GVector.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAEdisp.hpp"
#include "GCTASupport.hpp"
#include "GModelSpatialDiffuseMap.hpp"
#include "GModelSpatialDiffuseCube.hpp"
#include "GModelSpatialDiffuseConst.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PSF_RADIAL_KERNS_PHI      "cta_psf_radial_kerns_phi::eval(double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_OBSDIR_FOR_AEFF    //!< Use event offset for Aeff computation

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_INTEGRAL                             //!< Debug integration
//#define G_DEBUG_MODEL_ZERO           //!< Debug check for zero model values

/* __ Constants __________________________________________________________ */
const double g_kludge_radius = 1.0e-12;        //!< Tiny angle (radians)
const double g_ellipse_kludge_radius = 1.0e-6; //!< About 0.2 arc seconds
                                               //   A larger radius is used
                                               //   as the ellipse is defined
                                               //   in cartesian coordinates
                                               //   while distances are
                                               //   computed in spherical
                                               //   trigonometry


/*==========================================================================
 =                                                                         =
 =                              Helper functions                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Limit omega interval
 *
 * @param[in] min Interval minimum (radians).
 * @param[in] max Interval maximum (radians).
 * @param[in] domega Half length of interval (radians).
 * @return Vector of intervals.
 *
 * Limits an omega interval [@p min,@p max] to the interval specified by
 * [-@p domega,@p domega]. This may lead to a split of [@p min,@p max] in
 * several intervals. The method thus returns a vector of intervals that
 * overlap with [-@p domega,@p domega]. If there is no overlap with the
 * interval, the method returns an empty vector.
 *
 * The method takes care of wrap arounds. It is assumed that on input
 * [@p min,@p max] is contained within [-2pi,+2pi].
 ***************************************************************************/
cta_omega_intervals gammalib::limit_omega(const double& min,
                                          const double& max,
                                          const double& domega)
{
    // Allocate intervals
    cta_omega_intervals intervals;

    // Continue only if domega is smaller than pi
    if (domega < gammalib::pi) {

        // Set limiting intervals. To take care of a possible wrap around in
        // omega we consider also intervals that are shifted by +/-2pi
        double omega_min       = -domega;
        double omega_max       = +domega;
        double omega_min_plus  = omega_min + gammalib::twopi;
        double omega_max_plus  = omega_max + gammalib::twopi;
        double omega_min_minus = omega_min - gammalib::twopi;
        double omega_max_minus = omega_max - gammalib::twopi;

        // If the [min,max] interval overlaps with the unshifted
        // [-domega,domega] interval then constrain the interval
        if (max > omega_min && min < omega_max) {
            double interval_min = min;
            double interval_max = max;
            if (interval_min < omega_min) {
                interval_min = omega_min;
            }
            if (interval_max > omega_max) {
                interval_max = omega_max;
            }
            intervals.push_back(std::make_pair(interval_min,interval_max));
        }

        // If the [min,max] interval overlaps with the [-domega,domega]
        // interval shifted by +2pi then constrain the interval using the
        // shifted interval
        if (max > omega_min_plus && min < omega_max_plus) {
            double interval_min = min;
            double interval_max = max;
            if (interval_min < omega_min_plus) {
                interval_min = omega_min_plus;
            }
            if (interval_max > omega_max_plus) {
                interval_max = omega_max_plus;
            }
            intervals.push_back(std::make_pair(interval_min,interval_max));
        }

        // If the [min,max] interval overlaps with the [-domega,domega]
        // interval shifted by -2pi then constrain the interval using the
        // shifted interval
        if (max > omega_min_minus && min < omega_max_minus) {
            double interval_min = min;
            double interval_max = max;
            if (interval_min < omega_min_minus) {
                interval_min = omega_min_minus;
            }
            if (interval_max > omega_max_minus) {
                interval_max = omega_max_minus;
            }
            intervals.push_back(std::make_pair(interval_min,interval_max));
        }

    } // endif: interval was not the full circle

    // ... otherwise append the provided interval
    else {
        intervals.push_back(std::make_pair(min,max));
    }

    // Return intervals
    return intervals;
}


/***********************************************************************//**
 * @brief Determine resolution of spatial model
 *
 * @param[in] model Pointer to spatial model.
 * @return Resolution of spatial model (radians).
 *
 * Determine the resolution of a spatial model. So far the method only works
 * for a spatial map or cube model holding a WCS projection. If a constant
 * spatial model is encountered a resolution of 180 deg is returned.
 *
 * If the resolution of the model could not be determined, the method returns
 * a resolution of 0.01 deg.
 ***************************************************************************/
double gammalib::resolution(const GModelSpatial* model)
{
    // Initialise resolution to default resolution
    double resolution = 0.01 * gammalib::deg2rad;

    // Extract pointer to spatial map. The pointer will be NULL if no spatial
    // map exists.
    const GSkyMap* map = NULL;
    const GModelSpatialDiffuseMap* pmap =
          dynamic_cast<const GModelSpatialDiffuseMap*>(model);
    if (pmap != NULL) {
        map = &(pmap->map());
    }
    else {
        const GModelSpatialDiffuseCube* pcube =
              dynamic_cast<const GModelSpatialDiffuseCube*>(model);
        if (pcube != NULL) {
            map = &(pcube->cube());
        }
        else {
            const GModelSpatialDiffuseConst* pconst =
                  dynamic_cast<const GModelSpatialDiffuseConst*>(model);
            if (pconst != NULL) {
                resolution = gammalib::pi;
            }
        }
    }

    // If a spatial map exists then get it's resolution. This so far only
    // works for WCS maps.
    if (map != NULL) {
        const GWcs* wcs = dynamic_cast<const GWcs*>(map->projection());
        if (wcs != NULL) {
            double dx = std::abs(wcs->cdelt(0));
            double dy = std::abs(wcs->cdelt(1));
            resolution = ((dx < dy) ? dx : dy) * gammalib::deg2rad;
        }
    }

    // Return resolution
    return resolution;
}


/*==========================================================================
 =                                                                         =
 =              Helper class methods for response computation              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Integration kernel for npsf() method
 *
 * @param[in] delta Distance from PSF centre (radians).
 * @return Azimuthally integrated PSF.
 *
 * Computes
 *
 * \f[
 *    \int_{0}^{\phi} PSF(\delta) d\phi
 * \f]
 * 
 * for a given offset angle \f$\delta\f$.
 *
 * The azimuthal integration is performed over an arclength given by
 * \f$\phi\f$. The method actually assumes that the PSF is azimuthally
 * symmetric, hence it just multiplies the PSF value by the arclength times
 * the sinus of the offset angle.
 ***************************************************************************/
double cta_npsf_kern_rad_azsym::eval(const double& delta)
{
    // Initialise PSF value
    double value = 0.0;

    // Get arclength for given radius in radians
    double phi = gammalib::roi_arclength(delta,
                                         m_psf,
                                         m_cospsf,
                                         m_sinpsf,
                                         m_roi,
                                         m_cosroi);

    // If arclength is positive then compute the PSF value
    if (phi > 0) {

        // Compute PSF value
        value = m_rsp->psf(delta, m_theta, m_phi, m_zenith, m_azimuth, m_logE) * 
                           phi * std::sin(delta);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: cta_npsf_kern_rad_azsym::eval";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", delta=" << delta;
            std::cout << ", phi=" << phi << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: arclength was positive

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Integration kernel for GCTAResponseIrf::nroi method
 *
 * @param[in] etrue True photon energy in MeV.
 * @return Nroi.
 ***************************************************************************/
double cta_nroi_kern::eval(const double& etrue)
{
    // Set true energy
    GEnergy srcEng;
    srcEng.MeV(etrue);

    // Compute response components
    double nroi_spatial  = m_rsp->nroi(*m_model, srcEng, m_srcTime, m_obsEng, m_obsTime, *m_obs);
    double nroi_spectral = m_model->spectral()->eval(srcEng, m_srcTime);
    double nroi_temporal = m_model->temporal()->eval(m_srcTime);

    // Compute response
    double nroi = nroi_spatial * nroi_spectral * nroi_temporal;

    // Return response
    return nroi;
}


/*==========================================================================
 =                                                                         =
 =                  Helper class methods for radial models                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for radial model zenith angle integration of Irf
 *
 * @param[in] rho Zenith angle with respect to model centre [radians].
 *
 * Computes the kernel 
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     IRF(\rho, \omega) d\omega
 * \f]
 *
 * for the zenith angle integration of radial models.
 ***************************************************************************/
double cta_irf_radial_kern_rho::eval(const double& rho)
{
    // Initialise result
    double irf = 0.0;

    // Continue only if rho is positive (otherwise the integral will be
    // zero)
    if (rho > 0.0) {

        // Compute half length of arc that lies within PSF validity circle
        // (in radians)
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_zeta,
                                                      m_cos_zeta,
                                                      m_sin_zeta,
                                                      m_delta_max,
                                                      m_cos_delta_max);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Compute omega integration range
            double omega_min = -domega;
            double omega_max = +domega;

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model
            double rho_kludge = rho - g_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Evaluate sky model
            double model = m_model->eval(rho_kludge, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_irf_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model=" << (rho-m_model.theta_max());
                std::cout << " radians" << std::endl;
            }
            #endif

            // Continue only if model is positive
            if (model > 0.0) {

                // Precompute cosine and sine terms for azimuthal
                // integration
                double cos_rho = std::cos(rho);
                double sin_rho = std::sin(rho);
                double cos_psf = cos_rho*m_cos_zeta;
                double sin_psf = sin_rho*m_sin_zeta;
                double cos_ph  = cos_rho*m_cos_lambda;
                double sin_ph  = sin_rho*m_sin_lambda;

                // Setup integration kernel
                cta_irf_radial_kern_omega integrand(m_rsp,
                                                    m_zenith,
                                                    m_azimuth,
                                                    m_srcEng,
                                                    m_obsEng,
                                                    m_zeta,
                                                    m_lambda,
                                                    m_omega0,
                                                    rho,
                                                    cos_psf,
                                                    sin_psf,
                                                    cos_ph,
                                                    sin_ph);

                // Integrate over phi
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);
                irf = integral.romberg(omega_min, omega_max, m_iter) *
                      model * sin_rho;

                // Compile option: Check for NaN/Inf
                #if defined(G_NAN_CHECK)
                if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                    std::cout << "*** ERROR: cta_irf_radial_kern_rho";
                    std::cout << "(rho=" << rho << "):";
                    std::cout << " NaN/Inf encountered";
                    std::cout << " (irf=" << irf;
                    std::cout << ", domega=" << domega;
                    std::cout << ", model=" << model;
                    std::cout << ", sin_rho=" << sin_rho << ")";
                    std::cout << std::endl;
                }
                #endif

            } // endif: model was positive

        } // endif: arclength was positive

    } // endif: rho was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for radial model azimuth angle IRF integration
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * Computes the kernel 
 *
 * \f[
 *    IRF(\rho,\omega)
 * \f]
 *
 * for the azimuth angle integration of radial models.
 *
 * From the model coordinates \f$(\rho,\omega)\f$, the method computes the
 * angle between the true (\f$\vec{p}\f$) and observed (\f$\vec{p'}\f$) 
 * photon arrival direction using
 * 
 * \f[\delta = \arccos(\cos \rho \cos \zeta + 
 *                     \sin \rho \sin \zeta \cos \omega)\f]
 *
 * where
 * \f$\zeta\f$ is the angular distance between the observed photon arrival
 * direction \f$\vec{p'}\f$ and the model centre \f$\vec{m}\f$. This angle
 * \f$\delta\f$ is used to compute the \f$PSF(\delta)\f$ value.
 *
 * The method computes also the angle \f$\theta\f$ between the observed
 * photon arrival direction \f$\vec{p'}\f$ and the camera pointing
 * \f$\vec{d}\f$ using
 *
 * \f[\theta = \arccos(\cos \rho \cos \lambda + 
 *                     \sin \rho \sin \lambda \cos \omega_0 - \omega)\f]
 *
 * where
 * \f$\lambda\f$ is the angular distance between the model centre 
 * \f$\vec{m}\f$ and the camera pointing direction \f$\vec{d}\f$.
 * The angle \f$\theta\f$ is used in the computation of the IRF (no
 * azimuthal dependence is so far implemented for the IRF computation).
 ***************************************************************************/
double cta_irf_radial_kern_omega::eval(const double& omega)
{
    // Get log10(E/TeV) of true photon energy
    double srcLogEng = m_srcEng.log10TeV();

    // Compute PSF offset angle [radians]
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));

    // Compute true photon offset angle in camera system [radians]
    double offset = std::acos(m_cos_ph + m_sin_ph * std::cos(m_omega0 - omega));

    //TODO: Compute true photon azimuth angle in camera system [radians]
    double azimuth = 0.0;

    // Evaluate IRF
    double irf = m_rsp->aeff(offset, azimuth, m_zenith, m_azimuth, srcLogEng) *
                 m_rsp->psf(delta, offset, azimuth, m_zenith, m_azimuth, srcLogEng);

    // Optionally take energy dispersion into account
    if (m_rsp->use_edisp() && irf > 0.0) {
        irf *= m_rsp->edisp(m_obsEng, m_srcEng, offset, azimuth, m_zenith, m_azimuth);
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: cta_irf_radial_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", delta=" << delta;
        std::cout << ", offset=" << offset;
        std::cout << ", azimuth=" << azimuth << ")";
        std::cout << std::endl;
    }
    #endif

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for zenith angle Nroi integration or radial model
 *
 * @param[in] rho Radial model zenith angle (radians).
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     N_{\rm pred}(\rho,\omega) d\omega
 * \f]
 * 
 * The azimuth angle integration range 
 * \f$[\omega_{\rm min}, \omega_{\rm max}]\f$
 * is limited to an arc around the vector connecting the model centre to
 * the ROI centre. This limitation assures that the integration converges
 * properly.
 ***************************************************************************/
double cta_nroi_radial_kern_rho::eval(const double& rho)
{
    // Initialise Nroi value
    double nroi = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of arc that lies within ROI+PSF radius (radians)
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_dist,
                                                      m_cos_dist,
                                                      m_sin_dist,
                                                      m_radius,
                                                      m_cos_radius);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Compute omega integration range
            double omega_min = m_omega0 - domega;
            double omega_max = m_omega0 + domega;

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model
            double rho_kludge = rho - g_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Get radial model value
            double model = m_model->eval(rho_kludge, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_nroi_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model=" << (rho-m_model.theta_max());
                std::cout << " radians" << std::endl;
            }
            #endif

            // Continue only if we have a positive model value
            if (model > 0.0) {

                // Compute sine and cosine of offset angle
                double sin_rho = std::sin(rho);
                double cos_rho = std::cos(rho);

                // Setup phi integration kernel
                cta_nroi_radial_kern_omega integrand(m_rsp,
                                                     m_obs,
                                                     m_rot,
                                                     m_srcEng,
                                                     m_srcTime,
                                                     m_obsEng,
                                                     m_obsTime,
                                                     sin_rho,
                                                     cos_rho);

                // Integrate over phi
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);
                nroi = integral.romberg(omega_min, omega_max, m_iter) *
                       sin_rho * model;

                // Debug: Check for NaN
                #if defined(G_NAN_CHECK)
                if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
                    std::cout << "*** ERROR: cta_nroi_radial_kern_rho::eval";
                    std::cout << "(rho=" << rho << "):";
                    std::cout << " NaN/Inf encountered";
                    std::cout << " (nroi=" << nroi;
                    std::cout << ", model=" << model;
                    std::cout << ", omega=[" << omega_min << "," << omega_max << "]";
                    std::cout << ", sin_rho=" << sin_rho;
                    std::cout << ")" << std::endl;
                }
                #endif

            } // endif: model was positive

        } // endif: arc length was positive

    } // endif: rho was positive

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Nroi integration of radial model
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double cta_nroi_radial_kern_omega::eval(const double& omega)
{
    // Compute sky direction vector in native coordinates
    double  cos_omega = std::cos(omega);
    double  sin_omega = std::sin(omega);
    GVector native(-cos_omega*m_sin_rho, sin_omega*m_sin_rho, m_cos_rho);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set Photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute point source Nroi for this sky direction
    double nroi = m_rsp->nirf(photon, m_obsEng, m_obsTime, *m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
        std::cout << "*** ERROR: cta_nroi_radial_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (nroi=" << nroi;
        std::cout << ", cos_omega=" << cos_omega;
        std::cout << ", sin_omega=" << sin_omega;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Nroi
    return nroi;
}


/*==========================================================================
 =                                                                         =
 =                Helper class methods for elliptical models               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for elliptical model integration over model's zenith angle
 *
 * @param[in] rho Radial distance from model centre [radians].
 * @return Radial IRF integration kernel.
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     S_{\rm p}(\rho, \omega | E, t) \, IRF(\rho, \omega)
 *                     d\omega
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical source model
 * for a given true photon energy and photon arrival time,
 * \f$IRF(\rho, \omega)\f$ is the instrument response function.
 *
 * The method performs the required coordinate transformations from the
 * model system, spanned by \f$(\rho, \omega)\f$, to the system needed for
 * IRF computations. Furthermore, the method limits the integration range
 * to area where the ellipse intersects the IRF.
 ***************************************************************************/
double cta_irf_elliptical_kern_rho::eval(const double& rho)
{
    // Initialise result
    double irf = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of the arc (in radians) from a circle with
        // radius rho that intersects with the point spread function, defined
        // as a circle with maximum radius m_delta_max
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_rho_obs,
                                                      m_cos_rho_obs,
                                                      m_sin_rho_obs,
                                                      m_delta_max,
                                                      m_cos_delta_max);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Precompute cosine and sine terms for azimuthal integration
            double cos_rho = std::cos(rho);
            double sin_rho = std::sin(rho);
            double cos_psf = cos_rho * m_cos_rho_obs;
            double sin_psf = sin_rho * m_sin_rho_obs;
            double cos_ph  = cos_rho * m_cos_rho_pnt;
            double sin_ph  = sin_rho * m_sin_rho_pnt;

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model (e.g. an elliptical
            // disk model)
            double rho_kludge = rho - g_ellipse_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Setup integration kernel
            cta_irf_elliptical_kern_omega integrand(m_rsp,
                                                    m_model,
                                                    m_zenith,
                                                    m_azimuth,
                                                    m_srcEng,
                                                    m_srcTime,
                                                    m_obsEng,
                                                    m_posangle_obs,
                                                    m_omega_pnt,
                                                    rho_kludge,
                                                    cos_psf,
                                                    sin_psf,
                                                    cos_ph,
                                                    sin_ph);

            // Setup integrator
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);

            // If the radius rho is not larger than the semiminor axis
            // boundary, the circle with that radius is fully contained in
            // the ellipse and we can just integrate over the relevant arc
            if (rho <= m_semiminor) {

                // Compute omega integration range
                double omega_min = -domega;
                double omega_max = +domega;

                // Integrate over omega
                irf = integral.romberg(omega_min, omega_max, m_iter) *
                      sin_rho;

            } // endif: circle comprised in ellipse

            // ... otherwise there are arcs that intersect with the Psf circle
            else {

                // Compute half the arc length (in radians) of a circle of
                // radius rho, centred on the model, that intersects with
                // the ellipse boundary
                double arg1 = 1.0 - (m_semiminor*m_semiminor) / (rho*rho);
                double arg2 = 1.0 - (m_semiminor*m_semiminor) /
                                    (m_semimajor*m_semimajor);
                double omega_width = std::acos(std::sqrt(arg1/arg2));

                // Continue only if the arclength is positive
                if (omega_width > 0.0) {

                    // Compute azimuth angle difference between ellipse
                    // position angle and position angle of observed
                    // photon in the model system. This angle will define
                    // the reference point around which the circle arc. Make
                    // sure that omega_0 is within [-pi,pi] thus that the omega
                    // intervals are between [-2pi,2pi]
                    double omega_0 = m_posangle - m_posangle_obs;
                    if (omega_0 > gammalib::pi) {
                        omega_0 -= gammalib::pi;
                    }
                    else if (omega_0 < -gammalib::pi) {
                        omega_0 += gammalib::pi;
                    }

                    // Compute azimuth angle intervals
                    double omega1_min = omega_0    - omega_width;
                    double omega1_max = omega_0    + omega_width;
                    double omega2_min = omega1_min + gammalib::pi;
                    double omega2_max = omega1_max + gammalib::pi;

                    // Limit intervals to the intersection of the ellipse with
                    // the Psf circle. This may lead to a split of intervals,
                    // and we gather all these intervals in a special interval
                    // pair containers
                    cta_omega_intervals intervals1 = 
                        gammalib::limit_omega(omega1_min, omega1_max, domega);
                    cta_omega_intervals intervals2 = 
                        gammalib::limit_omega(omega2_min, omega2_max, domega);

                    // Integrate over all intervals for omega1
                    for (int i = 0; i < intervals1.size(); ++i) {
                        double min = intervals1[i].first;
                        double max = intervals1[i].second;
                        irf       += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                    // Integrate over all intervals for omega2
                    for (int i = 0; i < intervals2.size(); ++i) {
                        double min = intervals2[i].first;
                        double max = intervals2[i].second;
                        irf       += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                } // endif: arc length was positive

            } // endelse: circle was not comprised in ellipse

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                std::cout << "*** ERROR: cta_irf_elliptical_kern_rho";
                std::cout << "(rho=" << rho << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (irf=" << irf;
                std::cout << ", domega=" << domega;
                std::cout << ", sin_rho=" << sin_rho << ")";
                std::cout << std::endl;
            }
            #endif

        } // endif: arc length was positive

    } // endif: rho was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for elliptical model integration over model's azimuth angle
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * Computes
 *
 * \f[
 *    S_{\rm p}(\omega | \rho, E, t) \, IRF(\omega | \rho)
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 *
 * From the coordinates \f$(\rho,\omega)\f$ in the model system, the method
 * computes the angle between the true (\f$\vec{p}\f$) and observed
 * (\f$\vec{p'}\f$) photon arrival direction using
 * 
 * \f[\delta = \arccos(\cos \rho \cos \zeta + 
 *                     \sin \rho \sin \zeta \cos \omega)\f]
 *
 * where
 * \f$\zeta\f$ is the angular distance between the observed photon arrival
 * direction \f$\vec{p'}\f$ and the model centre \f$\vec{m}\f$. This angle
 * \f$\delta\f$ is used to compute the \f$PSF(\delta)\f$ value.
 *
 * The method computes also the angle \f$\theta\f$ between the observed
 * photon arrival direction \f$\vec{p'}\f$ and the camera pointing
 * \f$\vec{d}\f$ using
 *
 * \f[\theta = \arccos(\cos \rho \cos \lambda + 
 *                     \sin \rho \sin \lambda \cos \omega_0 - \omega)\f]
 *
 * where
 * \f$\lambda\f$ is the angular distance between the model centre 
 * \f$\vec{m}\f$ and the camera pointing direction \f$\vec{d}\f$.
 * The angle \f$\theta\f$ is used in the computation of the IRF (no
 * azimuthal dependence is so far implemented for the IRF computation).
 ***************************************************************************/
double cta_irf_elliptical_kern_omega::eval(const double& omega)
{
    // Initialise IRF value
    double irf = 0.0;

    // Compute azimuth angle in model coordinate system (radians)
    double omega_model = omega + m_posangle_obs;

    // Evaluate sky model
    double model = m_model->eval(m_rho, omega_model, m_srcEng, m_srcTime);

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (model <= 0.0) {
        double m_semiminor_rad = m_model->semiminor() * gammalib::deg2rad;
        double m_semimajor_rad = m_model->semimajor() * gammalib::deg2rad;
        double diff_angle      = omega_model - m_model->posangle() * gammalib::deg2rad;
        double cosinus         = std::cos(diff_angle);
        double sinus           = std::sin(diff_angle);
        double arg1            = m_semiminor_rad * cosinus;
        double arg2            = m_semimajor_rad * sinus;
        double r_ellipse       = m_semiminor_rad * m_semimajor_rad /
                                 std::sqrt(arg1*arg1 + arg2*arg2);
        std::cout << "*** WARNING: cta_irf_elliptical_kern_omega::eval:";
        std::cout << " zero model for (rho,omega)=(";
        std::cout << m_rho*gammalib::rad2deg << ",";
        std::cout << omega*gammalib::rad2deg << ")";
        std::cout << " semiminor=" << m_semiminor_rad;
        std::cout << " semimajor=" << m_semimajor_rad;
        std::cout << " posang=" << m_model->posangle() * gammalib::deg2rad;
        std::cout << " rel_posang=" << diff_angle * gammalib::deg2rad;
        std::cout << " rho/r_ellipse=" << (m_rho/r_ellipse);
        std::cout << std::endl;
    }
    #endif

    // Continue only if model is positive
    if (model > 0.0) {

        // Get log10(E/TeV) of true photon energy
        double srcLogEng = m_srcEng.log10TeV();

        // Compute Psf offset angle [radians]
        double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));

        // Compute true photon offset and azimuth angle in camera system
        // [radians]
        double theta = std::acos(m_cos_ph + m_sin_ph * std::cos(m_omega_pnt - omega));
        double phi   = 0.0; //TODO: Implement IRF Phi dependence

        // Evaluate IRF * model
        irf = m_rsp->aeff(theta, phi, m_zenith, m_azimuth, srcLogEng) *
              m_rsp->psf(delta, theta, phi, m_zenith, m_azimuth, srcLogEng) *
              model;

        // Optionally take energy dispersion into account
        if (m_rsp->use_edisp() && irf > 0.0) {
            irf *= m_rsp->edisp(m_obsEng, m_srcEng, theta, phi, m_zenith, m_azimuth);
        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: cta_irf_elliptical_kern_omega::eval";
            std::cout << "(omega=" << omega << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", model=" << model;
            std::cout << ", delta=" << delta;
            std::cout << ", theta=" << theta;
            std::cout << ", phi=" << phi << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: model is positive

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for zenith angle Nroi integration of elliptical model
 *
 * @param[in] rho Elliptical model offset angle (radians).
 *
 * Computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     S_{\rm p}(\rho,\omega | E, t) \,
 *                     N_{\rm pred}(\rho,\omega) d\omega
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the ROI centre, and
 * \f$\rho\f$ is the radial distance from the model centre.
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical source model
 * for a given true photon energy and photon arrival time,
 * \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral over the
 * response function.
 *
 * The method performs the required coordinate transformations from the
 * model system, spanned by \f$(\rho, \omega)\f$, to the system needed for
 * Npred computations. Furthermore, the method limits the integration range
 * to area where the ellipse intersects the ROI.
 ***************************************************************************/
double cta_nroi_elliptical_kern_rho::eval(const double& rho)
{
    // Initialise Nroi value
    double nroi = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of arc that lies within ROI+PSF radius (radians)
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_rho_roi,
                                                      m_cos_rho_roi,
                                                      m_sin_rho_roi,
                                                      m_radius_roi,
                                                      m_cos_radius_roi);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Compute sine and cosine terms for azimuthal integration
            double sin_rho = std::sin(rho);
            double cos_rho = std::cos(rho);

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model (e.g. an elliptical
            // disk model)
            double rho_kludge = rho - g_ellipse_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Setup phi integration kernel
            cta_nroi_elliptical_kern_omega integrand(m_rsp,
                                                     m_obs,
                                                     m_model,
                                                     m_rot,
                                                     m_srcEng,
                                                     m_srcTime,
                                                     m_obsEng,
                                                     m_obsTime,
                                                     rho_kludge,
                                                     sin_rho,
                                                     cos_rho,
                                                     m_posangle_roi);

            // Setup integrator
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);

            // If the radius rho is not larger than the semiminor axis
            // boundary, the circle with that radius is fully contained in
            // the ellipse and we can just integrate over the relevant arc
            if (rho <= m_semiminor) {

                // Compute omega integration range
                double omega_min = -domega;
                double omega_max = +domega;

                // Integrate over omega
                nroi = integral.romberg(omega_min, omega_max, m_iter) *
                       sin_rho;

            } // endif: circle comprised in ellipse

            // ... otherwise there are arcs that intersect with the ROI
            // circle
            else {

                // Compute half the arc length (in radians) of a circle of
                // radius rho, centred on the model, that intersects with
                // the ellipse boundary
                double arg1 = 1.0 - (m_semiminor*m_semiminor) / (rho*rho);
                double arg2 = 1.0 - (m_semiminor*m_semiminor) /
                                    (m_semimajor*m_semimajor);
                double omega_width = std::acos(std::sqrt(arg1/arg2));

                // Continue only if the arclength is positive
                if (omega_width > 0.0) {

                    // Compute azimuth angle difference between ellipse
                    // position angle and position angle of ROI in the model
                    // system. This angle will define the reference point for
                    // the circle arcs. Make sure that omega_0 is within [-pi,pi]
                    // thus that the omega intervals are between [-2pi,2pi]
                    double omega_0 = m_posangle - m_posangle_roi;
                    if (omega_0 > gammalib::pi) {
                        omega_0 -= gammalib::pi;
                    }
                    else if (omega_0 < -gammalib::pi) {
                        omega_0 += gammalib::pi;
                    }

                    // Compute azimuth angle intervals
                    double omega1_min = omega_0    - omega_width;
                    double omega1_max = omega_0    + omega_width;
                    double omega2_min = omega1_min + gammalib::pi;
                    double omega2_max = omega1_max + gammalib::pi;

                    // Limit intervals to the intersection of the ellipse with
                    // the Psf circle. This may lead to a split of intervals,
                    // and we gather all these intervals in a special interval
                    // pair containers
                    cta_omega_intervals intervals1 = 
                        gammalib::limit_omega(omega1_min, omega1_max, domega);
                    cta_omega_intervals intervals2 = 
                        gammalib::limit_omega(omega2_min, omega2_max, domega);

                    // Integrate over all intervals for omega1
                    for (int i = 0; i < intervals1.size(); ++i) {
                        double min = intervals1[i].first;
                        double max = intervals1[i].second;
                        nroi      += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                    // Integrate over all intervals for omega2
                    for (int i = 0; i < intervals2.size(); ++i) {
                        double min = intervals2[i].first;
                        double max = intervals2[i].second;
                        nroi      += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                } // endif: arc length was positive

            } // endelse: circle was not comprised in ellipse

            // Debug: Check for NaN
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
                std::cout << "*** ERROR: cta_nroi_elliptical_kern_rho::eval";
                std::cout << "(rho=" << rho << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (nroi=" << nroi;
                std::cout << ", sin_rho=" << sin_rho;
                std::cout << ", cos_rho=" << cos_rho;
                std::cout << ")" << std::endl;
            }
            #endif

        } // endif: arc length was positive

    } // endif: rho was positive

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Nroi integration of elliptical model
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * Computes
 *
 * \f[
 *    S_{\rm p}(\omega | \rho, E, t) \, N_{\rm pred}(\omega | \rho)
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the centre of the Region of Interest (ROI), and
 * \f$\rho\f$ is the radial distance from the model centre.
 *
 * @todo Npred computation goes over sky coordinates. This can maybe be
 *       optimized to reduce the number of coordinate transformations.
 * @todo Check whether the Npred omega argument is the right one.
 ***************************************************************************/
double cta_nroi_elliptical_kern_omega::eval(const double& omega)
{
    // Initialise Nroi value
    double nroi = 0.0;

    // Compute azimuth angle in model coordinate system (radians)
    double omega_model = omega + m_posangle_roi;

    // Evaluate sky model
    double model = m_model->eval(m_rho, omega_model, m_srcEng, m_srcTime);

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (model <= 0.0) {
        double m_semiminor_rad = m_model->semiminor() * gammalib::deg2rad;
        double m_semimajor_rad = m_model->semimajor() * gammalib::deg2rad;
        double diff_angle      = omega_model - m_model->posangle() * gammalib::deg2rad;
        double cosinus         = std::cos(diff_angle);
        double sinus           = std::sin(diff_angle);
        double arg1            = m_semiminor_rad * cosinus;
        double arg2            = m_semimajor_rad * sinus;
        double r_ellipse       = m_semiminor_rad * m_semimajor_rad /
                                 std::sqrt(arg1*arg1 + arg2*arg2);
        std::cout << "*** WARNING: cta_nroi_elliptical_kern_omega::eval";
        std::cout << " zero model for (rho,omega)=(";
        std::cout << m_rho*gammalib::rad2deg << ",";
        std::cout << omega*gammalib::rad2deg << ")";
        std::cout << " rho-r_ellipse=" << (m_rho-r_ellipse) << " radians";
        std::cout << std::endl;
    }
    #endif

    // Continue only if model is positive
    if (model > 0.0) {

        // Compute sky direction vector in native coordinates
        double  cos_omega = std::cos(omega_model);
        double  sin_omega = std::sin(omega_model);
        GVector native(-cos_omega*m_sin_rho, sin_omega*m_sin_rho, m_cos_rho);

        // Rotate from native into celestial system
        GVector cel = *m_rot * native;

        // Set sky direction
        GSkyDir srcDir;
        srcDir.celvector(cel);

        // Set Photon
        GPhoton photon(srcDir, m_srcEng, m_srcTime);

        // Compute Nroi for this sky direction
        nroi = m_rsp->nirf(photon, m_obsEng, m_obsTime, *m_obs) * model;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
            std::cout << "*** ERROR: cta_nroi_elliptical_kern_omega::eval";
            std::cout << "(omega=" << omega << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (nroi=" << nroi;
            std::cout << ", model=" << model;
            std::cout << ", cos_omega=" << cos_omega;
            std::cout << ", sin_omega=" << sin_omega;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: sky intensity was positive

    // Return Nroi
    return nroi;
}


/*==========================================================================
 =                                                                         =
 =                  Helper class methods for diffuse models                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for IRF offest angle integration of the diffuse source model
 *
 * @param[in] theta Offset angle with respect to observed photon direction
 *                  (radians).
 *
 * Computes
 *
 * \f[
 *    K(\theta | E, t) = \sin \theta \times PSF(\theta)
 *                       \int_{0}^{2\pi}
 *                       S_{\rm p}(\theta, \phi | E, t) \,
 *                       Aeff(\theta, \phi) \,
 *                       Edisp(\theta, \phi) d\phi
 * \f]
 *
 * The PSF is assumed to be azimuthally symmetric, hence the PSF is computed
 * outside the azimuthal integration.
 *
 * Note that the integration is only performed for \f$\theta>0\f$. Otherwise
 * zero is returned.
 ***************************************************************************/
double cta_irf_diffuse_kern_theta::eval(const double& theta)
{
    // Initialise result
    double irf = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Get PSF value. We can do this externally to the azimuthal
        // integration as the PSF is so far azimuthally symmetric. Once
        // we introduce asymmetries, we have to move this done into the
        // Phi kernel method/
        double psf = m_rsp->psf(theta, m_theta, m_phi, m_zenith, m_azimuth, m_srcLogEng);

        // Continue only if PSF is positive
        if (psf > 0.0) {

            // Precompute terms needed for the computation of the angular
            // distance between the true photon direction and the pointing
            // direction (i.e. the camera centre)
            double sin_theta = std::sin(theta);
            double cos_theta = std::cos(theta);
            double sin_ph    = sin_theta * m_sin_eta;
            double cos_ph    = cos_theta * m_cos_eta;

            // Setup kernel for azimuthal integration
            cta_irf_diffuse_kern_phi integrand(m_rsp,
                                               m_model,
                                               m_rot,
                                               m_zenith,
                                               m_azimuth,
                                               m_srcEng,
                                               m_srcTime,
                                               m_srcLogEng,
                                               m_obsEng,
                                               sin_theta,
                                               cos_theta,
                                               sin_ph,
                                               cos_ph);

            // Set number of azimuthal integration iterations
            int iter = gammalib::iter_phi(theta, m_resolution,
                                          m_min_iter, m_max_iter);

            // Integrate over phi
            GIntegral integral(&integrand);
            integral.fixed_iter(iter);
            irf = integral.romberg(0.0, gammalib::twopi) * psf * sin_theta;

            // Compile option: Debugging
            #if defined(G_DEBUG_INTEGRAL)
            if (!integral.isvalid()) {
                std::cout << "cta_irf_diffuse_kern_theta(theta=";
                std::cout << theta*gammalib::rad2deg << " deg) psf=";
                std::cout << psf << ":" << std::endl;
                std::cout << integral.print() << std::endl;
            }
            #endif

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                std::cout << "*** ERROR: cta_irf_diffuse_kern_theta";
                std::cout << "(theta=" << theta << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (irf=" << irf;
                std::cout << ", psf=" << psf;
                std::cout << ", sin_theta=" << sin_theta;
                std::cout << ", cos_theta=" << cos_theta;
                std::cout << ", sin_ph=" << sin_ph;
                std::cout << ", cos_ph=" << cos_ph;
                std::cout << ")";
                std::cout << std::endl;
            }
            #endif

        } // endif: PSF was positive

    } // endif: offset angle was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for IRF azimuth angle integration of the diffuse source model
 *
 * @param[in] phi Azimuth angle around observed photon direction (radians).
 *
 * Computes
 *
 * \f[
 *    S_{\rm p}(\theta, \phi | E, t) \,
 *    Aeff(\theta, \phi) \,
 *    Edisp(\theta, \phi)
 * \f]
 *
 * As the coordinates \f$(\theta, \phi)\f$ are given in the reference frame
 * of the observed photon direction, some coordinate transformations have
 * to be performed.
 *
 * First, \f$(\theta, \phi)\f$ are transformed into the celestial reference
 * frame using the rotation matrix.
 *
 * Then, the offset angle of the true photon direction is computed in the
 * camera system (so far we do not compute the azimuth angle as we assume an
 * azimuthally symmetric response).
 *
 * @todo Optimize computation of sky direction in native coordinates
 * @todo Implement azimuth angle computation of true photon in camera
 * @todo Replace (theta,phi) by (delta,alpha)
 ***************************************************************************/
double cta_irf_diffuse_kern_phi::eval(const double& phi)
{
    // Initialise result
    double irf = 0.0;

    // Compute sine and cosine of azimuth angle
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    // Compute sky direction vector in native coordinates
    m_native[0] = -cos_phi * m_sin_theta;
    m_native[1] =  sin_phi * m_sin_theta;
    m_native[2] = m_cos_theta;

    // Rotate from native into celestial system
    GVector cel = *m_rot * m_native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    m_photon.dir(srcDir);
    m_photon.energy(m_srcEng);
    m_photon.time(m_srcTime);

    // Get sky intensity for photon
    double intensity = m_model->eval(m_photon);

    // Continue only if sky intensity is positive
    if (intensity > 0.0) {

        // Compute true photon offset angle in camera system [radians]
        double offset = std::acos(m_cos_ph + m_sin_ph * cos_phi);

        //TODO: Compute true photon azimuth angle in camera system [radians]
        double azimuth = 0.0;

        // Evaluate model times the effective area
        irf = intensity *
              m_rsp->aeff(offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);

        // Optionally take energy dispersion into account
        if (m_rsp->use_edisp() && irf > 0.0) {
            irf *= m_rsp->edisp(m_obsEng, m_srcEng, offset, azimuth,
                                m_zenith, m_azimuth);
        }

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: cta_irf_diffuse_kern_phi::eval";
            std::cout << "(phi=" << phi << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", intensity=" << intensity;
            std::cout << ", offset=" << offset;
            std::cout << ", azimuth=" << azimuth;
            std::cout << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: sky intensity was positive

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for Nroi offest angle integration of diffuse model
 *
 * @param[in] theta Offset angle with respect to ROI centre (radians).
 *
 * Computes
 *
 * \f[
 *    K(\theta | E, t) = \sin \theta \times
 *                       \int_{0}^{2\pi}
 *                       S_{\rm p}(\theta, \phi | E, t) \,
 *                       N_{\rm pred}(\theta, \phi) d\phi
 * \f]
 * 
 * This method integrates a diffuse model for a given offset angle with
 * respect to the ROI centre over all azimuth angles (from 0 to 2pi). The
 * integration is only performed for positive offset angles, otherwise 0 is
 * returned.
 *
 * Integration is done using the Romberg algorithm. The integration kernel
 * is defined by the helper class cta_npred_diffuse_kern_phi.
 *
 * Note that the integration precision was adjusted trading-off between
 * computation time and computation precision. A value of 1e-4 was judged
 * appropriate.
 ***************************************************************************/
double cta_nroi_diffuse_kern_theta::eval(const double& theta)
{
    // Initialise Nroi value
    double nroi = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Compute sine of offset angle
        double sin_theta = std::sin(theta);

        // Setup phi integration kernel
        cta_nroi_diffuse_kern_phi integrand(m_rsp,
                                            m_obs,
                                            m_model,
                                            m_rot,
                                            m_srcEng,
                                            m_srcTime,
                                            m_obsEng,
                                            m_obsTime,
                                            theta,
                                            sin_theta);

        // Integrate over phi
        GIntegral integral(&integrand);
        integral.fixed_iter(m_iter);
        nroi = integral.romberg(0.0, gammalib::twopi, m_iter) * sin_theta;
        #if defined(G_DEBUG_INTEGRAL)
        if (!integral.isvalid()) {
            std::cout << "cta_nroi_diffuse_kern_theta(theta=";
            std::cout << theta*gammalib::rad2deg << " deg):" << std::endl;
            std::cout << integral.print() << std::endl;
        }
        #endif

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
            std::cout << "*** ERROR: cta_nroi_radial_kern_theta::eval";
            std::cout << "(theta=" << theta << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (nroi=" << nroi;
            std::cout << ", sin_theta=" << sin_theta;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: offset angle was positive

    // Return Nroi
    return nroi;
}


/***********************************************************************//**
 * @brief Kernel for Nroi azimuth angle integration of diffuse model
 *
 * @param[in] phi Azimuth angle with respect to ROI centre (radians).
 *
 * Computes
 *
 * \f[
 *    S_{\rm p}(\theta, \phi | E, t) \, N_{\rm pred}(\theta, \phi)
 * \f]
 *
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double cta_nroi_diffuse_kern_phi::eval(const double& phi)
{
    // Initialise Nroi value
    double nroi = 0.0;

    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set Photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Get sky intensity for this photon
    double intensity = m_model->eval(photon);

    // Continue only if sky intensity is positive
    if (intensity > 0.0) {

        // Compute Nroi for this sky direction
        nroi = m_rsp->nirf(photon, m_obsEng, m_obsTime, *m_obs) * intensity;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(nroi) || gammalib::is_infinite(nroi)) {
            std::cout << "*** ERROR: cta_nroi_diffuse_kern_phi::eval";
            std::cout << "(phi=" << phi << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (nroi=" << nroi;
            std::cout << ", intensity=" << intensity;
            std::cout << ", cos_phi=" << cos_phi;
            std::cout << ", sin_phi=" << sin_phi;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: sky intensity was positive

    // Return Nroi
    return nroi;
}


/*==========================================================================
 =                                                                         =
 =                Helper class methods for stacked analysis                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for radial model integration over zenith angle
 *
 * @param[in] rho Radial distance from model centre (radians).
 * @return Integration kernel.
 *
 * Computes the integration kernel
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega} PSF(\rho, \omega) d\omega
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 * \f$S_{\rm p}(\rho | E, t)\f$ is the radial source model
 * for a given true photon energy and photon arrival time,
 * \f$PSF(\rho, \omega)\f$ is the point spread function.
 ***************************************************************************/
double cta_psf_radial_kern_rho::eval(const double& rho)
{
    // Initialise result
    double irf = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of the arc (in radians) from a circle with
        // radius rho that intersects with the point spread function, defined
        // as a circle with maximum radius m_delta_max
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_rho_obs,
                                                      m_cos_rho_obs,
                                                      m_sin_rho_obs,
                                                      m_delta_max,
                                                      m_cos_delta_max);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model
            double rho_kludge = rho - g_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Evaluate sky model
            double model = m_model->eval(rho_kludge, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_psf_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model=" << (rho-m_model->theta_max());
                std::cout << " radians" << std::endl;
            }
            #endif

            // Continue only if model is positive
            if (model > 0.0) {

                // Compute omega integration range
                double omega_min = -domega;
                double omega_max = +domega;

                // Precompute cosine and sine terms for azimuthal integration
                double cos_rho = std::cos(rho);
                double sin_rho = std::sin(rho);
                double cos_psf = cos_rho * m_cos_rho_obs;
                double sin_psf = sin_rho * m_sin_rho_obs;

                // Setup integration kernel
                cta_psf_radial_kern_omega integrand(m_rsp,
                                                    m_model,
                                                    m_srcDir,
                                                    m_srcEng,
                                                    m_srcTime,
                                                    cos_psf,
                                                    sin_psf);

                // Setup integrator
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);

                // Integrate over omega
                irf = integral.romberg(omega_min, omega_max, m_iter) *
                      model * sin_rho;

                // Compile option: Check for NaN/Inf
                #if defined(G_NAN_CHECK)
                if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                    std::cout << "*** ERROR: cta_psf_radial_kern_rho";
                    std::cout << "(rho=" << rho << "):";
                    std::cout << " NaN/Inf encountered";
                    std::cout << " (irf=" << irf;
                    std::cout << ", domega=" << domega;
                    std::cout << ", model=" << model;
                    std::cout << ", sin_rho=" << sin_rho << ")";
                    std::cout << std::endl;
                }
                #endif

            } // endif: model was positive

        } // endif: arclength was positive

    } // endif: rho was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for radial model integration over azimuth angle
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * Computes
 *
 * \f[
 *    K(\omega | \rho, E, t) = PSF(\omega | \rho)
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 *
 * From the coordinates \f$(\rho,\omega)\f$ in the model system, the method
 * computes the angle between the true (\f$\vec{p}\f$) and observed
 * (\f$\vec{p'}\f$) photon arrival direction using
 * 
 * \f[\delta = \arccos(\cos \rho \cos \rho_{\rm obs} + 
 *                     \sin \rho \sin \rho_{\rm obs} \cos \omega)\f]
 *
 * where
 * \f$\rho_{\rm obs}\f$ is the angular distance between the observed photon
 * arrival direction \f$\vec{p'}\f$ and the model centre \f$\vec{m}\f$.
 * \f$\delta\f$ is used to compute the value of the point spread function.
 ***************************************************************************/
double cta_psf_radial_kern_omega::eval(const double& omega)
{
    // Compute Psf offset angle (radians)
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));

    // Evaluate Psf * model for this delta
    double irf = m_rsp->psf()(m_srcDir, delta, m_srcEng);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: cta_psf_radial_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", delta=" << delta;
        std::cout << ", omega=" << omega << ")";
        std::cout << std::endl;
    }
    #endif

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for PSF integration of radial model
 *
 * @param[in] delta PSF offset angle (radians).
 * @return Azimuthally integrated product between PSF and radial model.
 *
 * Computes the azimuthally integrated product of point spread function and
 * the radial model intensity. As the PSF is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated model. The method returns thus
 *
 * \f[
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm M}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm M}(\delta, \phi)\f$ is the radial model in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 ***************************************************************************/
double cta_psf_radial_kern_delta::eval(const double& delta)
{
    // Initialise value
    double value = 0.0;

    // If we're at the Psf peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas.
    if (delta > 0.0) {

        // Get Psf for this delta
        double psf = m_rsp->psf()(m_srcDir, delta, m_srcEng);

        // Continue only if Psf is positive
        if (psf > 0.0) {

            // Compute half length of the arc (in radians) from a circle with
            // radius delta that intersects with the model, defined as a circle
            // with maximum radius m_theta_max
            double dphi = 0.5 * gammalib::roi_arclength(delta,
                                                        m_delta_mod,
                                                        m_cos_delta_mod,
                                                        m_sin_delta_mod,
                                                        m_theta_max,
                                                        m_cos_theta_max);

            // Continue only if arc length is positive
            if (dphi > 0.0) {

                // Compute phi integration range
                double phi_min = -dphi;
                double phi_max = +dphi;

                // Precompute cosine and sine terms for azimuthal integration
                double sin_delta = std::sin(delta);
                double cos_delta = std::cos(delta);
                double sin_fact  = sin_delta * m_sin_delta_mod;
                double cos_fact  = cos_delta * m_cos_delta_mod;

                // Setup kernel for azimuthal integration of the spatial model
                cta_psf_radial_kern_phi integrand(m_model,
                                                  m_srcEng,
                                                  m_srcTime,
                                                  sin_fact,
                                                  cos_fact);

                // Setup integrator
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);

                // Integrate over azimuth
                value = integral.romberg(phi_min, phi_max, m_iter) *
                        psf * sin_delta;

                // Debug: Check for NaN
                #if defined(G_NAN_CHECK)
                if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
                    std::cout << "*** ERROR: cta_psf_radial_kern_delta::eval";
                    std::cout << "(delta=" << delta << "):";
                    std::cout << " NaN/Inf encountered";
                    std::cout << " (value=" << value;
                    std::cout << ")" << std::endl;
                }
                #endif

            } // endif: Psf value was positive

        } // endif: arc length was positive

    } // endif: delta was positive

    // Return kernel value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for azimuthal radial model integration
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Radial model value.
 *
 * Computes the value of the radial model at the position \f$(\delta,\phi)\f$
 * given in point spread function coordinates. The \f$\theta\f$ angle of the
 * radial model is computed using
 *
 * \f[
 *    \theta = \arccos \left( \cos \delta \cos \zeta +
 *                            \sin \delta \sin \zeta \cos \phi \right)
 * \f]
 *
 * where \f$\delta\f$ is the angle between true and measured photon
 * direction, \f$\zeta\f$ is the angle between model centre and measured
 * photon direction, and \f$\phi\f$ is the azimuth angle with respect to the
 * measured photon direction, where \f$\phi=0\f$ corresponds to the 
 * connecting line between model centre and measured photon direction.
 ***************************************************************************/
double cta_psf_radial_kern_phi::eval(const double& phi)
{
    // Compute radial model theta angle
    double theta = std::acos(m_cos_fact + m_sin_fact * std::cos(phi));

    // Reduce theta by an infinite amount to avoid rounding errors at the
    // boundary of a sharp edged model
    double theta_kluge = theta - 1.0e-12;
    if (theta_kluge < 0.0) {
        theta_kluge = 0.0;
    }

    // Get radial model value
    double value = m_model->eval(theta_kluge, m_srcEng, m_srcTime); 

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (value <= 0.0) {
        std::cout << "*** WARNING: cta_psf_radial_kern_phi::eval";
        std::cout << " zero model for (phi)=(";
        std::cout << phi*gammalib::rad2deg << ")";
        std::cout << " theta-r_model=" << (theta-m_model->theta_max());
        std::cout << " radians" << std::endl;
    }
    #endif

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_psf_radial_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
}


/***********************************************************************//**
 * @brief Return size of vector kernel for PSF integration of radial model
 *
 * @return Size of vector kernel.
 *
 * The size of the vector kernel depends on whether gradient computation is
 * requested or not. With gradient computation, the size is given by the
 * product of the number of true energies and the number of model parameters
 * plus one, without gradient computation the size is simply the number of
 * true energies.
 ***************************************************************************/
int cta_psf_radial_kerns_delta::size(void) const
{
    // Set size
    int nengs = m_srcEngs.size();
    int npars = m_model->size();
    int size  = (m_grad) ? nengs * (npars+1) : nengs;

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Kernel constructor for PSF integration of radial model
 *
 * @param[in] delta PSF offset angle (radians).
 ***************************************************************************/
cta_psf_radial_kerns_delta::cta_psf_radial_kerns_delta(const GCTAResponseCube*    rsp,
                                                       const GModelSpatialRadial* model,
                                                       const GSkyDir&             obsDir,
                                                       const GEnergies&           srcEngs,
                                                       const double&              zeta,
                                                       const double&              theta_max,
                                                       const int&                 iter,
                                                       const bool&                grad)
{
    // Store constructor arguments
    m_rsp       = rsp;
    m_model     = model;
    m_obsDir    = obsDir;
    m_srcEngs   = srcEngs;
    m_zeta      = zeta;
    m_theta_max = theta_max;
    m_iter      = iter;
    m_grad      = grad;

    // Pre-compute trigonometric functions
    m_cos_zeta      = std::cos(zeta);
    m_sin_zeta      = std::sin(zeta);
    m_cos_theta_max = std::cos(theta_max);

    // Store references to Right Ascension and Declination parameters
    m_par_ra  = &((*(const_cast<GModelSpatialRadial*>(model)))["RA"]);
    m_par_dec = &((*(const_cast<GModelSpatialRadial*>(model)))["DEC"]);

    // Pre-compute constants needed for gradient computation. The following
    // exceptions need to be handled
    // * beta_0 = +/-90 deg
    // * m_sin_zeta = 0
    // * denom = 0
    if (grad) {
        double alpha_0        = model->ra()  * gammalib::deg2rad;
        double beta_0         = model->dec() * gammalib::deg2rad;
        double alpha_reco     = obsDir.ra();
        double beta_reco      = obsDir.dec();
        double sin_beta_0     = std::sin(beta_0);
        double cos_beta_0     = std::cos(beta_0);
        double tan_beta_0     = std::tan(beta_0);  // Exception: beta_0 = 90 deg
        double sin_beta_reco  = std::sin(beta_reco);
        double cos_beta_reco  = std::cos(beta_reco);
        double sin_dalpha     = std::sin(alpha_0 - alpha_reco);
        double cos_dalpha     = std::cos(alpha_0 - alpha_reco);
        double arg            = cos_beta_reco * tan_beta_0 -
                                sin_beta_reco * cos_dalpha;
        double denom          = sin_dalpha * sin_dalpha + arg * arg;
        m_dzeta_dalpha_0      = cos_beta_0 * cos_beta_reco * sin_dalpha / m_sin_zeta;
        m_dzeta_dbeta_0       = (sin_beta_0 * cos_beta_reco * cos_dalpha -
                                 cos_beta_0 * sin_beta_reco) / m_sin_zeta;
        m_dphi_dalpha_0       = (sin_beta_reco - cos_dalpha * cos_beta_reco * tan_beta_0) /
                                denom;
        m_dphi_dbeta_0        = (sin_dalpha * cos_beta_reco * (1.0 + tan_beta_0 * tan_beta_0)) /
                                denom;
    }
    else {
        m_dzeta_dalpha_0 = 0.0;
        m_dzeta_dbeta_0  = 0.0;
        m_dphi_dalpha_0  = 0.0;
        m_dphi_dbeta_0   = 0.0;
    }

    // Return
    return;
}

/***********************************************************************//**
 * @brief Kernel for PSF integration of radial model
 *
 * @param[in] delta PSF offset angle (radians).
 * @return Azimuthally integrated product between PSF and radial model
 *         values for all energies.
 *
 * Computes the azimuthally integrated product of point spread function and
 * the radial model intensity. As the PSF is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated model. The method returns thus
 *
 * \f[
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm M}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm M}(\delta, \phi)\f$ is the radial model in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 ***************************************************************************/
GVector cta_psf_radial_kerns_delta::eval(const double& delta)
{
    // Determine size of the result vector
    int npars = m_model->size();
    int nengs = m_srcEngs.size();

    // Initialise result vector
    GVector values(size());

    // If we're at the Psf peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas.
    if (delta > 0.0) {

        // Compute half length of the arc (in radians) from a circle with
        // radius delta that intersects with the model, defined as a circle
        // with maximum radius m_theta_max
        double dphi = 0.5 * gammalib::roi_arclength(delta,
                                                    m_zeta,
                                                    m_cos_zeta,
                                                    m_sin_zeta,
                                                    m_theta_max,
                                                    m_cos_theta_max);

        // Continue only if arc length is positive
        if (dphi > 0.0) {

            // Compute phi integration range
            double phi_min = -dphi;
            double phi_max = +dphi;

            // Precompute cosine and sine terms for azimuthal integration
            double sin_delta = std::sin(delta);
            double cos_delta = std::cos(delta);
            double sin_delta_sin_zeta = sin_delta * m_sin_zeta;
            double sin_delta_cos_zeta = sin_delta * m_cos_zeta;
            double cos_delta_sin_zeta = cos_delta * m_sin_zeta;
            double cos_delta_cos_zeta = cos_delta * m_cos_zeta;

            // Setup kernel for azimuthal integration of the spatial model
            cta_psf_radial_kerns_phi integrand(this,
                                               sin_delta_sin_zeta,
                                               sin_delta_cos_zeta,
                                               cos_delta_sin_zeta,
                                               cos_delta_cos_zeta);

            // Setup integrator
            GIntegrals integral(&integrand);
            integral.fixed_iter(m_iter);

            // Integrate over azimuth
            GVector irf = integral.romberg(phi_min, phi_max, m_iter) * sin_delta;

            // Extract value and gradients
            double value = irf[0];

            // Multiply in energy dependent Psf
            for (int i = 0; i < nengs; ++i) {

                // Get Psf for this energy. We approximate here the true sky
                // direction by the reconstructed sky direction.
                double psf = m_rsp->psf()(m_obsDir, delta, m_srcEngs[i]);

                // Continue only if Psf is positive
                if (psf > 0.0) {

                    // Compute value
                    values[i] = value * psf;

                    // If gradient computation is requested the multiply also
                    // all factor gradients with the Psf value
                    if (m_grad) {
                        for (int k = 1, ig = i+nengs; k <= npars; ++k, ig += nengs) {
                            values[ig] = irf[k] * psf;
                        }
                    }

                } // endif: Psf was positive

            } // endfor: looped over energies

        } // endif: arc length was positive

    } // endif: delta was positive

    // Return kernel values
    return values;
}


/***********************************************************************//**
 * @brief Kernel for azimuthal radial model integration
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Vector of radial model value and parameter factor gradients.
 *
 * Computes the value and parameter factor gradients of the radial model at
 * the position \f$(\delta,\phi)\f$ given in point spread function
 * coordinates.
 *
 * The angle \f$\theta\f$ of the radial model is computed using
 *
 * \f[
 *    \theta = \arccos \left( \cos \delta \cos \zeta +
 *                            \sin \delta \sin \zeta \cos \phi \right)
 * \f]
 *
 * where \f$\delta\f$ is the angle between true and measured photon
 * direction, \f$\zeta\f$ is the angle between model centre and measured
 * photon direction, and \f$\phi\f$ is the azimuth angle with respect to the
 * measured photon direction, where \f$\phi=0\f$ corresponds to the
 * connecting line between model centre and measured photon direction.
 *
 * If gradient computation is requested, the method also returns the gradients
 *
 * \f[
 *    \frac{\partial f}{\partial \alpha_0}
 * \f]
 *
 * and
 *
 * \f[
 *    \frac{\partial f}{\partial \beta_0}
 * \f]
 *
 * in the second and third slot of the returned vector. The following slots
 * contains gradients with respect to other model parameters.
 ***************************************************************************/
GVector cta_psf_radial_kerns_phi::eval(const double& phi)
{
    // Allocate dummy energy and time to satisfy interface
    static const GEnergy srcEng;
    static const GTime   srcTime;

    // Compute sin(phi) and cos(phi). The following code allows usage of the
    // sincos() method in case that gradients are requested.
    double sin_phi;
    double cos_phi;
    if (m_outer->m_grad) {
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);
    }
    else {
        cos_phi = std::cos(phi);
    }

    // Compute radial model theta angle
    double cos_theta = m_cos_delta_cos_zeta + m_sin_delta_sin_zeta * cos_phi;
    double theta     = std::acos(cos_theta);

    // Reduce theta by an infinite amount to avoid rounding errors at the
    // boundary of a sharp edged model
    double theta_kluge = theta - 1.0e-12;
    if (theta_kluge < 0.0) {
        theta_kluge = 0.0;
    }

    // Initialise kernel values
    //GVector values(m_size);

    // If gradients are requested then compute partial derivatives of theta
    // and phi with respect to Right Ascension and Declination and store them
    // into the respective factor gradients so that the
    // GModelSpatialRadial::eval() method can use the information.
    if (m_outer->m_grad) {
        double g_ra  = 0.0;
        double g_dec = 0.0;
        if (theta > 0.0) {

            // Compute sin(theta) by avoiding to use a call to sine
            //double sin_theta = std::sin(theta);
            double sin_theta2 = 1.0 - cos_theta * cos_theta;
            double sin_theta  = (sin_theta2 > 0.0) ? std::sqrt(sin_theta2) : 0.0;

            // Continue only is sine is non-zero
            if (sin_theta != 0.0) {
                double norm         = 1.0 / sin_theta;
                double dtheta_dzeta = (m_cos_delta_sin_zeta -
                                       m_sin_delta_cos_zeta * cos_phi) *
                                       norm;
                double dtheta_dphi  = (m_sin_delta_sin_zeta * sin_phi) *
                                       norm;
                g_ra  = dtheta_dzeta * m_outer->m_dzeta_dalpha_0 +
                        dtheta_dphi  * m_outer->m_dphi_dalpha_0;
                g_dec = dtheta_dzeta * m_outer->m_dzeta_dbeta_0  +
                        dtheta_dphi  * m_outer->m_dphi_dbeta_0;
            }
        }
        m_outer->m_par_ra->factor_gradient(g_ra);
        m_outer->m_par_dec->factor_gradient(g_dec);
    }

    // Compute model value and optionally model parameter gradients
    //values[0] = m_outer->m_model->eval(theta_kluge, srcEng, srcTime,
    //                                   m_outer->m_grad);
    m_values[0] = m_outer->m_model->eval(theta_kluge, srcEng, srcTime,
                                         m_outer->m_grad);

    // If gradients are requested, extract now the fully computed model
    // parameter gradients into the kernel values vector
    if (m_outer->m_grad) {
        for (int i = 1, k = 0; i < m_size; ++i, ++k) {
            //values[i] = (*(m_outer->m_model))[k].factor_gradient();
            m_values[i] = (*(m_outer->m_model))[k].factor_gradient();
        }
    }

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
//    if (gammalib::is_notanumber(values[0]) || gammalib::is_infinite(values[0])) {
    if (gammalib::is_notanumber(m_values[0]) || gammalib::is_infinite(m_values[0])) {
        std::string msg = "NaN/Inf encountered for phi="+gammalib::str(phi);
        gammalib::warning(G_PSF_RADIAL_KERNS_PHI, msg);
    }
    #endif

    // Return kernel values
    //return values;
    return m_values;
}


/***********************************************************************//**
 * @brief Kernel for elliptical model integration over zenith angle
 *
 * @param[in] rho Radial distance from model centre (radians).
 * @return Integration kernel.
 *
 * Computes the integration kernel
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega} 
 *                     S_{\rm p}(\rho, \omega | E, t) \, PSF(\rho, \omega)
 *                     d\omega
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical source model
 * for a given true photon energy and photon arrival time,
 * \f$PSF(\rho, \omega)\f$ is the point spread function.
 * The integration is over the \f$\omega\f$ values that are comprised within
 * the ellipse boundaries.
 ***************************************************************************/
double cta_psf_elliptical_kern_rho::eval(const double& rho)
{
    // Initialise result
    double irf = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of the arc (in radians) from a circle with
        // radius rho that intersects with the point spread function, defined
        // as a circle with maximum radius m_delta_max
        double domega = 0.5 * gammalib::roi_arclength(rho,
                                                      m_rho_obs,
                                                      m_cos_rho_obs,
                                                      m_sin_rho_obs,
                                                      m_delta_max,
                                                      m_cos_delta_max);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Precompute cosine and sine terms for azimuthal integration
            double cos_rho = std::cos(rho);
            double sin_rho = std::sin(rho);
            double cos_psf = cos_rho * m_cos_rho_obs;
            double sin_psf = sin_rho * m_sin_rho_obs;

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model (e.g. an elliptical
            // disk model)
            double rho_kludge = rho - g_ellipse_kludge_radius;
            if (rho_kludge < 0.0) {
                rho_kludge = 0.0;
            }

            // Setup integration kernel
            cta_psf_elliptical_kern_omega integrand(m_rsp,
                                                    m_model,
                                                    m_srcDir,
                                                    m_srcEng,
                                                    m_srcTime,
                                                    m_posangle_obs,
                                                    rho_kludge,
                                                    cos_psf,
                                                    sin_psf);

            // Setup integrator
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);

            // If the radius rho is not larger than the semiminor axis
            // boundary, the circle with that radius is fully contained in
            // the ellipse and we can just integrate over the relevant arc
            if (rho <= m_semiminor) {

                // Compute omega integration range
                double omega_min = -domega;
                double omega_max = +domega;

                // Integrate over omega
                irf = integral.romberg(omega_min, omega_max, m_iter) *
                      sin_rho;

            } // endif: circle comprised in ellipse

            // ... otherwise there are arcs that intersect with the Psf circle
            else {

                // Compute half the arc length (in radians) of a circle of
                // radius rho, centred on the model, that intersects with
                // the ellipse boundary
                double arg1 = 1.0 - (m_semiminor*m_semiminor) / (rho*rho);
                double arg2 = 1.0 - (m_semiminor*m_semiminor) /
                                    (m_semimajor*m_semimajor);
                double omega_width = std::acos(std::sqrt(arg1/arg2));

                // Continue only if the arclength is positive
                if (omega_width > 0.0) {

                    // Compute azimuth angle difference between ellipse
                    // position angle and position angle of observed
                    // photon in the model system. This angle will define
                    // the reference point around which the circle arc. Make
                    // sure that omega_0 is within [-pi,pi] thus that the omega
                    // intervals are between [-2pi,2pi]
                    double omega_0 = m_posangle - m_posangle_obs;
                    if (omega_0 > gammalib::pi) {
                        omega_0 -= gammalib::pi;
                    }
                    else if (omega_0 < -gammalib::pi) {
                        omega_0 += gammalib::pi;
                    }

                    // Compute azimuth angle intervals
                    double omega1_min = omega_0    - omega_width;
                    double omega1_max = omega_0    + omega_width;
                    double omega2_min = omega1_min + gammalib::pi;
                    double omega2_max = omega1_max + gammalib::pi;

                    // Limit intervals to the intersection of the ellipse with
                    // the Psf circle. This may lead to a split of intervals,
                    // and we gather all these intervals in a special interval
                    // pair containers
                    cta_omega_intervals intervals1 = 
                        gammalib::limit_omega(omega1_min, omega1_max, domega);
                    cta_omega_intervals intervals2 = 
                        gammalib::limit_omega(omega2_min, omega2_max, domega);

                    // Integrate over all intervals for omega1
                    for (int i = 0; i < intervals1.size(); ++i) {
                        double min = intervals1[i].first;
                        double max = intervals1[i].second;
                        irf       += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                    // Integrate over all intervals for omega2
                    for (int i = 0; i < intervals2.size(); ++i) {
                        double min = intervals2[i].first;
                        double max = intervals2[i].second;
                        irf       += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                } // endif: arc length was positive

            } // endelse: circle was not comprised in ellipse

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                std::cout << "*** ERROR: cta_psf_elliptical_kern_rho";
                std::cout << "(rho=" << rho << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (irf=" << irf;
                std::cout << ", domega=" << domega;
                std::cout << ", sin_rho=" << sin_rho << ")";
                std::cout << std::endl;
            }
            #endif

        } // endif: arc length was positive

    } // endif: rho was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for elliptical model integration over azimuth angle
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * Computes
 *
 * \f[
 *    K(\omega | \rho, E, t) = S_{\rm p}(\omega | \rho, E, t) \,
 *                             PSF(\omega | \rho)
 * \f]
 *
 * where
 * \f$\omega\f$ is the azimuth angle with respect to the model centre,
 * counted counterclockwise from the vector connecting the model centre
 * to the observed photon direction, and
 * \f$\rho\f$ is the radial distance from the model centre.
 *
 * From the coordinates \f$(\rho,\omega)\f$ in the model system, the method
 * computes the angle between the true (\f$\vec{p}\f$) and observed
 * (\f$\vec{p'}\f$) photon arrival direction using
 * 
 * \f[\delta = \arccos(\cos \rho \cos \rho_{\rm obs} + 
 *                     \sin \rho \sin \rho_{\rm obs} \cos \omega)\f]
 *
 * where
 * \f$\rho_{\rm obs}\f$ is the angular distance between the observed photon
 * arrival direction \f$\vec{p'}\f$ and the model centre \f$\vec{m}\f$.
 * \f$\delta\f$ is used to compute the value of the point spread function.
 ***************************************************************************/
double cta_psf_elliptical_kern_omega::eval(const double& omega)
{
    // Initialise Irf value
    double irf = 0.0;

    // Compute azimuth angle in model coordinate system (radians)
    double omega_model = omega + m_posangle_obs;

    // Evaluate sky model
    double model = m_model->eval(m_rho, omega_model, m_srcEng, m_srcTime);

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (model <= 0.0) {
        double semiminor_rad = m_model->semiminor() * gammalib::deg2rad;
        double semimajor_rad = m_model->semimajor() * gammalib::deg2rad;
        double diff_angle    = omega_model - m_model->posangle() * gammalib::deg2rad;
        double cosinus       = std::cos(diff_angle);
        double sinus         = std::sin(diff_angle);
        double arg1          = semiminor_rad * cosinus;
        double arg2          = semimajor_rad * sinus;
        double r_ellipse     = semiminor_rad * semimajor_rad /
                               std::sqrt(arg1*arg1 + arg2*arg2);
        std::cout << "*** WARNING: cta_psf_elliptical_kern_omega::eval";
        std::cout << " zero model for (rho,omega)=";
        std::cout << m_rho*gammalib::rad2deg << ",";
        std::cout << omega*gammalib::rad2deg << ")";
        std::cout << " rho-r_ellipse=" << (m_rho-r_ellipse) << " radians";
        std::cout << std::endl;
    }
    #endif

    // Continue only if model is positive
    if (model > 0.0) {

        // Compute Psf offset angle (radians)
        double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));

        // Evaluate Psf * model for this delta
        irf = m_rsp->psf()(m_srcDir, delta, m_srcEng) * model;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
            std::cout << "*** ERROR: cta_psf_elliptical_kern_omega::eval";
            std::cout << "(omega=" << omega << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (irf=" << irf;
            std::cout << ", model=" << model;
            std::cout << ", delta=" << delta;
            std::cout << ", rho=" << m_rho;
            std::cout << ", omega=" << omega << ")";
            std::cout << std::endl;
        }
        #endif

    } // endif: model is positive

    // Return
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for PSF integration of spatial model
 *
 * @param[in] delta PSF offset angle (radians).
 * @return Azimuthally integrated product between PSF and model.
 *
 * Computes the azimuthally integrated product of point spread function and
 * the spatial model intensity. As the PSF is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated model. The method returns thus
 *
 * \f[
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm M}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm M}(\delta, \phi)\f$ is the spatial model in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 ***************************************************************************/
double cta_psf_diffuse_kern_delta::eval(const double& delta)
{
    // Initialise value
    double value = 0.0;

    // If we're at the PSF peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas.
    if (delta > 0.0) {

        // Get PSF for this delta
        value = m_rsp->psf()(m_srcDir, delta, m_srcEng);

        // Continue only if PSF is positive
        if (value > 0.0) {

            // Compute sine and cosine of delta
            double sin_delta = std::sin(delta);
            double cos_delta = std::cos(delta);

            // Setup kernel for azimuthal integration of the diffuse model
            cta_psf_diffuse_kern_phi integrand(m_model, m_rot,
                                               m_srcEng, m_srcTime,
                                               sin_delta, cos_delta);

            // Set number of azimuthal integration iterations
            int iter = gammalib::iter_phi(delta, m_resolution,
                                          m_min_iter, m_max_iter);

            // Setup untegration
            GIntegral integral(&integrand);

            // Set fixed number of iterations
            integral.fixed_iter(iter);

            // Azimuthally integrate model
            value *= integral.romberg(0.0, gammalib::twopi) * sin_delta;

        } // endif: PSF value was positive

    } // endif: delta was positive

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_psf_diffuse_kern_delta::eval";
        std::cout << "(delta=" << delta << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
}


/***********************************************************************//**
 * @brief Kernel for map integration of spatial model
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Spatial model value.
 *
 * Computes the value of the spatial model at the position (delta,phi) given
 * in point spread function coordinates. The transformation from point
 * spread function coordinates into sky coordinates is done using a rotation
 * matrix that is pre-computed on entry.
 ***************************************************************************/
double cta_psf_diffuse_kern_phi::eval(const double& phi)
{
    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_delta, sin_phi*m_sin_delta, m_cos_delta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute map value for this sky direction
    double value = m_model->eval(photon); 

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_psf_diffuse_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
}
