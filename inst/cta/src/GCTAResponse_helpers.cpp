/***************************************************************************
 *         GCTAResponse_helpers.cpp - CTA response helper classes          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @brief CTA response hepler classes implementation
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
#include "GVector.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAEdisp.hpp"
#include "GCTASupport.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_USE_OBSDIR_FOR_AEFF    //!< Use event offset for Aeff computation

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_INTEGRAL                             //!< Debug integration
#define G_DEBUG_MODEL_ZERO           //!< Debug check for zero model values

/* __ Constants __________________________________________________________ */


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
    double phi = gammalib::cta_roi_arclength(delta,
                                             m_psf,
                                             m_cospsf,
                                             m_sinpsf,
                                             m_roi,
                                             m_cosroi);

    // If arclength is positive then compute the PSF value
    if (phi > 0) {
    
        // Compute PSF value
        value = m_rsp.psf(delta, m_theta, m_phi, m_zenith, m_azimuth, m_logE) * 
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
 * @brief Integration kernel for nedisp() method
 *
 * @param[in] logEobs Observed event energy.
 * @return Energy dispersion PDF value.
 ***************************************************************************/
double cta_nedisp_kern::eval(const double& logEobs)
{
    // Get value
    double value = (*m_rsp.edisp())(logEobs,
                                    m_logEsrc,
                                    m_theta,
                                    m_phi,
                                    m_zenith,
                                    m_azimuth);

    // Correct for variable substitution
    //value *= std::exp(logEobs);

    // Return
    return value;
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
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
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

            // Evaluate sky model
            double model = m_model.eval(rho, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_irf_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model" << (rho-m_model.theta_max());
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
                                                    m_srcLogEng,
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
    // Compute PSF offset angle [radians]
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));
    
    // Compute true photon offset angle in camera system [radians]
    double offset = std::acos(m_cos_ph + m_sin_ph * std::cos(m_omega0 - omega));
    
    //TODO: Compute true photon azimuth angle in camera system [radians]
    double azimuth = 0.0;

    // Evaluate IRF
    double irf = m_rsp.aeff(offset, azimuth, m_zenith, m_azimuth, m_srcLogEng) *
                 m_rsp.psf(delta, offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);

    // Optionally take energy dispersion into account
    if (m_rsp.use_edisp() && irf > 0.0) {
        irf *= m_rsp.edisp(m_obsEng, offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);
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
 * @brief Kernel for zenith angle Npred integration or radial model
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
double cta_npred_radial_kern_rho::eval(const double& rho)
{
    // Initialise Npred value
    double npred = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of arc that lies within ROI+PSF radius (radians)
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
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

            // Get radial model value
            double model = m_model.eval(rho, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_npred_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model" << (rho-m_model.theta_max());
                std::cout << " radians" << std::endl;
            }
            #endif

            // Continue only if we have a positive model value
            if (model > 0.0) {

                // Compute sine and cosine of offset angle
                double sin_rho = std::sin(rho);
                double cos_rho = std::cos(rho);

                // Setup phi integration kernel
                cta_npred_radial_kern_omega integrand(m_rsp,
                                                      m_srcEng,
                                                      m_srcTime,
                                                      m_obs,
                                                      m_rot,
                                                      sin_rho,
                                                      cos_rho);

                // Integrate over phi
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);
                npred = integral.romberg(omega_min, omega_max, m_iter) *
                        sin_rho * model;

                // Debug: Check for NaN
                #if defined(G_NAN_CHECK)
                if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
                    std::cout << "*** ERROR: cta_npred_radial_kern_rho::eval";
                    std::cout << "(rho=" << rho << "):";
                    std::cout << " NaN/Inf encountered";
                    std::cout << " (npred=" << npred;
                    std::cout << ", model=" << model;
                    std::cout << ", omega=[" << omega_min << "," << omega_max << "]";
                    std::cout << ", sin_rho=" << sin_rho;
                    std::cout << ")" << std::endl;
                }
                #endif

            } // endif: model was positive

        } // endif: arc length was positive

    } // endif: rho was positive

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of radial model
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double cta_npred_radial_kern_omega::eval(const double& omega)
{
    // Compute sky direction vector in native coordinates
    double  cos_omega = std::cos(omega);
    double  sin_omega = std::sin(omega);
    GVector native(-cos_omega*m_sin_rho, sin_omega*m_sin_rho, m_cos_rho);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set Photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute point source Npred for this sky direction
    double npred = m_rsp.npred(photon, m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
        std::cout << "*** ERROR: cta_npred_radial_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", cos_omega=" << cos_omega;
        std::cout << ", sin_omega=" << sin_omega;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
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
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
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
            double rho_kluge = rho - 1.0e-12;
            if (rho_kluge < 0.0) {
                rho_kluge = 0.0;
            }

            // Setup integration kernel
            cta_irf_elliptical_kern_omega integrand(m_rsp,
                                                    m_model,
                                                    m_zenith,
                                                    m_azimuth,
                                                    m_srcEng,
                                                    m_srcTime,
                                                    m_srcLogEng,
                                                    m_obsEng,
                                                    m_posangle_obs,
                                                    m_omega_pnt,
                                                    rho_kluge,
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
            if (rho < m_semiminor) {

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
    double model = m_model.eval(m_rho, omega_model, m_srcEng, m_srcTime);

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (model <= 0.0) {
        double m_semiminor_rad = m_model.semiminor() * gammalib::deg2rad;
        double m_semimajor_rad = m_model.semimajor() * gammalib::deg2rad;
        double diff_angle      = omega_model - m_model.posangle() * gammalib::deg2rad;
        double cosinus         = std::cos(diff_angle);
        double sinus           = std::sin(diff_angle);
        double arg1            = m_semiminor_rad * cosinus;
        double arg2            = m_semimajor_rad * sinus;
        double r_ellipse       = m_semiminor_rad * m_semimajor_rad /
                                 std::sqrt(arg1*arg1 + arg2*arg2);
        std::cout << "*** WARNING: cta_irf_elliptical_kern_omega::eval";
        std::cout << " zero model for (rho,omega)=(";
        std::cout << m_rho*gammalib::rad2deg << ",";
        std::cout << omega*gammalib::rad2deg << ")";
        std::cout << " rho-r_ellipse=" << (m_rho-r_ellipse) << " radians";
        std::cout << std::endl;
    }
    #endif

    // Continue only if model is positive
    if (model > 0.0) {

        // Compute Psf offset angle [radians]
        double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));
    
        // Compute true photon offset and azimuth angle in camera system
        // [radians]
        double theta = std::acos(m_cos_ph + m_sin_ph * std::cos(m_omega_pnt - omega));
        double phi   = 0.0; //TODO: Implement IRF Phi dependence

        // Evaluate IRF * model
        irf = m_rsp.aeff(theta, phi, m_zenith, m_azimuth, m_srcLogEng) *
              m_rsp.psf(delta, theta, phi, m_zenith, m_azimuth, m_srcLogEng) *
              model;

        // Optionally take energy dispersion into account
        if (m_rsp.use_edisp() && irf > 0.0) {
            irf *= m_rsp.edisp(m_obsEng, theta, phi, 
                               m_zenith, m_azimuth, m_srcLogEng);
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
 * @brief Kernel for zenith angle Npred integration of elliptical model
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
double cta_npred_elliptical_kern_rho::eval(const double& rho)
{
    // Initialise Npred value
    double npred = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute half length of arc that lies within ROI+PSF radius (radians)
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
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
            double rho_kluge = rho - 1.0e-12;
            if (rho_kluge < 0.0) {
                rho_kluge = 0.0;
            }

            // Setup phi integration kernel
            cta_npred_elliptical_kern_omega integrand(m_rsp,
                                                      m_model,
                                                      m_srcEng,
                                                      m_srcTime,
                                                      m_obs,
                                                      m_rot,
                                                      rho_kluge,
                                                      sin_rho,
                                                      cos_rho,
                                                      m_posangle_roi);

            // Setup integrator
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);

            // If the radius rho is not larger than the semiminor axis
            // boundary, the circle with that radius is fully contained in
            // the ellipse and we can just integrate over the relevant arc
            if (rho < m_semiminor) {

                // Compute omega integration range
                double omega_min = -domega;
                double omega_max = +domega;

                // Integrate over omega
                npred = integral.romberg(omega_min, omega_max, m_iter) *
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
                        npred     += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                    // Integrate over all intervals for omega2
                    for (int i = 0; i < intervals2.size(); ++i) {
                        double min = intervals2[i].first;
                        double max = intervals2[i].second;
                        npred     += integral.romberg(min, max, m_iter) * sin_rho;
                    }

                } // endif: arc length was positive

            } // endelse: circle was not comprised in ellipse

            // Debug: Check for NaN
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
                std::cout << "*** ERROR: cta_npred_elliptical_kern_rho::eval";
                std::cout << "(rho=" << rho << "):";
                std::cout << " NaN/Inf encountered";
                std::cout << " (npred=" << npred;
                std::cout << ", sin_rho=" << sin_rho;
                std::cout << ", cos_rho=" << cos_rho;
                std::cout << ")" << std::endl;
            }
            #endif

        } // endif: arc length was positive
    
    } // endif: rho was positive

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of elliptical model
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
double cta_npred_elliptical_kern_omega::eval(const double& omega)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute azimuth angle in model coordinate system (radians)
    double omega_model = omega + m_posangle_roi;

    // Evaluate sky model
    double model = m_model.eval(m_rho, omega_model, m_srcEng, m_srcTime);

    // Debug: test if model is non positive
    #if defined(G_DEBUG_MODEL_ZERO)
    if (model <= 0.0) {
        double m_semiminor_rad = m_model.semiminor() * gammalib::deg2rad;
        double m_semimajor_rad = m_model.semimajor() * gammalib::deg2rad;
        double diff_angle      = omega_model - m_model.posangle() * gammalib::deg2rad;
        double cosinus         = std::cos(diff_angle);
        double sinus           = std::sin(diff_angle);
        double arg1            = m_semiminor_rad * cosinus;
        double arg2            = m_semimajor_rad * sinus;
        double r_ellipse       = m_semiminor_rad * m_semimajor_rad /
                                 std::sqrt(arg1*arg1 + arg2*arg2);
        std::cout << "*** WARNING: cta_npred_elliptical_kern_omega::eval";
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
        GVector cel = m_rot * native;

        // Set sky direction
        GSkyDir srcDir;
        srcDir.celvector(cel);

        // Set Photon
        GPhoton photon(srcDir, m_srcEng, m_srcTime);

        // Compute Npred for this sky direction
        npred = m_rsp.npred(photon, m_obs) * model;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::cout << "*** ERROR: cta_npred_elliptical_kern_omega::eval";
            std::cout << "(omega=" << omega << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", model=" << model;
            std::cout << ", cos_omega=" << cos_omega;
            std::cout << ", sin_omega=" << sin_omega;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: sky intensity was positive
    
    // Return Npred
    return npred;
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
        double psf = m_rsp.psf(theta, m_theta, m_phi, m_zenith, m_azimuth, m_srcLogEng);

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
                                               m_zenith,
                                               m_azimuth,
                                               m_srcEng,
                                               m_srcTime,
                                               m_srcLogEng,
                                               m_obsEng,
                                               m_rot,
                                               sin_theta,
                                               cos_theta,
                                               sin_ph,
                                               cos_ph);
            // Integrate over phi
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);
            irf = integral.romberg(0.0, gammalib::twopi) * psf * sin_theta;
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
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Get sky intensity for this sky direction
    double intensity = m_model.eval(GPhoton(srcDir, m_srcEng, m_srcTime));

    // Continue only if sky intensity is positive
    if (intensity > 0.0) {

        // Compute true photon offset angle in camera system [radians]
        double offset = std::acos(m_cos_ph + m_sin_ph * cos_phi);

        //TODO: Compute true photon azimuth angle in camera system [radians]
        double azimuth = 0.0;

        // Evaluate model times the effective area
        irf = intensity *
              m_rsp.aeff(offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);

        // Optionally take energy dispersion into account
        if (m_rsp.use_edisp() && irf > 0.0) {
            irf *= m_rsp.edisp(m_obsEng, offset, azimuth,
                               m_zenith, m_azimuth, m_srcLogEng);
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
 * @brief Kernel for Npred offest angle integration of diffuse model
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
double cta_npred_diffuse_kern_theta::eval(const double& theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Continue only if offset angle is positive
    if (theta > 0.0) {

        // Compute sine of offset angle
        double sin_theta = std::sin(theta);

        // Setup phi integration kernel
        cta_npred_diffuse_kern_phi integrand(m_rsp,
                                             m_model,
                                             m_srcEng,
                                             m_srcTime,
                                             m_obs,
                                             m_rot,
                                             theta,
                                             sin_theta);

        // Integrate over phi
        GIntegral integral(&integrand);
        integral.fixed_iter(m_iter);
        npred = integral.romberg(0.0, gammalib::twopi) * sin_theta;
        #if defined(G_DEBUG_INTEGRAL)
        if (!integral.isvalid()) {
            std::cout << "cta_npred_diffuse_kern_theta(theta=";
            std::cout << theta*gammalib::rad2deg << " deg):" << std::endl;
            std::cout << integral.print() << std::endl;
        }
        #endif

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::cout << "*** ERROR: cta_npred_radial_kern_theta::eval";
            std::cout << "(theta=" << theta << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", sin_theta=" << sin_theta;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: offset angle was positive

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for Npred azimuth angle integration of diffuse model
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
double cta_npred_diffuse_kern_phi::eval(const double& phi)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set Photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Get sky intensity for this photon
    double intensity = m_model.eval(photon);

    // Continue only if sky intensity is positive
    if (intensity > 0.0) {

        // Compute Npred for this sky direction
        npred = m_rsp.npred(photon, m_obs) * intensity;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(npred) || gammalib::is_infinite(npred)) {
            std::cout << "*** ERROR: cta_npred_diffuse_kern_phi::eval";
            std::cout << "(phi=" << phi << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", intensity=" << intensity;
            std::cout << ", cos_phi=" << cos_phi;
            std::cout << ", sin_phi=" << sin_phi;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: sky intensity was positive

    // Return Npred
    return npred;
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
            cta_psf_diffuse_kern_phi integrand(m_model, m_srcEng, m_srcTime,
                                               m_rot,
                                               sin_delta, cos_delta);

            // Azimuthally integrate model
            GIntegral integral(&integrand);
            integral.eps(m_eps);
            value *= integral.romberg(0.0, gammalib::twopi, m_order) * sin_delta;

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
    GVector cel = m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Set photon
    GPhoton photon(srcDir, m_srcEng, m_srcTime);

    // Compute map value this sky direction
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

    // If we're at the PSF peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas.
    if (delta > 0.0) {

        // Get PSF for this delta
        value = m_rsp->psf()(m_srcDir, delta, m_srcEng);

        // Continue only if PSF is positive
        if (value > 0.0) {

            // Compute sin(delta), sin(delta)*sin(zeta), cos(delta)*cos(zeta)
            double sin_delta = std::sin(delta);
            double sin_fact  = sin_delta * m_sin_zeta;
            double cos_fact  = std::cos(delta)* m_cos_zeta;

            // Setup kernel for azimuthal integration of the spatial model
            cta_psf_radial_kern_phi integrand(m_model, m_srcEng, m_srcTime,
                                              sin_fact, cos_fact);

            // Azimuthally integrate model
            GIntegral integral(&integrand);
            integral.eps(m_eps);
            value *= integral.romberg(0.0, gammalib::twopi, m_order) * sin_delta;

        } // endif: PSF value was positive

    } // endif: delta was positive

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
 *                            \sin \delta \sin \zeta \cos \phi
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

    // Get radial model value
    double value = m_model->eval(theta, m_srcEng, m_srcTime); 

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
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
                                                          m_rho_obs,
                                                          m_cos_rho_obs,
                                                          m_sin_rho_obs,
                                                          m_delta_max,
                                                          m_cos_delta_max);

        // Continue only if arc length is positive
        if (domega > 0.0) {

            // Reduce rho by an infinite amount to avoid rounding errors
            // at the boundary of a sharp edged model (e.g. an elliptical
            // disk model)
            double rho_kluge = rho - 1.0e-12;
            if (rho_kluge < 0.0) {
                rho_kluge = 0.0;
            }

            // Evaluate sky model
            double model = m_model->eval(rho_kluge, m_srcEng, m_srcTime);

            // Debug: test if model is non positive
            #if defined(G_DEBUG_MODEL_ZERO)
            if (model <= 0.0) {
                std::cout << "*** WARNING: cta_psf_radial_kern_rho::eval";
                std::cout << " zero model for (rho)=(";
                std::cout << rho*gammalib::rad2deg << ")";
                std::cout << " rho-r_model" << (rho-m_model->theta_max());
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
    // Initialise Irf value
    double irf = 0.0;

    // Compute Psf offset angle (radians)
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));

    // Evaluate Psf * model for this delta
    irf = m_rsp->psf()(m_srcDir, delta, m_srcEng);

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
 * @brief Kernel for PSF integration of elliptical model
 *
 * @param[in] delta PSF offset angle (radians).
 * @return Azimuthally integrated product between PSF and elliptical model.
 *
 * Computes the azimuthally integrated product of point spread function and
 * the elliptical model intensity. As the PSF is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated model. The method returns thus
 *
 * \f[
 *    {\rm PSF}(\delta) \times
 *    \int_0^{2\pi} {\rm M}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm M}(\delta, \phi)\f$ is the elliptical model in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 ***************************************************************************/
double cta_psf_elliptical_kern_delta::eval(const double& delta)
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

            // Compute sin(delta), sin(delta)*sin(zeta), cos(delta)*cos(zeta)
            double sin_delta = std::sin(delta);
            double sin_fact  = sin_delta * m_sin_zeta;
            double cos_fact  = std::cos(delta)* m_cos_zeta;

            // Setup kernel for azimuthal integration of the elliptical model
            cta_psf_elliptical_kern_phi integrand(m_model, m_srcEng, m_srcTime,
                                                  m_omega,
                                                  sin_delta, sin_fact, cos_fact);

            // Azimuthally integrate model
            GIntegral integral(&integrand);
            integral.eps(m_eps);
            value *= integral.romberg(0.0, gammalib::twopi, m_order) * sin_delta;

        } // endif: PSF value was positive

    } // endif: delta was positive

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_psf_elliptical_kern_delta::eval";
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
 * @brief Kernel for azimuthal elliptical model integration
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Radial model value.
 *
 * Computes the value of the elliptical model at the position \f$(\delta,\phi)\f$
 * given in point spread function coordinates. The \f$\theta\f$ angle of the
 * elliptical model is computed using
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
 * The position angle of the elliptical model is computed using
 *
 * \f[
 *    {\rm posangle} = \arcsin \left( \frac{\sin \phi \sin \delta}
 *                                         {\sin \theta} \right)
 * \f]
 ***************************************************************************/
double cta_psf_elliptical_kern_phi::eval(const double& phi)
{
    // Compute radial model theta angle
    double theta = std::acos(m_cos_fact + m_sin_fact * std::cos(phi));

    // Compute radial model position angle
    double posangle = std::asin((std::sin(phi) * m_sin_delta)/std::sin(theta)) +
                      m_omega;

    // Get radial model value
    double value = m_model->eval(theta, posangle, m_srcEng, m_srcTime);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_psf_elliptical_kern_phi::eval";
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
        double domega = 0.5 * gammalib::cta_roi_arclength(rho,
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
            double rho_kluge = rho - 1.0e-12;
            if (rho_kluge < 0.0) {
                rho_kluge = 0.0;
            }

            // Setup integration kernel
            cta_psf_elliptical_kern_omega integrand(m_rsp,
                                                    m_model,
                                                    m_srcDir,
                                                    m_srcEng,
                                                    m_srcTime,
                                                    m_posangle_obs,
                                                    rho_kluge,
                                                    cos_psf,
                                                    sin_psf);

            // Setup integrator
            GIntegral integral(&integrand);
            integral.fixed_iter(m_iter);

            // If the radius rho is not larger than the semiminor axis
            // boundary, the circle with that radius is fully contained in
            // the ellipse and we can just integrate over the relevant arc
            if (rho < m_semiminor) {

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


/*==========================================================================
 =                                                                         =
 =          Helper class methods for response computation testing          =
 =                                                                         =
 = WARNING: These helper classes are for testing purposes only and not for =
 =          production. These classes may be removed at any moment without =
 =          any warning.                                                   =
 ==========================================================================*/


/***********************************************************************//**
 * @brief Kernel for integration of radial model in Psf system
 *
 * @param[in] delta Psf offset angle (radians).
 * @return Azimuthally integrated product between Psf and radial model.
 *
 * Computes the azimuthally integrated product of point spread function and
 * the radial model intensity. As the Psf is azimuthally symmetric, it is
 * not included in the azimuthally integration, but just multiplied on the
 * azimuthally integrated model. The method returns thus
 *
 * \f[
 *    {\rm Psf}(\delta) \times
 *    \int_0^{2\pi} {\rm M}(\delta, \phi) \sin \delta {\rm d}\phi
 * \f]
 *
 * where \f${\rm M}(\delta, \phi)\f$ is the radial model in the coordinate
 * system of the point spread function, defined by the angle \f$\delta\f$
 * between the true and the measured photon direction and the azimuth angle
 * \f$\phi\f$ around the measured photon direction.
 ***************************************************************************/
double cta_irf_radial_kern_delta::eval(const double& delta)
{
    // Initialise value
    double value = 0.0;

    // If we're at the Psf peak the model is zero (due to the sin(delta)
    // term. We thus only integrate for positive deltas.
    if (delta > 0.0) {

        // Get Psf for this delta. We use here the event offset instead of
        // the true photon offset because the Psf varies little over the
        // dimension of the Psf, hence it should basically be invariant
        // over such small scales
        value = m_rsp->psf(delta, m_obsOffset, 0.0, m_pnt.zenith(), m_pnt.azimuth(), m_srcLogEng);

        // Compile option: Get Aeff now using the event offset instead of
        // computing Aeff later in the azimuthal integration kernel where
        // the true photon offset is available. This is an approximation and
        // assumes that Aeff varies little over the dimension of the Psf.
        // This option is only for testing purposes and should not be used
        // for production.
        #if defined(G_USE_OBSDIR_FOR_AEFF)
        value *= m_rsp->aeff(m_obsOffset, 0.0, m_pnt.zenith(), m_pnt.azimuth(), m_srcLogEng);
        #endif

        // Continue only if Psf is positive
        if (value > 0.0) {

            // Compute half length of arc that lies within model (radians)
            double dphi = 0.5 * gammalib::cta_roi_arclength(delta,
                                                            m_zeta,
                                                            m_cos_zeta,
                                                            m_sin_zeta,
                                                            m_theta_max,
                                                            m_cos_theta_max);

            // Continue only if arc length is positive
            if (dphi > 0.0) {

                // Compute azimuth integration range [phi_min,phi_max].
                double phi_min = -dphi;
                double phi_max = +dphi;

                // Pre-compute terms for angular distance computation in
                // azimuth integration kernel.
                double sin_delta      = std::sin(delta);
                double cos_delta      = std::cos(delta);
                double sin_fact_model = sin_delta * m_sin_zeta;
                double cos_fact_model = cos_delta * m_cos_zeta;
                double sin_fact_inst  = sin_delta * m_sin_obs_offset;
                double cos_fact_inst  = cos_delta * m_cos_obs_offset;

                // Setup kernel for azimuthal integration of the spatial model
                cta_irf_radial_kern_phi integrand(m_rsp, m_model, m_pnt,
                                                  m_srcEng, m_srcLogEng, m_srcTime,
                                                  m_obsEng,
                                                  m_phi0,
                                                  sin_fact_model, cos_fact_model,
                                                  sin_fact_inst,  cos_fact_inst);

                // Integrate kernel
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);
                value *= integral.romberg(phi_min, phi_max, m_iter) * sin_delta;

            } // endif: arc length was positive

        } // endif: Psf value was positive

    } // endif: delta was positive

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_irf_radial_kern_delta::eval";
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
 *                            \sin \delta \sin \zeta \cos \phi
 * \f]
 *
 * where \f$\delta\f$ is the angle between true and measured photon
 * direction, \f$\zeta\f$ is the angle between model centre and measured
 * photon direction, and \f$\phi\f$ is the azimuth angle with respect to the
 * measured photon direction, where \f$\phi=0\f$ corresponds to the 
 * connecting line between model centre and measured photon direction.
 ***************************************************************************/
double cta_irf_radial_kern_phi::eval(const double& phi)
{
    // Compute radial model theta angle
    double theta = gammalib::acos(m_cos_fact_model +
                                  m_sin_fact_model * std::cos(phi));

    // Get radial model value
    double value = m_model->eval(theta, m_srcEng, m_srcTime);

    // For debugging: Check boundary violation
    if (value == 0.0) {
        std::cout << "Model 0 at theta=" << theta*gammalib::rad2deg;
        std::cout << " deg; phi=";
        std::cout << phi*gammalib::rad2deg << " deg";
        std::cout << std::endl;
        throw;
    }

    // If model is positive the multiply it with the effective area and
    // optionally fold in energy dispersion
    #if defined(G_USE_OBSDIR_FOR_AEFF)
    if (m_rsp->use_edisp() && value > 0.0) {

        // Compute offset angle of true photon in camera system (radians)
        double offset = gammalib::acos(m_cos_fact_inst +
                                       m_sin_fact_inst * std::cos(m_phi0 - phi));

        // Take energy dispersion into account
        value *= m_rsp->edisp(m_obsEng, offset, 0.0, m_pnt.zenith(), m_pnt.azimuth(), m_srcLogEng);

    } // endif: model value was positive
    #else
    if (value > 0.0) {

        // Compute offset angle of true photon in camera system (radians)
        double offset = gammalib::acos(m_cos_fact_inst +
                                       m_sin_fact_inst * std::cos(m_phi0 - phi));

        // Multiply with effective area
        value *= m_rsp->aeff(offset, 0.0, m_pnt.zenith(), m_pnt.azimuth(), m_srcLogEng);

        // Optionally take energy dispersion into account
        if (m_rsp->use_edisp() && value > 0.0) {
            value *= m_rsp->edisp(m_obsEng, offset, 0.0, m_pnt.zenith(), m_pnt.azimuth(), m_srcLogEng);
        }

    } // endif: model value was positive
    #endif
    
    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: cta_irf_radial_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return kernel value
    return value;
}
