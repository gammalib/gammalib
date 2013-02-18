/***************************************************************************
 *         GCTAResponse_helpers.cpp - CTA response helper classes          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
#include "GCTAResponse_helpers.hpp"
#include "GCTASupport.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GVector.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =              Helper class methods for response computation              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Integration kernel for npsf() method
 *
 * @param[in] delta Distance from PSF centre (radians).
 *
 * This method implements the integration kernel needed for the npsf()
 * method. It performs an azimuthal integration of the PSF for a given
 * offset angle delta (an offset angle of 0 corresponds to the centre of
 * the PSF):
 * \f[\int_{0}^{\phi} PSF(\delta) d\phi\f]
 * 
 * The azimuthal integration is performed over an arclength given by
 * \f$\phi\f$. The method actually assumes that the PSF is azimuthally
 * symmetric, hence it just multiplies the PSF value by the arclength times
 * the sinus of the offset angle.
 ***************************************************************************/
double cta_npsf_kern_rad_azsym::eval(double delta)
{
    // Initialise PSF value
    double value = 0.0;
    
    // Get arclength for given radius in radians
    double phi = cta_roi_arclength(delta,
                                   m_psf,
                                   m_cospsf,
                                   m_sinpsf,
                                   m_roi,
                                   m_cosroi);

    // If arclength is positive then compute the PSF value
    if (phi > 0) {
    
        // Compute PSF value
        value = m_rsp->psf(delta, m_theta, m_phi, m_zenith, m_azimuth, m_logE) * phi * std::sin(delta);

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(value) || isinfinite(value)) {
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

/*==========================================================================
 =                                                                         =
 =                  Helper class methods for extended sources              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for radial model zenith angle integration of IRF
 *
 * @param[in] rho Zenith angle with respect to model centre [radians].
 *
 * This method evaluates the kernel \f$K(\rho)\f$ for the zenith angle
 * integration
 * \f[\int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho) d\rho\f]
 * of the product between model and IRF, where
 * \f[K(\rho) = \int_{\omega_{\rm min}}^{\omega_{\rm max}} M(\rho)
 *              IRF(\rho, \omega) d\omega\f],
 * \f$M(\rho)\f$ is the azimuthally symmetric source model, and
 * \f$IRF(\rho, \omega)\f$ is the instrument response function.
 ***************************************************************************/
double cta_irf_radial_kern_rho::eval(double rho)
{
    // Compute half length of arc that lies within PSF validity circle
    // (in radians)
    double domega = 0.5 * cta_roi_arclength(rho,
                                            m_zeta,
                                            m_cos_zeta,
                                            m_sin_zeta,
                                            m_delta_max,
                                            m_cos_delta_max);

    // Initialise result
    double irf = 0.0;

    // Continue only if arc length is positive
    if (domega > 0.0) {

        // Compute omega integration range
        double omega_min = -domega;
        double omega_max = +domega;

        // Evaluate sky model M(rho)
        double model = m_model->eval(rho);

        // Precompute cosine and sine terms for azimuthal integration
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
                                            m_obsLogEng,
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
        integral.eps(m_rsp->eps());
        irf = integral.romb(omega_min, omega_max) * model * sin_rho;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
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

    } // endif: arc length was positive

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for radial model azimuth angle IRF integration
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * This method evaluates the instrument response function
 * \f$IRF(\rho,\omega)\f$ for the azimuth angle integration of the IRF.
 *
 * From the model coordinates \f$(\rho,\omega)\f$ it computes the PSF
 * offset angle \f$\delta\f$, defined as the angle between true 
 * (\f$\vec{p}\f$) and observed (\f$\vec{p'}\f$) photon arrival direction,
 * using
 * \f[\delta = \arccos(\cos \rho \cos \zeta + 
 *                     \sin \rho \sin \zeta \cos \omega)\f]
 * where
 * \f$\zeta\f$ is the angular distance between the observed photon direction
 * \f$\vec{p}\f$ and the model centre \f$\vec{m}\f$.
 *
 * Furthermore, it computes the observed photon offset angle \f$\theta\f$,
 * defined as the angle between observed photon direction and camera pointing,
 * using
 * \f[\theta = \arccos(\cos \rho \cos \lambda + 
 *                     \sin \rho \sin \lambda \cos \omega_0 - \omega)\f]
 * where
 * \f$\lambda\f$ is the angular distance between the model centre and the
 * camera pointing direction.
 ***************************************************************************/
double cta_irf_radial_kern_omega::eval(double omega)
{
    // Compute PSF offset angle [radians]
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));
    
    // Compute observed photon offset angle in camera system [radians]
    double offset = std::acos(m_cos_ph + m_sin_ph * std::cos(m_omega0 - omega));
    
    //TODO: Compute true photon azimuth angle in camera system [radians]
    double azimuth = 0.0;

    // Evaluate IRF
    double irf = m_rsp->aeff(offset, azimuth, m_zenith, m_azimuth, m_srcLogEng) *
                 m_rsp->psf(delta, offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);

    // Optionally take energy dispersion into account
    if (m_rsp->hasedisp() && irf > 0.0) {
        irf *= m_rsp->edisp(m_obsLogEng, offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);
    }
    
    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(irf) || isinfinite(irf)) {
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
 * @param[in] theta Radial model zenith angle (radians).
 *
 * This method integrates a radial model for a given zenith angle theta over
 * all azimuth angles that fall within the ROI+PSF radius. The limitation to
 * an arc assures that the integration converges properly.
 ***************************************************************************/
double cta_npred_radial_kern_theta::eval(double theta)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute half length of arc that lies within ROI+PSF radius (radians)
    double dphi = 0.5 * cta_roi_arclength(theta,
                                          m_dist,
                                          m_cos_dist,
                                          m_sin_dist,
                                          m_radius,
                                          m_cos_radius);

    // Continue only if arc length is positive
    if (dphi > 0.0) {

        // Compute phi integration range
        double phi_min = m_phi - dphi;
        double phi_max = m_phi + dphi;

        // Get radial model value
        double model = m_model->eval(theta);

        // Compute sine of offset angle
        double sin_theta = std::sin(theta);

        // Setup phi integration kernel
        cta_npred_radial_kern_phi integrand(m_rsp,
                                            m_srcEng,
                                            m_srcTime,
                                            m_obs,
                                            m_rot,
                                            theta,
                                            sin_theta);

        // Integrate over phi
        GIntegral integral(&integrand);
        npred = integral.romb(phi_min, phi_max) * sin_theta * model;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (isnotanumber(npred) || isinfinite(npred)) {
            std::cout << "*** ERROR: cta_npred_radial_kern_theta::eval";
            std::cout << "(theta=" << theta << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (npred=" << npred;
            std::cout << ", model=" << model;
            std::cout << ", phi=[" << phi_min << "," << phi_max << "]";
            std::cout << ", sin_theta=" << sin_theta;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: arc length was positive

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Kernel for azimuth angle Npred integration of radial model
 *
 * @param[in] phi Azimuth angle (radians).
 *
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double cta_npred_radial_kern_phi::eval(double phi)
{
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
    GPhoton photon(srcDir, *m_srcEng, *m_srcTime);

    // Compute Npred for this sky direction
    double npred = m_rsp->npred(photon, *m_obs);

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: cta_npred_radial_kern_phi::eval";
        std::cout << "(phi=" << phi << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (npred=" << npred;
        std::cout << ", cos_phi=" << cos_phi;
        std::cout << ", sin_phi=" << sin_phi;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/*==========================================================================
 =                                                                         =
 =                Helper class methods for elliptical sources              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for elliptical model integration over model's zenith angle
 *
 * @param[in] rho Zenith angle with respect to model centre [radians].
 *
 * This method evaluates the kernel \f$K(\rho)\f$ for the zenith angle
 * integration
 * \f[\int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho) d\rho\f]
 * of the product between model and IRF, where
 * \f[K(\rho) = \int_{\omega_{\rm min}}^{\omega_{\rm max}} M(\rho, \omega)
 *              IRF(\rho, \omega) d\omega\f],
 * \f$M(\rho)\f$ is the azimuthally symmetric source model, and
 * \f$IRF(\rho, \omega)\f$ is the instrument response function.
 ***************************************************************************/
double cta_irf_elliptical_kern_rho::eval(double rho)
{
    // Compute half length of arc that lies within PSF validity circle
    // (in radians)
    double domega = 0.5 * cta_roi_arclength(rho,
                                            m_zeta,
                                            m_cos_zeta,
                                            m_sin_zeta,
                                            m_delta_max,
                                            m_cos_delta_max);

    // Initialise result
    double irf = 0.0;

    // Continue only if arc length is positive
    if (domega > 0.0) {

        // Compute omega integration range
        double omega_min = -domega;
        double omega_max = +domega;

        // Precompute cosine and sine terms for azimuthal integration
        double cos_rho = std::cos(rho);
        double sin_rho = std::sin(rho);
        double cos_psf = cos_rho*m_cos_zeta;
        double sin_psf = sin_rho*m_sin_zeta;
        double cos_ph  = cos_rho*m_cos_lambda;
        double sin_ph  = sin_rho*m_sin_lambda;

        // Setup integration kernel
        cta_irf_elliptical_kern_omega integrand(this,
                                                rho,
                                                cos_psf,
                                                sin_psf,
                                                cos_ph,
                                                sin_ph);

        // Integrate over phi
        GIntegral integral(&integrand);
        integral.eps(m_rsp->eps());
        irf = integral.romb(omega_min, omega_max) * sin_rho;

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
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

    // Return result
    return irf;
}


/***********************************************************************//**
 * @brief Kernel for elliptical model integration over model's azimuth angle
 *
 * @param[in] omega Azimuth angle (radians).
 *
 * This method evaluates the product
 *
 * \f[IRF(\rho,\omega) \times M(\rho, \omega)\f
 *
 * where 
 
 instrument response function
 * \f$IRF(\rho,\omega)\f$ for the azimuth angle integration of the IRF.
 *
 * From the model coordinates \f$(\rho,\omega)\f$ it computes the PSF
 * offset angle \f$\delta\f$, defined as the angle between true 
 * (\f$\vec{p}\f$) and observed (\f$\vec{p'}\f$) photon arrival direction,
 * using
 * \f[\delta = \arccos(\cos \rho \cos \zeta + 
 *                     \sin \rho \sin \zeta \cos \omega)\f]
 * where
 * \f$\zeta\f$ is the angular distance between the observed photon direction
 * \f$\vec{p}\f$ and the model centre \f$\vec{m}\f$.
 *
 * Furthermore, it computes the observed photon offset angle \f$\theta\f$,
 * defined as the angle between observed photon direction and camera pointing,
 * using
 * \f[\theta = \arccos(\cos \rho \cos \lambda + 
 *                     \sin \rho \sin \lambda \cos \omega_0 - \omega)\f]
 * where
 * \f$\lambda\f$ is the angular distance between the model centre and the
 * camera pointing direction.
 ***************************************************************************/
double cta_irf_elliptical_kern_omega::eval(double omega)
{
    // Compute PSF offset angle [radians]
    double delta = std::acos(m_cos_psf + m_sin_psf * std::cos(omega));
    
    // Compute observed photon offset angle in camera system [radians]
    double offset = std::acos(m_cos_ph + m_sin_ph * std::cos(m_kernel->m_omega0 - omega));
    
    // TODO: Compute true photon azimuth angle in camera system [radians]
    double azimuth = 0.0;

    // Evaluate sky model M(rho, omega)
    double model = m_kernel->m_model->eval(m_rho, omega);

    // Evaluate IRF
    double irf = m_kernel->m_rsp->aeff(offset, 
                                       azimuth,
                                       m_kernel->m_zenith, 
                                       m_kernel->m_azimuth, 
                                       m_kernel->m_srcLogEng) *
                 m_kernel->m_rsp->psf(delta, 
                                      offset,
                                      azimuth,
                                      m_kernel->m_zenith,
                                      m_kernel->m_azimuth,
                                      m_kernel->m_srcLogEng);

    // Optionally take energy dispersion into account
    if (m_kernel->m_rsp->hasedisp() && irf > 0.0) {
        irf *= m_kernel->m_rsp->edisp(m_kernel->m_obsLogEng,
                                      offset,
                                      azimuth, 
                                      m_kernel->m_zenith,
                                      m_kernel->m_azimuth,
                                      m_kernel->m_srcLogEng);
    }

    // Multiply IRF * model
    irf *= model;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(irf) || isinfinite(irf)) {
        std::cout << "*** ERROR: cta_irf_elliptical_kern_omega::eval";
        std::cout << "(omega=" << omega << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", model=" << model;
        std::cout << ", delta=" << delta;
        std::cout << ", offset=" << offset;
        std::cout << ", azimuth=" << azimuth << ")";
        std::cout << std::endl;
    }
    #endif

    // Return
    return irf;
}


/*==========================================================================
 =                                                                         =
 =                  Helper class methods for diffuse sources               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Kernel for IRF offest angle integration of the diffuse source model
 *
 * @param[in] theta Offset angle with respect to observed photon direction
 *                  (radians).
 *
 * This method provides the kernel for the IRF offset angle integration of
 * the diffuse source model.
 *
 * It computes
 *
 * \f[irf(\theta) = \sin \theta \times
 *    \int_{0}^{2\pi} M(\theta, \phi) IRF(\theta, \phi) d\phi\f]
 *
 * where \f$M(\theta, \phi)\f$ is the spatial component of the diffuse
 * source model and \f$IRF(\theta, \phi)\f$ is the response function.
 *
 * As we assume so far an azimuthally symmetric point spread function, the
 * PSF computation is done outside the azimuthal integration kernel?
 *
 * Note that the integration is only performed for \f$\theta>0\f$. Otherwise
 * zero is returned.
 ***************************************************************************/
double cta_irf_diffuse_kern_theta::eval(double theta)
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
                                               m_zenith,
                                               m_azimuth,
                                               m_srcLogEng,
                                               m_obsLogEng,
                                               m_rot,
                                               sin_theta,
                                               cos_theta,
                                               sin_ph,
                                               cos_ph);

            // Integrate over phi
            GIntegral integral(&integrand);
            integral.eps(1.0e-2);
            irf = integral.romb(0.0, twopi) * psf * sin_theta;

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (isnotanumber(irf) || isinfinite(irf)) {
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
 * This method provides the kernel for the IRF azimuth angle integration of
 * the diffuse source model.
 *
 * It computes
 *
 * \f[irf = M(\theta, \phi) Aeff(\theta, \phi) Edisp(\theta, \phi)\f]
 *
 * where \f$M(\theta, \phi)\f$ is the spatial component of the diffuse
 * source model, \f$Aeff(\theta, \phi)\f$ is the effective area, and
 * \f$Edisp(\theta, \phi)\f$ is the energy dispersion.
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
double cta_irf_diffuse_kern_phi::eval(double phi)
{
    // Initialise result
    double irf = 0.0;

    // Compute sine and cosine of azimuth angle
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    // Compute sky direction vector in native coordinates
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Get sky intensity for this sky direction
    double intensity = m_model->eval(srcDir);

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
        if (m_rsp->hasedisp() && irf > 0.0) {
            irf *= m_rsp->edisp(m_obsLogEng, offset, azimuth, m_zenith, m_azimuth, m_srcLogEng);
        }
    
        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (isnotanumber(irf) || isinfinite(irf)) {
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
double cta_npred_diffuse_kern_theta::eval(double theta)
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
        integral.eps(1.0e-4);
        npred = integral.romb(0.0, twopi) * sin_theta;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (isnotanumber(npred) || isinfinite(npred)) {
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
 * @todo Re-consider formula for possible simplification (dumb matrix
 *       multiplication is definitely not the fastest way to do that
 *       computation).
 ***************************************************************************/
double cta_npred_diffuse_kern_phi::eval(double phi)
{
    // Initialise Npred value
    double npred = 0.0;

    // Compute sky direction vector in native coordinates
    double  cos_phi = std::cos(phi);
    double  sin_phi = std::sin(phi);
    GVector native(-cos_phi*m_sin_theta, sin_phi*m_sin_theta, m_cos_theta);

    // Rotate from native into celestial system
    GVector cel = *m_rot * native;

    // Set sky direction
    GSkyDir srcDir;
    srcDir.celvector(cel);

    // Get sky intensity for this sky direction
    double intensity = m_model->eval(srcDir);

    // Continue only if sky intensity is positive
    if (intensity > 0.0) {

        // Set Photon
        GPhoton photon(srcDir, *m_srcEng, *m_srcTime);

        // Compute Npred for this sky direction
        npred = m_rsp->npred(photon, *m_obs) * intensity;

        // Debug: Check for NaN
        #if defined(G_NAN_CHECK)
        if (isnotanumber(npred) || isinfinite(npred)) {
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
