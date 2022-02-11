/***************************************************************************
 *              CTA helper classes for stacked vector response             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2022 by Juergen Knoedlseder                         *
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
 * @file cta_helpers_response_stacked_vector.hpp
 * @brief Implementation of CTA helper classes for stacked vector response
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "cta_helpers_response_stacked_vector.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GModelPar.hpp"
#include "GModelSpatialRadial.hpp"
#include "GCTAResponseCube.hpp"
#include "GIntegrals.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PSF_RADIAL_KERNS_PHI      "cta_psf_radial_kerns_phi::eval(double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                             Helper methods                              =
 =                                                                         =
 ==========================================================================*/

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
 * @param[in] rsp Pointer to CTA response cube.
 * @param[in] model Pointer to radial spatial model.
 * @param[in] obsDir Observed sky direction.
 * @param[in] srcEngs True energies.
 * @param[in] zeta Angular distance between event direction of model centre (radians).
 * @param[in] theta_max Maximum model radius (radians).
 * @param[in] iter Number of iterations for azimuth integration.
 * @param[in] grad Compute parameter gradients?
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

    // Set coordinate system
    m_par_cel = (model->coordsys() == "CEL");

    // Store references to longitude and latitude parameters
    if (m_par_cel) {
        m_par_lon = &((*(const_cast<GModelSpatialRadial*>(model)))["RA"]);
        m_par_lat = &((*(const_cast<GModelSpatialRadial*>(model)))["DEC"]);
    }
    else {
        m_par_lon = &((*(const_cast<GModelSpatialRadial*>(model)))["GLON"]);
        m_par_lat = &((*(const_cast<GModelSpatialRadial*>(model)))["GLAT"]);
    }

    // Pre-compute constants needed for gradient computation. The following
    // exceptions need to be handled
    // * beta_0 = +/-90 deg
    // * m_sin_zeta = 0
    // * denom = 0
    if (grad) {
        double alpha_0        = m_par_cel ? model->dir().ra()  : model->dir().l();
        double beta_0         = m_par_cel ? model->dir().dec() : model->dir().b();
        double alpha_reco     = m_par_cel ? obsDir.ra()        : obsDir.l();
        double beta_reco      = m_par_cel ? obsDir.dec()       : obsDir.b();
        double sin_beta_0     = std::sin(beta_0);
        double cos_beta_0     = std::cos(beta_0);
        double tan_beta_0     = std::tan(beta_0); // Exception: beta_0 = 90 deg
        double sin_beta_reco  = std::sin(beta_reco);
        double cos_beta_reco  = std::cos(beta_reco);
        double sin_dalpha     = std::sin(alpha_0 - alpha_reco);
        double cos_dalpha     = std::cos(alpha_0 - alpha_reco);
        double arg            = cos_beta_reco * tan_beta_0 -
                                sin_beta_reco * cos_dalpha;
        double denom          = sin_dalpha * sin_dalpha + arg * arg;
        if (m_sin_zeta != 0.0) {
            m_dzeta_dalpha_0  = cos_beta_0 * cos_beta_reco * sin_dalpha / m_sin_zeta;
            m_dzeta_dbeta_0   = (sin_beta_0 * cos_beta_reco * cos_dalpha -
                                 cos_beta_0 * sin_beta_reco) / m_sin_zeta;
        }
        else {
            m_dzeta_dalpha_0 = cos_beta_0 * cos_beta_reco;
            m_dzeta_dbeta_0  = 1.0;
        }
        if (denom != 0.0) {
            m_dphi_dalpha_0  = (sin_beta_reco - cos_dalpha * cos_beta_reco * tan_beta_0) /
                                denom;
            m_dphi_dbeta_0   = (sin_dalpha * cos_beta_reco * (1.0 + tan_beta_0 * tan_beta_0)) /
                                denom;
        }
        else {
            m_dphi_dalpha_0  = 0.0;
            m_dphi_dbeta_0   = 0.0;
        }
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

            // If gradients are requested then use vector integrator
            if (m_grad) {

                // Precompute terms for gradient computation
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

                // Multiply in energy dependent PSF
                for (int i = 0; i < nengs; ++i) {

                    // Get Psf for this energy. We approximate here the true sky
                    // direction by the reconstructed sky direction.
                    double psf = m_rsp->psf()(m_obsDir, delta, m_srcEngs[i]);

                    // Continue only if PSF is positive
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

                    } // endif: PSF was positive

                } // endfor: looped over energies

            } // endif: gradient computation was requested

            // ... otherwise use scalar integration
            else {
 
                // Allocate dummy energy and time to satisfy interface
                static const GEnergy srcEng;
                static const GTime   srcTime;

                // Setup kernel for azimuthal integration of the spatial model
                cta_psf_radial_kern_phi integrand(m_model,
                                                  srcEng,
                                                  srcTime,
                                                  sin_delta_sin_zeta,
                                                  sin_delta_cos_zeta);

                // Setup integrator
                GIntegral integral(&integrand);
                integral.fixed_iter(m_iter);

                // Integrate over azimuth
                double value = integral.romberg(phi_min, phi_max, m_iter) * sin_delta;

                // Multiply in energy dependent PSF
                for (int i = 0; i < nengs; ++i) {

                    // Get Psf for this energy. We approximate here the true sky
                    // direction by the reconstructed sky direction.
                    double psf = m_rsp->psf()(m_obsDir, delta, m_srcEngs[i]);

                    // If PSF is positive then compute corresponding value
                    if (psf > 0.0) {
                        values[i] = value * psf;
                    }

                } // endfor: looped over energies

            } // endelse: scalar integration requested

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

    // If gradients are requested then compute partial derivatives of theta
    // and phi with respect to Right Ascension and Declination and store them
    // into the respective factor gradients so that the
    // GModelSpatialRadial::eval() method can use the information.
    if (m_outer->m_grad) {
        double g_lon = 0.0;
        double g_lat = 0.0;
        if (theta > 0.0) {

            // Compute sin(theta) by avoiding to use a call to sine
            double sin_theta2 = 1.0 - cos_theta * cos_theta;
            double sin_theta  = (sin_theta2 > 0.0) ? std::sqrt(sin_theta2) : 0.0;

            // Continue only if sine is non-zero
            if (sin_theta != 0.0) {
                double norm         = 1.0 / sin_theta;
                double dtheta_dzeta = (m_cos_delta_sin_zeta -
                                       m_sin_delta_cos_zeta * cos_phi) *
                                       norm;
                double dtheta_dphi  = (m_sin_delta_sin_zeta * sin_phi) *
                                       norm;
                if (m_outer->m_zeta != 0.0) {
                    g_lon = dtheta_dzeta * m_outer->m_dzeta_dalpha_0 +
                            dtheta_dphi  * m_outer->m_dphi_dalpha_0;
                    g_lat = dtheta_dzeta * m_outer->m_dzeta_dbeta_0  +
                            dtheta_dphi  * m_outer->m_dphi_dbeta_0;
                }
                else {
                    g_lon =  0.0;    //TODO: Set correct value
                    g_lat = -cos_phi;
                }
            }
        }
        m_outer->m_par_lon->factor_gradient(g_lon);
        m_outer->m_par_lat->factor_gradient(g_lat);
    }

    // Compute model value and optionally model parameter gradients
    m_values[0] = m_outer->m_model->eval(theta_kluge, srcEng, srcTime,
                                         m_outer->m_grad);

    // If gradients are requested, extract now the fully computed model
    // parameter gradients into the kernel values vector
    if (m_outer->m_grad) {
        for (int i = 1, k = 0; i < m_size; ++i, ++k) {
            m_values[i] = (*(m_outer->m_model))[k].factor_gradient();
        }
    }

    // Debug: Check for NaN
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(m_values[0]) ||
        gammalib::is_infinite(m_values[0])) {
        std::string msg = "NaN/Inf encountered for phi="+gammalib::str(phi);
        gammalib::warning(G_PSF_RADIAL_KERNS_PHI, msg);
    }
    #endif

    // Return kernel values
    return m_values;
}
