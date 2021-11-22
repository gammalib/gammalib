/***************************************************************************
 *                COMPTEL helper classes for vector response               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file com_helpers_response_vector.cpp
 * @brief Implementation of COMPTEL helper classes for vector response
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "com_helpers_response_vector.hpp"
#include "GCOMResponse.hpp"
#include "GIntegral.hpp"
#include "GIntegrals.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * @brief Kernel for radial integration of radial models
 *
 * @param[in] rho Rho angle (radians).
 * @return Azimuthally integrated radial model.
 ***************************************************************************/
GVector com_radial_kerns_rho::eval(const double& rho)
{
    // Initialise kernel values
    m_irfs = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Compute radial model value
        double model = m_model.eval(rho, m_bin->energy(), m_bin->time());

        // Continue only if model is positive
        if (model > 0.0) {

            // Precompute cosine and sine terms for azimuthal integration
            double sin_rho = std::sin(rho);
            double cos_rho = std::cos(rho);

            // Setup azimuthal integration kernel
            com_radial_kerns_omega integrands(m_iaq,
                                              m_irfs,
                                              m_bin,
                                              m_rot,
                                              m_drx,
                                              m_phigeo_bin_size,
                                              m_phigeo_bins,
                                              m_phibar_bins,
                                              sin_rho,
                                              cos_rho);

            // Setup integrator
            GIntegrals integral(&integrands);
            integral.fixed_iter(m_iter);

            // Integrate over Omega angle
            m_irfs = integral.romberg(0.0, gammalib::twopi, m_iter) *
                     model * sin_rho;

        } // endif: radial model was positive

    } // endif: phigeo was positive

    // Return kernel values
    return m_irfs;
}


/***********************************************************************//**
 * @brief Kernel for azimuthal integration of radial models
 *
 * @param[in] omega Omega angle (radians).
 * @return Kernel value for radial model.
 ***************************************************************************/
GVector com_radial_kerns_omega::eval(const double& omega)
{
    // Initialise kernel values
    m_irfs = 0.0;

    // Compute sine and cosine of azimuth angle
    double sin_omega = std::sin(omega);
    double cos_omega = std::cos(omega);

    // Get sky direction
    GVector native(-cos_omega*m_sin_rho, sin_omega*m_sin_rho, m_cos_rho);
    GVector dir = m_rot * native;
    GSkyDir skyDir;
    skyDir.celvector(dir);

    // Compute Phigeo of current position
    double phigeo = m_bin->dir().dir().dist(skyDir);

    // Precompute interpolated IAQ vector for Phigeo angle
    double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
    int    iphigeo = int(phirat);                // index into which Phigeo falls
    double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre [-0.5, 0.5[

    // Continue only if Phigeo index is valid
    if (iphigeo < m_phigeo_bins) {

        // Get DRX value
        double intensity = m_drx(skyDir);

        // Multiply result with IAQ for all Phibar layers
        for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {

            // Initialise IAQ
            double iaq = 0.0;

            // Get IAQ index
            int i = iphibar * m_phigeo_bins + iphigeo;

            // Get interpolated IAQ value
            if (eps < 0.0) { // interpolate towards left
                if (iphigeo > 0) {
                    iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                }
                else {
                    iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                }
            }
            else {           // interpolate towards right
                if (iphigeo < m_phigeo_bins-1) {
                    iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                }
                else {
                    iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                }
            }

            // If interpolated IAQ value is positive then compute IRF value
            if (iaq > 0.0) {
                m_irfs[iphibar] = intensity * iaq;
            }

        } // endfor: looped over Phibar

    } // endif: Phigeo index was valid

    // Return kernel values
    return m_irfs;
}
/***********************************************************************//**
 * @brief Kernel for radial integration of elliptical models
 *
 * @param[in] rho Rho angle (radians).
 * @return Azimuthally integrated elliptical model.
 ***************************************************************************/
GVector com_elliptical_kerns_rho::eval(const double& rho)
{
    // Initialise kernel values
    m_irfs = 0.0;

    // Continue only if rho is positive
    if (rho > 0.0) {

        // Precompute cosine and sine terms for azimuthal integration
        double sin_rho = std::sin(rho);
        double cos_rho = std::cos(rho);

        // Setup azimuthal integration kernel
        com_elliptical_kerns_omega integrands(m_iaq,
                                              m_model,
                                              m_irfs,
                                              m_bin,
                                              m_rot,
                                              m_drx,
                                              m_phigeo_bin_size,
                                              m_phigeo_bins,
                                              m_phibar_bins,
                                              sin_rho,
                                              cos_rho);

        // Setup integrator
        GIntegrals integral(&integrands);
        integral.fixed_iter(m_iter);

        // Integrate over Omega angle
        m_irfs = integral.romberg(0.0, gammalib::twopi, m_iter) * sin_rho;

    } // endif: phigeo was positive

    // Return kernel values
    return m_irfs;
}


/***********************************************************************//**
 * @brief Kernel for azimuthal integration of elliptical models
 *
 * @param[in] omega Omega angle (radians).
 * @return Kernel value for elliptical model.
 ***************************************************************************/
GVector com_elliptical_kerns_omega::eval(const double& omega)
{
    // Initialise kernel values
    m_irfs = 0.0;

    // Compute sine and cosine of azimuth angle
    double sin_omega = std::sin(omega);
    double cos_omega = std::cos(omega);

    // Get sky direction
    GVector native(-cos_omega*m_sin_rho, sin_omega*m_sin_rho, m_cos_rho);
    GVector dir = m_rot * native;
    GSkyDir skyDir;
    skyDir.celvector(dir);

    // Compute Phigeo of current position
    double phigeo = m_bin->dir().dir().dist(skyDir);

    // Precompute interpolated IAQ vector for Phigeo angle
    double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
    int    iphigeo = int(phirat);                // index into which Phigeo falls
    double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre [-0.5, 0.5[

    // Continue only if Phigeo index is valid
    if (iphigeo < m_phigeo_bins) {

        // Set photon
        GPhoton photon(skyDir, m_bin->energy(), m_bin->time());

        // Get model sky intensity for photon (unit: sr^-1)
        double intensity = m_model.spatial()->eval(photon);

        // Continue only if intensity is positive
        if (intensity > 0.0) {

            // Multiply-in DRX value
            intensity *= m_drx(skyDir);

            // Multiply result with IAQ for all Phibar layers
            for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {

                // Initialise IAQ
                double iaq = 0.0;

                // Get IAQ index
                int i = iphibar * m_phigeo_bins + iphigeo;

                // Get interpolated IAQ value
                if (eps < 0.0) { // interpolate towards left
                    if (iphigeo > 0) {
                        iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                    }
                    else {
                        iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                    }
                }
                else {           // interpolate towards right
                    if (iphigeo < m_phigeo_bins-1) {
                        iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                    }
                    else {
                        iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                    }
                }

                // If interpolated IAQ value is positive then compute IRF value
                if (iaq > 0.0) {
                    m_irfs[iphibar] = intensity * iaq;
                }

            } // endfor: looped over Phibar

        } // endif: intensity was valid

    } // endif: Phigeo index was valid

    // Return kernel values
    return m_irfs;
}
