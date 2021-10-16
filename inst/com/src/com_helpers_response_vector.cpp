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
 * @brief Kernel for radial integration of extended model
 *
 * @param[in] phigeo Phigeo angle (radians).
 * @return Azimuthally integrated extended model.
 ***************************************************************************/
GVector com_extended_kerns_phigeo::eval(const double& phigeo)
{
    // Continue only if phigeo is positive
    if (phigeo > 0.0) {

        // Compute half length of the arc (in radians) from a circle with
        // radius phigeo that intersects with the model, defined as a
        // circle with maximum radius m_theta_max
        double dphi = 0.5 * gammalib::roi_arclength(phigeo,
                                                    m_zeta,
                                                    m_cos_zeta,
                                                    m_sin_zeta,
                                                    m_theta_max,
                                                    m_cos_theta_max);

        // Continue only if arc length is positive
        if (dphi > 0.0) {

            // Compute sine and cosine of Phigeo
            double sin_phigeo = std::sin(phigeo);
            double cos_phigeo = std::cos(phigeo);

            // Compute phi integration range
            double phi_min = m_phi0 - dphi;
            double phi_max = m_phi0 + dphi;

            // Precompute interpolated IAQ vector for Phigeo angle
            double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
            int    iphigeo = int(phirat);                // index into which Phigeo falls
            double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre [-0.5, 0.5[

            // Continue only if Phigeo index is valid
            if (iphigeo < m_phigeo_bins) {

                // Initialise IAQ vector to hold precomputation results
                GVector iaq(m_phibar_bins);

                // Loop over Phibar
                for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {

                    // Get IAQ index
                    int i = iphibar * m_phigeo_bins + iphigeo;

                    // Get interpolated IAQ value
                    if (eps < 0.0) { // interpolate towards left
                        if (iphigeo > 0) {
                            iaq[iphibar] = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                        }
                        else {
                            iaq[iphibar] = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                        }
                    }
                    else {           // interpolate towards right
                        if (iphigeo < m_phigeo_bins-1) {
                            iaq[iphibar] = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                        }
                        else {
                            iaq[iphibar] = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                        }
                    }

                    // Normalise IAQ value
                    iaq[iphibar] *= m_iaq_norm;

                } // endfor: looped over Phibar

                // Setup kernel for azimuthal integration
                com_extended_kerns_phi integrands(m_model,
                                                  m_irfs,
                                                  m_srcEng,
                                                  m_srcTime,
                                                  m_rot,
                                                  m_drx,
                                                  iaq,
                                                  sin_phigeo,
                                                  cos_phigeo);

                // Setup integrator
                GIntegrals integral(&integrands);
                integral.fixed_iter(m_iter);

                // Integrate over azimuth
                m_irfs = integral.romberg(phi_min, phi_max, m_iter) * sin_phigeo;

            } // endif: Phigeo bin was valid

        } // endif: arc length was positive

    } // endif: phigeo was positive

    // Return kernel values
    return m_irfs;
}


/***********************************************************************//**
 * @brief Kernel for azimuthal integration of extended model
 *
 * @param[in] phi Azimuth angle (radians).
 * @return Vector of azimuthally integrated model.
 ***************************************************************************/
GVector com_extended_kerns_phi::eval(const double& phi)
{
    // Initialise kernel values
    m_irfs = 0.0;

    // Compute sine and cosine of azimuth angle
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    // Get sky direction
    GVector native(-cos_phi*m_sin_phigeo, sin_phi*m_sin_phigeo, m_cos_phigeo);
    GVector dir = m_rot * native;
    GSkyDir skyDir;
    skyDir.celvector(dir);

    // Set photon
    GPhoton photon(skyDir, m_srcEng, m_srcTime);

    // Get model sky intensity for photon (unit: sr^-1)
    double intensity = m_model.spatial()->eval(photon);

    // Continue only if intensity is positive
    if (intensity > 0.0) {

        // Multiply-in DRX value if pointer is valid
        if (m_drx != NULL) {
            intensity *= (*m_drx)(skyDir);
        }

        // Continue only if intensity is still positive
        if (intensity > 0.0) {

            // Get number of Phibar layers
            int nphibar = m_iaq.size();

            // Loop over Phibar
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                // Continue only if IAQ value is positive
                if (m_iaq[iphibar] > 0.0) {

                    // Compute IRF value (unit: cm^2)
                    double irf = m_iaq[iphibar] * intensity;

                    // Set IRF value
                    m_irfs[iphibar] = irf;

                } // endif: IAQ value was valid

            } // endfor: looped over Phibar

        } // endif: intensity was valid

    } // endif: intensity was valid

    // Return kernel values
    return m_irfs;
}
