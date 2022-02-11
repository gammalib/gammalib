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
 * @brief Defintion of CTA helper classes for stacked vector response
 * @author Juergen Knoedlseder
 */

#ifndef CTA_HELPERS_RESPONSE_STACKED_VECTOR_HPP
#define CTA_HELPERS_RESPONSE_STACKED_VECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include "GSkyDir.hpp"
#include "GEnergies.hpp"
#include "GVector.hpp"
#include "GFunctions.hpp"
#include "GModelSpatialRadial.hpp"

/* __ Forward declarations _______________________________________________ */
class GModelPar;
class GCTAResponseCube;


/***********************************************************************//**
 * @class cta_psf_radial_kerns_delta
 *
 * @brief Kernel for radial spatial model PSF delta angle integration
 *
 * This class provides the kernel for the radial spatial model integration
 * of the delta angle of the PSF system.
 ***************************************************************************/
class cta_psf_radial_kerns_delta : public GFunctions {
    friend class cta_psf_radial_kerns_phi;
public:
    cta_psf_radial_kerns_delta(const GCTAResponseCube*    rsp,
                               const GModelSpatialRadial* model,
                               const GSkyDir&             obsDir,
                               const GEnergies&           srcEngs,
                               const double&              zeta,
                               const double&              theta_max,
                               const int&                 iter,
                               const bool&                grad);
    int     size(void) const;
    GVector eval(const double& delta);
protected:
    const GCTAResponseCube*    m_rsp;            //!< Response cube
    const GModelSpatialRadial* m_model;          //!< Radial model
    GModelPar*                 m_par_lon;        //!< Longitude parameter
    GModelPar*                 m_par_lat;        //!< Latitude parameter
    bool                       m_par_cel;        //!< Celestial or galactic coordinates
    GSkyDir                    m_obsDir;         //!< Reconstructed event direction
    GEnergies                  m_srcEngs;        //!< True photon energies
    double                     m_zeta;           //!< Distance of model from Psf
    double                     m_cos_zeta;       //!< Cosine of m_zeta
    double                     m_sin_zeta;       //!< Sine of m_zeta
    double                     m_theta_max;      //!< Maximum model radius
    double                     m_cos_theta_max;  //!< Cosine of m_theta_max
    double                     m_dzeta_dalpha_0; //!< d(zeta)/d(alpha0)
    double                     m_dzeta_dbeta_0;  //!< d(zeta)/d(beta0)
    double                     m_dphi_dalpha_0;  //!< d(phi)/d(alpha0)
    double                     m_dphi_dbeta_0;   //!< d(phi)/d(beta0)
    int                        m_iter;           //!< Integration iterations
    bool                       m_grad;           //!< Compute gradients
};


/***********************************************************************//**
 * @class cta_psf_radial_kerns_phi
 *
 * @brief Kernel for radial spatial model PSF phi angle integration
 *
 * This class provides the kernel for the radial spatial model integration
 * of the azimuth angle of the PSF system.
 *
 * The eval() method of the kernel returns a vector that contains for a
 * given azimuth angle phi the model value and all model parameter gradients.
 ***************************************************************************/
class cta_psf_radial_kerns_phi : public GFunctions {
public:
    cta_psf_radial_kerns_phi(cta_psf_radial_kerns_delta* outer,
                             const double&               sin_delta_sin_zeta,
                             const double&               sin_delta_cos_zeta,
                             const double&               cos_delta_sin_zeta,
                             const double&               cos_delta_cos_zeta) :
                             m_outer(outer),
                             m_size(outer->m_model->size()+1),
                             m_sin_delta_sin_zeta(sin_delta_sin_zeta),
                             m_sin_delta_cos_zeta(sin_delta_cos_zeta),
                             m_cos_delta_sin_zeta(cos_delta_sin_zeta),
                             m_cos_delta_cos_zeta(cos_delta_cos_zeta),
                             m_values(GVector(outer->m_model->size()+1)) { }
    int     size(void) const { return m_size; }
    GVector eval(const double& phi);
protected:
    cta_psf_radial_kerns_delta* m_outer;              //!< Pointer to outer integr.
    int                         m_size;               //!< Result size
    double                      m_sin_delta_sin_zeta; //!< sin(delta) * sin(zeta)
    double                      m_sin_delta_cos_zeta; //!< sin(delta) * cos(zeta)
    double                      m_cos_delta_sin_zeta; //!< cos(delta) * sin(zeta)
    double                      m_cos_delta_cos_zeta; //!< cos(delta) * cos(zeta)
    GVector                     m_values;             //!< Return value
};

#endif /* CTA_HELPERS_RESPONSE_STACKED_VECTOR_HPP */
