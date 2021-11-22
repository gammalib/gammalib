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
 * @file com_helpers_response_vector.hpp
 * @brief Defintion of COMPTEL helper classes for vector response
 * @author Juergen Knoedlseder
 */

#ifndef COM_HELPERS_RESPONSE_VECTOR_HPP
#define COM_HELPERS_RESPONSE_VECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GFunction.hpp"
#include "GFunctions.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialRadial.hpp"
#include "GSkyMap.hpp"
#include "GCOMResponse.hpp"
#include "GCOMEventBin.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class com_radial_kerns_rho
 *
 * @brief Kernel for rho angle integration of radial models
 *
 * This class provides the kernel for the rho angle integration of radial
 * sky models.
 ***************************************************************************/
class com_radial_kerns_rho : public GFunctions {
public:
    com_radial_kerns_rho(const std::vector<double>& iaq,
                         const GModelSpatialRadial& model,
                         GVector&                   irfs,
                         const GCOMEventBin*        bin,
                         const GMatrix&             rot,
                         const GSkyMap&             drx,
                         const double&              phigeo_bin_size,
                         const int&                 phigeo_bins,
                         const int&                 phibar_bins,
                         const int&                 iter) :
                         m_iaq(iaq),
                         m_model(model),
                         m_irfs(irfs),
                         m_bin(bin),
                         m_rot(rot),
                         m_drx(drx),
                         m_phigeo_bin_size(phigeo_bin_size),
                         m_phigeo_bins(phigeo_bins),
                         m_phibar_bins(phibar_bins),
                         m_iter(iter) { }
    int     size(void) const { return m_irfs.size(); }
    GVector eval(const double& phigeo);
protected:
    const std::vector<double>& m_iaq;             //!< IAQ vector
    const GModelSpatialRadial& m_model;           //!< Radial spatial model
    GVector&                   m_irfs;            //!< IRF vector to update
    const GCOMEventBin*        m_bin;             //!< Event bin
    const GMatrix&             m_rot;             //!< Rotation matrix
    const GSkyMap&             m_drx;             //!< DRX
    const double&              m_phigeo_bin_size; //!< Phigeo bin size
    const int&                 m_phigeo_bins;     //!< Number of phigeo bins
    const int&                 m_phibar_bins;     //!< Number of phibar bins
    const int&                 m_iter;            //!< Number of omega iterations
};


/***********************************************************************//**
 * @class com_radial_kerns_omega
 *
 * @brief Kernel for azimuth angle integration of radial models
 *
 * This class provides the kernel for the azimuth angle integration of
 * radial sky models.
 ***************************************************************************/
class com_radial_kerns_omega : public GFunctions {
public:
    com_radial_kerns_omega(const std::vector<double>& iaq,
                           GVector&                   irfs,
                           const GCOMEventBin*        bin,
                           const GMatrix&             rot,
                           const GSkyMap&             drx,
                           const double&              phigeo_bin_size,
                           const int&                 phigeo_bins,
                           const int&                 phibar_bins,
                           const double&              sin_rho,
                           const double&              cos_rho) :
                           m_iaq(iaq),
                           m_irfs(irfs),
                           m_bin(bin),
                           m_rot(rot),
                           m_drx(drx),
                           m_phigeo_bin_size(phigeo_bin_size),
                           m_phigeo_bins(phigeo_bins),
                           m_phibar_bins(phibar_bins),
                           m_sin_rho(sin_rho),
                           m_cos_rho(cos_rho) { }
    int     size(void) const { return m_irfs.size(); }
    GVector eval(const double& omega);
protected:
    const std::vector<double>& m_iaq;             //!< IAQ vector
    GVector&                   m_irfs;            //!< IRF vector to update
    const GCOMEventBin*        m_bin;             //!< Event bin
    const GMatrix&             m_rot;             //!< Rotation matrix
    const GSkyMap&             m_drx;             //!< DRX
    const double&              m_phigeo_bin_size; //!< Phigeo bin size
    const int&                 m_phigeo_bins;     //!< Number of phigeo bins
    const int&                 m_phibar_bins;     //!< Number of phibar bins
    const double&              m_sin_rho;         //!< Sine of Rho
    const double&              m_cos_rho;         //!< Cosine of Rho
};


/***********************************************************************//**
 * @class com_elliptical_kerns_rho
 *
 * @brief Kernel for rho angle integration of elliptical models
 *
 * This class provides the kernel for the rho angle integration of elliptical
 * sky models.
 ***************************************************************************/
class com_elliptical_kerns_rho : public GFunctions {
public:
    com_elliptical_kerns_rho(const std::vector<double>& iaq,
                             const GModelSky&           model,
                             GVector&                   irfs,
                             const GCOMEventBin*        bin,
                             const GMatrix&             rot,
                             const GSkyMap&             drx,
                             const double&              phigeo_bin_size,
                             const int&                 phigeo_bins,
                             const int&                 phibar_bins,
                             const int&                 iter) :
                             m_iaq(iaq),
                             m_model(model),
                             m_irfs(irfs),
                             m_bin(bin),
                             m_rot(rot),
                             m_drx(drx),
                             m_phigeo_bin_size(phigeo_bin_size),
                             m_phigeo_bins(phigeo_bins),
                             m_phibar_bins(phibar_bins),
                             m_iter(iter) { }
    int     size(void) const { return m_irfs.size(); }
    GVector eval(const double& phigeo);
protected:
    const std::vector<double>& m_iaq;             //!< IAQ vector
    const GModelSky&           m_model;           //!< Sky model
    GVector&                   m_irfs;            //!< IRF vector to update
    const GCOMEventBin*        m_bin;             //!< Event bin
    const GMatrix&             m_rot;             //!< Rotation matrix
    const GSkyMap&             m_drx;             //!< DRX
    const double&              m_phigeo_bin_size; //!< Phigeo bin size
    const int&                 m_phigeo_bins;     //!< Number of phigeo bins
    const int&                 m_phibar_bins;     //!< Number of phibar bins
    const int&                 m_iter;            //!< Number of omega iterations
};


/***********************************************************************//**
 * @class com_elliptical_kerns_omega
 *
 * @brief Kernel for azimuth angle integration of elliptical models
 *
 * This class provides the kernel for the azimuth angle integration of
 * elliptical sky models.
 ***************************************************************************/
class com_elliptical_kerns_omega : public GFunctions {
public:
    com_elliptical_kerns_omega(const std::vector<double>& iaq,
                               const GModelSky&           model,
                               GVector&                   irfs,
                               const GCOMEventBin*        bin,
                               const GMatrix&             rot,
                               const GSkyMap&             drx,
                               const double&              phigeo_bin_size,
                               const int&                 phigeo_bins,
                               const int&                 phibar_bins,
                               const double&              sin_rho,
                               const double&              cos_rho) :
                               m_iaq(iaq),
                               m_model(model),
                               m_irfs(irfs),
                               m_bin(bin),
                               m_rot(rot),
                               m_drx(drx),
                               m_phigeo_bin_size(phigeo_bin_size),
                               m_phigeo_bins(phigeo_bins),
                               m_phibar_bins(phibar_bins),
                               m_sin_rho(sin_rho),
                               m_cos_rho(cos_rho) { }
    int     size(void) const { return m_irfs.size(); }
    GVector eval(const double& omega);
protected:
    const std::vector<double>& m_iaq;             //!< IAQ vector
    const GModelSky&           m_model;           //!< Sky model
    GVector&                   m_irfs;            //!< IRF vector to update
    const GCOMEventBin*        m_bin;             //!< Event bin
    const GMatrix&             m_rot;             //!< Rotation matrix
    const GSkyMap&             m_drx;             //!< DRX
    const double&              m_phigeo_bin_size; //!< Phigeo bin size
    const int&                 m_phigeo_bins;     //!< Number of phigeo bins
    const int&                 m_phibar_bins;     //!< Number of phibar bins
    const double&              m_sin_rho;         //!< Sine of Rho
    const double&              m_cos_rho;         //!< Cosine of Rho
};

#endif /* COM_HELPERS_RESPONSE_VECTOR_HPP */
