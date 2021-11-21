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
 * @class com_extended_kerns_phigeo
 *
 * @brief Kernel for Phigeo angle integration of extended model
 *
 * This class provides the kernel for the Phigeo angle integration of
 * extended sky models.
 ***************************************************************************/
class com_extended_kerns_phigeo : public GFunctions {
public:
    com_extended_kerns_phigeo(const std::vector<double>& iaq,
                              const GModelSky&           model,
                              GVector&                   irfs,
                              const GEnergy&             srcEng,
                              const GTime&               srcTime,
                              const GMatrix&             rot,
                              const GSkyMap*             drx,
                              const double&              phigeo_bin_size,
                              const int&                 phigeo_bins,
                              const int&                 phibar_bins,
                              const int&                 ipix,
                              const int&                 npix,
                              const double&              iaq_norm,
                              const double&              zeta,
                              const double&              phi0,
                              const double&              theta_max,
                              const int&                 iter) :
                              m_iaq(iaq),
                              m_model(model),
                              m_irfs(irfs),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_rot(rot),
                              m_drx(drx),
                              m_phigeo_bin_size(phigeo_bin_size),
                              m_phigeo_bins(phigeo_bins),
                              m_phibar_bins(phibar_bins),
                              m_ipix(ipix),
                              m_npix(npix),
                              m_iaq_norm(iaq_norm),
                              m_zeta(zeta),
                              m_sin_zeta(std::sin(zeta)),
                              m_cos_zeta(std::cos(zeta)),
                              m_phi0(phi0),
                              m_theta_max(theta_max),
                              m_cos_theta_max(std::cos(theta_max)),
                              m_iter(iter) { }
    int     size(void) const { return m_irfs.size(); }
    GVector eval(const double& phigeo);
protected:
    const std::vector<double>& m_iaq;             //!< IAQ vector
    const GModelSky&           m_model;           //!< Sky model
    GVector&                   m_irfs;            //!< IRF vector to update
    const GEnergy&             m_srcEng;          //!< Source energy
    const GTime&               m_srcTime;         //!< Source time
    const GMatrix&             m_rot;             //!< Rotation matrix
    const GSkyMap*             m_drx;             //!< DRX
    const double&              m_phigeo_bin_size; //!< Phigeo bin size
    const int&                 m_phigeo_bins;     //!< Number of phigeo bins
    const int&                 m_phibar_bins;     //!< Number of phibar bins
    const int&                 m_ipix;            //!< DRI pixel index
    const int&                 m_npix;            //!< Number of DRI pixels
    const double&              m_iaq_norm;        //!< IAQ normalisation
    const double&              m_zeta;            //!< Zeta angle
    double                     m_sin_zeta;        //!< Sine of zeta angle
    double                     m_cos_zeta;        //!< Cosine of zeta angle
    const double&              m_phi0;            //!< Model position angle
    const double&              m_theta_max;       //!< Maximum model radius
    double                     m_cos_theta_max;   //!< Cosine of max. mod. rad.
    const int&                 m_iter;            //!< Number of phi iterations
};


/***********************************************************************//**
 * @class com_extended_kerns_phi
 *
 * @brief Kernel for azimuth angle integration of extended model
 *
 * This class provides the kernel for the azimuth angle integration of
 * extended sky models.
 ***************************************************************************/
class com_extended_kerns_phi : public GFunction {
public:
    com_extended_kerns_phi(const GModelSky&           model,
                           const GEnergy&             srcEng,
                           const GTime&               srcTime,
                           const GMatrix&             rot,
                           const GSkyMap*             drx,
                           const double&              sin_phigeo,
                           const double&              cos_phigeo) :
                           m_model(model),
                           m_srcEng(srcEng),
                           m_srcTime(srcTime),
                           m_rot(rot),
                           m_drx(drx),
                           m_sin_phigeo(sin_phigeo),
                           m_cos_phigeo(cos_phigeo) { }
    double eval(const double& phi);
protected:
    const GModelSky& m_model;       //!< Sky model
    const GEnergy&   m_srcEng;      //!< Source energy
    const GTime&     m_srcTime;     //!< Source time
    const GMatrix&   m_rot;         //!< Rotation matrix
    const GSkyMap*   m_drx;         //!< DRX
    const double&    m_sin_phigeo;  //!< Sine of Phigeo
    const double&    m_cos_phigeo;  //!< Cosine of Phigeo
};


/***********************************************************************//**
 * @class com_extended_kerns_rho
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
 * @class com_extended_kerns_omega
 *
 * @brief Kernel for azimuth angle integration of extended models
 *
 * This class provides the kernel for the azimuth angle integration of
 * extended sky models.
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
 * @class com_radial_kerns_rho
 *
 * @brief Kernel for rho angle integration of extended models
 *
 * This class provides the kernel for the rho angle integration of extended
 * sky models.
 ***************************************************************************/
class com_extended_kerns_rho : public GFunctions {
public:
    com_extended_kerns_rho(const std::vector<double>& iaq,
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
 * @class com_radial_kerns_omega
 *
 * @brief Kernel for azimuth angle integration of extended models
 *
 * This class provides the kernel for the azimuth angle integration of
 * extended sky models.
 ***************************************************************************/
class com_extended_kerns_omega : public GFunctions {
public:
    com_extended_kerns_omega(const std::vector<double>& iaq,
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
