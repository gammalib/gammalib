/***************************************************************************
 *         GCTAResponse_helpers.hpp - CTA response helper classes          *
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
 * @file GCTAResponse_helpers.hpp
 * @brief CTA response hepler classes definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSE_HELPERS_HPP
#define GCTARESPONSE_HELPERS_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GCTAResponse.hpp"
#include "GCTAObservation.hpp"
#include "GMatrix.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GModelSpatialRadial.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GFunction.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */


/***********************************************************************//**
 * @class cta_npsf_kern_rad_azsym
 *
 * @brief Integration kernel for npsf() method
 *
 * This class implements the integration kernel needed for the npsf() method.
 * method. It performs an azimuthal integration of the Point Spread Function
 * (PSF) for a given offset angle \f$delta\f$ using
 *
 * \f[
 *    K(\delta) = \int_{0}^{\phi} PSF(\delta) d\phi
 * \f]
 * (an offset angle of \f$\delta=0\f$ corresponds to the centre of the PSF).
 ***************************************************************************/
class cta_npsf_kern_rad_azsym : public GFunction {
public:
    cta_npsf_kern_rad_azsym(const GCTAResponse& rsp,
                            const double&       roi,
                            const double&       psf,
                            const double&       logE,
                            const double&       theta,
                            const double&       phi,
                            const double&       zenith,
                            const double&       azimuth) :
                            m_rsp(rsp),
                            m_roi(roi),
                            m_cosroi(std::cos(roi)),
                            m_psf(psf),
                            m_cospsf(std::cos(psf)),
                            m_sinpsf(std::sin(psf)),
                            m_logE(logE),
                            m_theta(theta),
                            m_phi(phi),
                            m_zenith(zenith),
                            m_azimuth(azimuth) { }
    double eval(double delta);
protected:
    const GCTAResponse& m_rsp;     //!< CTA response function
    const double&       m_roi;     //!< ROI radius in radians
    double              m_cosroi;  //!< Cosine of ROI radius
    const double&       m_psf;     //!< PSF-ROI centre distance in radians
    double              m_cospsf;  //!< Cosine of PSF-ROI centre distance
    double              m_sinpsf;  //!< Sine of PSF-ROI centre distance
    const double&       m_logE;    //!< Log10 of true photon energy (E/TeV).
    const double&       m_theta;   //!< Offset angle of source in camera system
    const double&       m_phi;     //!< Azimuth angle of source in camera system
    const double&       m_zenith;  //!< Zenith angle of source in Earth system
    const double&       m_azimuth; //!< Azimuth angle of source in Earth system
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_rho
 *
 * @brief Kernel for radial model zenith angle integration of IRF
 *
 * This class implements the integration kernel \f$K(\rho)\f$ for the
 * integration
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho | E, t) d\rho
 * \f]
 *
 * of radial spatial models. The eval() method computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     IRF(\rho, \omega) d\omega
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho | E, t)\f$ is the radial model,
 * - \f$IRF(\rho, \omega)\f$ is the IRF
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_irf_radial_kern_rho : public GFunction {
public:
    cta_irf_radial_kern_rho(const GCTAResponse&        rsp,
                            const GModelSpatialRadial& model,
                            const double&              zenith,
                            const double&              azimuth,
                            const GEnergy&             srcEng,
                            const GTime&               srcTime,
                            const double&              srcLogEng,
                            const double&              obsLogEng,
                            const double&              zeta,
                            const double&              lambda,
                            const double&              omega0,
                            const double&              delta_max) :
                            m_rsp(rsp),
                            m_model(model),
                            m_zenith(zenith),
                            m_azimuth(azimuth),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime),
                            m_srcLogEng(srcLogEng),
                            m_obsLogEng(obsLogEng),
                            m_zeta(zeta),
                            m_cos_zeta(std::cos(zeta)),
                            m_sin_zeta(std::sin(zeta)),
                            m_lambda(lambda),
                            m_cos_lambda(std::cos(lambda)),
                            m_sin_lambda(std::sin(lambda)),
                            m_omega0(omega0),
                            m_delta_max(delta_max),
                            m_cos_delta_max(std::cos(delta_max)) { }
    double eval(double rho);
protected:
    const GCTAResponse&        m_rsp;           //!< CTA response
    const GModelSpatialRadial& m_model;         //!< Radial spatial model
    const double&              m_zenith;        //!< Zenith angle
    const double&              m_azimuth;       //!< Azimuth angle
    const GEnergy&             m_srcEng;        //!< True photon energy
    const GTime&               m_srcTime;       //!< True photon time
    const double&              m_srcLogEng;     //!< True photon log10 energy
    const double&              m_obsLogEng;     //!< Measured photon energy
    const double&              m_zeta;          //!< Distance model centre - measured photon
    double                     m_cos_zeta;      //!< Cosine of zeta
    double                     m_sin_zeta;      //!< Sine of zeta
    const double&              m_lambda;        //!< Distance model centre - pointing
    double                     m_cos_lambda;    //!< Cosine of lambda
    double                     m_sin_lambda;    //!< Sine of lambda
    const double&              m_omega0;        //!< Azimuth of pointing in model system
    const double&              m_delta_max;     //!< Maximum PSF radius
    double                     m_cos_delta_max; //!< Cosine of maximum PSF radius
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_omega
 *
 * @brief Kernel for radial model azimuth angle IRF integration
 *
 * This class implements the computation of the IRF in the reference frame
 * of the radial source model. It computes
 *
 * \f[
 *    IRF(\rho, \omega)
 * \f]
 *
 * where
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_irf_radial_kern_omega : public GFunction {
public:
    cta_irf_radial_kern_omega(const GCTAResponse& rsp,
                              const double&       zenith,
                              const double&       azimuth,
                              const double&       srcLogEng,
                              const double&       obsLogEng,
                              const double&       zeta,
                              const double&       lambda,
                              const double&       omega0,
                              const double&       rho,
                              const double&       cos_psf,
                              const double&       sin_psf,
                              const double&       cos_ph,
                              const double&       sin_ph) :
                              m_rsp(rsp),
                              m_zenith(zenith),
                              m_azimuth(azimuth),
                              m_srcLogEng(srcLogEng),
                              m_obsLogEng(obsLogEng),
                              m_zeta(zeta),
                              m_lambda(lambda),
                              m_omega0(omega0),
                              m_rho(rho),
                              m_cos_psf(cos_psf),
                              m_sin_psf(sin_psf),
                              m_cos_ph(cos_ph),
                              m_sin_ph(sin_ph) { }
    double eval(double omega);
protected:
    const GCTAResponse& m_rsp;           //!< CTA response
    const double&       m_zenith;        //!< Zenith angle
    const double&       m_azimuth;       //!< Azimuth angle
    const double&       m_srcLogEng;     //!< True photon energy
    const double&       m_obsLogEng;     //!< Measured photon energy
    const double&       m_zeta;          //!< Distance model centre - measured photon
    const double&       m_lambda;        //!< Distance model centre - pointing
    const double&       m_omega0;        //!< Azimuth of pointing in model system
    const double&       m_rho;           //!< ...
    const double&       m_cos_psf;       //!< Cosine term for PSF offset angle computation
    const double&       m_sin_psf;       //!< Sine term for PSF offset angle computation
    const double&       m_cos_ph;        //!< Cosine term for photon offset angle computation
    const double&       m_sin_ph;        //!< Sine term for photon offset angle computation
};


/***********************************************************************//**
 * @class cta_npred_radial_kern_rho
 *
 * @brief Kernel for zenith angle Npred integration of radial model
 *
 * This class implements the integration kernel \f$K(\rho)\f$ for the
 * integration
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho | E, t) d\rho
 * \f]
 *
 * of radial spatial models. The eval() method computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     N_{\rm pred}(\rho,\omega) d\omega
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho | E, t)\f$ is the radial model,
 * - \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest,
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_npred_radial_kern_rho : public GFunction {
public:
    cta_npred_radial_kern_rho(const GCTAResponse&        rsp,
                              const GModelSpatialRadial& model,
                              const GEnergy&             srcEng,
                              const GTime&               srcTime,
                              const GCTAObservation&     obs,
                              const GMatrix&             rot,
                              double                     dist,
                              double                     radius,
                              double                     omega0) :
                              m_rsp(rsp),
                              m_model(model),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_obs(obs),
                              m_rot(rot),
                              m_dist(dist),
                              m_cos_dist(std::cos(dist)),
                              m_sin_dist(std::sin(dist)),
                              m_radius(radius),
                              m_cos_radius(std::cos(radius)),
                              m_omega0(omega0) { }
    double eval(double rho);
protected:
    const GCTAResponse&        m_rsp;        //!< CTA response
    const GModelSpatialRadial& m_model;      //!< Radial spatial model
    const GEnergy&             m_srcEng;     //!< True photon energy
    const GTime&               m_srcTime;    //!< True photon arrival time
    const GCTAObservation&     m_obs;        //!< CTA observation
    const GMatrix&             m_rot;        //!< Rotation matrix
    double                     m_dist;       //!< Distance model-ROI centre
    double                     m_cos_dist;   //!< Cosine of distance model-ROI centre
    double                     m_sin_dist;   //!< Sine of distance model-ROI centre
    double                     m_radius;     //!< ROI+PSF radius
    double                     m_cos_radius; //!< Cosine of ROI+PSF radius
    double                     m_omega0;     //!< Position angle of ROI
};


/***********************************************************************//**
 * @class cta_npred_radial_kern_omega
 *
 * @brief Kernel for azimuth angle Npred integration of radial model
 *
 * This class implements the computation of the data space integral of the
 * Instrument Response Function for a point spread function over the Region
 * Of Interest in the reference frame of the radial source model. It
 * computes
 *
 * \f[
 *    N_{\rm pred}(\rho,\omega)
 * \f]
 *
 * where
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_npred_radial_kern_omega : public GFunction {
public:
    cta_npred_radial_kern_omega(const GCTAResponse&    rsp,
                                const GEnergy&         srcEng,
                                const GTime&           srcTime,
                                const GCTAObservation& obs,
                                const GMatrix&         rot,
                                double                 sin_rho,
                                double                 cos_rho) :
                                m_rsp(rsp),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime),
                                m_obs(obs),
                                m_rot(rot),
                                m_sin_rho(sin_rho),
                                m_cos_rho(cos_rho) { }
    double eval(double omega);
protected:
    const GCTAResponse&    m_rsp;        //!< CTA response
    const GEnergy&         m_srcEng;     //!< True photon energy
    const GTime&           m_srcTime;    //!< True photon arrival time
    const GCTAObservation& m_obs;        //!< CTA observation
    const GMatrix&         m_rot;        //!< Rotation matrix
    double                 m_cos_rho;    //!< Cosine of offset angle
    double                 m_sin_rho;    //!< Sine of offset angle
};


/***********************************************************************//**
 * @class cta_irf_elliptical_kern_rho
 *
 * @brief Kernel for elliptical model zenith angle integration of IRF
 *
 * This class implements the integration kernel \f$K(\rho)\f$ for the
 * integration
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho | E, t) d\rho
 * \f]
 *
 * of elliptical spatial models. The eval() method computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     S_{\rm p}(\rho, \omega | E, t) \, IRF(\rho, \omega)
 *                     d\omega
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical model,
 * - \f$IRF(\rho, \omega)\f$ is the IRF
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_irf_elliptical_kern_rho : public GFunction {
public:
    cta_irf_elliptical_kern_rho(const GCTAResponse&            rsp,
                                const GModelSpatialElliptical& model,
                                const double&                  zenith,
                                const double&                  azimuth,
                                const GEnergy&                 srcEng,
                                const GTime&                   srcTime,
                                const double&                  srcLogEng,
                                const double&                  obsLogEng,
                                const double&                  zeta,
                                const double&                  lambda,
                                const double&                  obsOmega,
                                const double&                  omega0,
                                const double&                  delta_max) :
                                m_rsp(rsp),
                                m_model(model),
                                m_zenith(zenith),
                                m_azimuth(azimuth),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime),
                                m_srcLogEng(srcLogEng),
                                m_obsLogEng(obsLogEng),
                                m_zeta(zeta),
                                m_cos_zeta(std::cos(zeta)),
                                m_sin_zeta(std::sin(zeta)),
                                m_lambda(lambda),
                                m_cos_lambda(std::cos(lambda)),
                                m_sin_lambda(std::sin(lambda)),
                                m_obsOmega(obsOmega),
                                m_omega0(omega0),
                                m_delta_max(delta_max),
                                m_cos_delta_max(std::cos(delta_max)) { }
    double eval(double rho);
public:
    const GCTAResponse&            m_rsp;           //!< CTA response
    const GModelSpatialElliptical& m_model;         //!< Spatial model
    const double&                  m_zenith;        //!< Zenith angle
    const double&                  m_azimuth;       //!< Azimuth angle
    const GEnergy&                 m_srcEng;        //!< True photon energy
    const GTime&                   m_srcTime;       //!< True photon time
    const double&                  m_srcLogEng;     //!< True photon log energy
    const double&                  m_obsLogEng;     //!< Measured photon energy
    const double&                  m_zeta;          //!< Distance model centre - measured photon
    double                         m_cos_zeta;      //!< Cosine of zeta
    double                         m_sin_zeta;      //!< Sine of zeta
    const double&                  m_lambda;        //!< Distance model centre - pointing
    double                         m_cos_lambda;    //!< Cosine of lambda
    double                         m_sin_lambda;    //!< Sine of lambda
    const double&                  m_obsOmega;      //!< Measured photon position angle from model centre
    const double&                  m_omega0;        //!< Azimuth of pointing in model system
    const double&                  m_delta_max;     //!< Maximum PSF radius
    double                         m_cos_delta_max; //!< Cosine of maximum PSF radius
};


/***********************************************************************//**
 * @class cta_irf_elliptical_kern_omega
 *
 * @brief Kernel for ellitpical model azimuth angle IRF integration
 *
 * This class implements the computation of
 *
 * \f[
 *    S_{\rm p}(\rho, \omega | E, t) \, IRF(\rho, \omega)
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical model,
 * - \f$IRF(\rho, \omega)\f$ is the IRF
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_irf_elliptical_kern_omega : public GFunction {
public:
    cta_irf_elliptical_kern_omega(const GCTAResponse&            rsp,
                                  const GModelSpatialElliptical& model,
                                  const double&                  zenith,
                                  const double&                  azimuth,
                                  const GEnergy&                 srcEng,
                                  const GTime&                   srcTime,
                                  const double&                  srcLogEng,
                                  const double&                  obsLogEng,
                                  const double&                  obsOmega,
                                  const double&                  omega0,
                                  const double&                  rho,
                                  const double&                  cos_psf,
                                  const double&                  sin_psf,
                                  const double&                  cos_ph,
                                  const double&                  sin_ph) :
                                  m_rsp(rsp),
                                  m_model(model),
                                  m_zenith(zenith),
                                  m_azimuth(azimuth),
                                  m_srcEng(srcEng),
                                  m_srcTime(srcTime),
                                  m_srcLogEng(srcLogEng),
                                  m_obsLogEng(obsLogEng),
                                  m_obsOmega(obsOmega),
                                  m_omega0(omega0),
                                  m_rho(rho),
                                  m_cos_psf(cos_psf),
                                  m_sin_psf(sin_psf),
                                  m_cos_ph(cos_ph),
                                  m_sin_ph(sin_ph) { }
    double eval(double omega);
public:
    const GCTAResponse&            m_rsp;        //!< CTA response
    const GModelSpatialElliptical& m_model;      //!< Spatial model
    const double&                  m_zenith;     //!< Zenith angle
    const double&                  m_azimuth;    //!< Azimuth angle
    const GEnergy&                 m_srcEng;     //!< True photon energy
    const GTime&                   m_srcTime;    //!< True photon time
    const double&                  m_srcLogEng;  //!< True photon log energy
    const double&                  m_obsLogEng;  //!< Measured photon energy
    const double&                  m_obsOmega;   //!< Measured photon position angle from model centre
    const double&                  m_omega0;     //!< Azimuth of pointing in model system
    const double&                  m_rho;        //!< Model zenith angle
    const double&                  m_cos_psf;    //!< Cosine term for PSF offset angle computation
    const double&                  m_sin_psf;    //!< Sine term for PSF offset angle computation
    const double&                  m_cos_ph;     //!< Cosine term for photon offset angle computation
    const double&                  m_sin_ph;     //!< Sine term for photon offset angle computation
};


/***********************************************************************//**
 * @class cta_npred_elliptical_kern_rho
 *
 * @brief Kernel for zenith angle Npred integration of elliptical model
 *
 * This class implements the integration kernel \f$K(\rho)\f$ for the
 * integration
 *
 * \f[
 *    \int_{\rho_{\rm min}}^{\rho_{\rm max}} K(\rho | E, t) d\rho
 * \f]
 *
 * of elliptical elliptical models. The eval() method computes
 *
 * \f[
 *    K(\rho | E, t) = \sin \rho \times
 *                     \int_{\omega_{\rm min}}^{\omega_{\rm max}} 
 *                     S_{\rm p}(\rho,\omega | E, t) \,
 *                     N_{\rm pred}(\rho,\omega) d\omega
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho,\omega | E, t)\f$ is the elliptical model,
 * - \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest,
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_npred_elliptical_kern_rho : public GFunction {
public:
    cta_npred_elliptical_kern_rho(const GCTAResponse&            rsp,
                                  const GModelSpatialElliptical& model,
                                  const GEnergy&                 srcEng,
                                  const GTime&                   srcTime,
                                  const GCTAObservation&         obs,
                                  const GMatrix&                 rot,
                                  const double&                  dist,
                                  const double&                  radius,
                                  const double&                  omega0) :
                                  m_rsp(rsp),
                                  m_model(model),
                                  m_srcEng(srcEng),
                                  m_srcTime(srcTime),
                                  m_obs(obs),
                                  m_rot(rot),
                                  m_dist(dist),
                                  m_cos_dist(std::cos(dist)),
                                  m_sin_dist(std::sin(dist)),
                                  m_radius(radius),
                                  m_cos_radius(std::cos(radius)),
                                  m_omega0(omega0) { }
    double eval(double rho);
protected:
    const GCTAResponse&            m_rsp;        //!< CTA response
    const GModelSpatialElliptical& m_model;      //!< Spatial model
    const GEnergy&                 m_srcEng;     //!< True photon energy
    const GTime&                   m_srcTime;    //!< True photon arrival time
    const GCTAObservation&         m_obs;        //!< CTA observation
    const GMatrix&                 m_rot;        //!< Rotation matrix
    const double&                  m_dist;       //!< Distance model-ROI centre
    double                         m_cos_dist;   //!< Cosine of distance model-ROI centre
    double                         m_sin_dist;   //!< Sine of distance model-ROI centre
    const double&                  m_radius;     //!< ROI+PSF radius
    double                         m_cos_radius; //!< Cosine of ROI+PSF radius
    const double&                  m_omega0;     //!< Position angle of ROI
};


/***********************************************************************//**
 * @class cta_npred_elliptical_kern_omega
 *
 * @brief Kernel for azimuth angle Npred integration of elliptical model
 *
 * This class implements the computation of
 *
 * \f[
 *    S_{\rm p}(\rho,\omega | E, t) \, N_{\rm pred}(\rho,\omega)
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\rho,\omega | E, t)\f$ is the elliptical model,
 * - \f$N_{\rm pred}(\rho,\omega)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest in the reference frame of the elliptical source
 *   model
 * - \f$\rho\f$ is the distance from the model centre, and
 * - \f$\omega\f$ is the azimuth angle is the position angle with respect to
 *   the connecting line between the model centre and the observed photon
 *   arrival direction.
 ***************************************************************************/
class cta_npred_elliptical_kern_omega : public GFunction {
public:
    cta_npred_elliptical_kern_omega(const GCTAResponse&            rsp,
                                    const GModelSpatialElliptical& model,
                                    const GEnergy&                 srcEng,
                                    const GTime&                   srcTime,
                                    const GCTAObservation&         obs,
                                    const GMatrix&                 rot,
                                    const double&                  sin_rho,
                                    const double&                  cos_rho) :
                                    m_rsp(rsp),
                                    m_model(model),
                                    m_srcEng(srcEng),
                                    m_srcTime(srcTime),
                                    m_obs(obs),
                                    m_rot(rot),
                                    m_sin_rho(sin_rho),
                                    m_cos_rho(cos_rho) { }
    double eval(double omega);
protected:
    const GCTAResponse&            m_rsp;      //!< CTA response
    const GModelSpatialElliptical& m_model;    //!< Model
    const GEnergy&                 m_srcEng;   //!< True photon energy
    const GTime&                   m_srcTime;  //!< True photon arrival time
    const GCTAObservation&         m_obs;      //!< Pointer to observation
    const GMatrix&                 m_rot;      //!< Rotation matrix
    const double&                  m_sin_rho;  //!< Sine of offset angle
    const double&                  m_cos_rho;  //!< Cosine of offset angle
};


/***********************************************************************//**
 * @class cta_irf_diffuse_kern_theta
 *
 * @brief Kernel for IRF offest angle integration of the diffuse source model
 *
 * This class implements the integration kernel \f$K(\theta)\f$ for the
 * integration
 *
 * \f[
 *    \int_{0}^{\theta_{\rm max}} K(\theta | E, t) d\theta
 * \f]
 *
 * of diffuse models. The eval() method computes
 *
 * \f[
 *    K(\theta | E, t) = \sin \theta \times PSF(\theta)
 *                       \int_{0}^{2\pi}
 *                       S_{\rm p}(\theta, \phi | E, t) \,
 *                       Aeff(\theta, \phi) \,
 *                       Edisp(\theta, \phi) d\phi
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$PSF(\theta)\f$ is the azimuthally symmetric Point Spread Function,
 * - \f$Aeff(\theta, \phi)\f$ is the effective area,
 * - \f$Edisp(\theta, \phi)\f$ is the energy dispersion,
 * - \f$\theta\f$ is the distance from the PSF centre, and
 * - \f$\phi\f$ is the azimuth angle.
 ***************************************************************************/
class cta_irf_diffuse_kern_theta : public GFunction {
public:
    cta_irf_diffuse_kern_theta(const GCTAResponse&  rsp,
                               const GModelSpatial& model,
                               const double&        theta,
                               const double&        phi,
                               const double&        zenith,
                               const double&        azimuth,
                               const GEnergy&       srcEng,
                               const GTime&         srcTime,
                               const double&        srcLogEng,
                               const double&        obsLogEng,
                               const GMatrix&       rot,
                               const double&        eta) :
                               m_rsp(rsp),
                               m_model(model),
                               m_theta(theta),
                               m_phi(phi),
                               m_zenith(zenith),
                               m_azimuth(azimuth),
                               m_srcEng(srcEng),
                               m_srcTime(srcTime),
                               m_srcLogEng(srcLogEng),
                               m_obsLogEng(obsLogEng),
                               m_rot(rot),
                               m_sin_eta(std::sin(eta)),
                               m_cos_eta(std::cos(eta)) { }
    double eval(double theta);
protected:
    const GCTAResponse&  m_rsp;        //!< CTA response
    const GModelSpatial& m_model;      //!< Spatial model
    const double&        m_theta;      //!< Photon offset angle
    const double&        m_phi;        //!< Photon azimuth angle
    const double&        m_zenith;     //!< Pointing zenith angle
    const double&        m_azimuth;    //!< Pointing azimuth angle
    const GEnergy&       m_srcEng;     //!< True photon energy
    const GTime&         m_srcTime;    //!< True photon arrival time
    const double&        m_srcLogEng;  //!< True photon log energy
    const double&        m_obsLogEng;  //!< Measured photon energy
    const GMatrix&       m_rot;        //!< Rotation matrix
    double               m_sin_eta;    //!< Sine of angular distance between
                                       //   observed photon direction and
                                       //   camera centre
    double               m_cos_eta;    //!< Cosine of angular distance between
                                       //   observed photon direction and
                                       //   camera centre
};


/***********************************************************************//**
 * @class cta_irf_diffuse_kern_phi
 *
 * @brief Kernel for IRF azimuth angle integration of the diffuse source model
 *
 * This class implements the computation of
 *
 * \f[
 *    S_{\rm p}(\theta, \phi | E, t) \,
 *    Aeff(\theta, \phi) \,
 *    Edisp(\theta, \phi)
 * \f]
 *
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$Aeff(\theta, \phi)\f$ is the effective area,
 * - \f$Edisp(\theta, \phi)\f$ is the energy dispersion,
 * - \f$\theta\f$ is the distance from the PSF centre, and
 * - \f$\phi\f$ is the azimuth angle.
 ***************************************************************************/
class cta_irf_diffuse_kern_phi : public GFunction {
public:
    cta_irf_diffuse_kern_phi(const GCTAResponse&  rsp,
                             const GModelSpatial& model,
                             const double&        zenith,
                             const double&        azimuth,
                             const GEnergy&       srcEng,
                             const GTime&         srcTime,
                             const double&        srcLogEng,
                             const double&        obsLogEng,
                             const GMatrix&       rot,
                             const double&        sin_theta,
                             const double&        cos_theta,
                             const double&        sin_ph,
                             const double&        cos_ph) :
                             m_rsp(rsp),
                             m_model(model),
                             m_zenith(zenith),
                             m_azimuth(azimuth),
                             m_srcEng(srcEng),
                             m_srcTime(srcTime),
                             m_srcLogEng(srcLogEng),
                             m_obsLogEng(obsLogEng),
                             m_rot(rot),
                             m_sin_theta(sin_theta),
                             m_cos_theta(cos_theta),
                             m_sin_ph(sin_ph),
                             m_cos_ph(cos_ph) { }
    double eval(double phi);
protected:
    const GCTAResponse&  m_rsp;        //!< CTA response
    const GModelSpatial& m_model;      //!< Spatial model
    const double&        m_zenith;     //!< Zenith angle
    const double&        m_azimuth;    //!< Azimuth angle
    const GEnergy&       m_srcEng;     //!< True photon energy
    const GTime&         m_srcTime;    //!< True photon arrival time
    const double&        m_srcLogEng;  //!< True photon log energy
    const double&        m_obsLogEng;  //!< Measured photon energy
    const GMatrix&       m_rot;        //!< Rotation matrix
    const double&        m_sin_theta;  //!< Sine of offset angle
    const double&        m_cos_theta;  //!< Cosine of offset angle
    const double&        m_sin_ph;     //!< Sine term in angular distance equation
    const double&        m_cos_ph;     //!< Cosine term in angular distance equation    
};


/***********************************************************************//**
 * @class cta_npred_diffuse_kern_theta
 *
 * @brief Kernel for Npred offest angle integration of diffuse model
 *
 * This class implements the integration kernel \f$K(\theta)\f$ for the
 * integration
 *
 * \f[
 *    \int_{0}^{\theta_{\rm max}} K(\theta | E, t) d\theta
 * \f]
 *
 * of diffuse models. The eval() method computes
 *
 * \f[
 *    K(\theta | E, t) = \sin \theta \times
 *                       \int_{0}^{2\pi}
 *                       S_{\rm p}(\theta, \phi | E, t) \,
 *                       N_{\rm pred}(\theta, \phi) d\phi
 * \f]
 * 
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$N_{\rm pred}(\theta, \phi)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest in the reference frame of the diffuse source
 *   model
 * - \f$\theta\f$ is the distance from the ROI centre, and
 * - \f$\phi\f$ is the azimuth angle.
 ***************************************************************************/
class cta_npred_diffuse_kern_theta : public GFunction {
public:
    cta_npred_diffuse_kern_theta(const GCTAResponse&    rsp,
                                 const GModelSpatial&   model,
                                 const GEnergy&         srcEng,
                                 const GTime&           srcTime,
                                 const GCTAObservation& obs,
                                 const GMatrix&         rot) :
                                 m_rsp(rsp),
                                 m_model(model),
                                 m_srcEng(srcEng),
                                 m_srcTime(srcTime),
                                 m_obs(obs),
                                 m_rot(rot) { }
    double eval(double theta);
protected:
    const GCTAResponse&    m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GEnergy&         m_srcEng;     //!< True photon energy
    const GTime&           m_srcTime;    //!< True photon arrival time
    const GCTAObservation& m_obs;        //!< CTA observation
    const GMatrix&         m_rot;        //!< Rotation matrix
};


/***********************************************************************//**
 * @class cta_npred_diffuse_kern_phi
 *
 * @brief Kernel for Npred azimuth angle integration of diffuse model
 *
 * This class implements the computation of
 *
 * \f[
 *    S_{\rm p}(\theta, \phi | E, t) \, N_{\rm pred}(\theta, \phi)
 * \f]
 * 
 * where
 * - \f$S_{\rm p}(\theta, \phi | E, t)\f$ is the diffuse model,
 * - \f$N_{\rm pred}(\theta, \phi)\f$ is the data space integral of the
 *   Instrument Response Function for a point spread function over the
 *   Region Of Interest in the reference frame of the diffuse source
 *   model
 * - \f$\theta\f$ is the distance from the ROI centre, and
 * - \f$\phi\f$ is the azimuth angle.
 ***************************************************************************/
class cta_npred_diffuse_kern_phi : public GFunction {
public:
    cta_npred_diffuse_kern_phi(const GCTAResponse&    rsp,
                               const GModelSpatial&   model,
                               const GEnergy&         srcEng,
                               const GTime&           srcTime,
                               const GCTAObservation& obs,
                               const GMatrix&         rot,
                               const double&          theta,
                               const double&          sin_theta) :
                               m_rsp(rsp),
                               m_model(model),
                               m_srcEng(srcEng),
                               m_srcTime(srcTime),
                               m_obs(obs),
                               m_rot(rot),
                               m_theta(theta),
                               m_cos_theta(std::cos(theta)),
                               m_sin_theta(sin_theta) { }
    double eval(double phi);
protected:
    const GCTAResponse&    m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GEnergy&         m_srcEng;     //!< True photon energy
    const GTime&           m_srcTime;    //!< True photon arrival time
    const GCTAObservation& m_obs;        //!< CTA observation
    const GMatrix&         m_rot;        //!< Rotation matrix
    const double&          m_theta;      //!< Offset angle (radians)
    double                 m_cos_theta;  //!< Cosine of offset angle
    const double&          m_sin_theta;  //!< Sine of offset angle
};

#endif /* GCTARESPONSE_HELPERS_HPP */
