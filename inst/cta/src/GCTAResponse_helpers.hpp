/***************************************************************************
 *        GCTAResponse_helpers.hpp  -  CTA response helper classes         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @author J. Knoedlseder
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
#include "GModelRadial.hpp"
#include "GIntegrand.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */


/***********************************************************************//**
 * @class cta_irf_radial_kern_rho
 *
 * @brief Integration kernel for npsf() method
 ***************************************************************************/
class cta_npsf_kern_rad_azsym : public GIntegrand {
public:
    cta_npsf_kern_rad_azsym(const       GCTAResponse* rsp,
                            double      roi,
                            double      psf,
                            GCTAPsfPars pars) :
                            m_rsp(rsp),
                            m_roi(roi),
                            m_cosroi(std::cos(roi)),
                            m_psf(psf),
                            m_cospsf(std::cos(psf)),
                            m_sinpsf(std::sin(psf)),
                            m_pars(pars) { }
    double eval(double r);
protected:
    const GCTAResponse* m_rsp;     //!< Pointer to response function
    double              m_roi;     //!< ROI radius in radians
    double              m_cosroi;  //!< Cosine of ROI radius
    double              m_psf;     //!< PSF-ROI centre distance in radians
    double              m_cospsf;  //!< Cosine of PSF-ROI centre distance
    double              m_sinpsf;  //!< Sine of PSF-ROI centre distance
    GCTAPsfPars         m_pars;    //!< PSF parameters
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_rho
 *
 * @brief Kernel for radial model zenith angle integration of IRF
 ***************************************************************************/
class cta_irf_radial_kern_rho : public GIntegrand {
public:
    cta_irf_radial_kern_rho(const GCTAResponse* rsp,
                            const GModelRadial* model,
                            double              zenith,
                            double              azimuth,
                            double              srcLogEng,
                            double              obsLogEng,
                            GCTAPsfPars         pars,
                            double              zeta,
                            double              lambda,
                            double              omega0,
                            double              delta_max) :
                            m_rsp(rsp),
                            m_model(model),
                            m_zenith(zenith),
                            m_azimuth(azimuth),
                            m_srcLogEng(srcLogEng),
                            m_obsLogEng(obsLogEng),
                            m_pars(pars),
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
    const GCTAResponse* m_rsp;           //!< Pointer to CTA response
    const GModelRadial* m_model;         //!< Pointer to radial model
    double              m_zenith;        //!< Pointing zenith angle
    double              m_azimuth;       //!< Pointing azimuth angle
    double              m_srcLogEng;     //!< True photon energy
    double              m_obsLogEng;     //!< Measured photon energy
    GCTAPsfPars         m_pars;          //!< PSF parameters
    double              m_zeta;          //!< Distance model centre - measured photon
    double              m_cos_zeta;      //!< Cosine of zeta
    double              m_sin_zeta;      //!< Sine of zeta
    double              m_lambda;        //!< Distance model centre - pointing
    double              m_cos_lambda;    //!< Cosine of lambda
    double              m_sin_lambda;    //!< Sine of lambda
    double              m_omega0;        //!< Azimuth of pointing in model system
    double              m_delta_max;     //!< Maximum PSF radius
    double              m_cos_delta_max; //!< Cosine of maximum PSF radius
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_omega
 *
 * @brief Kernel for radial model azimuth angle IRF integration
 ***************************************************************************/
class cta_irf_radial_kern_omega : public GIntegrand {
public:
    cta_irf_radial_kern_omega(const GCTAResponse* rsp,
                              double              zenith,
                              double              azimuth,
                              double              srcLogEng,
                              double              obsLogEng,
                              GCTAPsfPars         pars,
                              double              zeta,
                              double              lambda,
                              double              omega0,
                              double              rho,
                              double              cos_psf,
                              double              sin_psf,
                              double              cos_ph,
                              double              sin_ph) :
                              m_rsp(rsp),
                              m_zenith(zenith),
                              m_azimuth(azimuth),
                              m_srcLogEng(srcLogEng),
                              m_obsLogEng(obsLogEng),
                              m_pars(pars),
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
    const GCTAResponse* m_rsp;           //!< Pointer to CTA response
    double              m_zenith;        //!< Pointing zenith angle
    double              m_azimuth;       //!< Pointing azimuth angle
    double              m_srcLogEng;     //!< True photon energy
    double              m_obsLogEng;     //!< Measured photon energy
    GCTAPsfPars         m_pars;          //!< PSF parameters
    double              m_zeta;          //!< Distance model centre - measured photon
    double              m_lambda;        //!< Distance model centre - pointing
    double              m_omega0;        //!< Azimuth of pointing in model system
    double              m_rho;           //!< ...
    double              m_cos_psf;       //!< Cosine term for PSF offset angle computation
    double              m_sin_psf;       //!< Sine term for PSF offset angle computation
    double              m_cos_ph;        //!< Cosine term for photon offset angle computation
    double              m_sin_ph;        //!< Sine term for photon offset angle computation
};


/***********************************************************************//**
 * @class cta_npred_radial_kern_theta
 *
 * @brief Kernel for zenith angle Npred integration of radial model
 ***************************************************************************/
class cta_npred_radial_kern_theta : public GIntegrand {
public:
    cta_npred_radial_kern_theta(const GCTAResponse*    rsp,
                                const GModelRadial*    model,
                                const GEnergy*         srcEng,
                                const GTime*           srcTime,
                                const GCTAObservation* obs,
                                const GMatrix*         rot,
                                double                 dist,
                                double                 radius,
                                double                 phi) :
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
                                m_phi(phi) { }
    double eval(double theta);
protected:
    const GCTAResponse*    m_rsp;        //!< Pointer to response
    const GModelRadial*    m_model;      //!< Pointer to radial model
    const GEnergy*         m_srcEng;     //!< Pointer to true photon energy
    const GTime*           m_srcTime;    //!< Pointer to true photon arrival time
    const GCTAObservation* m_obs;        //!< Pointer to observation
    const GMatrix*         m_rot;        //!< Rotation matrix
    double                 m_dist;       //!< Distance model-ROI centre
    double                 m_cos_dist;   //!< Cosine of distance model-ROI centre
    double                 m_sin_dist;   //!< Sine of distance model-ROI centre
    double                 m_radius;     //!< ROI+PSF radius
    double                 m_cos_radius; //!< Cosine of ROI+PSF radius
    double                 m_phi;        //!< Position angle of ROI
};


/***********************************************************************//**
 * @class cta_npred_radial_kern_phi
 *
 * @brief Kernel for azimuth angle Npred integration of radial model
 ***************************************************************************/
class cta_npred_radial_kern_phi : public GIntegrand {
public:
    cta_npred_radial_kern_phi(const GCTAResponse*    rsp,
                              const GEnergy*         srcEng,
                              const GTime*           srcTime,
                              const GCTAObservation* obs,
                              const GMatrix*         rot,
                              double                 theta,
                              double                 sin_theta) :
                              m_rsp(rsp),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_obs(obs),
                              m_rot(rot),
                              m_theta(theta),
                              m_cos_theta(std::cos(theta)),
                              m_sin_theta(sin_theta) { }
    double eval(double phi);
protected:
    const GCTAResponse*    m_rsp;        //!< Pointer to response
    const GEnergy*         m_srcEng;     //!< Pointer to true photon energy
    const GTime*           m_srcTime;    //!< Pointer to true photon arrival time
    const GCTAObservation* m_obs;        //!< Pointer to observation
    const GMatrix*         m_rot;        //!< Rotation matrix
    double                 m_theta;      //!< Offset angle (radians)
    double                 m_cos_theta;  //!< Cosine of offset angle
    double                 m_sin_theta;  //!< Sine of offset angle
};


/***********************************************************************//**
 * @class cta_irf_diffuse_kern_theta
 *
 * @brief Kernel for IRF offest angle integration of the diffuse source model
 ***************************************************************************/
class cta_irf_diffuse_kern_theta : public GIntegrand {
public:
    cta_irf_diffuse_kern_theta(const GCTAResponse*  rsp,
                               const GModelSpatial* model,
                               double               zenith,
                               double               azimuth,
                               double               srcLogEng,
                               double               obsLogEng,
                               const GMatrix*       rot,
                               GCTAPsfPars          pars,
                               double               eta) :
                               m_rsp(rsp),
                               m_model(model),
                               m_zenith(zenith),
                               m_azimuth(azimuth),
                               m_srcLogEng(srcLogEng),
                               m_obsLogEng(obsLogEng),
                               m_rot(rot),
                               m_pars(pars),
                               m_sin_eta(std::sin(eta)),
                               m_cos_eta(std::cos(eta)) { }
    double eval(double theta);
protected:
    const GCTAResponse*  m_rsp;        //!< Pointer to CTA response
    const GModelSpatial* m_model;      //!< Pointer to spatial model
    double               m_zenith;     //!< Pointing zenith angle
    double               m_azimuth;    //!< Pointing azimuth angle
    double               m_srcLogEng;  //!< True photon energy
    double               m_obsLogEng;  //!< Measured photon energy
    const GMatrix*       m_rot;        //!< Rotation matrix
    GCTAPsfPars          m_pars;       //!< PSF parameters
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
 ***************************************************************************/
class cta_irf_diffuse_kern_phi : public GIntegrand {
public:
    cta_irf_diffuse_kern_phi(const GCTAResponse*  rsp,
                             const GModelSpatial* model,
                             double               zenith,
                             double               azimuth,
                             double               srcLogEng,
                             double               obsLogEng,
                             const GMatrix*       rot,
                             double               sin_theta,
                             double               cos_theta,
                             double               sin_ph,
                             double               cos_ph) :
                             m_rsp(rsp),
                             m_model(model),
                             m_zenith(zenith),
                             m_azimuth(azimuth),
                             m_srcLogEng(srcLogEng),
                             m_obsLogEng(obsLogEng),
                             m_rot(rot),
                             m_sin_theta(sin_theta),
                             m_cos_theta(cos_theta),
                             m_sin_ph(sin_ph),
                             m_cos_ph(cos_ph) { }
    double eval(double phi);
protected:
    const GCTAResponse*  m_rsp;        //!< Pointer to CTA response
    const GModelSpatial* m_model;      //!< Pointer to spatial model
    double               m_zenith;     //!< Pointing zenith angle
    double               m_azimuth;    //!< Pointing azimuth angle
    double               m_srcLogEng;  //!< True photon energy
    double               m_obsLogEng;  //!< Measured photon energy
    const GMatrix*       m_rot;        //!< Rotation matrix
    double               m_sin_theta;  //!< Sine of offset angle
    double               m_cos_theta;  //!< Cosine of offset angle
    double               m_sin_ph;     //!< Sine term in angular distance equation
    double               m_cos_ph;     //!< Cosine term in angular distance equation    
};


/***********************************************************************//**
 * @class cta_npred_diffuse_kern_theta
 *
 * @brief Kernel for Npred offest angle integration of diffuse model
 ***************************************************************************/
class cta_npred_diffuse_kern_theta : public GIntegrand {
public:
    cta_npred_diffuse_kern_theta(const GCTAResponse*    rsp,
                                 const GModelSpatial*   model,
                                 const GEnergy*         srcEng,
                                 const GTime*           srcTime,
                                 const GCTAObservation* obs,
                                 const GMatrix*         rot) :
                                 m_rsp(rsp),
                                 m_model(model),
                                 m_srcEng(srcEng),
                                 m_srcTime(srcTime),
                                 m_obs(obs),
                                 m_rot(rot) { }
    double eval(double theta);
protected:
    const GCTAResponse*    m_rsp;        //!< Pointer to response
    const GModelSpatial*   m_model;      //!< Pointer to spatial model
    const GEnergy*         m_srcEng;     //!< Pointer to true photon energy
    const GTime*           m_srcTime;    //!< Pointer to true photon arrival time
    const GCTAObservation* m_obs;        //!< Pointer to observation
    const GMatrix*         m_rot;        //!< Rotation matrix
};


/***********************************************************************//**
 * @class cta_npred_diffuse_kern_phi
 *
 * @brief Kernel for Npred azimuth angle integration of diffuse model
 ***************************************************************************/
class cta_npred_diffuse_kern_phi : public GIntegrand {
public:
    cta_npred_diffuse_kern_phi(const GCTAResponse*    rsp,
                               const GModelSpatial*   model,
                               const GEnergy*         srcEng,
                               const GTime*           srcTime,
                               const GCTAObservation* obs,
                               const GMatrix*         rot,
                               double                 theta,
                               double                 sin_theta) :
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
    const GCTAResponse*    m_rsp;        //!< Pointer to response
    const GModelSpatial*   m_model;      //!< Pointer to  spatial model
    const GEnergy*         m_srcEng;     //!< Pointer to true photon energy
    const GTime*           m_srcTime;    //!< Pointer to true photon arrival time
    const GCTAObservation* m_obs;        //!< Pointer to observation
    const GMatrix*         m_rot;        //!< Rotation matrix
    double                 m_theta;      //!< Offset angle (radians)
    double                 m_cos_theta;  //!< Cosine of offset angle
    double                 m_sin_theta;  //!< Sine of offset angle
};

#endif /* GCTARESPONSE_HELPERS_HPP */
