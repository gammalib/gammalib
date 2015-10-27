/***************************************************************************
 *         GCTAResponse_helpers.hpp - CTA response helper classes          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
#include <vector>
#include <utility>
#include "GCTAResponseIrf.hpp"
#include "GCTAObservation.hpp"
#include "GMatrix.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GModelSky.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpatialRadial.hpp"
#include "GCTAResponseCube.hpp"
#include "GModelSpatialElliptical.hpp"
#include "GFunction.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */

/* __ Typedefs ___________________________________________________________ */
typedef std::vector<std::pair<double,double> > cta_omega_intervals;

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    cta_omega_intervals limit_omega(const double& min,
                                    const double& max,
                                    const double& domega);
}


/***********************************************************************//**
 * @class cta_npsf_kern_rad_azsym
 *
 * @brief Integration kernel for npsf() method
 *
 * This class implements the integration kernel needed for the npsf() method.
 * It performs an azimuthal integration of the Point Spread Function
 * (PSF) for a given offset angle \f$delta\f$ using
 *
 * \f[
 *    K(\delta) = \int_{0}^{\phi} PSF(\delta) d\phi
 * \f]
 * (an offset angle of \f$\delta=0\f$ corresponds to the centre of the PSF).
 ***************************************************************************/
class cta_npsf_kern_rad_azsym : public GFunction {
public:
    cta_npsf_kern_rad_azsym(const GCTAResponseIrf& rsp,
                            const double&          roi,
                            const double&          psf,
                            const double&          logE,
                            const double&          theta,
                            const double&          phi,
                            const double&          zenith,
                            const double&          azimuth) :
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
    double eval(const double& delta);
protected:
    const GCTAResponseIrf& m_rsp;     //!< CTA response function
    double                 m_roi;     //!< ROI radius in radians
    double                 m_cosroi;  //!< Cosine of ROI radius
    double                 m_psf;     //!< PSF-ROI centre distance in radians
    double                 m_cospsf;  //!< Cosine of PSF-ROI centre distance
    double                 m_sinpsf;  //!< Sine of PSF-ROI centre distance
    double                 m_logE;    //!< Log10 of true photon energy (E/TeV).
    double                 m_theta;   //!< Offset angle of source in camera system
    double                 m_phi;     //!< Azimuth angle of source in camera system
    double                 m_zenith;  //!< Zenith angle of source in Earth system
    double                 m_azimuth; //!< Azimuth angle of source in Earth system
};


/***********************************************************************//**
 * @class cta_nedisp_kern
 *
 * @brief Integration kernel for nedisp() method
 *
 * This class implements the integration kernel for the nedisp() method.
 * The cta_nedisp_kern::eval method evaluates the energy dispersion for a
 * given log10 of true photon energy as function of the log10 of observed
 * event energy.
 ***************************************************************************/
class cta_nroi_kern : public GFunction {
public:
    cta_nroi_kern(const GModelSky&       model,
                  const GCTAResponseIrf& rsp,
                  const GTime&           srcTime,
                  const GEnergy&         obsEng,
                  const GTime&           obsTime,
                  const GObservation&    obs) :
                  m_model(model),
                  m_rsp(rsp),
                  m_srcTime(srcTime),
                  m_obsEng(obsEng),
                  m_obsTime(obsTime),
                  m_obs(obs) {}
    double eval(const double& logEsrc);
protected:
    const GModelSky&       m_model;
    const GCTAResponseIrf& m_rsp;     //!< CTA response function
    const GObservation&    m_obs;
    GTime                  m_srcTime;
    GEnergy                m_obsEng;
    GTime                  m_obsTime;
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
    cta_irf_radial_kern_rho(const GCTAResponseIrf&     rsp,
                            const GModelSpatialRadial& model,
                            const double&              zenith,
                            const double&              azimuth,
                            const GEnergy&             srcEng,
                            const GTime&               srcTime,
                            const double&              srcLogEng,
                            const GEnergy&             obsEng,
                            const double&              zeta,
                            const double&              lambda,
                            const double&              omega0,
                            const double&              delta_max,
                            const int&                 iter) :
                            m_rsp(rsp),
                            m_model(model),
                            m_zenith(zenith),
                            m_azimuth(azimuth),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime),
                            m_srcLogEng(srcLogEng),
                            m_obsEng(obsEng),
                            m_zeta(zeta),
                            m_cos_zeta(std::cos(zeta)),
                            m_sin_zeta(std::sin(zeta)),
                            m_lambda(lambda),
                            m_cos_lambda(std::cos(lambda)),
                            m_sin_lambda(std::sin(lambda)),
                            m_omega0(omega0),
                            m_delta_max(delta_max),
                            m_cos_delta_max(std::cos(delta_max)),
                            m_iter(iter) { }
    double eval(const double& rho);
protected:
    const GCTAResponseIrf&     m_rsp;           //!< CTA response
    const GModelSpatialRadial& m_model;         //!< Radial spatial model
    double                     m_zenith;        //!< Zenith angle
    double                     m_azimuth;       //!< Azimuth angle
    GEnergy                    m_srcEng;        //!< True photon energy
    GTime                      m_srcTime;       //!< True photon time
    double                     m_srcLogEng;     //!< True photon log10 energy
    GEnergy                    m_obsEng;        //!< Measured event energy
    double                     m_zeta;          //!< Distance model centre - measured photon
    double                     m_cos_zeta;      //!< Cosine of zeta
    double                     m_sin_zeta;      //!< Sine of zeta
    double                     m_lambda;        //!< Distance model centre - pointing
    double                     m_cos_lambda;    //!< Cosine of lambda
    double                     m_sin_lambda;    //!< Sine of lambda
    double                     m_omega0;        //!< Azimuth of pointing in model system
    double                     m_delta_max;     //!< Maximum PSF radius
    double                     m_cos_delta_max; //!< Cosine of maximum PSF radius
    int                        m_iter;          //!< Integration iterations
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
    cta_irf_radial_kern_omega(const GCTAResponseIrf& rsp,
                              const double&          zenith,
                              const double&          azimuth,
                              const double&          srcLogEng,
                              const GEnergy&         obsEng,
                              const double&          zeta,
                              const double&          lambda,
                              const double&          omega0,
                              const double&          rho,
                              const double&          cos_psf,
                              const double&          sin_psf,
                              const double&          cos_ph,
                              const double&          sin_ph) :
                              m_rsp(rsp),
                              m_zenith(zenith),
                              m_azimuth(azimuth),
                              m_srcLogEng(srcLogEng),
                              m_obsEng(obsEng),
                              m_zeta(zeta),
                              m_lambda(lambda),
                              m_omega0(omega0),
                              m_rho(rho),
                              m_cos_psf(cos_psf),
                              m_sin_psf(sin_psf),
                              m_cos_ph(cos_ph),
                              m_sin_ph(sin_ph) { }
    double eval(const double& omega);
protected:
    const GCTAResponseIrf& m_rsp;       //!< CTA response
    double                 m_zenith;    //!< Zenith angle
    double                 m_azimuth;   //!< Azimuth angle
    double                 m_srcLogEng; //!< True photon energy
    GEnergy                m_obsEng;    //!< Measured event energy
    double                 m_zeta;      //!< Distance model centre - measured photon
    double                 m_lambda;    //!< Distance model centre - pointing
    double                 m_omega0;    //!< Azimuth of pointing in model system
    double                 m_rho;       //!< ...
    double                 m_cos_psf;   //!< Cosine term for PSF offset angle computation
    double                 m_sin_psf;   //!< Sine term for PSF offset angle computation
    double                 m_cos_ph;    //!< Cosine term for photon offset angle computation
    double                 m_sin_ph;    //!< Sine term for photon offset angle computation
};


/***********************************************************************//**
 * @class cta_nroi_radial_kern_rho
 *
 * @brief Kernel for zenith angle Nroi integration of radial model
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
class cta_nroi_radial_kern_rho : public GFunction {
public:
    cta_nroi_radial_kern_rho(const GCTAResponseIrf&     rsp,
                             const GModelSpatialRadial& model,
                             const GEnergy&             srcEng,
                             const GTime&               srcTime,
                             const GEnergy&             obsEng,
                             const GTime&               obsTime,
                             const GCTAObservation&     obs,
                             const GMatrix&             rot,
                             const double&              dist,
                             const double&              radius,
                             const double&              omega0,
                             const int&                 iter) :
                             m_rsp(rsp),
                             m_model(model),
                             m_srcEng(srcEng),
                             m_srcTime(srcTime),
                             m_obsEng(obsEng),
                             m_obsTime(obsTime),
                             m_obs(obs),
                             m_rot(rot),
                             m_dist(dist),
                             m_cos_dist(std::cos(dist)),
                             m_sin_dist(std::sin(dist)),
                             m_radius(radius),
                             m_cos_radius(std::cos(radius)),
                             m_omega0(omega0),
                             m_iter(iter) { }
    double eval(const double& rho);
protected:
    const GCTAResponseIrf&     m_rsp;        //!< CTA response
    const GModelSpatialRadial& m_model;      //!< Radial spatial model
    const GCTAObservation&     m_obs;        //!< CTA observation
    const GMatrix&             m_rot;        //!< Rotation matrix
    GEnergy                    m_srcEng;     //!< True photon energy
    GTime                      m_srcTime;    //!< True photon arrival time
    GEnergy                    m_obsEng;     //!< Observed photon energy
    GTime                      m_obsTime;    //!< Observed photon arrival time
    double                     m_dist;       //!< Distance model-ROI centre
    double                     m_cos_dist;   //!< Cosine of distance model-ROI centre
    double                     m_sin_dist;   //!< Sine of distance model-ROI centre
    double                     m_radius;     //!< ROI+PSF radius
    double                     m_cos_radius; //!< Cosine of ROI+PSF radius
    double                     m_omega0;     //!< Position angle of ROI
    int                        m_iter;       //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_roi_radial_kern_omega
 *
 * @brief Kernel for azimuth angle Nroi integration of radial model
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
class cta_nroi_radial_kern_omega : public GFunction {
public:
    cta_nroi_radial_kern_omega(const GCTAResponseIrf& rsp,
                               const GEnergy&         srcEng,
                               const GTime&           srcTime,
                               const GEnergy&         obsEng,
                               const GTime&           obsTime,
                               const GCTAObservation& obs,
                               const GMatrix&         rot,
                               double                 sin_rho,
                               double                 cos_rho) :
                               m_rsp(rsp),
                               m_srcEng(srcEng),
                               m_srcTime(srcTime),
                               m_obsEng(obsEng),
                               m_obsTime(obsTime),
                               m_obs(obs),
                               m_rot(rot),
                               m_cos_rho(cos_rho),
                               m_sin_rho(sin_rho) { }
    double eval(const double& omega);
protected:
    const GCTAResponseIrf& m_rsp;     //!< CTA response
    const GCTAObservation& m_obs;     //!< CTA observation
    GMatrix                m_rot;     //!< Rotation matrix
    GEnergy                m_srcEng;  //!< True photon energy
    GTime                  m_srcTime; //!< True photon arrival time
    GEnergy                m_obsEng;  //!< Observed photon energy
    GTime                  m_obsTime; //!< Observed photon arrival time
    double                 m_cos_rho; //!< Cosine of offset angle
    double                 m_sin_rho; //!< Sine of offset angle
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
    cta_irf_elliptical_kern_rho(const GCTAResponseIrf&         rsp,
                                const GModelSpatialElliptical& model,
                                const double&                  semimajor,
                                const double&                  semiminor,
                                const double&                  posangle,
                                const double&                  zenith,
                                const double&                  azimuth,
                                const GEnergy&                 srcEng,
                                const GTime&                   srcTime,
                                const double&                  srcLogEng,
                                const GEnergy&                 obsEng,
                                const double&                  rho_obs,
                                const double&                  posangle_obs,
                                const double&                  rho_pnt,
                                const double&                  omega_pnt,
                                const double&                  delta_max,
                                const int&                     iter) :
                                m_rsp(rsp),
                                m_model(model),
                                m_semimajor(semimajor),
                                m_semiminor(semiminor),
                                m_posangle(posangle),
                                m_zenith(zenith),
                                m_azimuth(azimuth),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime),
                                m_srcLogEng(srcLogEng),
                                m_obsEng(obsEng),
                                m_rho_obs(rho_obs),
                                m_cos_rho_obs(std::cos(rho_obs)),
                                m_sin_rho_obs(std::sin(rho_obs)),
                                m_posangle_obs(posangle_obs),
                                m_rho_pnt(rho_pnt),
                                m_cos_rho_pnt(std::cos(rho_pnt)),
                                m_sin_rho_pnt(std::sin(rho_pnt)),
                                m_omega_pnt(omega_pnt),
                                m_delta_max(delta_max),
                                m_cos_delta_max(std::cos(delta_max)),
                                m_iter(iter) { }
    double eval(const double& rho);
public:
    const GCTAResponseIrf&         m_rsp;           //!< CTA response
    const GModelSpatialElliptical& m_model;         //!< Elliptical model
    double                         m_semimajor;     //!< Ellipse boundary semimajor axis
    double                         m_semiminor;     //!< Ellipse boundary semiminor axis
    double                         m_posangle;      //!< Ellipse boundary position angle
    double                         m_zenith;        //!< Zenith angle
    double                         m_azimuth;       //!< Azimuth angle
    GEnergy                        m_srcEng;        //!< True photon energy
    GTime                          m_srcTime;       //!< True photon time
    double                         m_srcLogEng;     //!< True photon log energy
    GEnergy                        m_obsEng;        //!< Measured event energy
    double                         m_rho_obs;       //!< Distance of model centre from measured photon
    double                         m_cos_rho_obs;   //!< Cosine of m_rho_obs
    double                         m_sin_rho_obs;   //!< Sine of m_rho_obs
    double                         m_posangle_obs;  //!< Photon position angle measured from model centre
    double                         m_rho_pnt;       //!< Distance of model centre from pointing
    double                         m_cos_rho_pnt;   //!< Cosine of m_rho_pnt
    double                         m_sin_rho_pnt;   //!< Sine of m_rho_pnt
    double                         m_omega_pnt;     //!< Azimuth of pointing in model system
    double                         m_delta_max;     //!< Maximum PSF radius
    double                         m_cos_delta_max; //!< Cosine of maximum PSF radius
    int                            m_iter;          //!< Integration iterations
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
    cta_irf_elliptical_kern_omega(const GCTAResponseIrf&         rsp,
                                  const GModelSpatialElliptical& model,
                                  const double&                  zenith,
                                  const double&                  azimuth,
                                  const GEnergy&                 srcEng,
                                  const GTime&                   srcTime,
                                  const double&                  srcLogEng,
                                  const GEnergy&                 obsEng,
                                  const double&                  posangle_obs,
                                  const double&                  omega_pnt,
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
                                  m_obsEng(obsEng),
                                  m_posangle_obs(posangle_obs),
                                  m_omega_pnt(omega_pnt),
                                  m_rho(rho),
                                  m_cos_psf(cos_psf),
                                  m_sin_psf(sin_psf),
                                  m_cos_ph(cos_ph),
                                  m_sin_ph(sin_ph) { }
    double eval(const double& omega);
public:
    const GCTAResponseIrf&         m_rsp;          //!< CTA response
    const GModelSpatialElliptical& m_model;        //!< Spatial model
    double                         m_zenith;       //!< Zenith angle
    double                         m_azimuth;      //!< Azimuth angle
    GEnergy                        m_srcEng;       //!< True photon energy
    GTime                          m_srcTime;      //!< True photon time
    double                         m_srcLogEng;    //!< True photon log energy
    GEnergy                        m_obsEng;       //!< Measured event energy
    double                         m_posangle_obs; //!< Measured photon position angle from model centre
    double                         m_omega_pnt;    //!< Azimuth of pointing in model system
    double                         m_rho;          //!< Model zenith angle
    double                         m_cos_psf;      //!< Cosine term for PSF offset angle computation
    double                         m_sin_psf;      //!< Sine term for PSF offset angle computation
    double                         m_cos_ph;       //!< Cosine term for photon offset angle computation
    double                         m_sin_ph;       //!< Sine term for photon offset angle computation
};


/***********************************************************************//**
 * @class cta_nroi_elliptical_kern_rho
 *
 * @brief Kernel for zenith angle Nroi integration of elliptical model
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
class cta_nroi_elliptical_kern_rho : public GFunction {
public:
    cta_nroi_elliptical_kern_rho(const GCTAResponseIrf&         rsp,
                                 const GModelSpatialElliptical& model,
                                 const double&                  semimajor,
                                 const double&                  semiminor,
                                 const double&                  posangle,
                                 const GEnergy&                 srcEng,
                                 const GTime&                   srcTime,
                                 const GEnergy&                 obsEng,
                                 const GTime&                   obsTime,
                                 const GCTAObservation&         obs,
                                 const GMatrix&                 rot,
                                 const double&                  rho_roi,
                                 const double&                  posangle_roi,
                                 const double&                  radius_roi,
                                 const int&                     iter) :
                                 m_rsp(rsp),
                                 m_model(model),
                                 m_semimajor(semimajor),
                                 m_semiminor(semiminor),
                                 m_posangle(posangle),
                                 m_srcEng(srcEng),
                                 m_srcTime(srcTime),
                                 m_obsEng(obsEng),
                                 m_obsTime(obsTime),
                                 m_obs(obs),
                                 m_rot(rot),
                                 m_rho_roi(rho_roi),
                                 m_cos_rho_roi(std::cos(rho_roi)),
                                 m_sin_rho_roi(std::sin(rho_roi)),
                                 m_posangle_roi(posangle_roi),
                                 m_radius_roi(radius_roi),
                                 m_cos_radius_roi(std::cos(radius_roi)),
                                 m_iter(iter) { }
    double eval(const double& rho);
protected:
    const GCTAResponseIrf&         m_rsp;            //!< CTA response
    const GModelSpatialElliptical& m_model;          //!< Elliptical model
    const GCTAObservation&         m_obs;            //!< CTA observation
    const GMatrix&                 m_rot;            //!< Rotation matrix
    double                         m_semimajor;      //!< Ellipse boundary semimajor axis
    double                         m_semiminor;      //!< Ellipse boundary semiminor axis
    double                         m_posangle;       //!< Ellipse boundary position angle
    GEnergy                        m_srcEng;         //!< True photon energy
    GTime                          m_srcTime;        //!< True photon arrival time
    GEnergy                        m_obsEng;         //!< Observed photon energy
    GTime                          m_obsTime;        //!< Observed photon arrival time
    double                         m_rho_roi;        //!< Distance between model and ROI centre
    double                         m_cos_rho_roi;    //!< Cosine of m_rho_roi
    double                         m_sin_rho_roi;    //!< Sine of m_rho_roi
    double                         m_posangle_roi;   //!< Position angle of ROI
    double                         m_radius_roi;     //!< ROI+PSF radius
    double                         m_cos_radius_roi; //!< Cosine of m_radius_roi
    int                            m_iter;           //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_nroi_elliptical_kern_omega
 *
 * @brief Kernel for azimuth angle Nroi integration of elliptical model
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
class cta_nroi_elliptical_kern_omega : public GFunction {
public:
    cta_nroi_elliptical_kern_omega(const GCTAResponseIrf&         rsp,
                                   const GModelSpatialElliptical& model,
                                   const GEnergy&                 srcEng,
                                   const GTime&                   srcTime,
                                   const GEnergy&                 obsEng,
                                   const GTime&                   obsTime,
                                   const GCTAObservation&         obs,
                                   const GMatrix&                 rot,
                                   const double&                  rho,
                                   const double&                  sin_rho,
                                   const double&                  cos_rho,
                                   const double&                  posangle_roi) :
                                   m_rsp(rsp),
                                   m_model(model),
                                   m_srcEng(srcEng),
                                   m_srcTime(srcTime),
                                   m_obsEng(obsEng),
                                   m_obsTime(obsTime),
                                   m_obs(obs),
                                   m_rot(rot),
                                   m_rho(rho),
                                   m_sin_rho(sin_rho),
                                   m_cos_rho(cos_rho),
                                   m_posangle_roi(posangle_roi) { }
    double eval(const double& omega);
protected:
    const GCTAResponseIrf&         m_rsp;          //!< CTA response
    const GModelSpatialElliptical& m_model;        //!< Model
    const GCTAObservation&         m_obs;          //!< Pointer to observation
    const GMatrix&                 m_rot;          //!< Rotation matrix
    GEnergy                        m_srcEng;       //!< True photon energy
    GTime                          m_srcTime;      //!< True photon arrival time
    GEnergy                        m_obsEng;       //!< Observed photon energy
    GTime                          m_obsTime;      //!< Observed photon arrival time
    double                         m_rho;
    double                         m_sin_rho;      //!< Sine of offset angle
    double                         m_cos_rho;      //!< Cosine of offset angle
    double                         m_posangle_roi; //!< Position angle of ROI
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
    cta_irf_diffuse_kern_theta(const GCTAResponseIrf& rsp,
                               const GModelSpatial&   model,
                               const double&          theta,
                               const double&          phi,
                               const double&          zenith,
                               const double&          azimuth,
                               const GEnergy&         srcEng,
                               const GTime&           srcTime,
                               const double&          srcLogEng,
                               const GEnergy&         obsEng,
                               const GMatrix&         rot,
                               const double&          eta,
                               const int&             iter) :
                               m_rsp(rsp),
                               m_model(model),
                               m_theta(theta),
                               m_phi(phi),
                               m_zenith(zenith),
                               m_azimuth(azimuth),
                               m_srcEng(srcEng),
                               m_srcTime(srcTime),
                               m_srcLogEng(srcLogEng),
                               m_obsEng(obsEng),
                               m_rot(rot),
                               m_sin_eta(std::sin(eta)),
                               m_cos_eta(std::cos(eta)),
                               m_iter(iter) { }
    double eval(const double& theta);
protected:
    const GCTAResponseIrf& m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GMatrix&         m_rot;        //!< Rotation matrix
    double                 m_theta;      //!< Photon offset angle
    double                 m_phi;        //!< Photon azimuth angle
    double                 m_zenith;     //!< Pointing zenith angle
    double                 m_azimuth;    //!< Pointing azimuth angle
    GEnergy                m_srcEng;     //!< True photon energy
    GTime                  m_srcTime;    //!< True photon arrival time
    double                 m_srcLogEng;  //!< True photon log energy
    GEnergy                m_obsEng;     //!< Measured event energy
    double                 m_sin_eta;    //!< Sine of angular distance between
                                         //   observed photon direction and
                                         //   camera centre
    double                 m_cos_eta;    //!< Cosine of angular distance between
                                         //   observed photon direction and
                                         //   camera centre
    int                    m_iter;       // Integration iterations
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
    cta_irf_diffuse_kern_phi(const GCTAResponseIrf& rsp,
                             const GModelSpatial&   model,
                             const double&          zenith,
                             const double&          azimuth,
                             const GEnergy&         srcEng,
                             const GTime&           srcTime,
                             const double&          srcLogEng,
                             const GEnergy&         obsEng,
                             const GMatrix&         rot,
                             const double&          sin_theta,
                             const double&          cos_theta,
                             const double&          sin_ph,
                             const double&          cos_ph) :
                             m_rsp(rsp),
                             m_model(model),
                             m_zenith(zenith),
                             m_azimuth(azimuth),
                             m_srcEng(srcEng),
                             m_srcTime(srcTime),
                             m_srcLogEng(srcLogEng),
                             m_obsEng(obsEng),
                             m_rot(rot),
                             m_sin_theta(sin_theta),
                             m_cos_theta(cos_theta),
                             m_sin_ph(sin_ph),
                             m_cos_ph(cos_ph) { }
    double eval(const double& phi);
protected:
    const GCTAResponseIrf& m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GMatrix&         m_rot;        //!< Rotation matrix
    double                 m_zenith;     //!< Zenith angle
    double                 m_azimuth;    //!< Azimuth angle
    GEnergy                m_srcEng;     //!< True photon energy
    GTime                  m_srcTime;    //!< True photon arrival time
    double                 m_srcLogEng;  //!< True photon log energy
    GEnergy                m_obsEng;     //!< Measured event energy
    double                 m_sin_theta;  //!< Sine of offset angle
    double                 m_cos_theta;  //!< Cosine of offset angle
    double                 m_sin_ph;     //!< Sine term in angular distance equation
    double                 m_cos_ph;     //!< Cosine term in angular distance equation    
};


/***********************************************************************//**
 * @class cta_nroi_diffuse_kern_theta
 *
 * @brief Kernel for Nroi offest angle integration of diffuse model
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
class cta_nroi_diffuse_kern_theta : public GFunction {
public:
    cta_nroi_diffuse_kern_theta(const GCTAResponseIrf& rsp,
                                const GModelSpatial&   model,
                                const GEnergy&         srcEng,
                                const GTime&           srcTime,
                                const GEnergy&         obsEng,
                                const GTime&           obsTime,
                                const GCTAObservation& obs,
                                const GMatrix&         rot,
                                const int&             iter) :
                                m_rsp(rsp),
                                m_model(model),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime),
                                m_obsEng(obsEng),
                                m_obsTime(obsTime),
                                m_obs(obs),
                                m_rot(rot),
                                m_iter(iter) { }
    double eval(const double& theta);
protected:
    const GCTAResponseIrf& m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GCTAObservation& m_obs;        //!< CTA observation
    const GMatrix&         m_rot;        //!< Rotation matrix
    GEnergy                m_srcEng;     //!< True photon energy
    GTime                  m_srcTime;    //!< True photon arrival time
    GEnergy                m_obsEng;     //!< Observed photon energy
    GTime                  m_obsTime;    //!< Observed photon arrival time
    int                    m_iter;       //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_nroi_diffuse_kern_phi
 *
 * @brief Kernel for Nroi azimuth angle integration of diffuse model
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
class cta_nroi_diffuse_kern_phi : public GFunction {
public:
    cta_nroi_diffuse_kern_phi(const GCTAResponseIrf& rsp,
                              const GModelSpatial&   model,
                              const GEnergy&         srcEng,
                              const GTime&           srcTime,
                              const GEnergy&         obsEng,
                              const GTime&           obsTime,
                              const GCTAObservation& obs,
                              const GMatrix&         rot,
                              const double&          theta,
                              const double&          sin_theta) :
                              m_rsp(rsp),
                              m_model(model),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_obsEng(obsEng),
                              m_obsTime(obsTime),
                              m_obs(obs),
                              m_rot(rot),
                              m_theta(theta),
                              m_cos_theta(std::cos(theta)),
                              m_sin_theta(sin_theta) { }
    double eval(const double& phi);
protected:
    const GCTAResponseIrf& m_rsp;        //!< CTA response
    const GModelSpatial&   m_model;      //!< Spatial model
    const GCTAObservation& m_obs;        //!< CTA observation
    const GMatrix&         m_rot;        //!< Rotation matrix
    GEnergy                m_srcEng;     //!< True photon energy
    GTime                  m_srcTime;    //!< True photon arrival time
    GEnergy                m_obsEng;     //!< Observed photon energy
    GTime                  m_obsTime;    //!< Observed photon arrival time
    double                 m_theta;      //!< Offset angle (radians)
    double                 m_cos_theta;  //!< Cosine of offset angle
    double                 m_sin_theta;  //!< Sine of offset angle
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_rho
 *
 * @brief Kernel for radial model zenith angle integration
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
 *    K(\rho | E, t) = \sin \rho \times
 *                     S_{\rm p}(\rho | E, t) \times
 *                     \int_{\omega} PSF(\rho, \omega) d\omega
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho | E, t)\f$ is the radial model,
 * \f$PSF(\rho, \omega)\f$ is the point spread function,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 ***************************************************************************/
class cta_psf_radial_kern_rho : public GFunction {
public:
    cta_psf_radial_kern_rho(const GCTAResponseCube*    rsp,
                            const GModelSpatialRadial* model,
                            const GSkyDir&             srcDir,
                            const GEnergy&             srcEng,
                            const GTime&               srcTime,
                            const double&              rho_obs,
                            const double&              delta_max,
                            const int&                 iter) :
                            m_rsp(rsp),
                            m_model(model),
                            m_srcDir(srcDir),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime),
                            m_rho_obs(rho_obs),
                            m_cos_rho_obs(std::cos(rho_obs)),
                            m_sin_rho_obs(std::sin(rho_obs)),
                            m_delta_max(delta_max),
                            m_cos_delta_max(std::cos(delta_max)),
                            m_iter(iter) { }
    double eval(const double& rho);
public:
    const GCTAResponseCube*    m_rsp;           //!< CTA response
    const GModelSpatialRadial* m_model;         //!< Radial model
    GSkyDir                    m_srcDir;        //!< True photon arrival direction
    GEnergy                    m_srcEng;        //!< True photon energy
    GTime                      m_srcTime;       //!< True photon time
    double                     m_rho_obs;       //!< Distance of model centre from measured photon
    double                     m_cos_rho_obs;   //!< Cosine of m_rho_obs
    double                     m_sin_rho_obs;   //!< Sine of m_rho_obs
    double                     m_delta_max;     //!< Maximum PSF radius
    double                     m_cos_delta_max; //!< Cosine of maximum PSF radius
    int                        m_iter;          //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_irf_radial_kern_omega
 *
 * @brief Kernel for radial model azimuth angle integration
 *
 * This class implements the computation of
 *
 * \f[
 *    K(\omega | \rho, E, t) = PSF(\omega | \rho)
 * \f]
 *
 * where
 * \f$PSF(\omega | \rho)\f$ is the point spread function,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 ***************************************************************************/
class cta_psf_radial_kern_omega : public GFunction {
public:
    cta_psf_radial_kern_omega(const GCTAResponseCube*    rsp,
                              const GModelSpatialRadial* model,
                              const GSkyDir&             srcDir,
                              const GEnergy&             srcEng,
                              const GTime&               srcTime,
                              const double&              cos_psf,
                              const double&              sin_psf) :
                              m_rsp(rsp),
                              m_model(model),
                              m_srcDir(srcDir),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_cos_psf(cos_psf),
                              m_sin_psf(sin_psf) { }
    double eval(const double& omega);
public:
    const GCTAResponseCube*    m_rsp;     //!< CTA response
    const GModelSpatialRadial* m_model;   //!< Radial model
    GSkyDir                    m_srcDir;  //!< True photon sky direction
    GEnergy                    m_srcEng;  //!< True photon energy
    GTime                      m_srcTime; //!< True photon time
    double                     m_cos_psf; //!< Cosine term for PSF offset angle computation
    double                     m_sin_psf; //!< Sine term for PSF offset angle computation
};


/***********************************************************************//**
 * @class cta_psf_radial_kern_delta
 *
 * @brief Kernel for Psf delta angle integration used for stacked analysis
 ***************************************************************************/
class cta_psf_radial_kern_delta : public GFunction {
public:
    cta_psf_radial_kern_delta(const GCTAResponseCube*    rsp,
                              const GModelSpatialRadial* model,
                              const GSkyDir&             srcDir,
                              const GEnergy&             srcEng,
                              const GTime&               srcTime,
                              const double&              delta_mod,
                              const double&              theta_max,
                              const int&                 iter) :
                              m_rsp(rsp),
                              m_model(model),
                              m_srcDir(srcDir),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_delta_mod(delta_mod),
                              m_cos_delta_mod(std::cos(delta_mod)),
                              m_sin_delta_mod(std::sin(delta_mod)),
                              m_theta_max(theta_max),
                              m_cos_theta_max(std::cos(theta_max)),
                              m_iter(iter) { }
    double eval(const double& delta);
protected:
    const GCTAResponseCube*    m_rsp;           //!< Response cube
    const GModelSpatialRadial* m_model;         //!< Radial model
    GSkyDir                    m_srcDir;        //!< True photon arrival direction
    GEnergy                    m_srcEng;        //!< True photon energy
    GTime                      m_srcTime;       //!< True photon arrival time
    double                     m_delta_mod;     //!< Distance of model from Psf
    double                     m_cos_delta_mod; //!< Cosine of m_delta_mod
    double                     m_sin_delta_mod; //!< Sine of m_delta_mod
    double                     m_theta_max;     //!< Maximum model radius
    double                     m_cos_theta_max; //!< Cosine of m_theta_max
    int                        m_iter;          //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_psf_radial_kern_phi
 *
 * @brief Kernel for Psf phi angle integration used for stacked analysis
 ***************************************************************************/
class cta_psf_radial_kern_phi : public GFunction {
public:
    cta_psf_radial_kern_phi(const GModelSpatialRadial* model,
                            const GEnergy&             srcEng,
                            const GTime&               srcTime,
                            const double&              sin_fact,
                            const double&              cos_fact) :
                            m_model(model),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime),
                            m_sin_fact(sin_fact),
                            m_cos_fact(cos_fact) { }
    double eval(const double& phi);
protected:
    const GModelSpatialRadial* m_model;     //!< Radial model
    GEnergy                    m_srcEng;    //!< True photon energy
    GTime                      m_srcTime;   //!< True photon arrival time
    double                     m_sin_fact;  //!< sin(delta)*sin(delta_mod)
    double                     m_cos_fact;  //!< cos(delta)*cos(delta_mod)
};


/***********************************************************************//**
 * @class cta_irf_elliptical_kern_rho
 *
 * @brief Kernel for elliptical model zenith angle integration
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
 *                     \int_{\omega} 
 *                     S_{\rm p}(\rho, \omega | E, t) \, PSF(\rho, \omega)
 *                     d\omega
 * \f]
 *
 * where
 * \f$S_{\rm p}(\rho, \omega | E, t)\f$ is the elliptical model,
 * \f$PSF(\rho, \omega)\f$ is the point spread function,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 ***************************************************************************/
class cta_psf_elliptical_kern_rho : public GFunction {
public:
    cta_psf_elliptical_kern_rho(const GCTAResponseCube*        rsp,
                                const GModelSpatialElliptical* model,
                                const double&                  semimajor,
                                const double&                  semiminor,
                                const double&                  posangle,
                                const GSkyDir&                 srcDir,
                                const GEnergy&                 srcEng,
                                const GTime&                   srcTime,
                                const double&                  rho_obs,
                                const double&                  posangle_obs,
                                const double&                  delta_max,
                                const int&                     iter) :
                                m_rsp(rsp),
                                m_model(model),
                                m_semimajor(semimajor),
                                m_semiminor(semiminor),
                                m_posangle(posangle),
                                m_srcDir(srcDir),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime),
                                m_rho_obs(rho_obs),
                                m_cos_rho_obs(std::cos(rho_obs)),
                                m_sin_rho_obs(std::sin(rho_obs)),
                                m_posangle_obs(posangle_obs),
                                m_delta_max(delta_max),
                                m_cos_delta_max(std::cos(delta_max)),
                                m_iter(iter) { }
    double eval(const double& rho);
public:
    const GCTAResponseCube*        m_rsp;           //!< CTA response
    const GModelSpatialElliptical* m_model;         //!< Elliptical model
    double                         m_semimajor;     //!< Ellipse boundary semimajor axis
    double                         m_semiminor;     //!< Ellipse boundary semiminor axis
    double                         m_posangle;      //!< Ellipse boundary position angle
    GSkyDir                        m_srcDir;        //!< True photon arrival direction
    GEnergy                        m_srcEng;        //!< True photon energy
    GTime                          m_srcTime;       //!< True photon time
    double                         m_rho_obs;       //!< Distance of model centre from measured photon
    double                         m_cos_rho_obs;   //!< Cosine of m_rho_obs
    double                         m_sin_rho_obs;   //!< Sine of m_rho_obs
    double                         m_posangle_obs;  //!< Photon position angle measured from model centre
    double                         m_delta_max;     //!< Maximum PSF radius
    double                         m_cos_delta_max; //!< Cosine of maximum PSF radius
    int                            m_iter;          //!< Integration iterations
};


/***********************************************************************//**
 * @class cta_irf_elliptical_kern_omega
 *
 * @brief Kernel for elliptical model azimuth angle integration
 *
 * This class implements the computation of
 *
 * \f[
 *    S_{\rm p}(\omega | \rho, E, t) \, PSF(\omega | \rho)
 * \f]
 *
 * where
 * \f$S_{\rm p}(\omega | \rho, E, t)\f$ is the elliptical model,
 * \f$ PSF(\omega | \rho)\f$ is the point spread function,
 * \f$\rho\f$ is the distance from the model centre, and
 * \f$\omega\f$ is the position angle with respect to the connecting line
 * between the model centre and the observed photon arrival direction.
 ***************************************************************************/
class cta_psf_elliptical_kern_omega : public GFunction {
public:
    cta_psf_elliptical_kern_omega(const GCTAResponseCube*        rsp,
                                  const GModelSpatialElliptical* model,
                                  const GSkyDir&                 srcDir,
                                  const GEnergy&                 srcEng,
                                  const GTime&                   srcTime,
                                  const double&                  posangle_obs,
                                  const double&                  rho,
                                  const double&                  cos_psf,
                                  const double&                  sin_psf) :
                                  m_rsp(rsp),
                                  m_model(model),
                                  m_srcDir(srcDir),
                                  m_srcEng(srcEng),
                                  m_srcTime(srcTime),
                                  m_posangle_obs(posangle_obs),
                                  m_rho(rho),
                                  m_cos_psf(cos_psf),
                                  m_sin_psf(sin_psf) { }
    double eval(const double& omega);
public:
    const GCTAResponseCube*        m_rsp;          //!< CTA response
    const GModelSpatialElliptical* m_model;        //!< Spatial model
    GSkyDir                        m_srcDir;       //!< True photon sky direction
    GEnergy                        m_srcEng;       //!< True photon energy
    GTime                          m_srcTime;      //!< True photon time
    double                         m_posangle_obs; //!< Measured photon position angle from model centre
    double                         m_rho;          //!< Model zenith angle
    double                         m_cos_psf;      //!< Cosine term for PSF offset angle computation
    double                         m_sin_psf;      //!< Sine term for PSF offset angle computation
};


/***********************************************************************//**
 * @class cta_psf_diffuse_kern_delta
 *
 * @brief Kernel for Psf delta angle integration used for stacked analysis
 ***************************************************************************/
class cta_psf_diffuse_kern_delta : public GFunction {
public:
    cta_psf_diffuse_kern_delta(const GCTAResponseCube* rsp,
                               const GModelSpatial*    model,
                               const GSkyDir&          srcDir,
                               const GEnergy&          srcEng,
                               const GTime&            srcTime,
                               const GMatrix&          rot,
                               const int&              iter) :
                               m_rsp(rsp),
                               m_model(model),
                               m_srcDir(srcDir),
                               m_srcEng(srcEng),
                               m_srcTime(srcTime),
                               m_rot(rot),
                               m_iter (iter),
                               m_psf_max(rsp->psf()(srcDir, 0.0, srcEng)) { }
    double eval(const double& delta);
protected:
    const GCTAResponseCube* m_rsp;     //!< Response cube
    const GModelSpatial*    m_model;   //!< Spatial model
    const GMatrix&          m_rot;     //!< Rotation matrix
    GSkyDir                 m_srcDir;  //!< True photon arrival direction
    GEnergy                 m_srcEng;  //!< True photon energy
    GTime                   m_srcTime; //!< True photon arrival time
    int                     m_iter;    //!< Romberg iterations
    double                  m_psf_max; //!< Maximum PSF value
};


/***********************************************************************//**
 * @class cta_psf_diffuse_kern_phi
 *
 * @brief Kernel for Psf phi angle integration used for stacked analysis
 ***************************************************************************/
class cta_psf_diffuse_kern_phi : public GFunction {
public:
    cta_psf_diffuse_kern_phi(const GModelSpatial* model,
                             const GEnergy&       srcEng,
                             const GTime&         srcTime,
                             const GMatrix&       rot,
                             const double&        sin_delta,
                             const double&        cos_delta) :
                             m_model(model),
                             m_srcEng(srcEng),
                             m_srcTime(srcTime),
                             m_rot(rot),
                             m_sin_delta(sin_delta),
                             m_cos_delta(cos_delta) { }
    double eval(const double& phi);
protected:
    const GModelSpatial* m_model;     //!< Spatial model
    const GMatrix&       m_rot;       //!< Rotation matrix
    GEnergy              m_srcEng;    //!< True photon energy
    GTime                m_srcTime;   //!< True photon arrival time
    double               m_sin_delta; //!< sin(delta)
    double               m_cos_delta; //!< cos(delta)
};

#endif /* GCTARESPONSE_HELPERS_HPP */
