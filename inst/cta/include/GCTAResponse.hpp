/***************************************************************************
 *                  GCTAResponse.hpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAResponse.hpp
 * @brief CTA instrument response function class interface definition
 * @author J. Knodlseder
 */

#ifndef GCTARESPONSE_HPP
#define GCTARESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <vector>
#include "GEvent.hpp"
#include "GModelSky.hpp"
#include "GModelPointSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelRadial.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPhoton.hpp"
#include "GNodeArray.hpp"
#include "GIntegrand.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GCTARoi.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTADir.hpp"

/* __ Forward declaration ________________________________________________ */
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA instrument response function.
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    explicit GCTAResponse(const std::string& rspname, const std::string& caldb);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Operators
    virtual GCTAResponse& operator= (const GCTAResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCTAResponse* clone(void) const;
    virtual bool          hasedisp(void) const { return false; }
    virtual bool          hastdisp(void) const { return false; }
    virtual double        irf(const GInstDir&     obsDir,
                              const GEnergy&      obsEng,
                              const GTime&        obsTime,
                              const GSkyDir&      srcDir,
                              const GEnergy&      srcEng,
                              const GTime&        srcTime,
                              const GObservation& obs) const;
    virtual double        npred(const GSkyDir&      srcDir,
                                const GEnergy&      srcEng,
                                const GTime&        srcTime,
                                const GObservation& obs) const;
    virtual std::string   print(void) const;

    // Overload virtual base class methods
    virtual double irf_extended(const GInstDir&             obsDir,
                                const GEnergy&              obsEng,
                                const GTime&                obsTime,
                                const GModelExtendedSource& model,
                                const GEnergy&              srcEng,
                                const GTime&                srcTime,
                                const GObservation&         obs) const;
    virtual double irf_diffuse(const GInstDir&            obsDir,
                               const GEnergy&             obsEng,
                               const GTime&               obsTime,
                               const GModelDiffuseSource& model,
                               const GEnergy&             srcEng,
                               const GTime&               srcTime,
                               const GObservation&        obs) const;

    // Other Methods
    GCTAEventAtom* mc(const double& area, const GPhoton& photon,
                      const GPointing& pnt, GRan& ran) const;
    void           caldb(const std::string& caldb);
    std::string    caldb(void) const { return m_caldb; }
    void           load(const std::string& rspname);
    void           eps(const double& eps) { m_eps=eps; }
    const double&  eps(void) const { return m_eps; }

    // Low-level response methods
    double aeff(const double& theta,
                const double& phi,
                const double& zenith,
                const double& azimuth,
                const double& srcLogEng) const;
    double psf(const double& delta,
               const double& theta,
               const double& phi,
               const double& zenith,
               const double& azimuth,
               const double& srcLogEng) const;
    double psf_delta_max(const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const;
    double edisp(const double& obsLogEng,
                 const double& theta,
                 const double& phi,
                 const double& zenith,
                 const double& azimuth,
                 const double& srcLogEng) const;
    double npsf(const GSkyDir&      srcDir,
                const double&       srcLogEng,
                const GTime&        srcTime,
                const GCTAPointing& pnt,
                const GCTARoi&      roi) const;
    double nedisp(const GSkyDir&      srcDir,
                  const GEnergy&      srcEng,
                  const GTime&        srcTime,
                  const GCTAPointing& pnt,
                  const GEbounds&     ebds) const;

    // Analytical PSF implementation
    double psf_dummy(const double& delta, const double& sigma) const;
    double psf_dummy_sigma(const double& srcLogEng) const;
    double psf_dummy_max(const double& sigma) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCTAResponse& rsp);
    void free_members(void);
    void read_performance_table(const std::string& filename);

    // IRF theta integration kernel
    class irf_kern_theta : public GIntegrand {
    public:
        irf_kern_theta(const GCTAResponse* rsp,
                       const GCTAPointing* pnt,
                       const GModelRadial* radial,
                       const GCTADir*      obsCam,
                       const GCTADir*      srcCam,
                       double              delta_max,
                       double              src_max,
                       double              zenith,
                       double              azimuth,
                       double              srcLogEng,
                       double              obsLogEng,
                       double              sigma) :
                       m_rsp(rsp),
                       m_pnt(pnt),
                       m_radial(radial),
                       m_obs_cam(obsCam),
                       m_src_cam(srcCam),
                       m_delta_max(delta_max),
                       m_cos_delta_max(std::cos(delta_max)),
                       m_src_max(src_max),
                       m_cos_src_max(std::cos(src_max)),
                       m_zenith(zenith),
                       m_azimuth(azimuth),
                       m_srcLogEng(srcLogEng),
                       m_obsLogEng(obsLogEng),
                       m_sigma(sigma) { return; }
        double eval(double theta);
    protected:
        const GCTAResponse* m_rsp;           //!< Pointer to CTA response
        const GCTAPointing* m_pnt;           //!< Pointer to CTA pointing
        const GModelRadial* m_radial;        //!< Pointer to radial spatial model
        const GCTADir*      m_obs_cam;       //!< Pointer to photon direction in camera
        const GCTADir*      m_src_cam;       //!< Pointer to source location in camera
        double              m_delta_max;     //!< Maximum PSF angle delta_max
        double              m_cos_delta_max; //!< Cosine of delta_max
        double              m_src_max;       //!< Maximum source radius
        double              m_cos_src_max;   //!< Cosine of source radius
        double              m_zenith;        //!< Telescope zenith
        double              m_azimuth;       //!< Telescope azimuth
        double              m_srcLogEng;     //!< True photon energy
        double              m_obsLogEng;     //!< Measured photon energy
        double              m_sigma;         //!< Width of PSF in radians
    };

    // IRF phi integration kernel
    class irf_kern_phi : public GIntegrand {
    public:
        irf_kern_phi(const GCTAResponse* rsp,
                     const GCTAPointing* pnt,
                     const GModelRadial* radial,
                     const GCTADir*      obsCam,
                     const GCTADir*      srcCam,
                     double              theta,
                     double              zenith,
                     double              azimuth,
                     double              srcLogEng,
                     double              obsLogEng,
                     double              sigma,
                     double              cos_obs,
                     double              sin_obs,
                     double              cos_src,
                     double              sin_src) :
                     m_rsp(rsp),
                     m_pnt(pnt),
                     m_radial(radial),
                     m_obs_cam(obsCam),
                     m_src_cam(srcCam),
                     m_theta(theta),
                     m_zenith(zenith),
                     m_azimuth(azimuth),
                     m_srcLogEng(srcLogEng),
                     m_obsLogEng(obsLogEng),
                     m_sigma(sigma),
                     m_cos_obs(cos_obs),
                     m_sin_obs(sin_obs),
                     m_cos_src(cos_src),
                     m_sin_src(sin_src) { return; }
        double eval(double phi);
    protected:
        const GCTAResponse* m_rsp;           //!< Pointer to CTA response
        const GCTAPointing* m_pnt;           //!< Pointer to CTA pointing
        const GModelRadial* m_radial;        //!< Pointer to radial spatial model
        const GCTADir*      m_obs_cam;       //!< Pointer to photon direction in camera
        const GCTADir*      m_src_cam;       //!< Pointer to source location in camera
        double              m_theta;         //!< Radial offset in camera
        double              m_zenith;        //!< Telescope zenith
        double              m_azimuth;       //!< Telescope azimuth
        double              m_srcLogEng;     //!< True photon energy
        double              m_obsLogEng;     //!< Measured photon energy
        double              m_sigma;         //!< Width of PSF in radians
        double              m_cos_obs;       //!< Cosine term for photon distance
        double              m_sin_obs;       //!< Sine term for photon distance
        double              m_cos_src;       //!< Cosine term for source location distance
        double              m_sin_src;       //!< Sine term for source location distance
    };

    // Integration
    class npsf_kern_rad_azsym : public GIntegrand {
    public:
        npsf_kern_rad_azsym(const GCTAResponse* parent, double roi,
                            double psf, double sigma) :
                            m_parent(parent), m_roi(roi), m_cosroi(cos(roi)),
                            m_psf(psf), m_cospsf(cos(psf)), m_sinpsf(sin(psf)),
                            m_sigma(sigma) { return; }
        double eval(double r);
    protected:
        const GCTAResponse* m_parent;  //!< Pointer to parent
        GEnergy             m_srcEng;  //!< Source energy
        double              m_roi;     //!< ROI radius in radians
        double              m_cosroi;  //!< Cosine of ROI radius
        double              m_psf;     //!< PSF-ROI centre distance in radians
        double              m_cospsf;  //!< Cosine of PSF-ROI centre distance
        double              m_sinpsf;  //!< Sine of PSF-ROI centre distance
        double              m_sigma;   //!< Width of PSF in radians
    };

    // Private data members
    std::string         m_caldb;    //!< Name of or path to the calibration database
    std::string         m_rspname;  //!< Name of the instrument response
    GNodeArray          m_nodes;    //!< log(E) nodes for interpolation
    std::vector<double> m_logE;     //!< log(E) = log10(E/TeV) - bin centre
    std::vector<double> m_aeff;     //!< Effective area in square metres after all cuts
    std::vector<double> m_r68;      //!< 68% containment radius of PSF post cuts in degrees
    std::vector<double> m_r80;      //!< 80% containment radius of PSF post cuts in degrees
    double              m_eps;      //!< Integration precision

};

#endif /* GCTARESPONSE_HPP */
