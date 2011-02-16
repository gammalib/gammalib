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
 * @brief GCTAResponse class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTARESPONSE_HPP
#define GCTARESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include <vector>
#include "GEvent.hpp"
#include "GModel.hpp"
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
#include "GCTAEventAtom.hpp"
#include "GCTAEventBin.hpp"
#include "GCTAInstDir.hpp"

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
    virtual double        irf(const GEvent& event, const GModelSky& model,
                              const GEnergy& srcEng, const GTime& srcTime,
                              const GObservation& obs) const;
    virtual double        npred(const GModelSky& model, const GEnergy& srcEng,
                                const GTime& srcTime, const GObservation& obs) const;
    virtual std::string   print(void) const;

    // Other Methods
    GCTAEventAtom* mc(const double& area, const GPhoton& photon,
                      const GPointing& pnt, GRan& ran) const;
    void           caldb(const std::string& caldb);
    std::string    caldb(void) const { return m_caldb; }
    void           load(const std::string& rspname);

    // Other response methods
    double irf_ptsrc(const GCTAInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                     const GSkyDir& srcDir, const GEnergy& srcEng, const GTime& srcTime,
                     const GCTAObservation& obs) const;
    double irf_diffuse(const GCTAInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                       const GModelSky& model, const GEnergy& srcEng, const GTime& srcTime,
                       const GCTAObservation& obs) const;
    double aeff(const double& theta, const double& phi,
                const double& zenith, const double& azimuth,
                const double& srcLogEng) const;
    double psf(const double& delta,
               const double& theta, const double& phi,
               const double& zenith, const double& azimuth,
               const double& srcLogEng) const;
    double psf_delta_max(const double& theta, const double& phi,
                         const double& zenith, const double& azimuth,
                         const double& srcLogEng) const;
    double edisp(const double& obsLogEng,
                 const double& theta, const double& phi,
                 const double& zenith, const double& azimuth,
                 const double& srcLogEng) const;
    double npred(const GSkyDir& srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GObservation& obs) const;
    double npsf(const GSkyDir&  srcDir, const double& srcLogEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi) const;
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                  const GPointing& pnt, const GEbounds& ebds) const;
    double psf_dummy(const double& delta, const double& sigma) const;
    double psf_dummy_sigma(const double& srcLogEng) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCTAResponse& rsp);
    void free_members(void);
    void read_performance_table(const std::string& filename);

    // PSF theta integration kernel
    class psf_kern_theta : public GIntegrand {
    public:
        psf_kern_theta(const GCTAResponse* rsp,
                       const GCTAPointing* pnt,
                       const GModelSpatial* spatial,
                       double dist, double pa, double delta_max,
                       double zenith, double azimuth,
                       double srcLogEng, double obsLogEng, double sigma) :
                       m_rsp(rsp), m_pnt(pnt), m_spatial(spatial),
                       m_dist(dist),
                       m_cos_dist(std::cos(dist)), m_sin_dist(std::sin(dist)),
                       m_pa(pa), m_delta_max(delta_max),
                       m_cos_delta_max(std::cos(delta_max)),
                       m_zenith(zenith), m_azimuth(azimuth),
                       m_srcLogEng(srcLogEng), m_obsLogEng(obsLogEng),
                       m_sigma(sigma) { return; }
        double eval(double theta);
    protected:
        const GCTAResponse*  m_rsp;           //!< Pointer to CTA response
        const GCTAPointing*  m_pnt;           //!< Pointer to CTA pointing
        const GModelSpatial* m_spatial;       //!< Pointer to spatial model
        double               m_dist;          //!< PSF centre offset angle
        double               m_cos_dist;      //!< Cosine of PSF centre radial offset angle
        double               m_sin_dist;      //!< Sine of PSF centre radial offset angle
        double               m_pa;            //!< Position angle
        double               m_delta_max;     //!< Maximum PSF angle delta_max
        double               m_cos_delta_max; //!< Cosine of delta_max
        double               m_zenith;        //!< Telescope zenith
        double               m_azimuth;       //!< Telescope azimuth
        double               m_srcLogEng;     //!< True photon energy
        double               m_obsLogEng;     //!< Measured photon energy
        double               m_sigma;         //!< Width of PSF in radians
    };

    // PSF phi integration kernel
    class psf_kern_phi : public GIntegrand {
    public:
        psf_kern_phi(const GCTAResponse* rsp,
                     const GCTAPointing* pnt,
                     const GModelSpatial* spatial,
                     double theta, double pa, double zenith, double azimuth,
                     double srcLogEng, double obsLogEng, double sigma,
                     double cos_term, double sin_term) :
                     m_rsp(rsp), m_pnt(pnt), m_spatial(spatial),
                     m_theta(theta), m_pa(pa),
                     m_zenith(zenith), m_azimuth(azimuth),
                     m_srcLogEng(srcLogEng), m_obsLogEng(obsLogEng),
                     m_sigma(sigma),
                     m_cos_term(cos_term), m_sin_term(sin_term) { return; }
        double eval(double phi);
    protected:
        const GCTAResponse*  m_rsp;           //!< Pointer to CTA response
        const GCTAPointing*  m_pnt;           //!< Pointer to CTA pointing
        const GModelSpatial* m_spatial;       //!< Pointer to spatial model
        double               m_theta;         //!< Radial offset in camera
        double               m_pa;            //!< Position angle
        double               m_zenith;        //!< Telescope zenith
        double               m_azimuth;       //!< Telescope azimuth
        double               m_srcLogEng;     //!< True photon energy
        double               m_obsLogEng;     //!< Measured photon energy
        double               m_sigma;         //!< Width of PSF in radians
        double               m_cos_term;      //!< Cosine term
        double               m_sin_term;      //!< Sine term
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

};

#endif /* GCTARESPONSE_HPP */
