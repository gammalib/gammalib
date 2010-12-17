/***************************************************************************
 *                  GCTAResponse.hpp  -  CTA Response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
#include <vector>
#include "GCTAEventAtom.hpp"
#include "GCTAEventBin.hpp"
#include "GEvent.hpp"
#include "GModel.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GNodeArray.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA instrument response function.
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Operators
    GCTAResponse& operator= (const GCTAResponse & rsp);

    // Implement pure virtual base class methods
    void          clear(void);
    GCTAResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    double        irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                      const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                      const GPointing& pnt) const;
    double        irf(const GEvent& event, const GModel& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GPointing& pnt) const;
    double        nirf(const GSkyDir& srcDir, const GEnergy& srcEng, const GTime& srcTime,
                       const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                       const GGti& gti) const;
    std::string   print(void) const;

    // Other Methods
    double irf(const GCTAEventAtom& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double irf(const GCTAEventBin& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double live(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt) const;
    double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt) const;
    double psf(const GInstDir& obsDir,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double edisp(const GEnergy& obsEng,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const;
    double tdisp(const GTime& obsTime,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const;
    double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi) const;
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                  const GPointing& pnt, const GEbounds& ebds) const;
    double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                  const GPointing& pnt, const GGti& gti) const;

    // Other Methods
    double psf(const double& theta, const double& sigma) const;
    double psf_sigma(const GEnergy& srcEng) const;
    double npsf(const double& psf, const double& radroi, const double& sigma) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCTAResponse& rsp);
    void free_members(void);
    void read_performance_table(const std::string& filename);

    // Integration
    double npsf_kern_azsym(const double& rad,
                           const double& roi, const double& cosroi,
                           const double& psf, const double& cospsf,
                           const double& sinpsf) const;
    class npsf_kern_rad_azsym : public GIntegrand {
    public:
        npsf_kern_rad_azsym(const GCTAResponse* parent, double roi, double cosroi,
                            double psf, double cospsf, double sinpsf, double sigma) :
                            m_parent(parent), m_roi(roi), m_cosroi(cosroi),
                            m_psf(psf), m_cospsf(cospsf), m_sinpsf(sinpsf),
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
    GNodeArray          m_nodes; //!< log(E) nodes for interpolation
    std::vector<double> m_logE;  //!< log(E) = log10(E/TeV) - bin centre
    std::vector<double> m_aeff;  //!< Effective area in square metres after all cuts
    std::vector<double> m_r68;   //!< 68% containment radius of PSF post cuts in degrees
    std::vector<double> m_r80;   //!< 80% containment radius of PSF post cuts in degrees

};

#endif /* GCTARESPONSE_HPP */
