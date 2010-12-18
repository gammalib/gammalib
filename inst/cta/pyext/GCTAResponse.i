/***************************************************************************
 *          GCTAResponse.i  -  CTA Response class SWIG definition          *
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
 * @file GCTAResponse.i
 * @brief GCTAResponse class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponse.hpp"
%}
//%include stl.i


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief SWIG interface for the GLAST LAT instrument response function.
 ***************************************************************************/
class GCTAResponse : public GResponse {
public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Implement pure virtual base class methods
    void          clear(void);
    GCTAResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    double        irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                      const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                      const GObservation& obs) const;
    double        irf(const GEvent& event, const GModel& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GObservation& obs) const;
    double        nirf(const GSkyDir& srcDir, const GEnergy& srcEng, const GTime& srcTime,
                       const GObservation& obs) const;

    // Other Methods
    double irf(const GCTAEventAtom& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const;
    double irf(const GCTAEventBin& event, const GModel& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const;
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
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponse {
    GCTAResponse copy() {
        return (*self);
    }
};
