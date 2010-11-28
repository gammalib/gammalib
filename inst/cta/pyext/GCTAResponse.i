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

    // Reponse function computation methods
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

    // Implemented pure virtual base class methods
    GCTAResponse* clone(void) const;
    void          load(const std::string& irfname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }

    // Other Methods
    double psf(const double& theta, const double& sigma) const;
    double psf_sigma(const GEnergy& srcEng) const;
    double npsf(const double& psf, const double& radroi, const double& sigma) const;
};
