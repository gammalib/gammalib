/***************************************************************************
 *       GLATResponse.i  -  GLAST LAT Response class SWIG definition       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATResponse.i
 * @brief GLATResponse class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATResponse.hpp"
%}
%include stl.i


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief SWIG interface for the GLAST LAT instrument response function.
 ***************************************************************************/
class GLATResponse : public GResponse {
public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse(void);

    // Implemented pure virtual methods
    double live(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt) const;
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
    double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt, const GRoi& roi) const;
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                  const GTime& srcTime, const GPointing& pnt, const GEbounds& ebds) const;
    double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                  const GTime& srcTime, const GPointing& pnt, const GGti& gti) const;
    void   load(const std::string& rspname);
    bool   hasedisp(void) const { return false; }
    bool   hastdisp(void) const { return false; }
    GLATResponse* clone(void) const;

    // Other Methods
    double        aeff(const double& logE, const double& ctheta);
    void          aeff_ctheta_min(const double& ctheta);
    double        aeff_ctheta_min(void) const;
    double        psf(const double& delta, const double& logE, const double& ctheta);
    GVector       psf(const GVector& delta, const double& logE, const double& ctheta);
    void          save(const std::string& rspname) const;
};
