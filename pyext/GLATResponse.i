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
 * @brief SWIG interface for the GLAST LAT instrument response function classes.
 ***************************************************************************/
class GLATResponse : public GResponse {
public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse(void);

    // Implemented virtual base class methods
    double irf(GSkyDir& obsDir, const GEnergy& obsEng,
               GSkyDir& srcDir, const GEnergy& srcEng,
               const GPointing* pnt, const GTime& time);
    double aeff(GSkyDir& obsDir, const GEnergy& obsEng,
                GSkyDir& srcDir, const GEnergy& srcEng,
                const GPointing* pnt, const GTime& time);
    double  psf(GSkyDir& obsDir, const GEnergy& obsEng,
                GSkyDir& srcDir, const GEnergy& srcEng,
                const GPointing* pnt, const GTime& time);
    double edisp(GSkyDir& obsDir, const GEnergy& obsEng,
                 GSkyDir& srcDir, const GEnergy& srcEng,
                 const GPointing* pnt, const GTime& time);
    void   set_caldb(const std::string& caldb);

    // Other Methods
    double        aeff(const double& logE, const double& ctheta);
    void          aeff_ctheta_min(const double& ctheta);
    double        aeff_ctheta_min(void) const;
    double        psf(const double& delta, const double& logE, const double& ctheta);
    GVector       psf(const GVector& delta, const double& logE, const double& ctheta);
    void          load(const std::string& rspname, const std::string& rsptype);
    void          save(const std::string& rspname) const;
};
