/***************************************************************************
 *       GLATResponse.i  -  GLAST LAT Response class SWIG definition       *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
    GLATResponse();
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse();

    // IRF Methods
    double irf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);
    // Aeff Methods
    double aeff(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);
    double aeff(const double& logE, const double& ctheta);
    void   aeff_ctheta_min(const double& ctheta);
    double aeff_ctheta_min(void) const;

    // PSF Methods
    double  psf(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);
    double  psf(const double& delta, const double& logE, const double& ctheta);
    GVector psf(const GVector& delta, const double& logE, const double& ctheta);

    // Edisp Methods
    double edisp(const GSkyDir& obsDir, const double& obsEng,
                 const GSkyDir& srcDir, const double& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const double& time);

    // Other Methods
    void set_caldb(const std::string& caldb);
    void load(const std::string& rspname, const std::string& rsptype);
    void save(const std::string& rspname) const;
};
