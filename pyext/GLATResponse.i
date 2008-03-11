/***************************************************************************
 *      GLATResponse.hpp  -  GLAST LAT Response class SWIG definition      *
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

    // Methods
    double irf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);
    double psf(const GSkyDir& obsDir, const double& obsEng,
               const GSkyDir& srcDir, const double& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const double& time);
    double aeff(const GSkyDir& obsDir, const double& obsEng,
                const GSkyDir& srcDir, const double& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const double& time);
    double edisp(const GSkyDir& obsDir, const double& obsEng,
                 const GSkyDir& srcDir, const double& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const double& time);
    void   set_caldb(const std::string& caldb);
    void   load(const std::string& rspname, const std::string& rsptype);
    void   save(const std::string& rspname) const;
};
