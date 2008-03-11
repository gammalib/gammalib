/***************************************************************************
 *       GResponse.i  -  Response abstract base class SWIG interface       *
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
 * @file GResponse.i
 * @brief GResponse class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponse.hpp"
%}


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract SWIG interface for the instrument response function classes.
 ***************************************************************************/
class GResponse {
public:
    // Constructors and destructors
    GResponse();
    GResponse(const GResponse& rsp);
    virtual ~GResponse();

    // Methods
    virtual double irf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;
    virtual double psf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;
    virtual double aeff(const GSkyDir& obsDir, const double& obsEng,
                        const GSkyDir& srcDir, const double& srcEng,
                        const GSkyDir& instPntDir, const double& instPosAng,
                        const double& time) = 0;
    virtual double edisp(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time) = 0;
    virtual void set_caldb(const std::string& caldb) = 0;
};
