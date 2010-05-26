/***************************************************************************
 *       GResponse.i  -  Response abstract base class SWIG interface       *
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
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Pure virtual methods
    virtual double irf(GSkyDir& obsDir, const GEnergy& obsEng,
                       GSkyDir& srcDir, const GEnergy& srcEng,
                       const GPointing* pnt, const GTime& time) = 0;
    virtual double psf(GSkyDir& obsDir, const GEnergy& obsEng,
                       GSkyDir& srcDir, const GEnergy& srcEng,
                       const GPointing* pnt, const GTime& time) = 0;
    virtual double aeff(GSkyDir& obsDir, const GEnergy& obsEng,
                        GSkyDir& srcDir, const GEnergy& srcEng,
                        const GPointing* pnt, const GTime& time) = 0;
    virtual double edisp(GSkyDir& obsDir, const GEnergy& obsEng,
                         GSkyDir& srcDir, const GEnergy& srcEng,
                         const GPointing* pnt, const GTime& time) = 0;
    virtual void   set_caldb(const std::string& caldb) = 0;
};
