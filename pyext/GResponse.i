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

    // Virtual methods
    virtual double      irf(const GInstDir& obsDir, const GEnergy& obsEng,
                            const GTime& obsTime,
                            const GSkyDir&  srcDir, const GEnergy& srcEng,
                            const GTime& srcTime,
                            const GPointing& pnt);
    virtual double      aeff(const GInstDir& obsDir, const GEnergy& obsEng,
                             const GTime& obsTime,
                             const GSkyDir&  srcDir, const GEnergy& srcEng,
                             const GTime& srcTime,
                             const GPointing& pnt) = 0;
    virtual double      psf(const GInstDir& obsDir, const GEnergy& obsEng,
                            const GTime& obsTime,
                            const GSkyDir&  srcDir, const GEnergy& srcEng,
                            const GTime& srcTime,
                            const GPointing& pnt) = 0;
    virtual double      edisp(const GInstDir& obsDir, const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GSkyDir&  srcDir, const GEnergy& srcEng,
                              const GTime& srcTime,
                              const GPointing& pnt) = 0;
    virtual double      tdisp(const GInstDir& obsDir, const GEnergy& obsEng,
                              const GTime& obsTime,
                              const GSkyDir&  srcDir, const GEnergy& srcEng,
                              const GTime& srcTime,
                              const GPointing& pnt) = 0;
    virtual void        caldb(const std::string& caldb);
    virtual std::string caldb(void) const { return m_caldb; }
};
