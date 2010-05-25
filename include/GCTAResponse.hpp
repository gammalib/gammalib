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
//#include "GNodeArray.hpp"
//#include "GVector.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"
//#include "GLATResponseTable.hpp"
//#include "GFits.hpp"
//#include "GFitsHDU.hpp"


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA instrument response function classes.
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    ~GCTAResponse(void);

    // Operators
    GCTAResponse& operator= (const GCTAResponse & rsp);

    // Implemented virtual base class methods
    double irf(const GSkyDir& obsDir, const GEnergy& obsEng,
               const GSkyDir& srcDir, const GEnergy& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const GTime& time);
    double aeff(const GSkyDir& obsDir, const GEnergy& obsEng,
                const GSkyDir& srcDir, const GEnergy& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const GTime& time);
    double  psf(const GSkyDir& obsDir, const GEnergy& obsEng,
                const GSkyDir& srcDir, const GEnergy& srcEng,
                const GSkyDir& instPntDir, const double& instPosAng,
                const GTime& time);
    double edisp(const GSkyDir& obsDir, const GEnergy& obsEng,
                 const GSkyDir& srcDir, const GEnergy& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const GTime& time);

    // Other Methods
    void          set_caldb(const std::string& caldb);
    GCTAResponse* clone(void) const;

private:
    // Private methods
    void    init_members(void);
    void    copy_members(const GCTAResponse& rsp);
    void    free_members(void);
};

#endif /* GCTARESPONSE_HPP */
