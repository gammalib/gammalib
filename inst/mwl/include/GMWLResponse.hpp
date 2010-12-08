/***************************************************************************
 *          GMWLResponse.hpp  -  Multi-wavelength response class           *
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
 * @file GMWLResponse.hpp
 * @brief GMWLResponse class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLRESPONSE_HPP
#define GMWLRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GMWLResponse
 *
 * @brief Interface for the Multi-wavelength response class
 *
 * The Multi-wavelength response class is needed to evaluate the response
 * to a source model. Since the multi-wavelength instrument classes work
 * directly in photon space, the response is unity.
 ***************************************************************************/
class GMWLResponse : public GResponse {

public:
    // Constructors and destructors
    GMWLResponse(void);
    GMWLResponse(const GMWLResponse& rsp);
    virtual ~GMWLResponse(void);

    // Operators
    GMWLResponse& operator= (const GMWLResponse & rsp);

    // Reponse function computation methods
    double irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const { return 1.0; }
    double live(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt) const { return 1.0; }
    double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt) const { return 1.0; }
    double psf(const GInstDir& obsDir,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const { return 1.0; } 
    double edisp(const GEnergy& obsEng,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const { return 1.0; }
    double tdisp(const GTime& obsTime,
                 const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                 const GPointing& pnt) const { return 1.0; }
    double nirf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                const GGti& gti) const { return 1.0; }
    double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi) const { return 1.0; }
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                  const GPointing& pnt, const GEbounds& ebds) const { return 1.0; }
    double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                  const GPointing& pnt, const GGti& gti) const { return 1.0; }

    // Pure virtual base class methods
    void          clear(void);
    GMWLResponse* clone(void) const;
    void          load(const std::string& irfname) { return; }
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    std::string   print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLResponse& pnt);
    void free_members(void);
};

#endif /* GCTARESPONSE_HPP */
