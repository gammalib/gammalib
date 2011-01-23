/***************************************************************************
 *          GMWLResponse.hpp  -  Multi-wavelength response class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * directly in photon space, the response is by definition unity.
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
    double irf(const GEvent& event, const GModelSky& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const { return 1.0; }
    double npred(const GModelSky& model, const GEnergy& srcEng,
                 const GTime& srcTime,
                 const GObservation& obs) const { return 1.0; }

    // Pure virtual base class methods
    void          clear(void);
    GMWLResponse* clone(void) const;
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    std::string   print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLResponse& pnt);
    void free_members(void);
};

#endif /* GMWLRESPONSE_HPP */
