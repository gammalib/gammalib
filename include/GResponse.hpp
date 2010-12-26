/***************************************************************************
 *              GResponse.hpp  -  Abstract response base class             *
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
 * @file GResponse.hpp
 * @brief Abstract response base class definition
 * @author J. Knodlseder
 */

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GEvent.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GModel;
class GObservation;


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function
 *
 * The response function provides conversion between physical parameters
 * (such as source position, flux, ...) and the measured instrumental
 * parameters (such as measured energy, photon interaction, ...).
 * For a given observation, the irf method returns the instrument response
 * for a given event and source model as function of the true photon energy
 * and the true photon arrival time.
 * The npred method returns the integral of the instrument response function
 * over the dataspace. This method is only required for unbinned analysis.
 ***************************************************************************/
class GResponse {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GResponse& rsp);
    friend GLog&         operator<< (GLog& log, const GResponse& rsp);

public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Operators
    virtual GResponse& operator= (const GResponse& rsp);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GResponse*  clone(void) const = 0;
    virtual bool        hasedisp(void) const = 0;
    virtual bool        hastdisp(void) const = 0;
    virtual double      irf(const GEvent& event, const GModel& model,
                            const GEnergy& srcEng, const GTime& srcTime,
                            const GObservation& obs) const = 0;
    virtual double      npred(const GModel& model, const GEnergy& srcEng,
                              const GTime& srcTime,
                              const GObservation& obs) const = 0;
    virtual std::string print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GResponse& rsp);
    void free_members(void);
};

#endif /* GRESPONSE_HPP */
