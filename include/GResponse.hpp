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
 * @brief Abstract response base class definition.
 * @author J. Knodlseder
 */

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GEvent.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GModel;


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function
 *
 * The response function provides conversion methods between physical
 * parameters (such as source position, flux, ...) and the measured
 * instrumental parameters (such as measured energy, photon interaction, 
 * ...).
 *
 * @todo Replace (GInstDir,GEnergy,GTime) by (GEvent) in irf()?. This would
 *       make this method more symmetric to the diffuse diffrsp() method.
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
    virtual void        load(const std::string& irfname) = 0;
    virtual bool        hasedisp(void) const = 0;
    virtual bool        hastdisp(void) const = 0;
    virtual double      irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                            const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                            const GPointing& pnt) const = 0;
    virtual double      nirf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                             const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                             const GGti& gti) const = 0;
    virtual double      diffrsp(const GEvent& event, const GModel& model,
                                const GEnergy& srcEng, const GTime& srcTime,
                                const GPointing& pnt) const;
    virtual std::string print(void) const = 0;

    // Other methods
    virtual void        caldb(const std::string& caldb);
    virtual std::string caldb(void) const { return m_caldb; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GResponse& rsp);
    void free_members(void);
    
    // Protected data area 
    std::string m_caldb;    //!< Name of or path to the calibration database
    std::string m_rspname;  //!< Name of the instrument response
};

#endif /* GRESPONSE_HPP */
