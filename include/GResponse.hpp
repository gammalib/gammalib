/***************************************************************************
 *              GResponse.hpp  -  Response abstract base class             *
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
 * @brief GResponse class definition.
 * @author J. Knodlseder
 */

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GPointing.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function classes.
 ***************************************************************************/
class GResponse {

  // Friend classes
  friend class GObservation;

public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Operators
    virtual GResponse& operator= (const GResponse& rsp);

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

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GResponse& rsp);
    void    free_members(void);
    virtual GResponse* clone(void) const = 0;
    
    // Protected data area 
    std::string m_caldb;    //!< Name of or path to the calibration database
    std::string m_rspname;  //!< Name of the instrument response
};

#endif /* GRESPONSE_HPP */
