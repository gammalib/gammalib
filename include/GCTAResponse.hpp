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
#include <vector>
#include <iostream>
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GResponse.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA instrument response function classes.
 ***************************************************************************/
class GCTAResponse : public GResponse {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAResponse& rsp);

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
    double psf(const GSkyDir& obsDir, const GEnergy& obsEng,
               const GSkyDir& srcDir, const GEnergy& srcEng,
               const GSkyDir& instPntDir, const double& instPosAng,
               const GTime& time);
    double edisp(const GSkyDir& obsDir, const GEnergy& obsEng,
                 const GSkyDir& srcDir, const GEnergy& srcEng,
                 const GSkyDir& instPntDir, const double& instPosAng,
                 const GTime& time);

    // Other Methods
    void set_caldb(const std::string& caldb);
    void load(const std::string& irfname);

private:
    // Private methods
    void          init_members(void);
    void          copy_members(const GCTAResponse& rsp);
    void          free_members(void);
    GCTAResponse* clone(void) const;
    void          read_performance_table(const std::string& filename);
    
    // Private data members
    GNodeArray          m_nodes; //!< log(E) nodes for interpolation
    std::vector<double> m_logE;  //!< log(E) = log10(E/TeV) - bin centre
    std::vector<double> m_aeff;  //!< Effective area in square metres after all cuts
    std::vector<double> m_r68;   //!< 68% containment radius of PSF post cuts in degrees
    std::vector<double> m_r80;   //!< 80% containment radius of PSF post cuts in degrees

};

#endif /* GCTARESPONSE_HPP */
