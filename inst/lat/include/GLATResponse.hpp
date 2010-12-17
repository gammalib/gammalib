/***************************************************************************
 *               GLATResponse.hpp  -  GLAST LAT Response class             *
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
 * @file GLATResponse.hpp
 * @brief GLATResponse class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATRESPONSE_HPP
#define GLATRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATResponseTable.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GLATEdisp.hpp"
#include "GEvent.hpp"
#include "GModel.hpp"
#include "GResponse.hpp"
#include "GPointing.hpp"
#include "GInstDir.hpp"
#include "GRoi.hpp"
#include "GGti.hpp"
#include "GEbounds.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GVector.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Interface for the Fermi LAT instrument response function.
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Operators
    GLATResponse& operator= (const GLATResponse & rsp);

    // Implement pure virtual base class methods
    void          clear(void);
    GLATResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    double        irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                      const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                      const GPointing& pnt) const;
    double        irf(const GEvent& event, const GModel& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GPointing& pnt) const;
    double        nirf(const GSkyDir& srcDir, const GEnergy& srcEng, const GTime& srcTime,
                       const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                       const GGti& gti) const;
    std::string   print(void) const;

    // Other Methods
    int        size(void) const { return m_aeff.size(); }
    GLATAeff*  aeff(const int& index) const;
    GLATPsf*   psf(const int& index) const;
    GLATEdisp* edisp(const int& index) const;
    void       save(const std::string& rspname) const;
    double     irf(const GLATEventAtom& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GPointing& pnt) const;
    double     irf(const GLATEventBin& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GPointing& pnt) const;
    double     live(const GSkyDir& srcDir, const GEnergy& srcEng,
                    const GTime& srcTime, const GPointing& pnt) const;
    //double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //            const GPointing& pnt) const;
    //double psf(const GInstDir& obsDir,
    //           const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //           const GPointing& pnt) const;
    //double edisp(const GEnergy& obsEng,
    //             const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //             const GPointing& pnt) const;
    double     npsf(const GSkyDir&  srcDir, const GEnergy& srcEng,
                    const GTime& srcTime, const GPointing& pnt, const GRoi& roi) const;
    double     nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                      const GTime& srcTime, const GPointing& pnt, const GEbounds& ebds) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GLATResponse& rsp);
    void free_members(void);

    // Private members
    bool                    m_hasfront;    //!< Front IRF loaded
    bool                    m_hasback;     //!< Back IRF loaded
    std::vector<GLATAeff*>  m_aeff;        //!< Effective areas
    std::vector<GLATPsf*>   m_psf;         //!< Point spread functions
    std::vector<GLATEdisp*> m_edisp;       //!< Energy dispersions
};

#endif /* GLATRESPONSE_HPP */
