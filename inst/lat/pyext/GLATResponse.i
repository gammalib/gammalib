/***************************************************************************
 *               GLATResponse.i  -  Fermi LAT Response class               *
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
 * @file GLATResponse.i
 * @brief Fermi LAT Response class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATResponse.hpp"
%}


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Fermi LAT Response class python bindings
 ***************************************************************************/
class GLATResponse : public GResponse {
public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Implement pure virtual base class methods
    void          clear(void);
    GLATResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    double        irf(const GInstDir& obsDir, const GEnergy& obsEng,
                      const GTime& obsTime,
                      const GSkyDir&  srcDir, const GEnergy& srcEng,
                      const GTime& srcTime,
                      const GObservation& obs) const;
    double        irf(const GEvent& event, const GModel& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GObservation& obs) const;
    double        nirf(const GSkyDir& srcDir, const GEnergy& srcEng,
                       const GTime& srcTime,
                       const GObservation& obs) const;

    // Other Methods
    int        size(void) const { return m_aeff.size(); }
    GLATAeff*  aeff(const int& index) const;
    GLATPsf*   psf(const int& index) const;
    GLATEdisp* edisp(const int& index) const;
    void       save(const std::string& rspname) const;
    double     irf(const GLATEventAtom& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GObservation& obs) const;
    double     irf(const GLATEventBin& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GLATResponse class extension
 ***************************************************************************/
%extend GLATResponse {
    GLATResponse copy() {
        return (*self);
    }
};
