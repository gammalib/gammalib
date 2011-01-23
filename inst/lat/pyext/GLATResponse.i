/***************************************************************************
 *               GLATResponse.i  -  Fermi LAT Response class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
    bool          hasedisp(void) const;
    bool          hastdisp(void) const;
    double        irf(const GEvent& event, const GModelSky& model,
                      const GEnergy& srcEng, const GTime& srcTime,
                      const GObservation& obs) const;
    double        npred(const GModelSky& model, const GEnergy& srcEng,
                        const GTime& srcTime,
                        const GObservation& obs) const;

    // Other Methods
    void        caldb(const std::string& caldb);
    std::string caldb(void) const;
    void        load(const std::string& rspname);
    int         size(void) const;
    GLATAeff*   aeff(const int& index) const;
    GLATPsf*    psf(const int& index) const;
    GLATEdisp*  edisp(const int& index) const;
    void        save(const std::string& rspname) const;
    bool        force_mean(void);
    void        force_mean(const bool& value);

    // Reponse methods
    double irf(const GLATEventAtom& event, const GModelSky& model,
               const GEnergy& srcEng, const GTime& srcTime,
               const GObservation& obs) const;
    double irf(const GLATEventBin& event, const GModelSky& model,
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


/***********************************************************************//**
 * @brief GLATResponse type casts
 ***************************************************************************/
%inline %{
    GLATResponse* cast_GLATResponse(GResponse* rsp) {
        return dynamic_cast<GLATResponse*>(rsp);
    }
%}
