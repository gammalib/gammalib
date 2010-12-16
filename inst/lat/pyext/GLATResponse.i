/***************************************************************************
 *       GLATResponse.i  -  Fermi LAT Response class SWIG definition       *
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
 * @brief GLATResponse class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATResponse.hpp"
%}


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief SWIG interface for the GLAST LAT instrument response function.
 ***************************************************************************/
class GLATResponse : public GResponse {
public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    ~GLATResponse(void);

    // Implement pure virtual base class methods
    void          clear(void);
    GLATResponse* clone(void) const;
    void          load(const std::string& rspname);
    bool          hasedisp(void) const { return false; }
    bool          hastdisp(void) const { return false; }
    
    // Implemented response computation methods
    double irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
               const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
               const GPointing& pnt) const;
    double nirf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                const GGti& gti) const;
    double diffrsp(const GEvent& event, const GModel& model,
                   const GEnergy& srcEng, const GTime& srcTime,
                   const GPointing& pnt) const;
    double live(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt) const;
    //double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //            const GPointing& pnt) const;
    //double psf(const GInstDir& obsDir,
    //           const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //           const GPointing& pnt) const;
    //double edisp(const GEnergy& obsEng,
    //             const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //             const GPointing& pnt) const;
    //double tdisp(const GTime& obsTime,
    //             const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
    //             const GPointing& pnt) const;
    double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng,
                const GTime& srcTime, const GPointing& pnt, const GRoi& roi) const;
    double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
                  const GTime& srcTime, const GPointing& pnt, const GEbounds& ebds) const;
    //double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng,
    //              const GTime& srcTime, const GPointing& pnt, const GGti& gti) const;

    // Other Methods
    int        size(void) const { return m_aeff.size(); }
    GLATAeff*  aeff(const int& index) const;
    GLATPsf*   psf(const int& index) const;
    GLATEdisp* edisp(const int& index) const;
    void       save(const std::string& rspname) const;
};


/***********************************************************************//**
 * @brief GLATResponse class extension
 ***************************************************************************/
%extend GLATResponse {
    GLATResponse copy() {
        return (*self);
    }
};
