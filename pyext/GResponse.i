/***************************************************************************
 *       GResponse.i  -  Response abstract base class SWIG interface       *
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
 * @file GResponse.i
 * @brief GResponse class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponse.hpp"
%}


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract SWIG interface for the instrument response function.
 ***************************************************************************/
class GResponse {
public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GResponse*  clone(void) const = 0;
    virtual void        load(const std::string& irfname) = 0;
    virtual bool        hasedisp(void) const = 0;
    virtual bool        hastdisp(void) const = 0;

    // Reponse function computation methods
    virtual double irf(const GInstDir& obsDir, const GEnergy& obsEng, const GTime& obsTime,
                       const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                       const GPointing& pnt) const;
    virtual double diffrsp(const GEvent& event, const GModel& model,
                           const GEnergy& srcEng, const GTime& srcTime,
                           const GPointing& pnt) const;
    virtual double live(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                        const GPointing& pnt) const = 0;
    virtual double aeff(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                        const GPointing& pnt) const = 0;
    virtual double psf(const GInstDir& obsDir,
                       const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                       const GPointing& pnt) const = 0;
    virtual double edisp(const GEnergy& obsEng,
                         const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                         const GPointing& pnt) const = 0;
    virtual double tdisp(const GTime& obsTime,
                         const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                         const GPointing& pnt) const = 0;
    virtual double nirf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                        const GPointing& pnt, const GRoi& roi, const GEbounds& ebds,
                        const GGti& gti) const;
    virtual double npsf(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                        const GPointing& pnt, const GRoi& roi) const = 0;
    virtual double nedisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                          const GPointing& pnt, const GEbounds& ebds) const;
    virtual double ntdisp(const GSkyDir&  srcDir, const GEnergy& srcEng, const GTime& srcTime,
                          const GPointing& pnt, const GGti& gti) const;

    // Other methods
    virtual void        caldb(const std::string& caldb);
    virtual std::string caldb(void) const { return m_caldb; }
};


/***********************************************************************//**
 * @brief GResponse class extension
 ***************************************************************************/
%extend GResponse {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
};
