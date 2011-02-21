/***************************************************************************
 *        GCTAResponse.i  -  CTA instrument response function class        *
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
 * @file GCTAResponse.i
 * @brief CTA instrument response function Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponse.hpp"
%}


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief CTA instrument response function
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    explicit GCTAResponse(const std::string& rspname, const std::string& caldb);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCTAResponse* clone(void) const;
    virtual bool          hasedisp(void) const;
    virtual bool          hastdisp(void) const;
    virtual double        irf(const GInstDir&     obsDir,
                              const GEnergy&      obsEng,
                              const GTime&        obsTime,
                              const GSkyDir&      srcDir,
                              const GEnergy&      srcEng,
                              const GTime&        srcTime,
                              const GObservation& obs) const;
    virtual double        npred(const GSkyDir&      srcDir,
                                const GEnergy&      srcEng,
                                const GTime&        srcTime,
                                const GObservation& obs) const;

    // Overload virtual base class methods
    virtual double irf_extended(const GInstDir&             obsDir,
                                const GEnergy&              obsEng,
                                const GTime&                obsTime,
                                const GModelExtendedSource& model,
                                const GEnergy&              srcEng,
                                const GTime&                srcTime,
                                const GObservation&         obs) const;
    virtual double irf_diffuse(const GInstDir&            obsDir,
                               const GEnergy&             obsEng,
                               const GTime&               obsTime,
                               const GModelDiffuseSource& model,
                               const GEnergy&             srcEng,
                               const GTime&               srcTime,
                               const GObservation&        obs) const;

    // Other Methods
    GCTAEventAtom* mc(const double& area, const GPhoton& photon,
                      const GPointing& pnt, GRan& ran) const;
    void           caldb(const std::string& caldb);
    std::string    caldb(void) const;
    void           load(const std::string& rspname);
    void           eps(const double& eps);
    const double&  eps(void) const;

    // Low-level response methods
    double aeff(const double& theta,
                const double& phi,
                const double& zenith,
                const double& azimuth,
                const double& srcLogEng) const;
    double psf(const double& delta,
               const double& theta,
               const double& phi,
               const double& zenith,
               const double& azimuth,
               const double& srcLogEng) const;
    double psf_delta_max(const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth,
                         const double& srcLogEng) const;
    double edisp(const double& obsLogEng,
                 const double& theta,
                 const double& phi,
                 const double& zenith,
                 const double& azimuth,
                 const double& srcLogEng) const;
    double npsf(const GSkyDir&      srcDir,
                const double&       srcLogEng,
                const GTime&        srcTime,
                const GCTAPointing& pnt,
                const GCTARoi&      roi) const;
    double nedisp(const GSkyDir&      srcDir,
                  const GEnergy&      srcEng,
                  const GTime&        srcTime,
                  const GCTAPointing& pnt,
                  const GEbounds&     ebds) const;

    // Analytical PSF implementation
    double psf_dummy(const double& delta, const double& sigma) const;
    double psf_dummy_sigma(const double& srcLogEng) const;
    double psf_dummy_max(const double& sigma) const;
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponse {
    GCTAResponse copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GCTAResponse type casts
 ***************************************************************************/
%inline %{
    GCTAResponse* cast_GCTAResponse(GResponse* arg) {
        GCTAResponse* rsp = dynamic_cast<GCTAResponse*>(arg);
        if (rsp == NULL)
            throw GException::bad_type("cast_GCTAResponse(GResponse*)",
                                       "GResponse not of type GCTAResponse");
        return rsp;
    }
%}
