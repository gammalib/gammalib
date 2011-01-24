/***************************************************************************
 *       GModelSky.i  -  Abstract virtual sky model class python I/F       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSky.i
 * @brief GModelSky class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSky.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSky
 *
 * @brief Abstract virtual sky model class python interface defintion.
 ***************************************************************************/
class GModelSky : public GModel {
public:
    // Constructors and destructors
    GModelSky(void);
    explicit GModelSky(const GXmlElement& xml);
    explicit GModelSky(const GXmlElement& spatial, const GXmlElement& spectral);
    GModelSky(const GModelSky& model);
    virtual ~GModelSky(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModelSky*  clone(void) const = 0;
    virtual std::string type(void) const = 0;

    // Implemented pure virtual methods
    int         size(void) const;
    double      eval(const GEvent& event, const GObservation& obs);
    double      eval_gradients(const GEvent& event, const GObservation& obs);
    void        read(const GXmlElement& xml);
    void        write(GXmlElement& xml) const;

    // Other methods
    GModelSpatial*  spatial(void) const;
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
    double          value(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime);
    GVector         gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                              const GTime& srcTime);
    GPhotons        mc(const double& area, const GSkyDir& dir, const double& radius,
                       const GEnergy& emin, const GEnergy& emax,
                       const GTime& tmin, const GTime& tmax,
                       GRan& ran);
};


/***********************************************************************//**
 * @brief GModelSky class extension
 ***************************************************************************/
%extend GModelSky {
};


/***********************************************************************//**
 * @brief GModelSky type casts
 ***************************************************************************/
%inline %{
    GModelSky* cast_GModelSky(GModel* model) {
        return dynamic_cast<GModelSky*>(model);
    }
%};

