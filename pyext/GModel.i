/***************************************************************************
 *                 GModel.i  -  Model class python interface               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModel.i
 * @brief GModel class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModel.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class python interface defintion
 ***************************************************************************/
class GModel {
public:
    // Constructors and destructors
    GModel(void);
    GModel(const GModelSpatial& spatial, const GModelSpectral& spectral);
    explicit GModel(const GXmlElement& spatial, const GXmlElement& spectral);
    explicit GModel(const GModel& model);
    virtual ~GModel(void);

    // Methods
    void                 clear(void);
    GModel*              clone(void) const;
    int                  size(void) const { return m_npars; }
    std::string          name(void) const { return m_name; }
    void                 name(const std::string& name) { m_name=name; return; }
    void                 instruments(const std::string& instruments);
    GModelSpatial*       spatial(void) const { return m_spatial; }
    GModelSpectral*      spectral(void) const { return m_spectral; }
    GModelTemporal*      temporal(void) const { return m_temporal; }
    double               value(const GSkyDir& srcDir, const GEnergy& srcEng,
                               const GTime& srcTime);
    GVector              gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                                   const GTime& srcTime);
    double               eval(const GEvent& event, const GObservation& obs);
    double               eval_gradients(const GEvent& event, const GObservation& obs);
    std::vector<GPhoton> mc(const double& area, const GSkyDir& dir, const double& radius,
                            const GEnergy& emin, const GEnergy& emax,
                            const GTime& tmin, const GTime& tmax);
    bool                 isvalid(const std::string& name) const;
};


/***********************************************************************//**
 * @brief GModel class extension
 ***************************************************************************/
%extend GModel {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GModelPar& val) {
        if (index>=0 && index < self->size())
            (*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    GModel copy() {
        return (*self);
    }
};
