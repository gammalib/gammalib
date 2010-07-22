/***************************************************************************
 *                 GModel.i  -  Model class SWIG interface                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
 * @brief GModel class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModel.hpp"
%}


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class SWIG interface defintion.
 ***************************************************************************/
class GModel {
public:
    // Constructors and destructors
    GModel(void);
    GModel(const GModelSpatial& spatial, const GModelSpectral& spectral);
    GModel(const GModel& model);
    ~GModel(void);

    // Methods
    GModelPar* par(int index) const;
    double     value(const GSkyDir& srcDir, const GEnergy& srcEng,
                     const GTime& srcTime);
    GVector    gradients(const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GTime& srcTime);
    double     eval(const GInstDir& obsDir, const GEnergy& obsEng,
                    const GTime& obsTime, const GResponse& rsp,
                    const GPointing& pnt);
    double     eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng,
                              const GTime& obsTime, const GResponse& rsp,
                              const GPointing& pnt);
    bool       isvalid(const std::string& name);

    // Inline methods
    std::string           name(void) const { return m_name; }
    void                  name(const std::string& name) { m_name=name; return; }
    int                   npars(void) const { return m_npars; }
    const GModelSpatial*  spatial(void) const { return m_spatial; }
    const GModelSpectral* spectral(void) const { return m_spectral; }
    const GModelTemporal* temporal(void) const { return m_temporal; }
};


/***********************************************************************//**
 * @brief GModel class extension
 ***************************************************************************/
%extend GModel {
    char *__str__() {
        static char str_buffer[10001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 10001);
        str_buffer[10000] = '\0';
        return str_buffer;
    }
    GModel copy() {
        return (*self);
    }
};
