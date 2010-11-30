/***************************************************************************
 *            GModels.i  -  Model container class SWIG interface           *
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
 * @file GModels.i
 * @brief GModels class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModels.hpp"
%}


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class SWIG interface defintion.
 ***************************************************************************/
class GModels : public GOptimizerPars {
public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    GModels(const std::string& filename);
    ~GModels(void);
 
    // Methods
    void   clear(void);
    int    size(void) const { return m_elements; }
    void   append(const GModel& model);
    void   load(const std::string& filename);
    void   save(const std::string& filename) const;
    void   read(const GXml& xml);
    void   write(GXml& xml) const;
    double value(const GSkyDir& srcDir, const GEnergy& srcEng,
                 const GTime& srcTime);
    double eval(const GInstDir& obsDir, const GEnergy& obsEng,
                const GTime& obsTime, const GResponse& rsp,
                const GPointing& pnt);
    double eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng,
                          const GTime& obsTime, const GResponse& rsp,
                          const GPointing& pnt);
};


/***********************************************************************//**
 * @brief GModels class extension
 ***************************************************************************/
%extend GModels {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GModel* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, 
                                       self->size());
    }
    GModels copy() {
        return (*self);
    }
};
