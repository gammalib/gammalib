/***************************************************************************
 *                   GModels.i  -  Model container class                   *
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
 * @file GModels.i
 * @brief GModels class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModels.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class Python bindings.
 ***************************************************************************/
class GModels : public GOptimizerPars {
public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    explicit GModels(const std::string& filename);
    virtual ~GModels(void);
 
    // Methods
    void     clear(void);
    int      size(void) const { return m_elements; }
    void     append(const GModel& model);
    void     load(const std::string& filename);
    void     save(const std::string& filename) const;
    void     read(const GXml& xml);
    void     write(GXml& xml) const;
    /*
    double   value(const GSkyDir& srcDir, const GEnergy& srcEng,
                   const GTime& srcTime);
    */
    double   eval(const GEvent& event, const GObservation& obs);
    double   eval_gradients(const GEvent& event, const GObservation& obs);
};


/***********************************************************************//**
 * @brief GModels class extension
 ***************************************************************************/
%extend GModels {
    char *__str__() {
        return tochar(self->print());
    }
    GModel* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GModel& val) {
        if (index>=0 && index < self->size())
            *(*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    GModels copy() {
        return (*self);
    }
};
