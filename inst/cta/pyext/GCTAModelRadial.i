/***************************************************************************
 * GCTAModelRadial.i  -  Radial model abstract base class python interface *
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
 * @file GCTAModelRadial.i
 * @brief Abstract radial acceptance model class Python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadial.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadial
 *
 * @brief Abstract radial acceptance model class
 ***************************************************************************/
class GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadial(void);
    GCTAModelRadial(const GCTAModelRadial& model);
    virtual ~GCTAModelRadial(void);

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GCTAModelRadial* clone(void) const = 0;
    virtual std::string      type(void) const = 0;
    virtual double           eval(const double& offset) const = 0;
    virtual double           eval_gradients(const double& offset) const = 0;
    virtual GCTAInstDir      mc(const GCTAInstDir& dir, GRan& ran) const = 0;
    virtual double           omega(void) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;

    // Methods
    int size(void) const { return m_pars.size(); }
};


/***********************************************************************//**
 * @brief GCTAModelRadial class extension
 ***************************************************************************/
%extend GCTAModelRadial {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GModelPar __getitem__(const std::string& name) {
        return (*self)[name];
    }
    void __setitem__(const int& index, const GModelPar& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    void __setitem__(const std::string& name, const GModelPar& val) {
        (*self)[name] = val;
        return;
    }
};
