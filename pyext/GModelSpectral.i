/***************************************************************************
 *         GModelSpectral.i  -  Abstract spectral model base class         *
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
 * @file GModelSpectral.i
 * @brief Abstract spectral model base class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectral.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectral
 *
 * @brief Abstract spectral model class
 ***************************************************************************/
class GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectral(void);
    GModelSpectral(const GModelSpectral& model);
    virtual ~GModelSpectral(void);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GModelSpectral* clone(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual double          eval(const GEnergy& srcEng) const = 0;
    virtual double          eval_gradients(const GEnergy& srcEng) const = 0;
    virtual double          flux(const GEnergy& emin, const GEnergy& emax) const = 0;
    virtual GEnergy         mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;

    // Methods
    int size(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectral class extension
 ***************************************************************************/
%extend GModelSpectral {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar& __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GModelPar& __getitem__(const std::string& name) {
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
