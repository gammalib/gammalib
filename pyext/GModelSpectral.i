/***************************************************************************
 *        GModelSpectral.i  -  Spectral model class python interface       *
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
 * @brief GModelSpectral class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectral.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectral
 *
 * @brief Abstract python interface definition for the spectral model class
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
    virtual int             size(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual double          eval(const GEnergy& srcEng) = 0;
    virtual double          eval_gradients(const GEnergy& srcEng) = 0;
    virtual double          flux(const GEnergy& emin, const GEnergy& emax) const = 0;
    virtual GEnergy         mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;
};


/***********************************************************************//**
 * @brief GModelSpectral class extension
 ***************************************************************************/
%extend GModelSpectral {
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
};
