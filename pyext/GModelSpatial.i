/***************************************************************************
 *           GModelSpatial.i  -  Spatial model abstract base class         *
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
 * @file GModelSpatial.i
 * @brief GModelSpatial class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatial.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract SWIG interface definition for the spatial model class.
 ***************************************************************************/
class GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatial(void);
    GModelSpatial(const GModelSpatial& model);
    virtual ~GModelSpatial(void);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GModelSpatial* clone(void) const = 0;
    virtual int            size(void) const = 0;
    virtual std::string    type(void) const = 0;
    virtual double         eval(const GSkyDir& srcDir) = 0;
    virtual double         eval_gradients(const GSkyDir& srcDir) = 0;
    virtual void           read(const GXmlElement& xml) = 0;
    virtual void           write(GXmlElement& xml) const = 0;
    virtual bool           isptsource(void) const = 0;
};


/***********************************************************************//**
 * @brief GModelSpatial class extension
 ***************************************************************************/
%extend GModelSpatial {
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
