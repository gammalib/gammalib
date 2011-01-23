/***************************************************************************
 *     GModelDiffuseSource.i  -  Diffuse source model class python I/F     *
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
 * @file GModelDiffuseSource.i
 * @brief GModelDiffuseSource class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelDiffuseSource.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelDiffuseSource
 *
 * @brief Diffuse source model class python interface defintion
 ***************************************************************************/
class GModelDiffuseSource : public GModelSky {

public:
    // Constructors and destructors
    GModelDiffuseSource(void);
    explicit GModelDiffuseSource(const GXmlElement& xml);
    explicit GModelDiffuseSource(const GModelSpatial& spatial, const GModelSpectral& spectral);
    explicit GModelDiffuseSource(const GXmlElement& spatial, const GXmlElement& spectral);
    GModelDiffuseSource(const GModelDiffuseSource& model);
    virtual ~GModelDiffuseSource(void);

    // Implemented pure virtual methods
    void                 clear(void);
    GModelDiffuseSource* clone(void) const;
    std::string          type(void) const;
};


/***********************************************************************//**
 * @brief GModelDiffuseSource class extension
 ***************************************************************************/
%extend GModelDiffuseSource {
/*
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
*/
    GModelDiffuseSource copy() {
        return (*self);
    }
};
