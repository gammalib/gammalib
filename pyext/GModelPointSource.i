/***************************************************************************
 *       GModelPointSource.i  -  Point source model class python I/F       *
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
 * @file GModelPointSource.i
 * @brief GModelPointSource class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelPointSource.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelPointSource
 *
 * @brief Point source model class python interface defintion
 ***************************************************************************/
class GModelPointSource : public GModelSky {

public:
    // Constructors and destructors
    GModelPointSource(void);
    explicit GModelPointSource(const GXmlElement& xml);
    explicit GModelPointSource(const GModelSpatialPtsrc& ptsrc, const GModelSpectral& spectral);
    explicit GModelPointSource(const GXmlElement& ptsrc, const GXmlElement& spectral);
    GModelPointSource(const GModelPointSource& model);
    virtual ~GModelPointSource(void);

    // Implemented pure virtual methods
    void               clear(void);
    GModelPointSource* clone(void) const;
    std::string        type(void) const;
};


/***********************************************************************//**
 * @brief GModelPointSource class extension
 ***************************************************************************/
%extend GModelPointSource {
    GModelPointSource copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GModelPointSource type casts
 ***************************************************************************/
%inline %{
    GModelPointSource* cast_GModelPointSource(GModel* model) {
        return dynamic_cast<GModelPointSource*>(model);
    }
%};
