/***************************************************************************
 *            GModelPointSource.hpp  -  Point source model class           *
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
 * @file GModelPointSource.hpp
 * @brief GModelPointSource class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELPOINTSOURCE_HPP
#define GMODELPOINTSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSky.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpectral.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelPointSource
 *
 * @brief Point source model class interface defintion
 *
 * This class implements a factorised point source model.
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

    // Operators
    GModelPointSource& operator= (const GModelPointSource& model);

    // Implemented pure virtual methods
    void               clear(void);
    GModelPointSource* clone(void) const;
    std::string        type(void) const { return "PointSource"; }
    std::string        print(void) const;
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelPointSource& model);
    void free_members(void);
};

#endif /* GMODELPOINTSOURCE_HPP */
