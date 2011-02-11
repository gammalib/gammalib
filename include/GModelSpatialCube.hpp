/***************************************************************************
 *         GModelSpatialCube.hpp  -  Spatial map cube model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialCube.hpp
 * @brief Spatial map cube model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALCUBE_HPP
#define GMODELSPATIALCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialCube
 *
 * @brief Spatial map cube model
 *
 * This class implements the spatial component of the factorised source
 * model for a map cube. A map cube is a set of sky maps for different
 * energies.
 *
 * @todo Loading not yet implemented. Introduce a loading flag to load cube
 *       only when it is required. This avoid unnecessary memory allocation.
 * @todo eval() and eval_gradients() methods not yet implemented.
 ***************************************************************************/
class GModelSpatialCube : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialCube(void);
    explicit GModelSpatialCube(const GXmlElement& xml);
    GModelSpatialCube(const GModelSpatialCube& model);
    virtual ~GModelSpatialCube(void);

    // Operators
    virtual GModelSpatialCube& operator= (const GModelSpatialCube& model);

    // Implemented pure virtual methods
    virtual void               clear(void);
    virtual GModelSpatialCube* clone(void) const;
    virtual std::string        type(void) const { return "MapCubeFunction"; }
    virtual double             eval(const GSkyDir& srcDir) const;
    virtual double             eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir            mc(GRan& ran) const;
    virtual void               read(const GXmlElement& xml);
    virtual void               write(GXmlElement& xml) const;
    virtual std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialCube& model);
    void free_members(void);
    void load_cube(const std::string& filename);

    // Protected members
    GModelPar   m_value;           //!< Value
    std::string m_filename;        //!< Name of map cube
};

#endif /* GMODELSPATIALCUBE_HPP */
