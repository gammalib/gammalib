/***************************************************************************
 *        GModelSpatialConst.hpp  -  Spatial isotropic model class         *
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
 * @file GModelSpatialConst.hpp
 * @brief Isotropic spatial model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALCONST_HPP
#define GMODELSPATIALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialConst
 *
 * @brief Isotropic spatial model
 *
 * This class implements the spatial component of the factorised source
 * model for an isotropic source.
 ***************************************************************************/
class GModelSpatialConst : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialConst(void);
    explicit GModelSpatialConst(const GXmlElement& xml);
    GModelSpatialConst(const GModelSpatialConst& model);
    virtual ~GModelSpatialConst(void);

    // Operators
    GModelSpatialConst& operator= (const GModelSpatialConst& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialConst* clone(void) const;
    virtual std::string         type(void) const { return "ConstantValue"; }
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialConst& model);
    void free_members(void);

    // Protected members
    GModelPar m_value;         //!< Value
};

#endif /* GMODELSPATIALCONST_HPP */
