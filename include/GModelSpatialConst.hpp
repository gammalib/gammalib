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
 * @brief GModelSpatialConst class interface definition.
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
 * @brief Isotropic model interface definition.
 *
 * This class implements the spatial component of the factorised source
 * model for an isotropic source.
 ***************************************************************************/
class GModelSpatialConst  : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialConst(void);
    explicit GModelSpatialConst(const GXmlElement& xml);
    GModelSpatialConst(const GModelSpatialConst& model);
    virtual ~GModelSpatialConst(void);

    // Operators
    GModelSpatialConst& operator= (const GModelSpatialConst& model);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpatialConst* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "ConstantValue"; }
    double              eval(const GSkyDir& srcDir);
    double              eval_gradients(const GSkyDir& srcDir);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;
    bool                isptsource(void) const { return false; }

    // Other methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialConst& model);
    void free_members(void);

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[1];          //!< Pointers to parameters
    GModelPar  m_value;           //!< Value
};

#endif /* GMODELSPATIALCONST_HPP */
