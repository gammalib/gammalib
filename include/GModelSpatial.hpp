/***************************************************************************
 *         GModelSpatial.hpp  -  Abstract spatial model base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatial.hpp
 * @brief GModelSpatial abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIAL_HPP
#define GMODELSPATIAL_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelPar.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract interface definition for the spatial model class.
 *
 * This class implements the spatial component of the factorized gamma-ray
 * source model. Typical examples of spatial components are a point source
 * or an intensity map. The method isptsource() signals of the spatial
 * model is indeed a point source.
 ***************************************************************************/
class GModelSpatial {

    // Friend classes
    friend class GModel;

public:
    // Constructors and destructors
    explicit GModelSpatial(void);
    GModelSpatial(const GModelSpatial& model);
    virtual ~GModelSpatial();

    // Operators
    virtual GModelSpatial& operator= (const GModelSpatial& model);

    // Virtual methods
    virtual int        npars(void) const = 0;
    virtual GModelPar* par(int index) const = 0;
    virtual double     eval(const GSkyDir& srcDir) = 0;
    virtual double     eval_gradients(const GSkyDir& srcDir) = 0;
    virtual bool       isptsource(void) const { return false; }
    virtual void       read(const GXmlElement& xml) = 0;
    virtual void       write(GXmlElement& xml) const = 0;


protected:
    // Protected methods
    void                   init_members(void);
    void                   copy_members(const GModelSpatial& model);
    void                   free_members(void);
    virtual GModelSpatial* clone(void) const = 0;
};

#endif /* GMODELSPATIAL_HPP */
