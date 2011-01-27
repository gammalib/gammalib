/***************************************************************************
 *         GCTAModelRadial.hpp  -  Radial model abstract base class        *
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
 * @file GCTAModelRadial.hpp
 * @brief GCTAModelRadial abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIAL_HPP
#define GCTAMODELRADIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAModelRadial
 *
 * @brief Abstract interface definition for the radial CTA model class
 *
 * This class implements the radial component of the CTA radial acceptance
 * model.
 ***************************************************************************/
class GCTAModelRadial {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAModelRadial& model);
    friend GLog&         operator<< (GLog& log, const GCTAModelRadial& model);

public:
    // Constructors and destructors
    GCTAModelRadial(void);
    GCTAModelRadial(const GCTAModelRadial& model);
    virtual ~GCTAModelRadial(void);

    // Operators
    virtual GModelPar&       operator() (int index) = 0;
    virtual const GModelPar& operator() (int index) const = 0;
    virtual GCTAModelRadial& operator= (const GCTAModelRadial& model);

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GCTAModelRadial* clone(void) const = 0;
    virtual int              size(void) const = 0;
    virtual std::string      type(void) const = 0;
    virtual double           eval(const double& offset) = 0;
    virtual double           eval_gradients(const double& offset) = 0;
    virtual GCTAInstDir      mc(const GCTAInstDir& dir, GRan& ran) const = 0;
    virtual double           omega(void) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;
    virtual std::string      print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadial& model);
    void free_members(void);
};

#endif /* GCTAMODELRADIAL_HPP */
