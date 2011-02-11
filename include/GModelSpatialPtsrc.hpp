/***************************************************************************
 *      GModelSpatialPtsrc.hpp  -  Spatial point source model class        *
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
 * @file GModelSpatialPtsrc.hpp
 * @brief Point source spatial model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALPTSRC_HPP
#define GMODELSPATIALPTSRC_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialPtsrc
 *
 * @brief Point source spatial model
 *
 * This class implements the spatial component of the factorised source
 * model for a point source.
 ***************************************************************************/
class GModelSpatialPtsrc : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialPtsrc(void);
    explicit GModelSpatialPtsrc(const GSkyDir& dir);
    explicit GModelSpatialPtsrc(const GXmlElement& xml);
    GModelSpatialPtsrc(const GModelSpatialPtsrc& model);
    virtual ~GModelSpatialPtsrc(void);

    // Operators
    virtual GModelSpatialPtsrc& operator=(const GModelSpatialPtsrc& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialPtsrc* clone(void) const;
    virtual std::string         type(void) const { return "SkyDirFunction"; }
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double  ra(void) const { return m_ra.real_value(); }
    double  dec(void) const { return m_dec.real_value(); }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialPtsrc& model);
    void free_members(void);

    // Protected members
    GModelPar m_ra;          //!< Right Ascension (deg)
    GModelPar m_dec;         //!< Declination (deg)
};

#endif /* GMODELSPATIALPTSRC_HPP */
