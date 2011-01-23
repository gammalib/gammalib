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
 * @brief GModelSpatialPtsrc class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALPTSRC_HPP
#define GMODELSPATIALPTSRC_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialPtsrc
 *
 * @brief Point source model interface definition.
 *
 * This class implements the spatial component of the factorised source
 * model for a point source.
 ***************************************************************************/
class GModelSpatialPtsrc  : public GModelSpatial {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelSpatialPtsrc& model);
    friend GLog&         operator<< (GLog& log, const GModelSpatialPtsrc& model);

public:
    // Constructors and destructors
    GModelSpatialPtsrc(void);
    explicit GModelSpatialPtsrc(const GSkyDir& dir);
    explicit GModelSpatialPtsrc(const GXmlElement& xml);
    GModelSpatialPtsrc(const GModelSpatialPtsrc& model);
    virtual ~GModelSpatialPtsrc(void);

    // Operators
    GModelSpatialPtsrc& operator= (const GModelSpatialPtsrc& model);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpatialPtsrc* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "SkyDirFunction"; }
    double              eval(const GSkyDir& srcDir);
    double              eval_gradients(const GSkyDir& srcDir);
    GSkyDir             mc(GRan& ran) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;
    bool                isptsource(void) const { return true; }

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

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[2];          //!< Pointers to parameters
    GModelPar  m_ra;              //!< Right Ascension (deg)
    GModelPar  m_dec;             //!< Declination (deg)
};

#endif /* GMODELSPATIALPTSRC_HPP */
