/***************************************************************************
 *      GModelSpatialDisk.hpp  -  Spatial disk source model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Christoph Deil                                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialDisk.hpp
 * @brief Spatial disk model class interface definition
 * @author C. Deil
 */

#ifndef GMODELSPATIALDISK_HPP
#define GMODELSPATIALDISK_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelSpatialDisk : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialDisk(void);
    explicit GModelSpatialDisk(const GSkyDir& dir,
                                const double& radius);
    explicit GModelSpatialDisk(const GXmlElement& xml);
    GModelSpatialDisk(const GModelSpatialDisk& model);
    virtual ~GModelSpatialDisk(void);

    // Operators
    virtual GModelSpatialDisk& operator=(const GModelSpatialDisk& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialDisk*  clone(void) const;
    virtual std::string         type(void) const { return "DiskFunction"; }
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double  ra(void) const { return m_ra.real_value(); }
    double  dec(void) const { return m_dec.real_value(); }
    double  radius(void) const { return m_radius.real_value(); }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    radius(const double& radius) { m_radius.real_value(radius); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialDisk& model);
    void free_members(void);
    void update() const;

    // Protected members
    GModelPar      m_ra;            //!< Right Ascension of centre (deg)
    GModelPar      m_dec;           //!< Declination of centre (deg)
    GModelPar      m_radius;        //!< Disk radius (deg)

    // Cached members used for pre-computations
    mutable double m_last_radius;   //!< Last radius (deg)
    mutable double m_norm;          //!< Normalization
    mutable double m_gradient;      //!< Radius gradient
};

#endif /* GMODELSPATIALDISK_HPP */
