/***************************************************************************
 *       GModelRadialDisk.hpp  -  Radial disk source model class           *
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
 * @file GModelRadialDisk.hpp
 * @brief Radial disk model class interface definition
 * @author C. Deil
 */

#ifndef GMODELRADIALDISK_HPP
#define GMODELRADIALDISK_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelRadial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelRadialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelRadialDisk : public GModelRadial {

public:
    // Constructors and destructors
    GModelRadialDisk(void);
    explicit GModelRadialDisk(const GSkyDir& dir, const double& radius);
    explicit GModelRadialDisk(const GXmlElement& xml);
    GModelRadialDisk(const GModelRadialDisk& model);
    virtual ~GModelRadialDisk(void);

    // Operators
    virtual GModelRadialDisk& operator=(const GModelRadialDisk& model);

    // Implemented pure virtual methods
    virtual void              clear(void);
    virtual GModelRadialDisk* clone(void) const;
    virtual std::string       type(void) const { return "DiskFunction"; }
    virtual double            eval(const double& theta) const;
    virtual double            eval_gradients(const double& theta) const;
    virtual GSkyDir           mc(GRan& ran) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;
    virtual std::string       print(void) const;

    // Other methods
    double radius(void) const { return m_radius.real_value(); }
    void   radius(const double& radius) { m_radius.real_value(radius); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelRadialDisk& model);
    void free_members(void);
    void update(void) const;

    // Protected members
    GModelPar      m_radius;        //!< Disk radius (degrees)

    // Cached members used for pre-computations
    mutable double m_last_radius;   //!< Last disk radius
    mutable double m_radius_rad;    //!< Radius in radians
    mutable double m_norm;          //!< Normalization
};

#endif /* GMODELRADIALDISK_HPP */
