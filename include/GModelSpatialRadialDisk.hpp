/***************************************************************************
 *      GModelSpatialRadialDisk.hpp - Radial disk source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Christoph Deil                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialRadialDisk.hpp
 * @brief Radial disk model class interface definition
 * @author Christoph Deil
 */

#ifndef GMODELSPATIALRADIALDISK_HPP
#define GMODELSPATIALRADIALDISK_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialRadialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelSpatialRadialDisk : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialDisk(void);
    explicit GModelSpatialRadialDisk(const GSkyDir& dir, const double& radius);
    explicit GModelSpatialRadialDisk(const GXmlElement& xml);
    GModelSpatialRadialDisk(const GModelSpatialRadialDisk& model);
    virtual ~GModelSpatialRadialDisk(void);

    // Operators
    virtual GModelSpatialRadialDisk& operator=(const GModelSpatialRadialDisk& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialRadialDisk* clone(void) const;
    virtual std::string              type(void) const { return "DiskFunction"; }
    virtual double                   eval(const double& theta) const;
    virtual double                   eval_gradients(const double& theta) const;
    virtual GSkyDir                  mc(GRan& ran) const;
    virtual double                   theta_max(void) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(void) const;

    // Other methods
    double radius(void) const { return m_radius.real_value(); }
    void   radius(const double& radius) { m_radius.real_value(radius); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialRadialDisk& model);
    void free_members(void);
    void update(void) const;

    // Protected members
    GModelPar      m_radius;        //!< Disk radius (degrees)

    // Cached members used for pre-computations
    mutable double m_last_radius;   //!< Last disk radius
    mutable double m_radius_rad;    //!< Radius in radians
    mutable double m_norm;          //!< Normalization
};

#endif /* GMODELSPATIALRADIALDISK_HPP */
