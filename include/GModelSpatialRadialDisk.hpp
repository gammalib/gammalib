/***************************************************************************
 *      GModelSpatialRadialDisk.hpp - Radial disk source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Christoph Deil                              *
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

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GPhoton;
class GRan;
class GSkyDir;
class GSkyRegion;
class GXmlElement;


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
    GModelSpatialRadialDisk(const bool& dummy, const std::string& type);
    GModelSpatialRadialDisk(const GSkyDir&     dir,
                            const double&      radius,
                            const std::string& coordsys = "CEL");
    explicit GModelSpatialRadialDisk(const GXmlElement& xml);
    GModelSpatialRadialDisk(const GModelSpatialRadialDisk& model);
    virtual ~GModelSpatialRadialDisk(void);

    // Operators
    virtual GModelSpatialRadialDisk& operator=(const GModelSpatialRadialDisk& model);

    // Implemented pure virtual base class methods
    virtual void                     clear(void);
    virtual GModelSpatialRadialDisk* clone(void) const;
    virtual std::string              classname(void) const;
    virtual double                   eval(const double&  theta,
                                          const GEnergy& energy,
                                          const GTime&   time,
                                          const bool&    gradients = false) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime&   time,
                                        GRan&          ran) const;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const;
    virtual double                   theta_max(void) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double radius(void) const;
    void   radius(const double& radius);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialRadialDisk& model);
    void         free_members(void);
    void         update(void) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_radius; //!< Disk radius (degrees)

    // Cached members used for pre-computations
    mutable double m_last_radius;   //!< Last disk radius
    mutable double m_radius_rad;    //!< Radius in radians
    mutable double m_norm;          //!< Normalization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialDisk").
 ***************************************************************************/
inline
std::string GModelSpatialRadialDisk::classname(void) const
{
    return ("GModelSpatialRadialDisk");
}


/***********************************************************************//**
 * @brief Return disk radius
 *
 * @return Disk radius (degrees).
 *
 * Returns the radius of the disk in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialDisk::radius(void) const
{
    return (m_radius.value());
}


/***********************************************************************//**
 * @brief Set disk radius
 *
 * @param[in] radius Disk radius (degrees).
 *
 * Sets the radius of the disk in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialDisk::radius(const double& radius)
{
    m_radius.value(radius);
    return;
}

#endif /* GMODELSPATIALRADIALDISK_HPP */
