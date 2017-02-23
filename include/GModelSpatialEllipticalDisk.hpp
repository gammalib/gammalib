/***************************************************************************
 *   GModelSpatialEllipticalDisk.hpp - Elliptical disk source model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Michael Mayer                               *
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
 * @file GModelSpatialEllipticalDisk.hpp
 * @brief Elliptical disk model class interface definition
 * @author Michael Mayer
 */

#ifndef GMODELSPATIALELLIPTICALDISK_HPP
#define GMODELSPATIALELLIPTICALDISK_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialElliptical.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegionCircle.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialEllipticalDisk
 *
 * @brief Elliptical disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for an elliptical disk source, i.e. constant surface brightness
 * within ellipse and no emission outside.
 ***************************************************************************/
class GModelSpatialEllipticalDisk : public GModelSpatialElliptical {

public:
    // Constructors and destructors
    GModelSpatialEllipticalDisk(void);
    GModelSpatialEllipticalDisk(const GSkyDir& dir,
                                const double&  semimajor,
                                const double&  semiminor,
                                const double&  posangle);
    explicit GModelSpatialEllipticalDisk(const GXmlElement& xml);
    GModelSpatialEllipticalDisk(const GModelSpatialEllipticalDisk& model);
    virtual ~GModelSpatialEllipticalDisk(void);

    // Operators
    virtual GModelSpatialEllipticalDisk& operator=(const GModelSpatialEllipticalDisk& model);

    // Implemented pure virtual base class methods
    virtual void                         clear(void);
    virtual GModelSpatialEllipticalDisk* clone(void) const;
    virtual std::string                  classname(void) const;
    virtual std::string                  type(void) const;
    virtual double                       eval(const double&  theta,
                                              const double&  posangle,
                                              const GEnergy& energy,
                                              const GTime&   time,
                                              const bool&    gradients = false) const;
    virtual GSkyDir                      mc(const GEnergy& energy,
                                            const GTime& time,
                                            GRan& ran) const;
    virtual bool                         contains(const GSkyDir& dir,
                                                  const double&  margin = 0.0) const;
    virtual double                       theta_max(void) const;
    virtual GSkyRegion*                  region(void) const;
    virtual void                         read(const GXmlElement& xml);
    virtual void                         write(GXmlElement& xml) const;
    virtual std::string                  print(const GChatter& chatter = NORMAL) const;


protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialEllipticalDisk& model);
    void free_members(void);
    void update(void) const;
    void set_region(void) const;

    // Protected members
    mutable GSkyRegionCircle m_region; //!< Bounding circle

    // Cached members used for pre-computations
    mutable double m_last_semiminor;   //!< Last semi-minor axis
    mutable double m_last_semimajor;   //!< Last semi-major axis
    mutable double m_semiminor_rad;    //!< Radius in radians
    mutable double m_semimajor_rad;    //!< Radius in radians
    mutable double m_norm;             //!< Normalization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialEllipticalDisk").
 ***************************************************************************/
inline
std::string GModelSpatialEllipticalDisk::classname(void) const
{
    return ("GModelSpatialEllipticalDisk");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "EllipticalDisk".
 *
 * Returns the type of the elliptical disk model.
 ***************************************************************************/
inline
std::string GModelSpatialEllipticalDisk::type(void) const
{
    return "EllipticalDisk";
}


/***********************************************************************//**
 * @brief Return boundary sky region
 *
 * @return Boundary sky region.
 *
 * Returns a sky region that fully encloses the spatial model component.
 ***************************************************************************/
inline
GSkyRegion* GModelSpatialEllipticalDisk::region(void) const
{
    set_region();
    return (&m_region);
}

#endif /* GMODELSPATIALELLIPTICALDISK_HPP */
