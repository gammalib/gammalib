/***************************************************************************
 *      GModelSpatialRadialRing.hpp - Radial ring source model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Pierrick Martin                                  *
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
 * @file GModelSpatialRadialRing.hpp
 * @brief Radial ring model class interface definition
 * @author Pierrick Martin
 */

#ifndef GMODELSPATIALRADIALRING_HPP
#define GMODELSPATIALRADIALRING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GSkyRegionCircle.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialRadialRing
 *
 * @brief Ring source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a ring source, i.e. constant surface brightness between an
 * inner radius and an outer radius, and no emission elsewhere.
 ***************************************************************************/
class GModelSpatialRadialRing : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialRing(void);
    GModelSpatialRadialRing(const GSkyDir& dir, const double& radius, const double& width);
    explicit GModelSpatialRadialRing(const GXmlElement& xml);
    GModelSpatialRadialRing(const GModelSpatialRadialRing& model);
    virtual ~GModelSpatialRadialRing(void);

    // Operators
    virtual GModelSpatialRadialRing& operator=(const GModelSpatialRadialRing& model);

    // Implemented pure virtual base class methods
    virtual void                     clear(void);
    virtual GModelSpatialRadialRing* clone(void) const;
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
    virtual double                   theta_min(void) const;
    virtual double                   theta_max(void) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double radius(void) const;
    void   radius(const double& radius);
    double width(void) const;
    void   width(const double& width);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialRadialRing& model);
    void         free_members(void);
    void         update(void) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_radius;    //!< Ring inner radius (degrees)
    GModelPar m_width;     //!< Ring width (degrees)

    // Cached members used for pre-computations
    mutable double m_last_radius;          //!< Last ring radius
    mutable double m_last_width;           //!< Last ring width
    mutable double m_inner_radius_rad;     //!< Inner radius in radians
    mutable double m_outer_radius_rad;     //!< Outer radius in radians
    mutable double m_cos_inner_radius_rad; //!< Cosine of inner radius in radians
    mutable double m_cos_outer_radius_rad; //!< Cosine of outer radius in radians
    mutable double m_norm;                 //!< Normalization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialRing").
 ***************************************************************************/
inline
std::string GModelSpatialRadialRing::classname(void) const
{
    return ("GModelSpatialRadialRing");
}


/***********************************************************************//**
 * @brief Return ring inner radius
 *
 * @return Ring inner radius (degrees).
 *
 * Returns the inner radius of the ring in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialRing::radius(void) const
{
    return (m_radius.value());
}


/***********************************************************************//**
 * @brief Set ring inner radius
 *
 * @param[in] radius Ring inner radius (degrees).
 *
 * Sets the inner radius of the ring in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialRing::radius(const double& radius)
{
    m_radius.value(radius);
    return;
}


/***********************************************************************//**
 * @brief Return ring width
 *
 * @return Ring width (degrees).
 *
 * Returns the width of the ring in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialRing::width(void) const
{
    return (m_width.value());
}


/***********************************************************************//**
 * @brief Set ring width
 *
 * @param[in] width Ring width (degrees).
 *
 * Sets the width of the ring in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialRing::width(const double& width)
{
    m_width.value(width);
    return;
}

#endif /* GMODELSPATIALRADIALRING_HPP */
