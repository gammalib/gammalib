/***************************************************************************
 *  GModelSpatialElliptical.hpp - Abstract elliptical spatial model class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialElliptical.hpp
 * @brief Abstract elliptical spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALELLIPTICAL_HPP
#define GMODELSPATIALELLIPTICAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GPhoton;
class GRan;
class GSkyRegion;
class GXmlElement;


/***********************************************************************//**
 * @class GModelSpatialElliptical
 *
 * @brief Abstract elliptical spatial model base class
 *
 * This class defines the interface for an elliptical model as spatial
 * component of the factorized source model. Typical examples of elliptical
 * components are elliptical Disk, Gaussian or Shell shaped sources.
 ***************************************************************************/
class GModelSpatialElliptical : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialElliptical(void);
    GModelSpatialElliptical(const GModelSpatialElliptical& model);
    explicit GModelSpatialElliptical(const GXmlElement& xml);
    virtual ~GModelSpatialElliptical(void);

    // Operators
    virtual GModelSpatialElliptical& operator=(const GModelSpatialElliptical& model);

    // Pure virtual methods
    virtual void                     clear(void) = 0;
    virtual GModelSpatialElliptical* clone(void) const = 0;
    virtual std::string              classname(void) const = 0;
    virtual double                   eval(const double&  theta,
                                          const double&  posangle,
                                          const GEnergy& energy,
                                          const GTime&   time,
                                          const bool&    gradients = false) const = 0;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime&   time,
                                        GRan&          ran) const = 0;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const = 0;
    virtual double                   theta_max(void) const = 0;
    virtual std::string              print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented virtual base class methods
    virtual GClassCode code(void) const;
    virtual double     eval(const GPhoton& photon,
                            const bool&    gradients = false) const;
    virtual double     mc_norm(const GSkyDir& dir, const double&  radius) const;
    virtual void       read(const GXmlElement& xml);
    virtual void       write(GXmlElement& xml) const;

    // Other methods
    std::string    coordsys(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);
    double         posangle(void) const;
    void           posangle(const double& posangle);
    double         semiminor(void) const;
    double         semimajor(void) const;
    void           semiminor(const double& semiminor);
    void           semimajor(const double& semimajor);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialElliptical& model);
    void         free_members(void);
    bool         is_celestial(void) const;
    virtual void set_region(void) const = 0;

    // Protected members
    GModelPar m_lon;       //!< Right Ascension or Galactic longitude (deg)
    GModelPar m_lat;       //!< Declination or Galactic latitude (deg)
    GModelPar m_posangle;  //!< Position angle from North, counterclockwise (deg)
    GModelPar m_semiminor; //!< Semi-minor axis of ellipse (deg)
    GModelPar m_semimajor; //!< Semi-major axis of ellipse (deg)
    
    // Cached members for sky direction handling
    mutable GSkyDir m_dir;      //!< Sky direction representing parameters
    mutable double  m_last_lon; //!< Last longitude
    mutable double  m_last_lat; //!< Last latitude
};


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialElliptical.
 *
 * Returns the code GModelSpatialElliptical of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialElliptical::code(void) const
{
    return GMODEL_SPATIAL_ELLIPTICAL;
}


/***********************************************************************//**
 * @brief Return Position Angle of model
 *
 * @return Position Angle of model (degrees).
 *
 * Returns the Position Angle of model in degrees, measured counterclockwise
 * from celestial North.
 ***************************************************************************/
inline
double GModelSpatialElliptical::posangle(void) const
{
    return (m_posangle.value());
}


/***********************************************************************//**
 * @brief Set Position Angle of model
 *
 * @param[in] posangle Position Angle of model (degrees).
 *
 * Sets the Position Angle of model in degrees, measured counterclockwise
 * from celestial North.
 ***************************************************************************/
inline
void GModelSpatialElliptical::posangle(const double& posangle)
{
    m_posangle.value(posangle);
    return;
}


/***********************************************************************//**
 * @brief Return semi-minor axis of ellipse
 *
 * @return Semi-minor axis of ellipse (degrees).
 *
 * Returns the semi-minor axis of the ellipse in degrees.
 ***************************************************************************/
inline
double GModelSpatialElliptical::semiminor(void) const
{
    return (m_semiminor.value());
}


/***********************************************************************//**
 * @brief Set semi-minor axis of ellipse
 *
 * @param[in] semiminor Semi-minor axis of ellipse (degrees)
 *
 * Sets the semi-minor axis of the ellipse in degrees.
 ***************************************************************************/
inline
void GModelSpatialElliptical::semiminor(const double& semiminor)
{
    m_semiminor.value(semiminor);
    return;
}


/***********************************************************************//**
 * @brief Return semi-major axis of ellipse
 *
 * @return Semi-major axis of ellipse (degrees).
 *
 * Returns the semi-major axis of the ellipse in degrees.
 ***************************************************************************/
inline
double GModelSpatialElliptical::semimajor(void) const
{
    return (m_semimajor.value());
}


/***********************************************************************//**
 * @brief Set semi-major axis of ellipse
 *
 * @param[in] semimajor Semi-major axis of ellipse (degrees)
 *
 * Sets the semi-major axis of the ellipse in degrees.
 ***************************************************************************/
inline
void GModelSpatialElliptical::semimajor(const double& semimajor)
{
    m_semimajor.value(semimajor);
    return;
}


/***********************************************************************//**
 * @brief Return normalization of elliptical source for Monte Carlo
 *        simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (degrees).
 * @return Normalization.
 *
 * Returns the normalization for an elliptical source within a circular
 * region. The normalization is 1 if the elliptical source falls within
 * the circle define by @p dir and @p radius, 0 otherwise.
 ***************************************************************************/
inline
double GModelSpatialElliptical::mc_norm(const GSkyDir& dir,
                                        const double&  radius) const
{
    double norm = (dir.dist_deg(this->dir()) <= radius+theta_max()) ? 1.0 : 0.0;
    return (norm);
}


/***********************************************************************//**
 * @brief Return coordinate system
 *
 * @return Coordinate system of point source model.
 *
 * Returns "CEL" for a celestial coordinate system and "GAL" for a Galactic
 * coordinate system.
 ***************************************************************************/
inline
std::string GModelSpatialElliptical::coordsys(void) const
{
    return (is_celestial() ? "CEL" : "GAL");
}


/***********************************************************************//**
 * @brief Check if model holds celestial coordinates
 *
 * @return True if model holds celestial coordinates, false otherwise.
 ***************************************************************************/
inline
bool GModelSpatialElliptical::is_celestial(void) const
{
    return (m_lon.name() == "RA");
}

#endif /* GMODELSPATIALELLIPTICAL_HPP */
