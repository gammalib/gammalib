/***************************************************************************
 *     GModelSpatialPointSource.hpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialPointSource.hpp
 * @brief Point source spatial model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALPOINTSOURCE_HPP
#define GMODELSPATIALPOINTSOURCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
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
 * @class GModelSpatialPointSource
 *
 * @brief Point source spatial model
 *
 * This class implements a point source as the spatial component of the
 * factorised source model. The point source has two parameters: the
 * longitude and latitude of the point source location. The parameters may
 * either be specified in celestial or Galactic coordinates. For celestial
 * coordinates the parameter names are "RA" and "DEC", for Galactic
 * coordinates the parameter names are "GLON" and "GLAT".
 ***************************************************************************/
class GModelSpatialPointSource : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialPointSource(void);
    GModelSpatialPointSource(const bool& dummy, const std::string& type);
    GModelSpatialPointSource(const GSkyDir&     dir,
                             const std::string& coordsys = "CEL");
    GModelSpatialPointSource(const double&      lon,
                             const double&      lat,
                             const std::string& coordsys = "CEL");
    explicit GModelSpatialPointSource(const GXmlElement& xml);
    GModelSpatialPointSource(const GModelSpatialPointSource& model);
    virtual ~GModelSpatialPointSource(void);

    // Operators
    virtual GModelSpatialPointSource& operator=(const GModelSpatialPointSource& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialPointSource* clone(void) const;
    virtual std::string               classname(void) const;
    virtual GClassCode                code(void) const;
    virtual double                    eval(const GPhoton& photon,
                                           const bool& gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual double                    mc_norm(const GSkyDir& dir,
                                              const double&  radius) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Overloaded base class methods
    virtual double flux(const GSkyRegion& region,
                        const GEnergy&    srcEng  = GEnergy(),
                        const GTime&      srcTime = GTime()) const;

    // Other methods
    std::string    coordsys(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialPointSource& model);
    void         free_members(void);
    bool         is_celestial(void) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_lon;    //!< Right Ascension or Galactic longitude (deg)
    GModelPar m_lat;    //!< Declination or Galactic latitude (deg)

    // Cached members for sky direction handling
    mutable GSkyDir m_dir;      //!< Sky direction representing parameters
    mutable double  m_last_lon; //!< Last longitude
    mutable double  m_last_lat; //!< Last latitude
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialPointSource").
 ***************************************************************************/
inline
std::string GModelSpatialPointSource::classname(void) const
{
    return ("GModelSpatialPointSource");
}


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialPointSource.
 *
 * Returns the code GModelSpatialPointSource of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialPointSource::code(void) const
{
    return GMODEL_SPATIAL_POINT_SOURCE;
}


/***********************************************************************//**
 * @brief Return normalization of point source for Monte Carlo simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (degrees).
 * @return Normalization.
 *
 * Returns the normalization for a point source within a circular region.
 * The normalization is 1 if the point source falls within the circle
 * defined by @p dir and @p radius, 0 otherwise.
 ***************************************************************************/
inline
double GModelSpatialPointSource::mc_norm(const GSkyDir& dir,
                                         const double&  radius) const
{
    double norm = (dir.dist_deg(this->dir()) <= radius) ? 1.0 : 0.0;
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
std::string GModelSpatialPointSource::coordsys(void) const
{
    return (is_celestial() ? "CEL" : "GAL");
}


/***********************************************************************//**
 * @brief Check if model holds celestial coordinates
 *
 * @return True if model holds celestial coordinates, false otherwise.
 ***************************************************************************/
inline
bool GModelSpatialPointSource::is_celestial(void) const
{
    return (m_lon.name() == "RA");
}

#endif /* GMODELSPATIALPOINTSOURCE_HPP */
