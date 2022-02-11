/***************************************************************************
 *    GModelSpatialRadial.hpp - Abstract radial spatial model base class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadial.hpp
 * @brief Abstract radial spatial model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELSPATIALRADIAL_HPP
#define GMODELSPATIALRADIAL_HPP

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
 * @class GModelSpatialRadial
 *
 * @brief Abstract radial spatial model base class
 *
 * This class defines the interface for a radial model as spatial component
 * of the factorized source model. Typical examples of radial components are
 * axisymmetric Disk, Gaussian or Shell sources.
 ***************************************************************************/
class GModelSpatialRadial : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialRadial(void);
    GModelSpatialRadial(const GModelSpatialRadial& model);
    explicit GModelSpatialRadial(const GXmlElement& xml);
    virtual ~GModelSpatialRadial(void);

    // Operators
    virtual GModelSpatialRadial& operator=(const GModelSpatialRadial& model);

    // Pure virtual methods
    virtual void                 clear(void) = 0;
    virtual GModelSpatialRadial* clone(void) const = 0;
    virtual std::string          classname(void) const = 0;
    virtual double               eval(const double&  theta,
                                      const GEnergy& energy,
                                      const GTime&   time,
                                      const bool&    gradients = false) const = 0;
    virtual GSkyDir              mc(const GEnergy& energy,
                                    const GTime&   time,
                                    GRan&          ran) const = 0;
    virtual bool                 contains(const GSkyDir& dir,
                                          const double&  margin = 0.0) const = 0;
    virtual double               theta_max(void) const = 0;
    virtual std::string          print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented pure virtual base class methods
    virtual GClassCode code(void) const;
    virtual bool       is_energy_dependent(void) const;
    virtual bool       is_time_dependent(void) const;
    virtual double     eval(const GPhoton& photon,
                            const bool&    gradients = false) const;
    virtual double     mc_norm(const GSkyDir& dir, const double&  radius) const;
    virtual void       read(const GXmlElement& xml);
    virtual void       write(GXmlElement& xml) const;

    // Other methods
    std::string    coordsys(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialRadial& model);
    void         free_members(void);
    bool         is_celestial(void) const;
    virtual void set_region(void) const = 0;

    // Protected members
    GModelPar m_lon;    //!< Right Ascension or Galactic longitude (deg)
    GModelPar m_lat;    //!< Declination or Galactic latitude (deg)

    // Cached members for sky direction handling
    mutable GSkyDir m_dir;      //!< Sky direction representing parameters
    mutable double  m_last_lon; //!< Last longitude
    mutable double  m_last_lat; //!< Last latitude
};


/***********************************************************************//**
 * @brief Return class code
 *
 * @return GModelSpatialRadial.
 *
 * Returns the code GModelSpatialRadial of the class.
 ***************************************************************************/
inline
GClassCode GModelSpatialRadial::code(void) const
{
    return GMODEL_SPATIAL_RADIAL;
}


/***********************************************************************//**
 * @brief Signals whether radial model is energy dependent
 *
 * @return True if radial model is energy dependent, false otherwise.
 *
 * Signals whether the radial model is energy dependent. This method always
 * returns false.
 ***************************************************************************/
inline
bool GModelSpatialRadial::is_energy_dependent(void) const
{
    return (false);
}


/***********************************************************************//**
 * @brief Signals whether radial model is time dependent
 *
 * @return True if radial model is time dependent, false otherwise.
 *
 * Signals whether the radial model is time dependent. This method always
 * returns false.
 ***************************************************************************/
inline
bool GModelSpatialRadial::is_time_dependent(void) const
{
    return (false);
}


/***********************************************************************//**
 * @brief Return normalization of radial source for Monte Carlo simulations
 *
 * @param[in] dir Centre of simulation cone.
 * @param[in] radius Radius of simulation cone (degrees).
 * @return Normalization.
 *
 * Returns the normalization for a radial source within a circular region.
 * The normalization is 1 if the radial source falls within the circle
 * defined by @p dir and @p radius, 0 otherwise.
 ***************************************************************************/
inline
double GModelSpatialRadial::mc_norm(const GSkyDir& dir,
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
std::string GModelSpatialRadial::coordsys(void) const
{
    return (is_celestial() ? "CEL" : "GAL");
}


/***********************************************************************//**
 * @brief Check if model holds celestial coordinates
 *
 * @return True if model holds celestial coordinates, false otherwise.
 ***************************************************************************/
inline
bool GModelSpatialRadial::is_celestial(void) const
{
    return (m_lon.name() == "RA");
}

#endif /* GMODELSPATIALRADIAL_HPP */
