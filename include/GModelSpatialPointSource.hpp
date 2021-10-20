/***************************************************************************
 *     GModelSpatialPointSource.hpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2021 by Juergen Knoedlseder                         *
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
#include "GSkyRegionCircle.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialPointSource
 *
 * @brief Point source spatial model
 *
 * This class implements a point source as the spatial component of the
 * factorised source model. The point source has two parameters: the Right
 * Ascension and Declination of the point source location.
 ***************************************************************************/
class GModelSpatialPointSource : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialPointSource(void);
    GModelSpatialPointSource(const bool& dummy, const std::string& type);
    explicit GModelSpatialPointSource(const GSkyDir& dir);
    GModelSpatialPointSource(const double& ra, const double& dec);
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
    virtual double flux(const GSkyRegion&    region,
                        const GEnergy&       srcEng  = GEnergy(),
                        const GTime&         srcTime = GTime(),
                        const GPolarization& srcPol  = GPolarization()) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    void    ra(const double& ra);
    void    dec(const double& dec);
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialPointSource& model);
    void         free_members(void);
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_ra;     //!< Right Ascension (deg)
    GModelPar m_dec;    //!< Declination (deg)
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
 * @brief Return Right Ascencion of model centre
 *
 * @return Right Ascencion of model centre (degrees).
 *
 * Returns the Right Ascension of the model centre in degrees.
 ***************************************************************************/
inline
double GModelSpatialPointSource::ra(void) const
{
    return (m_ra.value());
}


/***********************************************************************//**
 * @brief Set Right Ascencion of model centre
 *
 * @param[in] ra Right Ascencion of model centre.
 *
 * Sets the Right Ascencion of model centre.
 ***************************************************************************/
inline
void GModelSpatialPointSource::ra(const double& ra)
{
    m_ra.value(ra);
    return;
}


/***********************************************************************//**
 * @brief Return Declination of model centre
 *
 * @return Declination of model centre (degrees).
 *
 * Returns the Declination of the model centre in degrees.
 ***************************************************************************/
inline
double GModelSpatialPointSource::dec(void) const
{
    return (m_dec.value());
}


/***********************************************************************//**
 * @brief Set Declination of model centre
 *
 * @param[in] dec Declination of model centre.
 *
 * Sets the Declination of model centre.
 ***************************************************************************/
inline
void GModelSpatialPointSource::dec(const double& dec)
{
    m_dec.value(dec);
    return;
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

#endif /* GMODELSPATIALPOINTSOURCE_HPP */
