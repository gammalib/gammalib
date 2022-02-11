/***************************************************************************
 *    GModelSpatialRadialGeneralGauss.hpp - Generalised radial Gaussian    *
 *                           source model class                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021-2022 by Luigi Tibaldo                               *
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
 * @file GModelSpatialRadialGeneralGauss.hpp
 * @brief Generalized radial Gaussian model class interface definition
 * @author Luigi Tibaldo
 */

#ifndef GMODELSPATIALRADIALGENERALGAUSS_HPP
#define GMODELSPATIALRADIALGENERALGAUSS_HPP

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


/***********************************************************************//**
 * @class GModelSpatialRadialGeneralGauss
 *
 * @brief Generalized radial Gaussian model class
 *
 * This class implements the spatial component of the factorised source
 * model for a generalised radial Gaussian source.
 ***************************************************************************/
class GModelSpatialRadialGeneralGauss : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialGeneralGauss(void);
    GModelSpatialRadialGeneralGauss(const GSkyDir&     dir,
                                    const double&      radius,
                                    const double&      ridx,
                                    const std::string& coordsys = "CEL");
    explicit GModelSpatialRadialGeneralGauss(const GXmlElement& xml);
    GModelSpatialRadialGeneralGauss(const GModelSpatialRadialGeneralGauss& model);
    virtual ~GModelSpatialRadialGeneralGauss(void);

    // Operators
    virtual GModelSpatialRadialGeneralGauss& operator=(const GModelSpatialRadialGeneralGauss& model);

    // Implemented pure virtual methods
    virtual void                             clear(void);
    virtual GModelSpatialRadialGeneralGauss* clone(void) const;
    virtual std::string                      classname(void) const;
    virtual double                           eval(const double&  theta,
                                                  const GEnergy& energy,
                                                  const GTime&   time,
                                                  const bool&    gradients = false) const;
    virtual GSkyDir                          mc(const GEnergy& energy,
                                                const GTime& time,
                                                GRan& ran) const;
    virtual bool                                contains(const GSkyDir& dir,
                                                         const double&  margin = 0.0) const;
    virtual double                           theta_max(void) const;
    virtual void                             read(const GXmlElement& xml);
    virtual void                             write(GXmlElement& xml) const;
    virtual std::string                      print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double radius(void) const;
    void   radius(const double& radius);
    double ridx(void) const;
    void   ridx(const double& ridx);

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModelSpatialRadialGeneralGauss& model);
    void         free_members(void);
    void         update(void) const;
    virtual void set_region(void) const;

    // Protected members
    GModelPar m_radius;                 //!< Gaussian width (deg)
    GModelPar m_ridx;                   //!< Reciprocal of exponent of the radial profile

    // Cached members used for pre-computations
    mutable double m_last_radius;       //!< Last radius
    mutable double m_inv_radius_rad;    //!< radius(rad)^-1
    mutable double m_last_ridx;         //!< Last reciprocal radial index
    mutable double m_inv_ridx;          //!< Spatial profile index
    mutable double m_value_norm;        //!< 1/(2pi radius(rad)^2 ridx Gamma(ridx))
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialGeneralGauss").
 ***************************************************************************/
inline
std::string GModelSpatialRadialGeneralGauss::classname(void) const
{
    return ("GModelSpatialRadialGeneralGauss");
}


/***********************************************************************//**
 * @brief Return radius
 *
 * @return Radius (deg).
 *
 * Returns the radius in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialGeneralGauss::radius(void) const
{
    return (m_radius.value());
}

/***********************************************************************//**
 * @brief Set radius
 *
 * @param[in] radius (deg).
 *
 * Sets the radius in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialGeneralGauss::radius(const double& radius)
{
    m_radius.value(radius);
    return;
}

/***********************************************************************//**
 * @brief Return ridx
 *
 * @return Reciprocal of the radial profile index.
 *
 * Returns the reciprocal of the radial profile index.
 ***************************************************************************/
inline
double GModelSpatialRadialGeneralGauss::ridx(void) const
{
    return (m_ridx.value());
}

/***********************************************************************//**
 * @brief Set reciprocal index
 *
 * @param[in] ridx Reciprocal of the radial profile index.
 *
 * Sets the reciprocal index of the radial profile.
 ***************************************************************************/
inline
void GModelSpatialRadialGeneralGauss::ridx(const double& ridx)
{
    m_ridx.value(ridx);
    return;
}

#endif /* GMODELSPATIALRADIALGENERALGAUSS_HPP */
