/***************************************************************************
 *   GModelSpatialRadialProfileDMBurkert.hpp - DM Burkert profile class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2020 by Nathan Kelley-Hoskins                       *
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
 * @file GModelSpatialRadialProfileDMBurkert.hpp
 * @brief Dark Matter Burkert profile model class interface definition
 * @author Nathan Kelley-Hoskins
 */

#ifndef GMODELSPATIALRADIALPROFILEDMBURKERT_HPP
#define GMODELSPATIALRADIALPROFILEDMBURKERT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadialProfile.hpp"
#include "GModelPar.hpp"
#include "GFunction.hpp"
#include "GIntegral.hpp"

/* __ Forward declaration ________________________________________________ */
class GXmlElement;


/**************************************************************************
 * @class GModelSpatialRadialProfileDMBurkert
 *
 * @brief Radial Dark Matter Burkert profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a radial Dark Matter profile, using an Burkert density halo.
 ***************************************************************************/
class GModelSpatialRadialProfileDMBurkert : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileDMBurkert(void);
    explicit GModelSpatialRadialProfileDMBurkert(const GXmlElement& xml);
    GModelSpatialRadialProfileDMBurkert(const GModelSpatialRadialProfileDMBurkert& model);
    virtual ~GModelSpatialRadialProfileDMBurkert(void);

    // Operators
    virtual GModelSpatialRadialProfileDMBurkert& operator=(const GModelSpatialRadialProfileDMBurkert& model);

    // Implemented pure virtual base class methods
    virtual void                                 clear(void);
    virtual GModelSpatialRadialProfileDMBurkert* clone(void) const;
    virtual std::string                          classname(void) const;
    virtual double                               theta_min(void) const;
    virtual double                               theta_max(void) const;
    virtual void                                 read(const GXmlElement& xml);
    virtual void                                 write(GXmlElement& xml) const;
    virtual std::string                          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double scale_radius(void) const;
    void   scale_radius(const double& scale_radius);
    double scale_density(void) const;
    void   scale_density(const double& scale_density);
    double halo_distance(void) const;
    void   halo_distance(const double& halo_distance);
    double mass_density(const double& radius ) const;
    double jfactor(const double& angle) const;

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GModelSpatialRadialProfileDMBurkert& model);
    void           free_members(void);
    virtual double profile_value(const double& theta) const;
    void           update(void) const;

    // Integration kernel for line-of-sight integral
    class halo_kernel_los : public GFunction {
    public :
        halo_kernel_los(const double& scale_radius,
                        const double& halo_distance,
                        const double& theta,
                        const double& core_radius) :
                        m_scale_radius(scale_radius),
                        m_halo_distance(halo_distance),
                        m_theta(theta),
                        m_core_radius(core_radius) {}
        double eval(const double& los);
    protected :
        double m_scale_radius;
        double m_halo_distance;
        double m_theta;
        double m_core_radius;
    };

    // Protected members
    GModelPar m_theta_min;     //!< Minimum theta angle
    GModelPar m_theta_max;     //!< Maximum theta angle
    GModelPar m_scale_radius;  //!< Scale radius of halo profile
    GModelPar m_scale_density; //!< Scale density of halo profile
    GModelPar m_halo_distance; //!< Distance from Earth to halo center
    GModelPar m_core_radius;   //!< Core radius
    
    // Cached members used for precomputation
    mutable double m_last_scale_radius;
    mutable double m_last_scale_density;
    mutable double m_mass_radius;
    mutable double m_scale_density_squared;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialProfileDMBurkert").
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileDMBurkert::classname(void) const
{
    return ("GModelSpatialRadialProfileDMBurkert");
}


/***********************************************************************//**
 * @brief Return scale radius 
 *
 * @return Scale radius (kpc).
 *
 * Returns the scale radius of the halo profile in kpc.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMBurkert::scale_radius(void) const
{
    return (m_scale_radius.value());
}


/***********************************************************************//**
 * @brief Set scale radius
 *
 * @param[in] radius Scale radius (kpc).
 *
 * Sets the scale radius of the halo profile in kpc.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMBurkert::scale_radius(const double& radius)
{
    m_scale_radius.value(radius);
    return;
}


/***********************************************************************//**
 * @brief Return scale density
 *
 * @return Scale density (GeV/cm^3).
 *
 * Returns the scale density (mass/volume density at the scale radius) of 
 * the halo profile in GeV/cm^3.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMBurkert::scale_density(void) const
{
    return (m_scale_density.value());
}


/***********************************************************************//**
 * @brief Set scale density
 *
 * @param[in] density Scale density (GeV/cm^3).
 *
 * Sets the scale density ( mass/volume density at the scale radius) of the 
 * halo profile in GeV/cm^3.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMBurkert::scale_density(const double& density)
{
    m_scale_density.value(density);
    return;
}


/***********************************************************************//**
 * @brief Return halo distance
 *
 * @return Halo distance (kpc).
 *
 * Returns the distance to the halo center in kpc.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMBurkert::halo_distance(void) const
{
    return (m_halo_distance.value());
}


/***********************************************************************//**
 * @brief Set halo distance
 *
 * @param[in] distance Halo distance (kpc).
 *
 * Sets the distance between the observer and the halo center in kpc.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMBurkert::halo_distance(const double& distance)
{
    m_halo_distance.value(distance);
    return;
}

#endif /* GMODELSPATIALRADIALPROFILEDMBURKERT_HPP */
