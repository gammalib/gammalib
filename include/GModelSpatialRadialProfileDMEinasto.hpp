/***************************************************************************
 *   GModelSpatialRadialProfileDMEinasto.hpp - DM Einasto profile class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Nathan Kelley-Hoskins                            *
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
 * @file GModelSpatialRadialProfileDMEinasto.hpp
 * @brief Dark Matter Einasto profile model class interface definition
 * @author Nathan Kelley-Hoskins
 */

#ifndef GMODELSPATIALRADIALPROFILEDMEINASTO_HPP
#define GMODELSPATIALRADIALPROFILEDMEINASTO_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadialProfile.hpp"
#include "GModelPar.hpp"
#include "GFunction.hpp"
#include "GIntegral.hpp"

/* __ Forward declaration ________________________________________________ */
class GXmlElement;


/**************************************************************************
 * @class GModelSpatialRadialProfileDMEinasto
 *
 * @brief Radial Dark Matter Einasto profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a radial Dark Matter profile, using an Einasto density halo.
 ***************************************************************************/
class GModelSpatialRadialProfileDMEinasto : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileDMEinasto(void);
    explicit GModelSpatialRadialProfileDMEinasto(const GXmlElement& xml);
    GModelSpatialRadialProfileDMEinasto(const GModelSpatialRadialProfileDMEinasto& model);
    virtual ~GModelSpatialRadialProfileDMEinasto(void);

    // Operators
    virtual GModelSpatialRadialProfileDMEinasto& operator=(const GModelSpatialRadialProfileDMEinasto& model);

    // Implemented pure virtual base class methods
    virtual void                                 clear(void);
    virtual GModelSpatialRadialProfileDMEinasto* clone(void) const;
    virtual std::string                          classname(void) const;
    virtual std::string                          type(void) const;
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
    double alpha(void) const;
    void   alpha(const double& alpha);
    double mass_density( const double& radius ) const ;
    double jfactor( const double& angle ) const ;

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GModelSpatialRadialProfileDMEinasto& model);
    void           free_members(void);
    virtual double profile_value(const double& theta) const;
    void           update(void) const;

    // Integration kernel for line-of-sight integral
    class halo_kernel_los : public GFunction {
    public :
        halo_kernel_los(const double& scale_radius,
                        const double& halo_distance,
                        const double& alpha,
                        const double& theta,
                        const double& core_radius ) :
                        m_scale_radius(scale_radius),
                        m_halo_distance(halo_distance),
                        m_alpha(alpha),
                        m_theta(theta),
                        m_core_radius(core_radius) {}
        double eval(const double& los);
    protected :
        double m_scale_radius;
        double m_halo_distance;
        double m_alpha;
        double m_theta;
        double m_core_radius;
    };

    // Protected members
    GModelPar m_theta_min;     //!< Minimum theta angle
    GModelPar m_theta_max;     //!< Maximum theta angle
    GModelPar m_scale_radius;  //!< Scale radius of halo profile
    GModelPar m_scale_density; //!< Scale density of halo profile
    GModelPar m_halo_distance; //!< Distance from earth to halo center
    GModelPar m_alpha;         //!< Einasto spatial power index
    GModelPar m_core_radius;   //!< Core radius
    
    // Cached members used for pre-computation
    mutable double m_last_scale_radius;
    mutable double m_last_scale_density;
    mutable double m_mass_radius;
    mutable double m_scale_density_squared;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialProfileDMEinasto").
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileDMEinasto::classname(void) const
{
    return ("GModelSpatialRadialProfileDMEinasto");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return "DMEinastoProfile".
 *
 * Returns the type of the radial profile model.
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileDMEinasto::type(void) const
{
    return "DMEinastoProfile";
}

/***********************************************************************//**
 * @brief Return scale radius
 *
 * @return scale radius (kpc).
 *
 * Returns the scale radius of the halo profile in kpc.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMEinasto::scale_radius(void) const
{
    return (m_scale_radius.value());
}

/***********************************************************************//**
 * @brief Set scale radius
 *  
 * @param[in] radius scale radius (kpc).
 *
 * Sets the scale radius of the halo profile in kpc.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMEinasto::scale_radius(const double& scale_radius)
{
    m_scale_radius.value(scale_radius);
    return;
}

/***********************************************************************//**
 * @brief Return scale density
 *
 * @return scale radius (GeV/cm^3).
 *
 * Returns the scale radius of the halo profile in kpc.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMEinasto::scale_density(void) const
{
    return (m_scale_density.value());
}

/***********************************************************************//**
 * @brief Set scale density 
 *  
 * @param[in] radius scale density (GeV/cm^3).
 *
 * Sets the scale density (mass/volume density at the scale radius) of 
 * the halo profile in GeV/cm^3.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMEinasto::scale_density(const double& scale_density)
{
    m_scale_density.value(scale_density);
    return;
}

/***********************************************************************//**
 * @brief Return halo distance
 *
 * @return halo distance (kpc).
 *
 * Returns the distance to the halo center in kpc.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMEinasto::halo_distance(void) const
{
    return (m_halo_distance.value());
}

/***********************************************************************//**
 * @brief Set halo distance
 *  
 * @param[in] halo distance (kpc).
 *
 * Sets the distance between the observer and the halo center in kpc.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMEinasto::halo_distance(const double& halo_distance)
{
    m_halo_distance.value(halo_distance);
    return;
}

/***********************************************************************//**
 * @brief Return Einasto alpha power index
 *
 * @return alpha (unitless).
 *
 * Returns the alpha power index in the Einasto halo density function.
 ***************************************************************************/
inline
double GModelSpatialRadialProfileDMEinasto::alpha(void) const
{
    return (m_alpha.value());
}

/***********************************************************************//**
 * @brief Set Einasto alpha power index
 *  
 * @param[in] alpha (unitless).
 *
 * Sets the Einasto profile power index, should be positive.
 ***************************************************************************/
inline
void GModelSpatialRadialProfileDMEinasto::alpha(const double& alpha)
{
    m_alpha.value(alpha);
    return;
}

#endif /* GMODELSPATIALRADIALPROFILEDMEINASTO_HPP */
