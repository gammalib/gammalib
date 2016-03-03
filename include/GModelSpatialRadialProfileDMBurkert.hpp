/***************************************************************************
 *   GModelSpatialRadialProfileDMBurkert.hpp - DM Einasto profile class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @brief Dark Matter Einasto profile model class interface definition
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
 * @brief Radial Dark Matter Einasto profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a radial Dark Matter profile, using an Einasto density halo.
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
    virtual void                             clear(void);
    virtual GModelSpatialRadialProfileDMBurkert* clone(void) const;
    virtual std::string                      classname(void) const;
    virtual std::string                      type(void) const;
    virtual double                           theta_max(void) const;
    virtual void                             read(const GXmlElement& xml);
    virtual void                             write(GXmlElement& xml) const;
    virtual std::string                      print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GModelSpatialRadialProfileDMBurkert& model);
    void           free_members(void);
    virtual double profile_value(const double& theta) const;

    // Integration kernel for line-of-sight integral
    class halo_kernel_los : public GFunction {
      public :
        halo_kernel_los(const double& scale_radius ,
                        const double& halo_distance,
                        const double& theta        ) :
                        m_scale_radius(scale_radius),
                        m_halo_distance(halo_distance),
                        m_theta(theta) {}
        double eval( const double& los ) ;
      protected :
        double m_scale_radius  ;
        double m_halo_distance ;
        double m_theta         ;
    } ;

    // Protected members
    GModelPar m_scale_radius  ; //!< scale radius of halo profile
    GModelPar m_halo_distance ; //!< distance from earth to halo center
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
 * @brief Return model type
 *
 * @return "DMBurkertProfile".
 *
 * Returns the type of the radial profile model.
 ***************************************************************************/
inline
std::string GModelSpatialRadialProfileDMBurkert::type(void) const
{
    return "DMBurkertProfile";
}

#endif /* GMODELSPATIALRADIALPROFILEDMBURKERT_HPP */
