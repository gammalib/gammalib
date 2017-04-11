/***************************************************************************
 *  GModelSpatialRadialProfileDMZhao.i - Zhao radial profile class         *
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
 * @file GModelSpatialRadialProfileDMZhao.hpp
 * @brief Radial Dark Matter halo for Zhao Density Profile
 * @author Nathan Kelley-Hoskins
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialProfileDMZhao.hpp"
%}


/**************************************************************************
 * @class GModelSpatialRadialProfileDMZhao
 *
 * @brief Radial DM Zhao profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Dark Matter Zhao halo radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfileDMZhao : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileDMZhao(void);
    explicit GModelSpatialRadialProfileDMZhao(const GXmlElement& xml);
    GModelSpatialRadialProfileDMZhao(const GModelSpatialRadialProfileDMZhao& model);
    virtual ~GModelSpatialRadialProfileDMZhao(void);

    // Implemented pure virtual base class methods
    virtual void                              clear(void);
    virtual GModelSpatialRadialProfileDMZhao* clone(void) const;
    virtual std::string                       classname(void) const;
    virtual std::string                       type(void) const;
    virtual double                            theta_min(void) const;
    virtual double                            theta_max(void) const;
    virtual void                              read(const GXmlElement& xml);
    virtual void                              write(GXmlElement& xml) const;

    // Other methods
    double scale_radius(void) const;
    void   scale_radius(const double& scale_radius);
    double scale_density(void) const;
    void   scale_density(const double& scale_density);
    double halo_distance(void) const ;
    void   halo_distance(const double& halo_distance);
    double alpha(void) const ;
    void   alpha(const double& alpha) ;
    double beta(void) const ;
    void   beta(const double& beta) ;
    double gamma(void) const ;
    void   gamma(const double& gamma) ;
    //double prof_val(const double& theta);
    double mass_density(const double& radius) const;
    double jfactor(const double& angle) const ;
};


/***********************************************************************//**
 * @brief GModelSpatialRadialGauss class extension
 ***************************************************************************/
%extend GModelSpatialRadialProfileDMZhao {
    GModelSpatialRadialProfileDMZhao copy() {
        return (*self);
    }
};
