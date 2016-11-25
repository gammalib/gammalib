/***************************************************************************
 *  GModelSpatialRadialProfileDMBurkert.i - Burkert radial profile class   *
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
 * @file GModelSpatialRadialProfileDMBurkert.hpp
 * @brief Radial Dark Matter halo for Burkert Density Profile
 * @author Nathan Kelley-Hoskins
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialProfileDMBurkert.hpp"
%}


/**************************************************************************
 * @class GModelSpatialRadialProfileDMBurkert
 *
 * @brief Radial DM Burkert profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Dark Matter Burkert halo radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfileDMBurkert : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileDMBurkert(void);
    explicit GModelSpatialRadialProfileDMBurkert(const GXmlElement& xml);
    GModelSpatialRadialProfileDMBurkert(const GModelSpatialRadialProfileDMBurkert& model);
    virtual ~GModelSpatialRadialProfileDMBurkert(void);

    // Implemented pure virtual base class methods
    virtual void                                 clear(void);
    virtual GModelSpatialRadialProfileDMBurkert* clone(void) const;
    virtual std::string                          classname(void) const;
    virtual std::string                          type(void) const;
    virtual double                               theta_min(void) const;
    virtual double                               theta_max(void) const;
    virtual void                                 read(const GXmlElement& xml);
    virtual void                                 write(GXmlElement& xml) const;

    // Other methods
    double scale_radius(void) const;
    void   scale_radius(const double& scale_radius);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialGauss class extension
 ***************************************************************************/
%extend GModelSpatialRadialProfileDMBurkert {
    GModelSpatialRadialProfileDMBurkert copy() {
        return (*self);
    }
};
