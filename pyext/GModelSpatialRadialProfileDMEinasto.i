/***************************************************************************
 *  GModelSpatialRadialProfileDMEinasto.i - Einasto radial profile class   *
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
 * @file GModelSpatialRadialProfileDMEinasto.hpp
 * @brief Radial Dark Matter halo for Einasto Density Profile
 * @author Nathan Kelley-Hoskins
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialProfileDMEinasto.hpp"
%}


/**************************************************************************
 * @class GModelSpatialRadialProfileDMEinasto
 *
 * @brief Radial DM Einasto profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Dark Matter Einasto halo radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfileDMEinasto : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileDMEinasto(void);
    explicit GModelSpatialRadialProfileDMEinasto(const GXmlElement& xml);
    GModelSpatialRadialProfileDMEinasto(const GModelSpatialRadialProfileDMEinasto& model);
    virtual ~GModelSpatialRadialProfileDMEinasto(void);

    // Implemented pure virtual base class methods
    virtual void                                 clear(void);
    virtual GModelSpatialRadialProfileDMEinasto* clone(void) const;
    virtual std::string                          classname(void) const;
    virtual double                               theta_min(void) const;
    virtual double                               theta_max(void) const;
    virtual void                                 read(const GXmlElement& xml);
    virtual void                                 write(GXmlElement& xml) const;

    // Other methods
    double scale_radius(void) const;
    void   scale_radius(const double& scale_radius);
    double scale_density(void) const;
    void   scale_density(const double& scale_density);
    double halo_distance(void) const;
    void   halo_distance(const double& halo_distance);
    double alpha(void) const;
    void   alpha(const double& alpha);
    double mass_density(const double& radius) const;
    double jfactor(const double& angle) const;
};


/***********************************************************************//**
 * @brief GModelSpatialRadialGauss class extension
 ***************************************************************************/
%extend GModelSpatialRadialProfileDMEinasto {
    GModelSpatialRadialProfileDMEinasto copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
