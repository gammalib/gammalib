/***************************************************************************
 *    GModelSpatialRadialProfileGauss.i - Gaussian radial profile class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadialProfileGauss.hpp
 * @brief Radial Gaussian profile model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialProfileGauss.hpp"
%}


/**************************************************************************
 * @class GModelSpatialRadialProfileGauss
 *
 * @brief Radial Gaussian profile source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian radial profile.
 ***************************************************************************/
class GModelSpatialRadialProfileGauss : public GModelSpatialRadialProfile {

public:
    // Constructors and destructors
    GModelSpatialRadialProfileGauss(void);
    explicit GModelSpatialRadialProfileGauss(const GXmlElement& xml);
    GModelSpatialRadialProfileGauss(const GSkyDir& dir, const double& sigma);
    GModelSpatialRadialProfileGauss(const GModelSpatialRadialProfileGauss& model);
    virtual ~GModelSpatialRadialProfileGauss(void);

    // Implemented pure virtual base class methods
    virtual void                             clear(void);
    virtual GModelSpatialRadialProfileGauss* clone(void) const;
    virtual std::string                      classname(void) const;
    virtual std::string                      type(void) const;
    virtual double                           theta_min(void) const;
    virtual double                           theta_max(void) const;
    virtual void                             read(const GXmlElement& xml);
    virtual void                             write(GXmlElement& xml) const;

    // Other methods
    double  sigma(void) const;
    void    sigma(const double& sigma);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialProfileGauss class extension
 ***************************************************************************/
%extend GModelSpatialRadialProfileGauss {
    GModelSpatialRadialProfileGauss copy() {
        return (*self);
    }
};
