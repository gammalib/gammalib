/***************************************************************************
 *        GModelSpatialComposite.i - Spatial point source model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Domenico Tiziani                            *
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
 * @file GModelSpatialComposite.i
 * @brief Spatial composite model class Python interface
 * @author Domenico Tiziani
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialComposite.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialComposite
 *
 * @brief Spatial composite model
 ***************************************************************************/
class GModelSpatialComposite  : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialComposite(void);
    explicit GModelSpatialComposite(const GXmlElement& xml);
    GModelSpatialComposite(const GModelSpatialComposite& model);
    virtual ~GModelSpatialComposite(void);

    // Implemented virtual methods
    virtual void                    clear(void);
    virtual GModelSpatialComposite* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const GPhoton& photon,
                                         const bool& gradients = false) const;
    virtual GSkyDir                 mc(const GEnergy& energy,
                                       const GTime& time,
                                       GRan& ran) const;
    virtual double                  mc_norm(const GSkyDir& dir,
                                            const double&  radius) const;
    virtual bool                    contains(const GSkyDir& dir,
                                             const double&  margin = 0.0) const;
    virtual GSkyRegion*             region(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    int                  components(void) const;
    void                 append(const GModelSpatial& component,
                                const std::string&   name = "",
                                const GModelPar&     par = GModelPar("", 1.0));
    const GModelSpatial* component(const int& index) const;
    const GModelSpatial* component(const std::string& name) const;
    double               scale(const int& index) const;
    double               sum_of_scales(void) const;
};


/***********************************************************************//**
 * @brief GModelSpatialComposite class extension
 ***************************************************************************/
%extend GModelSpatialComposite {
    GModelSpatialComposite copy() {
        return (*self);
    }
};
