/***************************************************************************
 *           GModelSpatialDiffuseMap.i - Spatial map model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialDiffuseMap.i
 * @brief Spatial map model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialDiffuseMap.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialDiffuseMap
 *
 * @brief Spatial map model
 ***************************************************************************/
class GModelSpatialDiffuseMap  : public GModelSpatialDiffuse {
public:
    // Constructors and destructors
    GModelSpatialDiffuseMap(void);
    explicit GModelSpatialDiffuseMap(const GXmlElement& xml);
    GModelSpatialDiffuseMap(const std::string& filename,
                            const double&      value = 1.0,
                            const bool&        normalize = true);
    GModelSpatialDiffuseMap(const GSkyMap& map,
                            const double&  value = 1.0,
                            const bool&    normalize = true);
    GModelSpatialDiffuseMap(const GModelSpatialDiffuseMap& model);
    virtual ~GModelSpatialDiffuseMap(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialDiffuseMap* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual double                   eval(const GPhoton& photon) const;
    virtual double                   eval_gradients(const GPhoton& photon) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime& time,
                                        GRan& ran) const;
    virtual double                   mc_norm(const GSkyDir& dir,
                                             const double&  radius) const;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;

    // Other methods
    double             value(void) const;
    void               value(const double& value);
    const std::string& filename(void) const;
    void               load(const std::string& filename);
    const GSkyMap&     map(void) const;
    void               map(const GSkyMap& map);
    bool               normalize(void) const;
    void               set_mc_cone(const GSkyDir& centre,
                                   const double&  radius) const;
};


/***********************************************************************//**
 * @brief GModelSpatialDiffuseMap class extension
 ***************************************************************************/
%extend GModelSpatialDiffuseMap {
    GModelSpatialDiffuseMap copy() {
        return (*self);
    }
};
