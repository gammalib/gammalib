/***************************************************************************
 *           GModelSpatialDiffuseMap.i - Spatial map model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
    GModelSpatialDiffuseMap(const GFilename& filename,
                            const double&    value = 1.0,
                            const bool&      normalize = true);
    GModelSpatialDiffuseMap(const GSkyMap& map,
                            const double&  value = 1.0,
                            const bool&    normalize = true);
    GModelSpatialDiffuseMap(const GModelSpatialDiffuseMap& model);
    virtual ~GModelSpatialDiffuseMap(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GModelSpatialDiffuseMap* clone(void) const;
    virtual std::string              classname(void) const;
    virtual double                   eval(const GPhoton& photon,
                                          const bool& gradients = false) const;
    virtual GSkyDir                  mc(const GEnergy& energy,
                                        const GTime& time,
                                        GRan& ran) const;
    virtual double                   mc_norm(const GSkyDir& dir,
                                             const double&  radius) const;
    virtual bool                     contains(const GSkyDir& dir,
                                              const double&  margin = 0.0) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;

    // Overloaded base class methods
    virtual double flux(const GSkyRegion& region,
                        const GEnergy&    srcEng  = GEnergy(),
                        const GTime&      srcTime = GTime()) const;

    // Other methods
    double                  value(void) const;
    void                    value(const double& value);
    const GFilename&        filename(void) const;
    void                    load(const GFilename& filename);
    const GSkyMap&          map(void) const;
    void                    map(const GSkyMap& map);
    bool                    normalize(void) const;
    void                    mc_cone(const GSkyRegionCircle& cone) const;
    const GSkyRegionCircle& mc_cone(void) const;
};


/***********************************************************************//**
 * @brief GModelSpatialDiffuseMap class extension
 ***************************************************************************/
%extend GModelSpatialDiffuseMap {
    GModelSpatialDiffuseMap copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml, self.mc_cone())
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
        self.mc_cone(state[1])
}
};
