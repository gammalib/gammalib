/***************************************************************************
 *       GModelSpatialPointSource.i - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2022 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialPointSource.i
 * @brief Point source spatial model class Python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialPointSource.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialPointSource
 *
 * @brief Point source spatial model
 ***************************************************************************/
class GModelSpatialPointSource  : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialPointSource(void);
    GModelSpatialPointSource(const GSkyDir&     dir,
                             const std::string& coordsys = "CEL");
    GModelSpatialPointSource(const double&      lon,
                             const double&      lat,
                             const std::string& coordsys = "CEL");
    explicit GModelSpatialPointSource(const GXmlElement& xml);
    GModelSpatialPointSource(const GModelSpatialPointSource& model);
    virtual ~GModelSpatialPointSource(void);

    // Implemented virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialPointSource* clone(void) const;
    virtual std::string               classname(void) const;
    virtual double                    eval(const GPhoton& photon,
                                           const bool& gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual double                    mc_norm(const GSkyDir& dir,
                                              const double&  radius) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Overloaded base class methods
    virtual double flux(const GSkyRegion& region,
                        const GEnergy&    srcEng  = GEnergy(),
                        const GTime&      srcTime = GTime()) const;

    // Other methods
    std::string    coordsys(void) const;
    const GSkyDir& dir(void) const;
    void           dir(const GSkyDir& dir);
};


/***********************************************************************//**
 * @brief GModelSpatialPointSource class extension
 ***************************************************************************/
%extend GModelSpatialPointSource {
    GModelSpatialPointSource copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = {'type': self.type(), 'dir': self.dir(), 'coordsys': self.coordsys()}
        return state
    def __setstate__(self, state):
        self.__init__(state['dir'],state['coordsys'])
        self.type(state['type'])
}
};
