/***************************************************************************
 *     GModelSpatialRadialGeneralGauss.i - Generalised radial Gaussian     *
 *                           source model class                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021-2022 by Luigi Tibaldo                               *
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
 * @file GModelSpatialRadialGeneralGauss.i
 * @brief Generalised radial Gaussian model class Python interface definition
 * @author Luigi Tibaldo
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialGeneralGauss.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialRadialGeneralGauss
 *
 * @brief Generalised radial Gaussian model class
 ***************************************************************************/
class GModelSpatialRadialGeneralGauss : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialGeneralGauss(void);
    GModelSpatialRadialGeneralGauss(const GSkyDir&     dir,
                                    const double&      radius,
                                    const double&      ridx,
                                    const std::string& coordsys = "CEL");
    explicit GModelSpatialRadialGeneralGauss(const GXmlElement& xml);
    GModelSpatialRadialGeneralGauss(const GModelSpatialRadialGeneralGauss& model);
    virtual ~GModelSpatialRadialGeneralGauss(void);

    // Implemented pure virtual methods
    virtual void                             clear(void);
    virtual GModelSpatialRadialGeneralGauss* clone(void) const;
    virtual std::string                      classname(void) const;
    virtual double                           eval(const double&  theta,
                                                  const GEnergy& energy,
                                                  const GTime&   time,
                                                  const bool&    gradients = false) const;
    virtual GSkyDir                          mc(const GEnergy& energy,
                                                const GTime& time,
                                                GRan& ran) const;
    virtual bool                             contains(const GSkyDir& dir,
                                                      const double&  margin = 0.0) const;
    virtual double                           theta_max(void) const;
    virtual void                             read(const GXmlElement& xml);
    virtual void                             write(GXmlElement& xml) const;

    // Other methods
    double radius(void) const;
    double ridx(void) const;
    void   radius(const double& radius);
    void   ridx(const double& ridx);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialGeneralGauss class extension
 *
 * The eval(GPhoton&, bool&) method need to be defined in the extension to
 * force swig to build also the interface for these methods that are
 * implemented in the base class only. It's not clear to me why these
 * methods are not inherited automatically. Maybe this could also be handled
 * by a %typemap(typecheck) construct.
 ***************************************************************************/
%extend GModelSpatialRadialGeneralGauss {
    GModelSpatialRadialGeneralGauss copy() {
        return (*self);
    }
    double eval(const GPhoton& photon) const {
        return self->GModelSpatialRadial::eval(photon);
    }
    double eval(const GPhoton& photon, const bool& gradients) const {
        return self->GModelSpatialRadial::eval(photon, gradients);
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
