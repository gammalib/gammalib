/***************************************************************************
 *     GModelSpatialRadialGauss.i - Radial Gaussian source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialRadialGauss.i
 * @brief Radial Gaussian model class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRadialGauss.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialRadialGauss
 *
 * @brief Radial Gaussian model class
 ***************************************************************************/
class GModelSpatialRadialGauss : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialGauss(void);
    GModelSpatialRadialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialRadialGauss(const GXmlElement& xml);
    GModelSpatialRadialGauss(const GModelSpatialRadialGauss& model);
    virtual ~GModelSpatialRadialGauss(void);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialRadialGauss* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const double&  theta,
                                           const GEnergy& energy,
                                           const GTime&   time,
                                           const bool&    gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual double                    theta_max(void) const;
    virtual GSkyRegion*               region(void) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    double  sigma(void) const;
    void    sigma(const double& sigma);
};


/***********************************************************************//**
 * @brief GModelSpatialRadialGauss class extension
 *
 * The eval(GPhoton&, bool&) method need to be defined in the extension to
 * force swig to build also the interface for these methods that are
 * implemented in the base class only. It's not clear to me why these
 * methods are not inherited automatically. Maybe this could also be handled
 * by a %typemap(typecheck) construct.
 ***************************************************************************/
%extend GModelSpatialRadialGauss {
    GModelSpatialRadialGauss copy() {
        return (*self);
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
