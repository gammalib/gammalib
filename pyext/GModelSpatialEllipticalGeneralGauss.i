/***************************************************************************
 *   GModelSpatialEllipticalGeneralGauss.i                                 *
 * - Elliptical gauss source model class                                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022- by Luigi Tibaldo                                   *
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
 * @file GModelSpatialEllipticalGeneralGauss.i
 * @brief Elliptical generalised gauss model class interface definition
 * @author Luigi Tibaldo
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialEllipticalGeneralGauss.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialEllipticalGeneralGauss
 *
 * @brief Elliptical generalised gauss model class
 ***************************************************************************/
class GModelSpatialEllipticalGeneralGauss : public GModelSpatialElliptical {

public:
    // Constructors and destructors
    GModelSpatialEllipticalGeneralGauss(void);
    GModelSpatialEllipticalGeneralGauss(const GSkyDir& dir,
					const double&  semiminor,
					const double&  semimajor,
					const double&  posangle,
					const double&  ridx);
    explicit GModelSpatialEllipticalGeneralGauss(const GXmlElement& xml);
    GModelSpatialEllipticalGeneralGauss(const GModelSpatialEllipticalGeneralGauss& model);
    virtual ~GModelSpatialEllipticalGeneralGauss(void);

    // Implemented pure virtual base class methods
    virtual void                          clear(void);
    virtual GModelSpatialEllipticalGauss* clone(void) const;
    virtual std::string                   classname(void) const;
    virtual double                        eval(const double&  theta,
                                               const double&  posangle,
                                               const GEnergy& energy,
                                               const GTime&   time,
                                               const bool&    gradients = false) const;
    virtual GSkyDir                       mc(const GEnergy& energy,
                                             const GTime& time,
                                             GRan& ran) const;
    virtual bool                          contains(const GSkyDir& dir,
                                                   const double&  margin = 0.0) const;
    virtual double                        theta_max(void) const;
    virtual void                          read(const GXmlElement& xml);
    virtual void                          write(GXmlElement& xml) const;

    // Other methods
    double  ridx(void) const;
    void    ridx(const double& ridx);
};


/***********************************************************************//**
 * @brief GModelSpatialEllipticalGeneralGauss class extension
 ***************************************************************************/
%extend GModelSpatialEllipticalGeneralGauss {
    GModelSpatialEllipticalGeneralGauss copy() {
        return (*self);
    }
    double eval(const GPhoton& photon) const {
        return self->GModelSpatialElliptical::eval(photon);
    }
    double eval(const GPhoton& photon, const bool& gradients) const {
        return self->GModelSpatialElliptical::eval(photon, gradients);
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
