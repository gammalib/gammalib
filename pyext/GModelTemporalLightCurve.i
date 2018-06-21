/***************************************************************************
 *      GModelTemporalLightCurve.i - Temporal light curve model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2018 by Juergen Knoedlseder                         *
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
 * @file GModelTemporalLightCurve.i
 * @brief Light curve model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalLightCurve.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalLightCurve
 *
 * @brief Light curve model class
 ***************************************************************************/
class GModelTemporalLightCurve : public GModelTemporal {

public:
    // Constructors and destructors
    GModelTemporalLightCurve(void);
    explicit GModelTemporalLightCurve(const GXmlElement& xml);
    GModelTemporalLightCurve(const GFilename& filename,
                             const double&    norm = 1.0);
    GModelTemporalLightCurve(const GModelTemporalLightCurve& model);
    virtual ~GModelTemporalLightCurve(void);

    // Implemented virtual base class methods
    virtual void                      clear(void);
    virtual GModelTemporalLightCurve* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GTime& srcTime,
                                           const bool& gradients = false) const;
    virtual GTimes                    mc(const double& rate, const GTime& tmin,
                                         const GTime& tmax, GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    const GFilename& filename(void) const;
    void             filename(const GFilename& filename);
    double           norm(void) const;
    void             norm(const double& norm);
};


/***********************************************************************//**
 * @brief GModelTemporalLightCurve class extension
 ***************************************************************************/
%extend GModelTemporalLightCurve {
    GModelTemporalLightCurve copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = xml,
        return state
    def __setstate__(self, state):
        self.__init__(state[0])
}
};
