/***************************************************************************
 *          GSPIModelDataSpace.i - INTEGRAL/SPI data space model           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIModelDataSpace.i
 * @brief INTEGRAL/SPI data space model interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSPIModelDataSpace.hpp"
%}


/***********************************************************************//**
 * @class GSPIModelDataSpace
 *
 * @brief INTEGRAL/SPI data space model
 *
 * This class implements an INTEGRAL/SPI data space model that is based on
 * model values that are found in the INTEGRAL/SPI event cube.
 ***************************************************************************/
class GSPIModelDataSpace : public GModelData {

public:
    // Constructors and destructors
    GSPIModelDataSpace(void);
    GSPIModelDataSpace(const GSPIObservation& obs,
                       const std::string&     name,
                       const std::string&     method,
                       const int&             index);
    explicit GSPIModelDataSpace(const GXmlElement& xml);
    GSPIModelDataSpace(const GSPIModelDataSpace& model);
    virtual ~GSPIModelDataSpace(void);

    // Operators
    virtual GSPIModelDataSpace& operator=(const GSPIModelDataSpace& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GSPIModelDataSpace* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         type(void) const;
    virtual bool                is_constant(void) const;
    virtual double              eval(const GEvent&       event,
                                     const GObservation& obs,
                                     const bool&         gradients = false) const;
    virtual double              npred(const GEnergy&      obsEng,
                                      const GTime&        obsTime,
                                      const GObservation& obs) const;
    virtual GSPIEventCube*      mc(const GObservation& obs, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GSPIModelDataSpace class extension
 ***************************************************************************/
%extend GSPIModelDataSpace {
    GSPIModelDataSpace copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        if state[0][0].elements('parameter') == 0:
            self.__init__()
        else:
            self.__init__(state[0][0])
}
};
